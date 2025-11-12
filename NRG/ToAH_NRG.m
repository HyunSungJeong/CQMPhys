function ToAH_NRG(parfn,varargin)

    try % in case of bug

        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});

        [PE, Nkeep, Lambda, Hyb, U, J, J_L, T, mu, JobName, getSusc, nz, N_SuscIter] = loadvar(partot, ...
        {'PE', 'Nkeep', 'Lambda', 'Hyb', 'U', 'J', 'J_L', 'T', 'mu', 'JobName', 'getSusc', 'nz', 'N_SuscIter'}, ...
            {[], [], [], [], [], [], [], [], [], [], [], [], []});

        %num_threads_SL(8);
        
        disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
        disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
        disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
        disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
        disp2(['LOG_FILE : ', getenv('LOG_FILE')]);

        disp2(sprintf('%d',nz));

        % Impurity Hamiltonian parameters
        % U : Hubbard interaction
        % J : spin Hund interaction
        % J_L : orbital-z Hund interaction
        % mu : chemical potential

        % T : Temperature
        D = 1;
        ozin = [-1;1]*D;
        RhoV2in = [1;1]*Hyb;

        % NRG parameters
        % nz : number of z-shifts
        N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
        Etrunc = [];
        ETRUNC = [];
        Nfit = round(-2*log(1e-8)/log(Lambda));
        emin = T;           % Minimum absolute value of frequency grid
        emax = 1e3;         % Maximum absolute value of frequency grid
        estep = 250;        % # of steps to increase frequency per decade(x10), used in getAdisc and getAcont

        % Logarithmic discretization of the hybridization function
        [ff, gg, dff, dgg] = doZLD(ozin, RhoV2in, Lambda, N, nz, 'Nfit', Nfit);

        %% Define local operators

        % local operators
        [FF,ZF,SF,IF] = getLocalSpace('FermionS','Acharge(:),SU2spin','NC',2);
        NF = quadOp(FF,FF,[]);

        [FF,ZF,SF,EF,NF] = setItag('s00','op',FF(:),ZF,SF(:),IF.E,NF);
        Ntot = sum(NF);
        L_z = (NF(1) - NF(2)) / 2;
        L_plus = quadOp(FF(1), FF(2), [+1,-1,0]);
        L_minus = quadOp(FF(2), FF(1), [-1,+1,0]);

        % local isometry and Hamiltonian
        A0 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);

        H0 = (U-3*J/2)*Ntot*(Ntot-EF)/2;
        H0 = H0 - J*contract(SF,'!2*',SF) + (3/4)*J*Ntot;
        H0 = H0 - J_L*contract(L_z,'!2*',L_z);
        H0 = H0 - mu*Ntot;

        H0 = contract(A0,'!2*',{H0,'!1',A0});   % H0 in the K00 space
        H0 = H0 + 1e-40*getIdentity(A0,2);
        
        STG = ['/data/',getenv('USER'),'/ToAH/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];


        %% NRG calculation
        for itz = 1:nz

            nrgdata{itz} = ['/tmp/hyunsung/NRGdata_nz=',sprintf('%d',itz)];
            NRG_SL(nrgdata{itz},H0,A0,Lambda,ff{itz},FF,ZF,gg{itz},NF,'Nkeep',Nkeep, ...
                                        'Etrunc',Etrunc,'ETRUNC',ETRUNC,'dff',dff{itz},'dgg',dgg{itz},'deps',1e-10);
            getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);         % calculating the full density matrix(FDM)
        end

        for it = 1:numel(FF)
            nu(it) = getEpVal(nrgdata{1}, NF(it));
        end
        save([STG,'/nu.mat'],'nu');
                
        [Etot,Qtot,Qdiff] = plotE(nrgdata{1},'Emax',10,'legmax',25);       % Data for Eflow diagram    
        save([STG,'/Etot.mat'],'Etot');
        save([STG,'/Qtot.mat'],'Qtot');
        

        %% Compute dynamic susceptibilities

        if getSusc      %% Compute dynamic susceptibilities

            SuscOps_all = [NF(:); Ntot; SF(:); L_z; L_plus];
            OpNames_all = {'ImpCharge_plus', 'ImpCharge_minus', 'ImpCharge_tot', 'ImpSp', 'ImpOrb_z', 'ImpOrb_plus'};

            N_calc = floor(numel(SuscOps_all)/N_SuscIter);
            N_rem = rem(numel(SuscOps_all), N_SuscIter);

            N_calc_cnt = 0;

            for itN = 1:N_SuscIter

                clear Adiscs_Susc
                clear Aconts_Susc

                if itN <= N_rem
                    SuscOps1 = SuscOps_all(N_calc_cnt+1 : N_calc_cnt+N_calc+1);
                    OpNames = OpNames_all(N_calc_cnt+1 : N_calc_cnt+N_calc+1);
                    N_calc_cnt = N_calc_cnt + N_calc + 1;
                else
                    SuscOps1 = SuscOps_all(N_calc_cnt+1 : N_calc_cnt+N_calc);
                    OpNames = OpNames_all(N_calc_cnt+1 : N_calc_cnt+N_calc);
                    N_calc_cnt = N_calc_cnt + N_calc;
                end
                
                SuscOps2 = SuscOps1;

                zflag = zeros(1,numel(SuscOps1));
                cflag = (zflag-0.5)*2;

                Adiscs_Susc = cell(numel(SuscOps1),nz); % discrete data
                Aconts_Susc = cell(1,size(Adiscs_Susc,1)); % continuous (i.e., broadened) spectral function
                
                for itz = (1:nz)
                    [odisc,Adiscs_Susc(:,itz),sigmak] = getAdisc(nrgdata{itz}, SuscOps1, SuscOps2, ZF, 'zflag', zflag, 'cflag', cflag, ...
                                                        'emin', emin, 'emax', emax, 'estep', 2*estep);
                end
                
                % file path to save Aconts
                SaveAconts = cellfun(@(x) [STG,'/NRG_Op=',OpNames{x},'.mat'], num2cell(1:numel(SuscOps1)), 'UniformOutput', false);
                SaveAdiscs = cellfun(@(x) [STG,'/Adiscs_NRG_Op=',OpNames{x},'_nz=',sprintf('%d',nz),'.mat'], num2cell(1:numel(SuscOps1)), 'UniformOutput', false);
                
                % Calculate dynamic susceptibilities
                for ita = (1:size(Adiscs_Susc,1))

                    Adisc_Susc = mean(cell2mat(reshape(Adiscs_Susc(ita,:),[1 1 nz])),3);
                    [ocont, Aconts_Susc{ita}] = getAcont(odisc, Adisc_Susc, sigmak,T/5, 'alphaz', 1/nz, ...
                                                            'emin', emin, 'emax', emax, 'estep', estep);

                    temp = Aconts_Susc{ita};
                    save(SaveAconts{ita},'temp');
                    save([STG,'/ocont.mat'],'ocont');
                end

            end % itN
        end % getSusc

        %% calculate impurity contribution to entropy

        [ff, ~] = doZLD(ozin,RhoV2in,Lambda,N+10,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));

        A02 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);
        H02 = contract(A02, '!2*', {1e7 * Ntot, '!1', A02}) +1e-40 * getIdentity(A02, 2);

        bathNRG_data = NRG_SL([], H02, A02, Lambda, ff{1}, FF, ZF, 'Nkeep', Nkeep, 'deps', 1e-10);

        beta = 0.5:0.1:2;
        EntData.beta = beta;
        EntData.Temps = cell(1,numel(beta));
        EntData.S_imp = cell(1,numel(beta));


        for itN = 1:numel(beta)

            [Temps,~,~,Sent_bath,~,~,~] = getTDconv(bathNRG_data,'useT','beta',beta(itN));
            [~,~,~,Sent,~,~,~] = getTDconv(nrgdata{1},'useT','beta',beta(itN));

            Sent = Sent(Temps > T);
            Sent_bath = Sent_bath(Temps > T);
            Sent_imp = Sent - Sent_bath;

            Temps = Temps(Temps > T);
            
            EntData.Temps{itN} = Temps;
            EntData.S_imp{itN} = Sent_imp;
        end

        save([STG,'/EntData.mat'],'EntData');
    
    catch Err
        disp2(getReport(Err));
        rethrow(Err);
    end

end