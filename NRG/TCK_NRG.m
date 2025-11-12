function TCK_NRG(parfn,varargin)

    try % in case of bug
        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [PE, Nkeep, Lambda, J0, T, JobName, getSusc, getCorr, nz] = loadvar(partot, ...
        {'PE', 'Nkeep', 'Lambda', 'J0', 'T', 'JobName', 'getSusc', 'getCorr', 'nz'}, ...
            {[], [], [], [], [], [], [], [], []});
    
        %num_threads_SL(8);
        
        disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
        disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
        disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
        disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
        disp2(['LOG_FILE : ', getenv('LOG_FILE')]);
    
        disp2(sprintf('%d',nz));
    
        % Impurity Hamiltonian parameters
        % J0 : spin-spin exchange coupling
    
        % T : Temperature
        emin = T;
        D = 1;
        Delta = pi;
        ozin = [-1;1]*D;
        RhoV2in = [1;1]*(Delta/pi);
    
        % NRG parameters
        N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
        Etrunc = [];
        ETRUNC = [];
    
        [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N,nz,'Nfit',round(-2*log(1e-8)/log(Lambda)));
    
        % impuity operators
        [~,Z_imp,S_imp,IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1);
        [Z_imp,S_imp,E_imp] = setItag('L00','op',Z_imp,S_imp,IF.E);

        S_imp = getsub(S_imp, find(all(S_imp.Q{1}(:,1) == 0, 2)));

        % bath operators
        [FF,ZF,J_sp,IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',2);
        [FF,ZF,J_sp,EF] = setItag('s00','op',FF(:),ZF,J_sp(:),IF.E);

        % local isometry and Hamiltonian
        A0 = getIdentity(E_imp,2,EF,2,'K00*',[1,3,2]);

        H0 = J0*contract(A0,'!2*',{J_sp,'!1',{S_imp,'!2*',A0}});
        H0 = H0 + 1e-40*contract(A0,'!2*',A0);
    
        % operators that define the two-point correlators
        Ops1 = [S_imp; J_sp];
        Ops2 = Ops1;
        OpNames = {'ImpSp','BathSp'};

        zflag = zeros(1,numel(Ops1));
        cflag = (zflag-0.5)*2;
    
        STG = ['/data/',getenv('USER'),'/TCK/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];

        Adiscs = cell(numel(Ops1),nz); % discrete data
        Aconts = cell(1,size(Adiscs,1)); % continuous (i.e., broadened) spectral function
    
        nrgdata = cell(1,nz);

        %% Perform NRG

        for itz = (1:nz)
            nrgdata{itz} = NRG_SL([],H0,A0,Lambda,ff{itz}(2:end),FF,ZF,'Nkeep',Nkeep,'deps',1e-10);

            if itz == 1
                [Etot,Qtot,Qdiff] = plotE(nrgdata{1},'Emax',10,'legmax',25);       % Data for Eflow diagram

                save([STG,'/Etot.mat'],'Etot');
                save([STG,'/Qtot.mat'],'Qtot');
                save([STG,'/Qdiff.mat'],'Qdiff');
            end

            nrgdata{itz} = getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);   % calculating the full density matrix(FDM)

            if getSusc
                [odisc,Adiscs(:,itz),sigmak] = getAdisc(nrgdata{itz},Ops1,Ops2,ZF,'Z_L00',Z_imp,'zflag',zflag,'cflag',cflag,'emin',emin);
            end
        end

        %% Compute dynamic susceptibilities

        if getSusc

            DiscData.odisc = odisc;
            DiscData.sigmak = sigmak;
            DiscData.Adiscs = Adiscs;
            DiscData.nz = nz;
            DiscData.emin = emin;
            save([STG,'/DiscData.mat'],'DiscData');

            % file path to save Aconts
            SaveAconts = cellfun(@(x) [STG,'/NRG_Op=',OpNames{x},'.mat'], num2cell(1:numel(Ops1)), 'UniformOutput', false);

            % Calculate dynamic susceptibilities
            for ita = (1:size(Adiscs,1))

                Adisc = mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3);

                [ocont, Aconts{ita}] = getAcont(odisc,Adisc,sigmak,T/5,'alphaz',1/nz,'emin',emin);
                % ocont : [numeric vector] Logarithimic frequency grid.
                % Aconts : [numeric vector] Smoothened spectral function.

                temp = Aconts{ita};
                save(SaveAconts{ita},'temp');
            end
            save([STG,'/ocont.mat'],'ocont');
        
        end % if getSusc

       %% calculate correlation functions between impurity and Wilson chain sites
    
        if getCorr

            Sp_corr = cell(N-1,1);        % cell array of spin-spin correlators between the impurity and Wilson chain sites

            Sp_corr_Left = S_imp;
            
            for site = 0:N-2

                disptime(['Calculating correlator between impurity and Wilson chain site #',sprintf('%02d',site)]);
            
                Wilson_sp = setItag(['s',sprintf('%02d',site)],'op',J_sp);      % spin operator on Wilson chain site of interest

                Sp_corrT = contract(Sp_corr_Left, '!1', nrgdata{1}.AT{site+1}, '!2');
                Sp_corrT = contract(Sp_corrT, conj(Wilson_sp));
                Sp_corrT = contract(nrgdata{1}.AT{site+1}, '!2*', Sp_corrT);
            
                Sp_corrK = contract(Sp_corr_Left, '!1', nrgdata{1}.AK{site+1}, '!2');
                Sp_corrK = contract(Sp_corrK, conj(Wilson_sp));
                Sp_corrK = contract(nrgdata{1}.AK{site+1}, '!2*', Sp_corrK);


                Sp_corr{site+1} = 0;
                corrT = contract(Sp_corrT, diag(nrgdata{1}.RhoT{site+1}));
                corrK = contract(Sp_corrK, nrgdata{1}.RhoK{site+1});

                if ~isempty(corrT)
                    Sp_corr{site+1} = Sp_corr{site+1} + corrT.data{1};
                end
                
                if ~isempty(corrK)
                    Sp_corr{site+1} = Sp_corr{site+1} + corrK.data{1};
                end
            

                Sp_corr_Left = contract(Sp_corr_Left,'!1',nrgdata{1}.AK{site+1});
                Sp_corr_Left = contract(nrgdata{1}.AK{site+1},'!2*',Sp_corr_Left,[1,3,2]);
                
            end

            Sp_corr = cell2mat(Sp_corr);

            save([STG,'/Spin_correlators.mat'],'Sp_corr');

        end % if getCorr


        %% calculate impurity contribution to entropy

        [ff, ~] = doZLD(ozin,RhoV2in,Lambda,N+10,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));

        A02 = getIdentity(setItag('L00',getvac(E_imp)),2,EF,2,'K00*',[1 3 2]);
        H02 = contract(A02,'!2*',A02) + 1e-40*getIdentity(A02,2);
        bathNRG_data = NRG_SL([],H02,A02,Lambda,ff{1}(2:end),FF,ZF,'Nkeep',Nkeep,'deps',1e-10);
            
        beta = 0.5:0.1:2;
        EntData.beta = beta;
        EntData.Temps = cell(1,numel(beta));
        EntData.S_imp = cell(1,numel(beta));

        for itN = 1:numel(beta)
            [Temps,~,~,Sent_bath,~,~,~] = getTDconv(bathNRG_data,'useT','beta',beta(itN));
            [~,~,~,Sent,~,~,~] = getTDconv(nrgdata{1},'useT','beta',beta(itN));
            
            Sent = Sent(Temps > T);
            Sent_bath = Sent_bath(Temps > T);
            Temps = Temps(Temps > T);

            Sent_imp = Sent - Sent_bath;
            
            EntData.Temps{itN} = Temps;
            EntData.S_imp{itN} = Sent_imp;
        end

        save([STG,'/EntData.mat'],'EntData');

    catch Err
        disp2(getReport(Err));
        rethrow(Err);
    end
    
end