function Kondo_Spin1_NRG(parfn,varargin)

    try % in case of bug
        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [PE, Nkeep, J0, T, JobName, nz] = loadvar(partot, ...
        {'PE', 'Nkeep', 'J0', 'T', 'JobName', 'nz'}, ...
            {[], [], [], [], [], []});
    
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
        D = 1;
        Delta = pi;
        ozin = [-1;1]*D;
        RhoV2in = [1;1]*(Delta/pi);
    
        % NRG parameters
        Lambda = 2.5;
        N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
        Etrunc = [];
        ETRUNC = [];
        emin = T;

        [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N,nz,'Nfit',round(-2*log(1e-8)/log(Lambda)));
    
        % local operators
        [S,IF] = getLocalSpace('Spin',1);
        [S_imp,E_imp] = setItag('L00','op',S,IF.E);

        S_imp = addSymmetry(S_imp, 'A', 'pos', 1);
        E_imp = addSymmetry(E_imp, 'A', 'pos', 1);

        [FF,ZF,J_sp,IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1);
        [FF,ZF,J_sp,EF] = setItag('s00','op',FF(:),ZF,J_sp(:),IF.E);
    
        % local isometry and Hamiltonian
        A0 = getIdentity(E_imp,2,EF,2,'K00*',[1,3,2]);
    
        H0 = J0*contract(A0,'!2*',{J_sp,'!1',{S_imp,'!2*',A0}});
        H0 = H0 + 1e-40*contract(A0,'!2*',A0);
    
        % operators that define the two-point correlators
        Ops1 = [S_imp,J_sp];
        Ops2 = Ops1;
        OpNames = {'ImpSp','BathSp'};
    
        STG = ['/data/',getenv('USER'),'/Kondo_Spin1/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep)];
    
        nrgdata = cell(1,nz);

        for itz = 1:nz
            nrgdata{itz} = NRG_SL([],H0,A0,Lambda,ff{itz}(2:end),FF,ZF,'Nkeep',Nkeep,'deps',1e-10);
            nrgdata{itz} = getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);   % calculating the full density matrix(FDM)
        end

        [Etot,Qtot,Qdiff] = plotE(nrgdata{1},'Emax',10,'legmax',25);       % Data for Eflow diagram    
        save([STG,'/Etot.mat'],'Etot');
        save([STG,'/Qtot.mat'],'Qtot');

        zflag = zeros(1,numel(Ops1));     
        cflag = (zflag-0.5)*2;

        Adiscs = cell(numel(Ops1),nz); % discrete data
        Aconts = cell(1,size(Adiscs,1)); % continuous (i.e., broadened) spectral function
        
        for itz = 1:nz
            [odisc,Adiscs(:,itz),sigmak] = getAdisc(nrgdata{itz},Ops1,Ops2,ZF,'Z_L00',E_imp,'zflag',zflag,'cflag',cflag,'emin',emin);
            % odisc: [Numeric vector] center values of the frequency grid
            % sigmak [Numeric vector] center values of the broadening width bins
            % Adiscs: [cell array of numeric matrix] spectral function binned along odisc
            % Adiscs(:,itz): cell vector of numeric matrix. 
            % nth matrix Adiscs(n,itz) corresponds to spectral function of Ops1(n) and Ops2(n) binned along odisc and sigmak
            % Length of cell vector Adisc(:,itz): numel(Ops1) = numel(Ops2)
        end

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

        % Calculate the impurity contribution to entropy
        [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N+10,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));

        A02 = getIdentity(setItag('L00',getvac(E_imp)),2,EF,2,'K00*',[1 3 2]);
        H02 = contract(A02,'!2*',A02) + 1e-40*getIdentity(A02,2);

        bathNRG_data = NRG_SL([],H02,A02,Lambda,ff{1}(2:end),FF,ZF,'Nkeep',Nkeep,'deps',1e-10);

        beta = 1:0.1:2;
        EntData.beta = beta;
        EntData.Temps = cell(1,numel(beta));
        EntData.S_imp = cell(1,numel(beta));
    
        for itN = 1:numel(beta)
            [Temps,~,~,Sent_bath,~,~,~] = getTDconv(bathNRG_data,'useT','beta',beta(itN));
            [Temps,~,~,Sent,~,~,~] = getTDconv(nrgdata{1},'useT','beta',beta(itN));

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