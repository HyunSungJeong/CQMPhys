function Anderson_4flav_NRG(parfn,varargin)

    try % in case of bug

        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [PE, Nkeep, Hyb, U, J, N0, T, JobName, nz] = loadvar(partot, ...
            {'PE', 'Nkeep', 'Hyb', 'U', 'J', 'N0', 'T', 'JobName', 'nz'}, ...
                {[], [], [], [], [], [], [], [], []});
    
        %num_threads_SL(8);
        
        disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
        disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
        disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
        disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
        disp2(['LOG_FILE : ', getenv('LOG_FILE')]);
    
        disp2(sprintf('%d',nz));
    
        
        % Impurity Hamiltonian parameters
        % U : Hubbard U
        % J : inter-valley exchange interaction
        % epsilon : on-site level
    
        % T : Temperature
        % Hyb : hybridization strength
        D = 1;                  % half-bandwidth
        ozin = [-1;1]*D;        % input frequency grid on which the hybridization function is evaluated    
        RhoV2in = [1;1]*Hyb;    % hybridization function
    
        % NRG parameters
        Lambda = 4;
        N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
        Etrunc = [];
        ETRUNC = [];

        % Logarithmic discretization of the hybridization function
        [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N,nz,'Nfit',round(-2*log(1e-8)/log(Lambda)));


        %% Define local operators
        [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge(:),SU2spin','NC',2);
        Fcell = {FF(1), FF(2), ZF, SF, IF.E; FF(1), FF(2), ZF, SF, IF.E};
        sym = {'A','A','SU2'};
        symnum = [3,4,6; 1,2,5];
        tagname = {'o1', 'o2'};

        for it1 = 1:size(Fcell, 2)          % F,Z,S,I, respectively
            for it2 = 1:size(symnum, 1)     % 2,4 to F1, 1,3 to F2
                for it3 = 1:numel(sym)      
                    Fcell{it2, it1} = addSymmetry(Fcell{it2, it1}, sym{it3}, 'pos', symnum(it2, it3));
                end

                Fcell{it2, it1} = setItag(tagname{it2}, 'op', Fcell{it2, it1});
            end
        end

        A = getIdentity(Fcell{1, 5}, 2, Fcell{2, 5}, 2, 's00');
        FF = QSpace(1, 4); ZF = QSpace(1, 2); SF = QSpace(1, 2); NF = QSpace(1, 2);
        FF(1) = contract(A, '!3*', {Fcell{1, 1}, '!1', {A, Fcell{2, 3}}}, [1 3 2]);
        FF(2) = contract(A, '!3*', {Fcell{1, 2}, '!1', {A, Fcell{2, 3}}}, [1 3 2]);
        FF(3) = contract(A, '!3*', {Fcell{2, 1}, '!1', A}, [1 3 2]);
        FF(4) = contract(A, '!3*', {Fcell{2, 2}, '!1', A}, [1 3 2]);
        % FF: orb1 valley+, orb1 valley-, orb2 valley+, orb2 valley-

        for it1 = 1:size(Fcell, 1)
            ZF(it1) = contract(A, '!3*', {Fcell{it1, 3}, '!1', A});
            SF(it1) = contract(A, '!3*', {Fcell{it1, 4}, '!1', A}, [1 3 2]);
        end
        Stot = SF(1) + SF(2);         % total spin

        for it1 = 1:numel(FF)
            NF(it1) = quadOp(FF(it1), FF(it1), []);
        end
        EF = contract(A,'!3*',A);

        % interacting term in the impurity Hamiltonian
        HU = (U/2)*sum(NF)*sum(NF);
        HU = HU + J_S*(contract(Stot,'!2*',Stot) - contract(SF(1),'!2*',SF(1)) - contract(SF(2),'!2*',SF(2)));
        HU = HU - J_L*(contract(L_z,'!2*',L_z) - sum(NF));

        % quadratic term in the impurity Hamiltonian
        Hmu = - U*N0*sum(NF);

        % isometry connecting the impurity and the first bath site
        A0 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);
        HU = contract(A0,'!2*',{HU,'!1',A0});           % HU in the K00 space
        Hmu = contract(A0,'!2*',{Hmu,'!1',A0});         % Hmu in the K00 space
        H0 = HU + Hmu;              % impurity Hamiltonian

        STG = ['/data/',getenv('USER'),'/Quartic/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep)];

        nrgdata = cell(1,nz);

        %% NRG & susceptibility calculations
        for itz = 1:nz
            nrgdata{itz} = NRG_SL([],H0,A0,Lambda,ff{itz},FF,ZF,gg{itz},NF,'Nkeep',Nkeep, ...
                                        'Etrunc',Etrunc,'ETRUNC',ETRUNC,'dff',dff{itz},'dgg',dgg{itz},'deps',1e-10);
            nrgdata{itz} = getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);         % calculating the full density matrix(FDM)

            if itz == nz
                plotE(nrgdata,'title',symstr{its});
            end

            %[ocont,SE,nrgdataz,ImSEle] = impSE (Hmu,HU,A0,FF,ZF,RhoV2in,T,Lambda,nz,Nkeep,nrgdata);

        end

    catch Err
        disp2(getReport(Err));
        rethrow(Err);
    end
    
end