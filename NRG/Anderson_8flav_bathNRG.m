function Anderson_8flav_bathNRG(parfn,varargin)

    try % in case of bug

        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [PE, Nkeep, Lambda, Hyb, T] = loadvar(partot, ...
            {'PE', 'Nkeep', 'Lambda', 'Hyb', 'T'}, ...
                {[], [], [], [], []});
    
        num_threads_SL(8);
        
        disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
        disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
        disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
        disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
        disp2(['LOG_FILE : ', getenv('LOG_FILE')]);
    
        % T : Temperature
        % Hyb : hybridization strength
        D = 1;                  % half-bandwidth
        ozin = [-1;1]*D;        % input frequency grid on which the hybridization function is evaluated    
        RhoV2in = [1;1]*Hyb;    % hybridization function
    
        % NRG parameters
        % Lambda
        N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
        Etrunc = [];
        ETRUNC = [];
        Nfit = round(-2*log(1e-8)/log(Lambda));
        emin = T;           % Minimum absolute value of frequency grid
        emax = 1e3;         % Maximum absolute value of frequency grid
        estep = 250;        % # of steps to increase frequency per decade(x10), used in getAdisc and getAcont

        %% Define local operators

        %% Construct operators for individual valley degrees of freedom
        [FF, ZF, SF, IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1);
        Fcell = {FF, ZF, SF, IF.E; FF, ZF, SF, IF.E};
        symnum = [2; 1];
        sym = {'A'};
        tagnames = {'v+','v-'};

        for it1 = 1:size(Fcell,2)           % F,Z,S,I, respectively
            for it2 = 1:size(symnum,1)      % 2,4 to F1, 1,3 to F2
                for it3 = 1:numel(sym)
                    Fcell{it2,it1} = addSymmetry(Fcell{it2,it1}, sym{it3}, 'pos', symnum(it2,it3));
                end

                Fcell{it2,it1} = setItag(tagnames{it2}, 'op', Fcell{it2, it1});
            end
        end

        A = getIdentity(Fcell{1,4},2,Fcell{2,4},2,'s00');
        FF = QSpace(1,2); SF = QSpace(1,2); NF = QSpace(1,2);
        
        FF(1) = contract(A,'!3*',{Fcell{1,1},'!1',{A,Fcell{2,2}}}, [1,3,2]);    % annihilation operator for valley index +. Note that the fermionic sign operator is contracted
        FF(2) = contract(A,'!3*',{Fcell{2,1},'!1',A}, [1,3,2]);                 % annihilation operator for valley index -. Fermionic sign operator is not contracted.

        for it1 = 1:size(Fcell, 1)
            SF(it1) = contract(A, '!3*', {Fcell{it1, 3}, '!1', A}, [1 3 2]);
        end
        ZF = contract(A, '!3*', {Fcell{2,2}, '!1', {Fcell{1,2}, '!1', A}});     % fermionic sign operator in the combined space of valley indices + and -
        EF = contract(A, '!3*', A);                                             % identity operator in the combined space of valley indices + and -


        %% Add orbital indices
        Fcell = {FF(1), FF(2), ZF, SF(1), SF(2), EF; ...
                    FF(1), FF(2), ZF, SF(1), SF(2), EF};
        sym = {'A','A','SU2'};
        symnum = [3,4,6; 1,2,5];
        tagname = {'o1', 'o2'};

        for it1 = 1:size(Fcell, 2)          % F,Z,S,I, respectively
            for it2 = 1:size(symnum, 1)     % 3,4,6 to F1, 1,2,5 to F2
                for it3 = 1:numel(sym)      
                    Fcell{it2, it1} = addSymmetry(Fcell{it2, it1}, sym{it3}, 'pos', symnum(it2, it3));
                end

                Fcell{it2, it1} = setItag(tagname{it2}, 'op', Fcell{it2, it1});
            end
        end

        A = getIdentity(Fcell{1, 6}, 2, Fcell{2, 6}, 2, 's00');
        FF = QSpace(1, 4); SF = QSpace(1, 4); NF = QSpace(1, 4);
        FF(1) = contract(A, '!3*', {Fcell{1, 1}, '!1', {A, Fcell{2, 3}}}, [1 3 2]);
        FF(2) = contract(A, '!3*', {Fcell{1, 2}, '!1', {A, Fcell{2, 3}}}, [1 3 2]);
        FF(3) = contract(A, '!3*', {Fcell{2, 1}, '!1', A}, [1 3 2]);
        FF(4) = contract(A, '!3*', {Fcell{2, 2}, '!1', A}, [1 3 2]);
        % FF: orb1(+) valley+, orb1(+) valley-, orb2(-) valley+, orb2(-) valley-

        for it1 = 1:size(Fcell, 1)
            for it2 = 1:2
                SF(2*(it1-1)+it2) = contract(A, '!3*', {Fcell{it1, 3+it2}, '!1', A}, [1 3 2]);
            end
        end

        ZF = contract(A, '!3*', {Fcell{2,3}, '!1', {Fcell{1,3}, '!1', A}});     % fermionic sign operator in the combined space of valley indices + and - and orbital indices 1 and 2
        EF = contract(A,'!3*',A);   % identity operator in the combined space of valley indices + and - and orbital indices 1 and 2

        for it1 = 1:numel(FF)
            NF(it1) = quadOp(FF(it1), FF(it1), []);
        end
        Ntot = sum(NF);             % total charge operator

        %% Define variables and paths for storing calculation results

        node_name = getenv('HOSTNAME');
        disp2(node_name);

        TempAvail = cell(1,14);
        for itN = 1:14
            TempAvail{itN} = ['b',sprintf('%02d',itN)];
        end

        if ismember(node_name,TempAvail)    % memory of compuational node enough to run this NRG calculation

            NRGdataSTG = '/tmp/hyunsung/NRG';
            STG_avail = true;
            if ~exist(NRGdataSTG,'dir')
                mkdir(NRGdataSTG);
            end
        else

            STG_avail = false;
            for itN = 1:6
                tmp = load(['/project/hyunsung/NRG',sprintf('%d',itN),'/Avail.mat']);
                Avail = tmp.Avail;

                if Avail
                    NRGdataSTG = ['/project/hyunsung/NRG',sprintf('%d',itN)];
                    Avail = false;
                    STG_avail = true;
                    save([NRGdataSTG,'/Avail.mat'], 'Avail');
                    break;
                end
            end
        end

        if STG_avail
            disp2(['NRGdata storage directory : ''', NRGdataSTG,'''']);
        else
            error('ERR: There is no available NRGdata storage directory!');
        end

        %% calculate bath contribution to entropy

        [ff, gg] = doZLD(ozin, RhoV2in, Lambda, N+6, 1, 'Nfit', Nfit);

        %A02 = getIdentity(setItag('L00',getvac(EF)), 2, setItag('s00',getvac(EF)), 2, 'K00*', [1 3 2]);
        %H02 = contract(A02,'!2*',A02) + 1e-40*getIdentity(A02,2);

        A02 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);
        H02 = contract(A02, '!2*', {1e7 * Ntot, '!1', A02}) +1e-40 * getIdentity(A02, 2);

        bathNRG_data = [NRGdataSTG,'/bathNRG_data'];
        NRG_SL(bathNRG_data, H02, A02, Lambda, ff{1}, FF, ZF, 'Nkeep', Nkeep, 'deps',1e-10);
        
        beta = 0.5:0.1:2;

        BathEntData.beta = beta;
        BathEntData.Temps = cell(1,numel(beta));
        BathEntData.S_bath = cell(1,numel(beta));

        for itN = 1:numel(beta)
            [Temps,~,~,Sent_bath,~,~,~] = getTDconv(bathNRG_data,'useT','beta',beta(itN));

            BathEntData.Temps{itN} = Temps;
            BathEntData.S_bath{itN} = Sent_bath;
        end

        BathEnt_folder = ['/data/hyunsung/8flav/BathEnt_Nkeep=', sprintf('%d',Nkeep), '_Lambda=', sprintf('%.15g',Lambda), ...
                                '_Hyb=', sprintf('%.15g',Hyb), '_T=', sprintf('%.15g',T)];

        if ~exist(BathEnt_folder, 'dir')
            mkdir(BathEnt_folder);
        end
        save([BathEnt_folder, '/BathEntData.mat'],'BathEntData');

        if ~ismember(node_name,TempAvail)
            Avail = true;
            save([NRGdataSTG,'/Avail.mat'], 'Avail');
        end


    catch Err
        if ~ismember(node_name,TempAvail) && STG_avail
            Avail = true;
            save([NRGdataSTG,'/Avail.mat'], 'Avail');
        end
        disp2(getReport(Err));
        rethrow(Err);
    end
    
end