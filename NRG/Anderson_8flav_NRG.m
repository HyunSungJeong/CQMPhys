function Anderson_8flav_NRG(parfn,varargin)

    try % in case of bug

        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [PE, Nkeep, Lambda, Hyb, U, J, N0, mu, T, JobName, nz, getSE, getSusc, N_SuscIter] = loadvar(partot, ...
            {'PE', 'Nkeep', 'Lambda', 'Hyb', 'U', 'J', 'N0', 'mu', 'T', 'JobName', 'nz', 'getSE', 'getSusc', 'N_SuscIter'}, ...
                {[], [], [], [], [], [], [], [], [], [], [], [], [], []});
    
        num_threads_SL(8);
        
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
        % Lambda
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

        %% Define impurity spin, orbital, and valley operators

        % spin operators
        S_o1 = SF(1) + SF(2);   % total spin operator for orbital 1(+)
        S_o2 = SF(3) + SF(4);   % total spin operator for orbital 2(-)
        Stot = sum(SF);         % total spin

        % orbital pseudospin operators
        Orb_z(1) = contract(FF(1),'!2*',FF(1)) / 2 - contract(FF(3),'!2*',FF(3)) / 2;       % total orbital pseudospin operator for valley index +
        Orb_z(2) = contract(FF(2),'!2*',FF(2)) / 2 - contract(FF(4),'!2*',FF(4)) / 2;       % total orbital pseudospin operator for valley index -
        Orb_tot = Orb_z(1) + Orb_z(2);      % total orbital pseudospin operator

        %Orb_plus(1) = quadOp(FF(1),FF(3),'*');      % orbital+ operator for valley index +
        %Orb_plus(2) = quadOp(FF(2),FF(4),'*');      % orbital+ operator for valley index -
        %Orb_plus_tot = Orb_plus(1) + Orb_plus(2);   % total orbital+ operator

        % valley operators
        V_z(1) = contract(FF(1),'!2*',FF(1)) / 2 - contract(FF(2),'!2*',FF(2)) / 2;     % total valley operator for orbital 1(+)
        V_z(2) = contract(FF(3),'!2*',FF(3)) / 2 - contract(FF(4),'!2*',FF(4)) / 2;     % total valley operator for orbital 2(-)
        Vtot_z = V_z(1) + V_z(2);   % total valley operator

        FF_quad = QSpace(1,2);
        Qnums = [-1,-1,0,0,0,0;
                 0,0,-1,-1,0,0];
        for it = 1:2
            FF_quad(it) = contract(setItag(FF(2*it-1),'s00','op1'), '!1', setItag(FF(2*it),'s00','op2'));
            A = getIdentity(FF_quad(it), 2, FF_quad(it), 4, 'op');
            FF_quad(it) = contract(FF_quad(it), '!1', A, '!3') / sqrt(2);
            FF_quad(it) = getsub(FF_quad(it), find(ismember(FF_quad(it).Q{3}, Qnums(it,:), 'rows')));
        end

        Lambda_p = 4*quadOp(FF_quad(1), FF_quad(2), '*');
        Lambda_p = getsub(Lambda_p, find(ismember(Lambda_p.Q{2}, [-1,-1,0,0,0,0], 'rows')));

        Lambda_z = 2*quadOp(FF_quad(1), FF_quad(1), '*') - 2*quadOp(FF_quad(2), FF_quad(2), '*');
        Lambda_z = getsub(Lambda_z, find(ismember(Lambda_z.Q{2}, Qnums, 'rows')));

        for it1 = 1:numel(FF)
            NF(it1) = quadOp(FF(it1), FF(it1), []);
        end
        Ntot = sum(NF);             % total charge operator

        % interacting term in the impurity Hamiltonian
        HU = (U/2)*(Ntot - N0*EF)*(Ntot - N0*EF);
        HU = HU + (J/2)*(contract(S_o1,'!2*',S_o1) - contract(SF(1),'!2*',SF(1)) - contract(SF(2),'!2*',SF(2)));
        HU = HU + (J/2)*(contract(S_o2,'!2*',S_o2) - contract(SF(3),'!2*',SF(3)) - contract(SF(4),'!2*',SF(4)));
        HU = HU - (J/4)*(NF(1)*NF(2) + NF(3)*NF(4));

        % quadratic term in the impurity Hamiltonian
        Hmu = - mu*Ntot;

        QF = QSpace(size(FF)); % QF = [FF,HU], to be used for the improved estimators of the self-energy
        for ito = (1:numel(FF))
            QF(ito) = contract(FF(ito),'!1',HU,[1 3 2])-contract(HU,'!1',FF(ito)); 
        end

        QFFF = QSpace(size(FF)); % commutator {QF, FF'}
        for ito = (1:numel(FF))
            QFFF(ito) = (contract(QF(ito),'!1',FF(ito),'*') + contract(FF(ito),'!2*',QF(ito)))/trace(getIdentity(FF(ito),3));
        end

        % isometry connecting the impurity and the first bath site
        A0 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);
        HU = contract(A0,'!2*',{HU,'!1',A0});           % HU in the K00 space
        Hmu = contract(A0,'!2*',{Hmu,'!1',A0});         % Hmu in the K00 space
        H0 = HU + Hmu + 1e-40*getIdentity(A0,2);        % impurity Hamiltonian

        STG = ['/data/',getenv('USER'),'/8flav/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];

        %% Define frequency grid

        ocont = getAcont(0, 0, 0, 0, 'emin', emin, 'emax', emax, 'estep', estep); 
        docon = zeros(numel(ocont), 1);
        docon(1) = 0.5 * (ocont(2) - ocont(1));
        docon(2:end - 1) = 0.5 * (ocont(3:end) - ocont(1:end - 2));
        docon(end) = 0.5 * (ocont(end) - ocont(end - 1));

        %% Define variables and paths for storing calculation results

        SE_H = zeros(numel(FF), nz);    % Hartree self-energy
        nu = zeros(numel(FF), 1);      % Filling of each orbital
        Adiscs = cell(36,nz);
        Aconts = cell(1,size(Adiscs,1));
        Acont1 = zeros(numel(ocont), numel(FF), numel(FF));
        Acont2 = Acont1;
        Acont3 = Acont1;
        save([STG,'/ocont.mat'],'ocont');

        nrgdata = cell(1,nz);

        node_name = getenv('HOSTNAME');
        disp2(node_name);

        TempAvail = cell(1,26);
        for itN = 1:14
            TempAvail{itN} = ['b',sprintf('%02d',itN)];
        end

        for itN = 1:12
            TempAvail{itN+14} = ['a',sprintf('%02d',itN+18)];
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

        %% NRG calculation
        for itz = 1:nz

            nrgdata{itz} = [NRGdataSTG,'/NRGdata_nz=',sprintf('%d',itz)];
            NRG_SL(nrgdata{itz},H0,A0,Lambda,ff{itz},FF,ZF,gg{itz},NF,'Nkeep',Nkeep, ...
                                        'Etrunc',Etrunc,'ETRUNC',ETRUNC,'dff',dff{itz},'dgg',dgg{itz},'deps',1e-10);
            getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);         % calculating the full density matrix(FDM)
        end

        for it1 = 1:numel(FF)
            nu(it1) = getEpVal(nrgdata{1}, NF(it1));
        end
        save([STG,'/nu.mat'],'nu');
                
        [Etot,Qtot,Qdiff] = plotE(nrgdata{1},'Emax',10,'legmax',25);       % Data for Eflow diagram    
        save([STG,'/Etot.mat'],'Etot');
        save([STG,'/Qtot.mat'],'Qtot');

        %% calculate impurity contribution to entropy

        BathEnt_path = ['/data/hyunsung/8flav/BathEnt_Nkeep=', sprintf('%d',Nkeep), '_Lambda=', sprintf('%.15g',Lambda), ...
                                '_Hyb=', sprintf('%.15g',Hyb), '_T=1e-12/BathEntData.mat'];

        BathEntData = load(BathEnt_path);   BathEntData = BathEntData.BathEntData;
        
        beta = 0.5:0.1:2;
        EntData.beta = beta;
        EntData.Temps = cell(1,numel(beta));
        EntData.S_imp = cell(1,numel(beta));


        for itN = 1:numel(beta)
            [Temps,~,~,Sent,~,~,~] = getTDconv(nrgdata{1},'useT','beta',beta(itN));
            Sent = Sent(Temps > T);
            Temps = Temps(Temps > T);

            idx = find(BathEntData.beta == beta(itN));
            Sent_bath = BathEntData.S_bath{idx};
            Temps_bath = BathEntData.Temps{idx};
            Sent_bath = Sent_bath(Temps_bath > T);

            Sent_imp = Sent - Sent_bath;
            
            EntData.Temps{itN} = Temps;
            EntData.S_imp{itN} = Sent_imp;
        end

        save([STG,'/EntData.mat'],'EntData');

        %% Compute self-energy

        if getSE    %% Compute self-energy

            for itz = 1:nz   
                % % Two-point correlators
                % Hartree self-energy
                QFz = QSpace(1, numel(FF)); % QFz = [FF,HU], to be used for the improved estimators of the self-energy
                for it1 = 1:numel(FF)
                    SE_H(it1,itz) = getEpVal(nrgdata{itz},QFFF(it1));
                    QFz(it1) = QF(it1) - SE_H(it1,itz)*FF(it1);
                end
                

                Ops1_FF = [repmat(FF(1), [1,4]), repmat(FF(2), [1,3]), repmat(FF(3), [1,2]), FF(4)];
                Ops2_FF = [FF(1:4), FF(2:4), FF(3:4), FF(4)];
                zflag_FF = ones(1, numel(Ops1_FF));
                cflag_FF = ones(1, numel(Ops1_FF));

                Ops1_QF = [repmat(FF(1), [1,4]), repmat(FF(2), [1,4]), repmat(FF(3), [1,4]), repmat(FF(4), [1,4])];
                Ops2_QF = repmat(QFz, [1,4]);
                zflag_QF = ones(1, numel(Ops1_QF));
                cflag_QF = ones(1, numel(Ops1_QF));

                Ops1_QQ = [repmat(QFz(1), [1,4]), repmat(QFz(2), [1,3]), repmat(QFz(3), [1,2]), QFz(4)];
                Ops2_QQ = [QFz(1:4), QFz(2:4), QFz(3:4), QFz(4)];
                zflag_QQ = ones(1, numel(Ops1_QQ));
                cflag_QQ = ones(1, numel(Ops1_QQ));

                Ops1 = [Ops1_FF, Ops1_QF, Ops1_QQ];
                Ops2 = [Ops2_FF, Ops2_QF, Ops2_QQ];
                zflag = [zflag_FF, zflag_QF, zflag_QQ];
                cflag = [cflag_FF, cflag_QF, cflag_QQ];

                [odisc,Adiscs(:,itz),sigmak] = getAdisc(nrgdata{itz}, Ops1, Ops2, ZF, 'zflag', zflag, 'cflag', cflag, ...
                                                            'emin', emin, 'emax', emax, 'estep', 2*estep);
            end % itz

            for ita = (1:size(Adiscs,1))
                [ocont,Aconts{ita}] = getAcont(odisc, mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3), sigmak, T/5, 'alphaz', 1/nz, ...
                                                    'emin', emin, 'emax', emax, 'estep', estep);
            end

            Gdep_FF = [1,2,3,4; -2,5,6,7; -3,-6,8,9; -4,-7,-9,10];
            Gdep_QF = 10 + [1,2,3,4; 5,6,7,8; 9,10,11,12; 13,14,15,16];
            Gdep_QQ = 26 + Gdep_FF;
            for it1 = 1:4
                for it2 = 1:4
                    Adisc2sum = mean(cellfun(@(x) sum(x(:)), Adiscs(Gdep_QF(it1,it2),:)));
                    Acont2(:,it1,it2) = Aconts{Gdep_QF(it1,it2)};
                    if it1 > it2
                        Acont1(:,it2,it1) = conj(Aconts{Gdep_FF(it2,it1)});
                        Acont3(:,it2,it1) = conj(Aconts{Gdep_QQ(it2,it1)});
                    else
                        Acont1(:,it1,it2) = Aconts{Gdep_FF(it1,it2)};
                        Acont3(:,it1,it2) = Aconts{Gdep_QQ(it1,it2)};
                    end
                end
            end

            for it1 = 1:numel(FF)
                Acont1(:, it1, it1) = Acont1(:, it1, it1) + mean(SE_H(it1,:));    % re-instate the Hartree self-energy
            end

            [SE,Aimp,~,ImSEle] = SEtrick(ocont, Acont1, Acont2, Adisc2sum, Acont3, ...
                                            'ozin', ozin, 'RhoV2in', RhoV2in, 'epsd', repmat(-U*N0,[1,4]));
            save([STG,'/SE.mat'],'SE');
            save([STG,'/Aimp.mat'],'Aimp');
            save([STG,'/ImSEle.mat'],'ImSEle');

        end % if getSE


        %% Compute dynamic susceptibilities

        if getSusc      %% Compute dynamic susceptibilities

            SuscOps_all = [Lambda_p; Lambda_z; SF(:); NF(:); Orb_z(:); V_z(:)];
            OpNames_all = {'Lambda_p', 'Lambda_z', 'ImpSp1p', 'ImpSp1m', 'ImpSp2p', 'ImpSp2m', 'Charge1p', 'Charge1m', 'Charge2p', 'Charge2m', 'Orbz_vp', 'Orbz_vm', 'Vz_o1', 'Vz_o2'};

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
                end

            end % itN
        end % getSusc

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