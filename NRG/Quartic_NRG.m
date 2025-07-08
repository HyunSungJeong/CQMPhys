function Quartic_NRG(parfn,varargin)

    try % in case of bug

        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [PE, Nkeep, Lambda, K_z, Q, T, JobName, nz] = loadvar(partot, ...
            {'PE', 'Nkeep', 'Lambda', 'K_z', 'Q', 'T', 'JobName', 'nz'}, ...
                {[], [], [], [], [], [], [], []});
    
        %num_threads_SL(8);
        
        disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
        disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
        disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
        disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
        disp2(['LOG_FILE : ', getenv('LOG_FILE')]);
    
        disp2(sprintf('%d',nz));
    
        % Impurity Hamiltonian parameters
        % K_z : spin-spin exchange coupling strength_z-component
        % Q : quartic-quartic (ddddcccc) interaction strength
    
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
    
        % local operators
        [FF,ZF,J_sp,IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1);
        [Fs,Zs,J_sp,Es] = setItag('s00','op',FF(:),ZF,J_sp(:),IF.E);

        Fcell = {Fs,Zs,J_sp,Es; Fs,Zs,J_sp,Es};
        symnum = [2,4; 1,3];
        sym = {'A', 'SU2'};
        tagnames = {'s1','s2'};

        for it1 = 1:size(Fcell,2)
            for it2 = 1:size(symnum,1)
                for it3 = 1:numel(sym)
                    Fcell{it2,it1} = addSymmetry(Fcell{it2,it1}, sym{it3}, 'pos', symnum(it2,it3));
                end

                Fcell{it2,it1} = setItag(tagnames{it2}, 'op', Fcell{it2, it1});
            end
        end

        A = getIdentity(Fcell{1,4},2,Fcell{2,4},2,'s00');
        Fs = QSpace(1,2); Zs = QSpace(1,2); J_sp = QSpace(1,2);
        
        Fs(1) = contract(A,'!3*',{Fcell{1,1},'!1',{A,Fcell{2,2}}}, [1,3,2]);    % annihilation operator for the first orbital. Note that the fermionic sign operator is contracted
        Fs(2) = contract(A,'!3*',{Fcell{2,1},'!1',A}, [1,3,2]);                 % annihilation operator for the second orbital. Fermionic sign operator is not contracted.

        for it = 1:size(Fcell,1)
            J_sp(it) = contract(A,'!3*',{Fcell{it,3},'!1',A}, [1,3,2]);     % spin operator for orbital # it(=1,2)
        end
        Zs = contract(A,'!3*',{Fcell{2,2},'!1',{Fcell{1,2},'!1',A}});       % fermionic sign operator
        J_sp_tot = J_sp(1) + J_sp(2);                                       % total spin(orb1+orb2) operator

        Es = contract(A,'!3*',{Fcell{2,4},{Fcell{1,4},'!1',A}});    % identity operator in the (2-orbital)local space
        
        % J_sp: bath spin op, Tr(S_a S_b) = 1/2 \delta_{ab}
    
        % bath orbital pseudospin op, Tr(T_a T_b) = 1/2 \delta_{ab}
        % Todo: J_orb_plus and J_orb_minus

        J_orb_z = (1/2)*( contract(Fs(1),'!2*',Fs(1),'!2') - contract(Fs(2),'!2*',Fs(2),'!2') );    % orbital-z operator
    
        %% project to the GS subspace to obtain impurity local operators
        % spin operators are zero in the GS(impurity) subspace
        % +- orbital pseudospin operators are zero in the GS(impurity) subspace

        Gsub = all(Zs.Q{1}(:,1) == -1, 2) & all(Zs.Q{1}(:,2) == 1, 2);
        Gsub = Gsub | all(Zs.Q{1}(:,1) == 1, 2) & all(Zs.Q{1}(:,2) == -1, 2);
        Z_imp = getsub(Zs, find(Gsub));         % fermionic sign operators at the impurity site

        Gsub = all(Es.Q{1}(:,1) == -1, 2) & all(Es.Q{1}(:,2) == 1, 2);
        Gsub = Gsub | all(Es.Q{1}(:,1) == 1, 2) & all(Es.Q{1}(:,2) == -1, 2);
        E_imp = getsub(Es, find(Gsub));                   % identity at the impurity site

        Gsub = all(J_orb_z.Q{1}(:,1) == -1, 2) & all(J_orb_z.Q{1}(:,2) == 1, 2);
        Gsub = Gsub | all(J_orb_z.Q{1}(:,1) == 1, 2) & all(J_orb_z.Q{1}(:,2) == -1, 2);
        T_z = getsub(J_orb_z, find(Gsub));      % orbital-z operator at the impurity site
    
        [Z_imp, E_imp, T_z] = setItag('L00','op',Z_imp, E_imp, T_z);

        %% Quartic operators
        CC = QSpace(1,2);
        for it = 1:2
            CC(it) = contract(setItag(Fs(it),'s00','op1'),'!1',setItag(Fs(it),'s00','op2'));
            A = getIdentity(CC(it),2,CC(it),4,['OP',sprintf('%d',it)]);
            CC(it) = -contract(CC(it),'!1',A,'!3')/sqrt(2);
        end

        Quart_c = quadOp(CC(2),CC(1),'*');
        Quart_c = setItag(Quart_c,'s00','op');
        Quart_imp = setItag(Quart_c,'L00','op');

        % local isometry and Hamiltonian
        A0 = getIdentity(E_imp,2,Es,2,'K00*',[1,3,2]);

        H0 = K_z*contract(A0,'!2*',{J_orb_z,'!2*',{T_z,A0}});               % orbz-orbz interaction
        H0 = H0 + Q*contract(A0,'!2*',{Quart_c,'!2*',{Quart_imp,'!1',A0}});     % quartic-quartic (ddddcccc) interaction
        H0 = H0 + Q*contract(A0,'!2*',{Quart_c,'!1',{Quart_imp,'!2*',A0}});
        H0 = H0 + 1e-40*contract(A0,'!2*',A0);
    
        STG = ['/data/',getenv('USER'),'/Quartic/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep)];
    
        nrgdata = cell(1,nz);

        for itz = 1:nz
            nrgdata{itz} = NRG_SL([],H0,A0,Lambda,ff{itz}(2:end),Fs,Zs,'Nkeep',Nkeep,'deps',1e-10);
            nrgdata{itz} = getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);   % calculating the full density matrix(FDM)
        end

        [Etot,Qtot,Qdiff] = plotE(nrgdata{1},'Emax',10,'legmax',25);       % Data for Eflow diagram    
        save([STG,'/Etot.mat'],'Etot');
        save([STG,'/Qtot.mat'],'Qtot');

        T_z_exp = getEpVal(nrgdata{1}, T_z);    % thermal expectation values
        J_orb_z_exp = getEpVal(nrgdata{1}, J_orb_z);
        disp2(['Thermal expectation value of T_z: ',sprintf('%.3g',T_z_exp)]);
        disp2(['Thermal expectation value of J_orb_z: ',sprintf('%.3g',J_orb_z_exp)]);
        save([STG,'/T_z_exp.mat'],'T_z_exp');
        save([STG,'/J_orb_z_exp.mat'],'J_orb_z_exp');

        % operators that define the two-point correlators
        Ops1 = [T_z-T_z_exp; J_sp(1); J_sp(2); J_sp_tot; J_orb_z-J_orb_z_exp];    % TODO: J_orb_plus/sqrt(2)];
        Ops2 = Ops1;
        OpNames = {'ImpOrb_z', 'BathSp1', 'BathSp2', 'BathSptot', 'BathOrb_z'};     % , 'BathOrb_plus'};

        % zflag and cflag are options in getAdisc(fdmNRG calc. of the spectral func. of the correlation function)
        % zflag: to use(1) or not to use(0) fermionic sign change for each operator.
        % 1 for fermionic operators, 0 for bosonic operators (e.g. if Ops1 = [fermionic, bosonic], zflag = [1,0]
        zflag = zeros(1,numel(Ops1));     
        % cflag: sign factors for each commutator(+ for fermionic, - for bosonic operators)
        % corresponds to the purple sign factor in lecture note 16.2 of tensor networks course
        cflag = (zflag-0.5)*2;

        Adiscs = cell(numel(Ops1),nz); % discrete data
        Aconts = cell(1,size(Adiscs,1)); % continuous (i.e., broadened) spectral function
        
        for itz = (1:nz)
            [odisc,Adiscs(:,itz),sigmak] = getAdisc(nrgdata{itz},Ops1,Ops2,ZF,'Z_L00',Z_imp,'zflag',zflag,'cflag',cflag,'emin',emin);
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
            %save(SaveAdiscs{ita},'Adisc');
        
            [ocont, Aconts{ita}] = getAcont(odisc,Adisc,sigmak,T/5,'alphaz',1/nz,'emin',emin);
            % ocont : [numeric vector] Logarithimic frequency grid.
            % Aconts : [numeric vector] Smoothened spectral function.
        
            temp = Aconts{ita};
            save(SaveAconts{ita},'temp');
        end
    
        save([STG,'/ocont.mat'],'ocont');
    
    
        
        %% calculate impurity contribution to entropy
        [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N+10,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));

        A02 = getIdentity(setItag('L00',getvac(E_imp)),2,Es,2,'K00*',[1 3 2]);
        H02 = contract(A02,'!2*',A02) + 1e-40*getIdentity(A02,2);
        bathNRG_data = NRG_SL([],H02,A02,Lambda,ff{1}(2:end),Fs,Zs,'Nkeep',Nkeep,'deps',1e-10);
        
        beta = 1.7;
        [Temps,~,~,Sent_bath,~,~,~] = getTDconv(bathNRG_data,'useT','beta',beta);
        [~,~,~,Sent,~,~,~] = getTDconv(nrgdata{1},'useT','beta',beta);
        Sent = Sent(Temps > T);
        Sent_bath = Sent_bath(Temps > T);
        Sent_imp = Sent - Sent_bath;
        Temps = Temps(Temps > T);
        save([STG,'/Temps.mat'],'Temps');
        save([STG,'/Sent_imp.mat'],'Sent_imp');
        
    catch Err
        disp2(getReport(Err));
        rethrow(Err);
    end
    
end