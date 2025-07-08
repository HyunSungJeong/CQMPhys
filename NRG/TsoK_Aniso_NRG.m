function TsoK_Aniso_NRG(parfn,varargin)

    setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

    partot = job_func_preamble(parfn, varargin{:});
  
    [PE, Nkeep, Lambda, J0, K_perp, K_z, I0, T, JobName, nz] = loadvar(partot, ...
      {'PE', 'Nkeep', 'Lambda', 'J0', 'K_perp', 'K_z', 'I0', 'T', 'JobName', 'nz'}, ...
        {[], [], [], [], [], [], [], [], [], []});

    getCorr = false;
  
    %num_threads_SL(8);
    
    disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
    disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
    disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
    disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
    disp2(['LOG_FILE : ', getenv('LOG_FILE')]);
  
    disp2(sprintf('%d',nz));
  
    % Impurity Hamiltonian parameters
    % J0 : spin-spin exchange coupling
    % K_perp : orbital-orbital exchange coupling_perpendicular component
    % K_z : spin-spin exchange coupling_z-component
    % I0 : spin-orbital exchange coupling
  
    % T : Temperature
    emin = T;
    D = 1;
    Delta = pi;
    ozin = [-1;1]*D;
    RhoV2in = [1;1]*(Delta/pi);
  
    % NRG parameters
    N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
    %nz = 2;    
    Etrunc = []; %9;
    ETRUNC = []; %inf(1,20);
  
    [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N,nz,'Nfit',round(-2*log(1e-8)/log(Lambda)));
  
    % local operators
    [FF,ZF,J_sp,IF] = getLocalSpace('FermionS','Acharge(:),SU2spin','NC',2);
    [Fs,Zs,J_sp,Es] = setItag('s00','op',FF(:),ZF,J_sp(:),IF.E);
  
    % J_sp: bath spin op, Tr(S_a S_b) = 1/2 \delta_{ab}
  
    % bath orbital pseudospin op, Tr(T_a T_b) = 1/2 \delta_{ab}
    F1F2all = quadOp(Fs(1),Fs(2),'*');
    J_orb_plus = -sqrt(2)*getsub(F1F2all, find(all(F1F2all.Q{3}(:,3) == 0, 2)));
    F2F1all = quadOp(Fs(2),Fs(1),'*');
    J_orb_minus = -sqrt(2)*getsub(F2F1all, find(all(F2F1all.Q{3}(:,3) == 0, 2)));
    J_orb_z = (1/2)*( contract(Fs(1),'!2*',Fs(1),'!2') - contract(Fs(2),'!2*',Fs(2),'!2') );
  
    % bath spin-orbital op
    J_sporb_plus = (1/sqrt(2))*getsub(F1F2all, find(all(F1F2all.Q{3}(:,3) == 2, 2)));
    J_sporb_minus = (1/sqrt(2))*getsub(F2F1all, find(all(F2F1all.Q{3}(:,3) == 2, 2)));
    F1F1all = quadOp(Fs(1),Fs(1),'*');
    F2F2all = quadOp(Fs(2),Fs(2),'*');
    J_sporb_z = 0.5*getsub(F1F1all, find(all(F1F1all.Q{3}(:,3) == 2, 2)));
    J_sporb_z = J_sporb_z - 0.5*getsub(F2F2all, find(all(F2F2all.Q{3}(:,3) == 2, 2)));
    J_sporb_z = (1/sqrt(2))*J_sporb_z;
  
    Z_imp = getsub(Zs, find(all(Zs.Q{1}(:,1)+Zs.Q{1}(:,2) == -1, 2)));
    E_imp = getsub(Es, find(all(Es.Q{1}(:,1)+Es.Q{1}(:,2) == -1, 2)));
    S_sp = getsub(J_sp, find(all(J_sp.Q{1}(:,1)+J_sp.Q{1}(:,2) == -1, 2)));
  
    S_orb_plus = getsub(J_orb_plus, find(all(J_orb_plus.Q{1}(:,1)+J_orb_plus.Q{1}(:,2) == -1, 2)));
    S_orb_minus = getsub(J_orb_minus, find(all(J_orb_minus.Q{1}(:,1)+J_orb_minus.Q{1}(:,2) == -1, 2)));
    S_orb_z = getsub(J_orb_z, find(all(J_orb_z.Q{1}(:,1)+J_orb_z.Q{1}(:,2) == -1, 2)));
  
    S_sporb_plus = getsub(J_sporb_plus, find(all(J_sporb_plus.Q{1}(:,1)+J_sporb_plus.Q{1}(:,2) == -1, 2)));
    S_sporb_minus = getsub(J_sporb_minus, find(all(J_sporb_minus.Q{1}(:,1)+J_sporb_minus.Q{1}(:,2) == -1, 2)));
    S_sporb_z = getsub(J_sporb_z, find(all(J_sporb_z.Q{1}(:,1)+J_sporb_z.Q{1}(:,2) == -1, 2)));
  
    [Z_imp,E_imp,S_sp,S_orb_plus, S_orb_minus, S_orb_z, S_sporb_plus, S_sporb_minus, S_sporb_z] = ...
        setItag('L00','op',Z_imp(:), E_imp,S_sp(:), S_orb_plus, S_orb_minus, S_orb_z, S_sporb_plus, S_sporb_minus, S_sporb_z);
  
    % local isometry and Hamiltonian
    A0 = getIdentity(E_imp,2,Es,2,'K00*',[1,3,2]);
  
    H0 = J0*contract(A0,'!2*',{J_sp,'!2*',{S_sp,A0}});              % spin-spin
    H0 = H0 + (K_perp/2)*contract(A0,'!2*',{J_orb_plus,'!2*',{S_orb_plus,A0}});       % orbital-orbital
    H0 = H0 + (K_perp/2)*contract(A0,'!2*',{J_orb_minus,'!2*',{S_orb_minus,A0}});
    H0 = H0 + K_z*contract(A0,'!2*',{J_orb_z,'!2*',{S_orb_z,A0}});
    A = getIdentity(S_sp,3,S_orb_plus,3,'op*');
    H0 = H0 + (I0/2)*contract(A0,'!2*',{contract(A,'1,2',contract(S_sp,'!3',S_orb_plus,'!2'),'2,4'),'!2',{A0,J_sporb_plus,'!2*'}});  % spin-orbital
    %H0 = H0 + (I0/2)*contract(A0,'!2*',{contract(A,'1,2',contract(S_sp,'!3*',S_orb_plus,'!1*'),'2,4'),'!2',{A0,J_sporb_plus,'!1'}});
    H0 = H0 + (I0/2)*contract(A0,'!2*',{contract(conj(A),'1,2',contract(S_sp,'!3*',S_orb_plus,'!1*'),'2,4'),'!2',{A0,J_sporb_plus,'!1'}});
    H0 = H0 + I0*contract(A0, '!2*', contract(S_sp,'!3',S_orb_z,'!2',[1,3,2]), {A0,J_sporb_z,'!2*'});
    H0 = H0 + 1e-40*contract(A0,'!2*',A0);
  
    % operators that define the two-point correlators
    Ops1 = [S_sp, S_orb_plus/sqrt(2), S_orb_z; ...
                J_sp, J_orb_plus/sqrt(2), J_orb_z];% J_sp; J_orb; J_sp_orb];
    Ops2 = Ops1;
    OpNames_All = {'ImpSp','ImpOrb_plus','ImpOrb_z';'BathSp','BathOrb_plus','BathOrb_z'};
  
    STG = ['/data/',getenv('USER'),'/TsoK_Aniso/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];
    %STG = ['/data/',getenv('USER'),'/TsoK/TEST_Nk=3000_',JobName];
  
    nrgdata = cell(1,nz);

    for itz = (1:nz)
        nrgdata{itz} = NRG_SL([],H0,A0,Lambda,ff{itz}(2:end),FF,ZF,'Nkeep',Nkeep,'deps',1e-10);
        nrgdata{itz} = getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);   % calculating the full density matrix(FDM)
    end

    [Etot,Qtot,Qdiff] = plotE(nrgdata{1},'Emax',10,'legmax',25);       % Data for Eflow diagram    
    save([STG,'/Etot.mat'],'Etot');
    save([STG,'/Qtot.mat'],'Qtot');


    for it = (1:2)
        OpNames = OpNames_All(it,:);

        % zflag and cflag are options in getAdisc(fdmNRG calc. of the spectral func. of the correlation function)
        % zflag: to use(1) or not to use(0) fermionic sign change for each operator.
        % 1 for fermionic operators, 0 for bosonic operators (e.g. if Ops1 = [fermionic, bosonic], zflag = [1,0]
        zflag = zeros(1,numel(Ops1(it,:)));     
        % cflag: sign factors for each commutator(+ for fermionic, - for bosonic operators)
        % corresponds to the purple sign factor in lecture note 16.2 of tensor networks course
        cflag = (zflag-0.5)*2;

        Adiscs = cell(numel(Ops1(it,:)),nz); % discrete data
        Aconts = cell(1,size(Adiscs,1)); % continuous (i.e., broadened) spectral function
    
        for itz = (1:nz)
            [odisc,Adiscs(:,itz),sigmak] = getAdisc(nrgdata{itz},Ops1(it,:),Ops2(it,:),ZF,'Z_L00',Z_imp,'zflag',zflag,'cflag',cflag,'emin',emin);
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
        SaveAconts = cellfun(@(x) [STG,'/NRG_Op=',OpNames{x},'.mat'], num2cell(1:numel(Ops1(it,:))), 'UniformOutput', false);
    
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
  
        if it == 1
            save([STG,'/ocont.mat'],'ocont');
        end
    end


    %% calculate correlation functions between impurity and Wilson chain sites

    if getCorr

        Sp_corr = cell(N,1);        % cell array of spin-spin correlators between the impurity and Wilson chain sites
        Orb_corr = cell(N,1);       % cell array of orbz-orbz correlators between the impurity and Wilson chain sites
        
        Sp_corr_Left = S_sp;
        Orb_corr_Left = S_orb_z;
        
        for site = (0:N-2)
            disptime(['Calculating correlator between impurity and Wilson chain site #',sprintf('%02d',site)]);
        
            Wilson_sp = setItag(['s',sprintf('%02d',site)],'op',J_sp);      % spin operator on Wilson chain site of interest
            Wilson_orb_z = setItag(['s',sprintf('%02d',site)],'op',J_orb_z);    % orbital psuedospin operator on Wilson chain site of interest
            
            if site > 0
                Sp_corr_Left = contract(Sp_corr_Left,'!1',nrgdata{1}.AK{site});
                Sp_corr_Left = contract(nrgdata{1}.AK{site},'!2*',Sp_corr_Left,[1,3,2]);

                Orb_corr_Left = contract(Orb_corr_Left,'!1',nrgdata{1}.AK{site});
                Orb_corr_Left = contract(nrgdata{1}.AK{site},'!2*',Orb_corr_Left,'!2');
            end

            Sp_corrT = contract(Sp_corr_Left,'!1',nrgdata{1}.AT{site+1},'!2');
            Sp_corrT = contract(Sp_corrT,conj(Wilson_sp));
            Sp_corrT = contract(nrgdata{1}.AT{site+1},'!2*',Sp_corrT);
        
            Sp_corrK = contract(Sp_corr_Left,'!1',nrgdata{1}.AK{site+1},'!2');
            Sp_corrK = contract(Sp_corrK,conj(Wilson_sp));
            Sp_corrK = contract(nrgdata{1}.AK{site+1},'!2*',Sp_corrK);

            Orb_corrT = contract(Orb_corr_Left,'!1',nrgdata{1}.AT{site+1},'!2');
            Orb_corrT = contract(Orb_corrT,Wilson_orb_z);
            Orb_corrT = contract(nrgdata{1}.AT{site+1},'!2*',Orb_corrT);
        
            Orb_corrK = contract(Orb_corr_Left,'!1',nrgdata{1}.AK{site+1},'!2');
            Orb_corrK = contract(Orb_corrK,Wilson_orb_z);
            Orb_corrK = contract(nrgdata{1}.AK{site+1},'!2*',Orb_corrK);


            Sp_corr{site+1} = 0;
            corrT = contract(Sp_corrT, diag(nrgdata{1}.RhoT{site+1}));
            corrK = contract(Sp_corrK, nrgdata{1}.RhoK{site+1});

            if ~isempty(corrT)
                Sp_corr{site+1} = Sp_corr{site+1} + corrT.data{1};
            end

            if ~isempty(corrK)
                Sp_corr{site+1} = Sp_corr{site+1} + corrK.data{1};
            end
            
            Orb_corr{site+1} = 0;
            corrT = contract(Orb_corrT, diag(nrgdata{1}.RhoT{site+1}));
            corrK = contract(Orb_corrK, nrgdata{1}.RhoK{site+1});

            if ~isempty(corrT)
                Orb_corr{site+1} = Orb_corr{site+1} + corrT.data{1};
            end

            if ~isempty(corrK)
                Orb_corr{site+1} = Orb_corr{site+1} + corrK.data{1};
            end

        end

        save([STG,'/spin_spin_correlators.mat'],'Sp_corr');
        save([STG,'/orbital_orbital_correlators.mat'],'Orb_corr');
        %%%
    end

    [ff, ~] = doZLD(ozin,RhoV2in,Lambda,N+10,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));

    A02 = getIdentity(setItag('L00',getvac(E_imp)),2,Es,2,'K00*',[1 3 2]);
    H02 = contract(A02,'!2*',A02) + 1e-40*getIdentity(A02,2);
    bathNRG_data = NRG_SL([],H02,A02,Lambda,ff{1}(2:end),FF,ZF,'Nkeep',Nkeep,'deps',1e-10);

    for beta = 1:0.1:2
        [Temps,~,~,Sent_bath,~,~,~] = getTDconv(bathNRG_data,'useT','beta',beta);
        [~,~,~,Sent,~,~,~] = getTDconv(nrgdata{1},'useT','beta',beta);

        Sent = Sent(Temps > T);
        Sent_bath = Sent_bath(Temps > T);
        Sent_imp = Sent - Sent_bath;
        Temps = Temps(Temps > T);
        save([STG,'/Temps.mat'],'Temps');
        save([STG,'/Sent_imp_beta=',sprintf('%.15g',beta),'.mat'],'Sent_imp');
    end
    
end