function ThsoK_NRG (parfn, varargin)

    setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

    partot = job_func_preamble(parfn, varargin{:});

    [PE, Nkeep, J0, K0, I0, T, JobName, nz] = loadvar(partot, ...
        {'PE', 'Nkeep', 'J0', 'K0', 'I0', 'T', 'JobName', 'nz'}, ...
            {[], [], [], [], [], [], [], []});

    num_threads_SL(8);

    disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
    disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
    disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
    disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
    disp2(['LOG_FILE : ', getenv('LOG_FILE')]);

    D = 1;
    Delta = pi;
    ozin = [-1;1]*D;
    RhoV2in = [1;1]*(Delta/pi);

    % NRG parameters
    Lambda = 4;
    N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
    Etrunc = [];
    ETRUNC = [];

    [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N,nz,'Nfit',round(-2*log(1e-8)/log(Lambda)));

    % Define Local Operators
    [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge, SU2spin, SU3channel');
    
    zf = getsub(ZF, find(ZF.Q{1}(:,3)==0));
    zf = getsub(zf, find(zf.Q{1}(:,4)==1));
    
    sf = getsub(SF,find(SF.Q{1}(:,3)==0)); % spin operators
    sf = getsub(sf,find(sf.Q{1}(:,4)==1));
    sf = getsub(sf,find(sf.Q{2}(:,3)==0));
    sf = getsub(sf,find(sf.Q{2}(:,4)==1));
         
    LF = quadOp(FF,FF,[0 0 1 1]); % orbital opertors;
    lf = getsub(LF,find(LF.Q{1}(:,3)==0));
    lf = getsub(lf,find(lf.Q{1}(:,4)==1));
    lf = getsub(lf,find(lf.Q{2}(:,3)==0));
    lf = getsub(lf,find(lf.Q{2}(:,4)==1));
    
    SLF = quadOp(FF,FF,[0 2 1 1]); % spin-orbital opertors;
    slf = getsub(SLF,find(SLF.Q{1}(:,3)==0));
    slf = getsub(slf,find(slf.Q{1}(:,4)==1));
    slf = getsub(slf,find(slf.Q{2}(:,3)==0));
    slf = getsub(slf,find(slf.Q{2}(:,4)==1));
    
    ef = getsub(IF.E,find(IF.E.Q{1}(:,3)==0));
    ef = getsub(ef,find(ef.Q{1}(:,4)==1));
    ef = getsub(ef,find(ef.Q{2}(:,3)==0));
    ef = getsub(ef,find(ef.Q{2}(:,4)==1));
    
    [F_imp, Z_imp, E_imp, S_sp, S_orb, S_sporb] = setItag('s00','op',ef,FF,zf,sf,lf,slf);
    [ZL,J_sp,J_orb,J_sporb,IL] = setItag('L00','op',ZF,SF,LF,SLF,IF.E);
    
    % local isometry and Hamiltonian
    A0 = getIdentity(E_imp,2,IL,2,'K00*',[1 3 2]); % isometry
    HSS = contract(J_sp, '!2*', S_sp, [2 1 3 4]);
    HLL = contract(J_orb, '!2*', S_orb, [2 1 3 4]);
    HSL = contract(J_sporb, '!2*', S_sporb, [2 1 3 4]);
   
    Htot = J0*HSS+K0*HLL+I0*HSL;    % Htot = HSS+HLL+HSL;
    
    H0 = contract(A0,'!2*',{Htot,'!13',A0}) + 1e-40*getIdentity(A0,2);
   
    % operators that define the two-point correlators
    Ops1 = [S_sp, S_orb, S_sporb; ...
                J_sp, J_orb, J_sporb];
    Ops2 = Ops1;
    OpNames_All = {'ImpSp','ImpOrb','ImpSpOrb';'BathSp','BathOrb','BathSpOrb'};
  
    STG = ['/data/',getenv('USER'),'/ThsoK/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep)];
    
    Adiscs = cell(numel(Ops1),nz); % discrete data
    Aconts = cell(1,size(Adiscs,1)); % continuous (i.e., broadened) spectral function

    nrgdata = cell(1,nz);

    for itz = (1:nz)
        nrgdata{itz} = NRG_SL([],H0,A0,Lambda,ff{itz}(2:end),FF,ZF,'Nkeep',Nkeep,'deps',1e-10);
        nrgdata{itz} = getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);   % calculating the full density matrix(FDM)
    end

    [Etot,Qtot,Qdiff] = plotE(nrgdata{1},'Emax',10,'legmax',25);       % Data for Eflow diagram    
    save([STG,'/Etot.mat'],'Etot');
    save([STG,'/Qtot.mat'],'Qtot');

    %{%}
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
            [odisc,Adiscs(:,itz),sigmak] = getAdisc(nrgdata{itz},Ops1(it,:),Ops2(it,:),ZF,'Z_L00',Z_imp,'zflag',zflag,'cflag',cflag,'emin',T);
            % odisc: [Numeric vector] center values of the frequency grid
            % sigmak [Numeric vector] center values of the broadening width bins
            % Adiscs: [cell array of numeric matrix] spectral function binned along odisc
            % Adiscs(:,itz): cell vector of numeric matrix. 
            % nth matrix Adiscs(n,itz) corresponds to spectral function of Ops1(n) and Ops2(n) binned along odisc and sigmak
            % Length of cell vector Adisc(:,itz): numel(Ops1) = numel(Ops2)
        end
    
        % file path to save Aconts
        SaveAconts = cellfun(@(x) [STG,'/NRG_Op=',OpNames{x},'.mat'], num2cell(1:numel(Ops1(it,:))), 'UniformOutput', false);
        SaveAdiscs = cellfun(@(x) [STG,'/Adiscs_NRG_Op=',OpNames{x},'_nz=',sprintf('%d',nz),'.mat'], num2cell(1:numel(Ops1(it,:))), 'UniformOutput', false);
    
        % Calculate dynamic susceptibilities
        for ita = (1:size(Adiscs,1))
  
            Adisc = mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3);
            %save(SaveAdiscs{ita},'Adisc');
    
            [ocont, Aconts{ita}] = getAcont(odisc,Adisc,sigmak,T/5,'alphaz',1/nz,'emin',T);
            % ocont : [numeric vector] Logarithimic frequency grid.
            % Aconts : [numeric vector] Smoothened spectral function.
    
            temp = Aconts{ita};
            save(SaveAconts{ita},'temp');
        end
  
        if it == 1
            save([STG,'/ocont.mat'],'ocont');
        end
    end
    %}

    %% calculate correlation functions between impurity and Wilson chain sites

    Sp_corr = cell(N,1);        % cell array of spin-spin correlators between the impurity and Wilson chain sites
    Orb_corr = cell(N,1);       % cell array of orbz-orbz correlators between the impurity and Wilson chain sites

    Sp_corr_Left = S_sp;
    Orb_corr_Left = S_orb;

    for site = (0:N-2)
        disptime(['Calculating correlator between impurity and Wilson chain site #',sprintf('%02d',site)]);
        
        Wilson_sp = setItag(['s',sprintf('%02d',site)],'op',J_sp);      % spin operator on Wilson chain site of interest
        Wilson_orb = setItag(['s',sprintf('%02d',site)],'op',J_orb);    % orbital psuedospin operator on Wilson chain site of interest
        
        if site > 0
            Sp_corr_Left = contract(Sp_corr_Left,'!1',nrgdata{1}.AK{site});
            Sp_corr_Left = contract(nrgdata{1}.AK{site},'!2*',Sp_corr_Left,[1,3,2]);
            
            Orb_corr_Left = contract(Orb_corr_Left,'!1',nrgdata{1}.AK{site});
            Orb_corr_Left = contract(nrgdata{1}.AK{site},'!2*',Orb_corr_Left,[1,3,2]);
            
        end
        
        Sp_corrT = contract(Sp_corr_Left,'!1',nrgdata{1}.AT{site+1},'!2');
        Sp_corrT = contract(Sp_corrT,conj(Wilson_sp));
        Sp_corrT = contract(nrgdata{1}.AT{site+1},'!2*',Sp_corrT);
        
        Sp_corrK = contract(Sp_corr_Left,'!1',nrgdata{1}.AK{site+1},'!2');
        Sp_corrK = contract(Sp_corrK,conj(Wilson_sp));
        Sp_corrK = contract(nrgdata{1}.AK{site+1},'!2*',Sp_corrK);
        
        Orb_corrT = contract(Orb_corr_Left,'!1',nrgdata{1}.AT{site+1},'!2');
        Orb_corrT = contract(Orb_corrT,conj(Wilson_orb));
        Orb_corrT = contract(nrgdata{1}.AT{site+1},'!2*',Orb_corrT);
        
        Orb_corrK = contract(Orb_corr_Left,'!1',nrgdata{1}.AK{site+1},'!2');
        Orb_corrK = contract(Orb_corrK,conj(Wilson_orb));
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