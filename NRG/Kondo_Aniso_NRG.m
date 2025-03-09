function Kondo_Aniso_NRG(parfn,varargin)

    setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

    partot = job_func_preamble(parfn, varargin{:});
  
    [PE, Nkeep, J_perp, J_z, T, h, JobName, nz] = loadvar(partot, ...
      {'PE', 'Nkeep', 'J_perp', 'J_z', 'T', 'h', 'JobName', 'nz'}, ...
        {[], [], [], [], [], [], [], []});
  
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
    D = 1;
    Delta = pi;
    ozin = [-1;1]*D;
    RhoV2in = [1;1]*(Delta/pi);
  
    % NRG parameters
    Lambda = 4;
    N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
    %nz = 2;    
    Etrunc = []; %9;
    ETRUNC = []; %inf(1,20);
  
    [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N,nz,'Nfit',round(-2*log(1e-8)/log(Lambda)));
  
    % local operators
    % J: bath spin op, Tr(S_a S_b) = 1/2 \delta_{ab}
    [FF,ZF,J_sp,IF] = getLocalSpace('FermionS','Acharge,Aspin','NC',1);
    [Fs,Zs,J_sp,Es] = setItag('s00','op',FF(:),ZF,J_sp(:),IF.E);
    J_sp_plus = -sqrt(2)*J_sp(1);
    J_sp_minus = sqrt(2)*J_sp(2);
    J_sp_z = J_sp(3);

    % Impurity operators
    Z_imp = getsub(Zs, find(all(Zs.Q{1}(:,1) == 0, 2)));
    E_imp = getsub(Es, find(all(Es.Q{1}(:,1) == 0, 2)));

    S_plus = getsub(J_sp_plus, find(all(J_sp_plus.Q{1}(:,1) == 0, 2)));
    S_minus = getsub(J_sp_minus, find(all(J_sp_minus.Q{1}(:,1) == 0, 2)));
    S_z = getsub(J_sp_z, find(all(J_sp_z.Q{1}(:,1) == 0, 2)));

    [Z_imp,E_imp,S_plus,S_minus,S_z] = setItag('L00','op',Z_imp(:), E_imp,S_plus,S_minus,S_z);
  
    % local isometry and Hamiltonian
    A0 = getIdentity(E_imp,2,Es,2,'K00*',[1,3,2]);
  
    H0 = (J_perp/2)*contract(A0,'!2*',{J_sp_plus,'!2*',{S_plus,A0}});       % spin-spin exchange interaction
    H0 = H0 + (J_perp/2)*contract(A0,'!2*',{J_sp_minus,'!2*',{S_minus,A0}});
    H0 = H0 + J_z*contract(A0,'!2*',{J_sp_z,'!2*',{S_z,A0}});
    H0 = H0 + h*contract(A0,'!2*',{S_z,A0});

    H0 = H0 + 1e-40*contract(A0,'!2*',A0);
  
    % operators that define the two-point correlators
    Ops1 = [S_plus/sqrt(2), S_z; J_sp_plus/sqrt(2), J_sp_z];
    Ops2 = Ops1;
    OpNames_All = {'ImpSp_plus','ImpSp_z'; 'BathSp_plus','BathSp_z'};
  
    STG = ['/data/',getenv('USER'),'/Kondo_Aniso/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep)];
  
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
  
  
    beta = 1.5;
    [ff, gg] = doZLD(ozin,RhoV2in,Lambda,90,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));
    A02 = getIdentity(setItag('L00',getvac(E_imp)),2,Es,2,'K00*',[1 3 2]);
    H02 = contract(A02,'!2*',A02) + 1e-40*getIdentity(A02,2);
    bathNRG_data = NRG_SL([],H02,A02,Lambda,ff{1}(2:end),FF,ZF,'Nkeep',Nkeep);

    [Temps,~,~,Sent_bath,~,~,~] = getTDconv(bathNRG_data,'useT','beta',beta);
    
    [Temps,~,~,Sent,~,~,~] = getTDconv(nrgdata{1},'useT','beta',beta);
    Sent = Sent(Temps > T);
    Sent_bath = Sent_bath(Temps > T);
    Sent_imp = Sent - Sent_bath;
    Temps = Temps(Temps > T);
    save([STG,'/Temps.mat'],'Temps');
    save([STG,'/Sent_imp.mat'],'Sent_imp');
    
end