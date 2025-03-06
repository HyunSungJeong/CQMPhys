function DMFT_Hubbard_1B(parfn,varargin)
    setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'
  
    partot = job_func_preamble(parfn, varargin{:});
  
    [PE, Nkeep, J0, K0, I0, T, JobName] = loadvar(partot, ...
      {'PE', 'Nkeep', 'U', 'epsd', 'T', 'JobName'}, ...
        {[], [], [], [], [], []});
    
    % Hybridization function parametrized by the frequency grid 'ozin' and the
    % function value 'RhoV2in' evaluated at 'ozin'. 
    % Box-shaped function as the initial hybridization function
    D = 1; % half-bandwidth
    Gamma = 1e-2; % hybridization strength
    ozin = [-D;D];
    RhoV2in = (Gamma/pi)*[1;1]; % values outside of the 'ozin' grid are assumed to be zero
    
    % NRG parameters
    Lambda = 4;
    N = max(ceil(-2*log(T/500)/log(Lambda)),20);
    nz = 1;
    Nkeep = 3000;
    Etrunc = []; %9;
    ETRUNC = []; %inf(1,20);

    [ff,gg,dff,dgg] = doZLD(ozin,RhoV2in,Lambda,N,nz,'Nfit',round(-2*log(1e-8)/log(Lambda)));

    [FF,ZF,SF,IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1);
    [FF,ZF,SF,EF] = setItag('s00','op',FF(:),ZF,SF(:),IF.E);
    NF = quadOp(FF,FF,[]);

    FHU = QSpace(size(FF));
    for ito = (1:numel(FF))
        FHU(ito) = contract(FF(ito),'!1',HU,[1 3 2])-contract(HU,'!1',FF(ito)); % commutator [FF,HU], to be used for the self-energy trick
    end

    A0 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]); % isometry
    H0 = (U/2)*NF*(NF-EF) - epsd*NF + 1e-40*getIdentity(A0,2); % add "infinitesimal" term to keep all symmetry sectors

    % % Two-point correlators
    % operators that define the two-point correlators
    Ops1 = [FF(1);FHU(1);FHU(1);sum(NF);SF]; % set of the first operators
    Ops2 = Ops1; Ops2(2) = FF(1); % set of the second operators
    zflag = [ones(3,1);zeros(numel(Ops1)-3,1)];
    cflag = (zflag-0.5)*2;

    % 1st corr.: bare correlator < FF || FF' >
    % 2nd corr.: auxiliary correlator < FHU || FF' >
    % 3rd corr.: auxiliary correlator < FHU || FHU' >
    % 4th corr.: charge susceptibility
    % 5th corr.: spin susceptibility

    STG = ['/data/',getenv('USER'),'/DMFT/',JobName];

    nrgdata = ['/tmp/',getenv('USER'),filesep,JobName,'/NRG/NRG'];

    for itz = (1:nz)
        NRG_SL(nrgdata,H0,A0,Lambda,ff{itz},FF,ZF,gg{itz},NF, ...
                        'Nkeep',Nkeep,'Etrunc',Etrunc,'ETRUNC',ETRUNC,'dff',dff{itz},'dgg',dgg{itz});

        if itz == nz
            [Etot,Qtot,Qdiff] = plotE(nrgdata,'Emax',10,'legmax',25);       % Data for Eflow diagram
            save([STG,'/Etot.mat'],'Etot');
            save([STG,'/Qtot.mat'],'Qtot');
        end

        getRhoFDM(nrgdata,T,'-v','Rdiag',true);

        [odisc,Adiscs(:,itz),sigmak] = getAdisc(nrgdata,Ops1,Ops2,ZF,'zflag',zflag,'cflag',cflag);
    end

    for ita = (1:size(Adiscs,1))
        if ita <= 3 % fermionic correlators
            [ocont,Aconts{ita}] = getAcont(odisc,mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3),sigmak,T/5,'alphaz',1/nz);
        else % bosonic correlators
            [ocont,Aconts{ita}] = getAcont(odisc,sum(mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3),2),log(Lambda),T/5,'alphaz',1/nz,'Hfun','SLG');
        end
    end

    Adisc2sum = mean(cellfun(@(x) sum(x(:)), Adiscs(2,:)));

    % Obtain improved estimates of the impurity spectral function, by using the equations of motion (EoM)
    % Kugler's improved estimator ([Kugler 2022])
    [SE,Aimp] = SEtrick(ocont,Aconts{1},Aconts{2},Adisc2sum,Aconts{3}, ...
                'ozin',ozin,'RhoV2in',RhoV2in);

    Imp_charge_suscep = Aconts{4};
    Imp_spin_suscep = Aconts{5};
    save([STG,'/ocont.mat'],'ocont');
    save([STG,'/Imp_charge_suscep.mat'],'Imp_charge_suscep'); 
    save([STG,'/Imp_spin_suscep.mat'],'Imp_Spin_suscep');
    save([STG,'/Aimp.mat'],'Aimp');
end