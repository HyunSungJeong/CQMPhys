function bathNRG(parfn,varargin)
    setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

    partot = job_func_preamble(parfn, varargin{:});

    [PE,T] = loadvar(partot, ...
                        {'PE', 'T'}, ...
                            {[], []});

    D = 1;
    Delta = pi;
    ozin = [-1;1]*D;
    RhoV2in = [1;1]*(Delta/pi);

    Lambda = 4;
    N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
    nz = 1;    
    Nkeep = 8000;
    
    if ~exist(['/data/',getenv('USER'),'/bathNRG/Lambda=',sprintf('%.15g',Lambda),'_Nkeep=',sprintf('%.15g',Nkeep)], 'dir')
        mkdir(['/data/',getenv('USER'),'/bathNRG/Lambda=',sprintf('%.15g',Lambda),'_Nkeep=',sprintf('%.15g',Nkeep)]);
    end

    [ff, gg] = doZLD(ozin,RhoV2in,Lambda,N,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));

    % local operators
    [FF,ZF,J_sp,IF] = getLocalSpace('FermionS','Acharge,SU2spin,SU2channel','NC',2);
    [Fs,Zs,J_sp,Es] = setItag('s00','op',FF(:),ZF,J_sp(:),IF.E);
    Ns = quadOp(Fs,Fs,[]);

    F_imp = getsub(FF,1);
    Z_imp = getsub(ZF,2);
    E_imp = getsub(IF.E,2);
    [F_imp,Z_imp,E_imp] = setItag('L00','op',F_imp(:),Z_imp(:),E_imp);

    A02 = getIdentity(setItag('L00',getvac(E_imp)),2,Es,2,'K00*',[1 3 2]);
    H02 = contract(A02,'!2*',A02) + 1e-40*getIdentity(A02,2);
    bathNRG_data = NRG_SL([],H02,A02,Lambda,ff{1}(2:end),FF,ZF,gg{1}(2:end),Ns,'Nkeep',Nkeep);
    

    %bathNRG_data = ['/data/',getenv('USER'),'/bathNRG/Lambda=',sprintf('%.15g',Lambda),'_Nkeep=',sprintf('%.15g',Nkeep), ...
    %                    '/T=',sprintf('%.15g',T),'/bathNRGdata'];

    beta = (0.5: 0.1: 2);

    for it = (1:numel(beta))
        [Temps,~,~,Sent,~,~,~] = getTDconv(bathNRG_data,'useT','beta',beta(it));
        disp2(['beta = ',sprintf('%.15g',beta(it)),'\n']);

        save(['/data/',getenv('USER'),'/bathNRG/Lambda=',sprintf('%.15g',Lambda),'_Nkeep=',sprintf('%.15g',Nkeep), ...
                '/Temps_MinT=',sprintf('%.15g',T),'_beta=',sprintf('%.2f',beta(it)),'.mat'],'Temps');
        
        save(['/data/',getenv('USER'),'/bathNRG/Lambda=',sprintf('%.15g',Lambda),'_Nkeep=',sprintf('%.15g',Nkeep), ...
                '/Sent_MinT=',sprintf('%.15g',T),'_beta=',sprintf('%.2f',beta(it)),'.mat'],'Sent');
    end

end