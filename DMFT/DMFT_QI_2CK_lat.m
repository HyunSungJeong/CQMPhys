function DMFT_QI_2CK_lat(parfn,varargin)

    try % in case of bug
    
        setenv('RC_STORE',go('rcs'));           % '/project/hyunsung/RCStore'
        partot = job_func_preamble(parfn, varargin{:});

        [PE, ndfix, T, Lambda, JobName, nz, sigmarat, Nkeep, D, U, V, t_0, phi_div_pi, initSE, emin, emax, estep] = loadvar(partot, ...
                {'PE', 'ndfix', 'T', 'Lambda', 'JobName', 'nz', 'sigmarat', 'Nkeep', 'D', 'U', 'V', 't_0', 'phi_div_pi', 'initSE', 'emin', 'emax', 'estep'}, ...
                    {[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []});

        parpool('IdleTimeout',3500);

        %% DMFT parameters
        cvgth = 1e-3;   % convergence threshold, option for updateHyb (Default: 1e-3)
        amix = 1;       % linear mixing parameter for updateHyb (0<=amix<=1, Default: 1)
        Bro = 0;        % Broyden parameter for updateHyb (# of last iterations to be used for Broyden scheme, Default: 0)
        Ndmft = 100;    % maximum number of DMFT iterations before breaking the loop
        ieta = 1i * 1e-4;   % convergence generating factor (infinitesimal pure imaginary on the upper plane)
        abstol = 1e-8;      % absolute error tolerance for DOS integration. Option for MATLAB function "integral"
        SigmaTol = 1e-3;    % maximum tolerance for the absolute value of off-block-diagonal element of self-energy
        
        %% NRG parameters
        % Lambda
        N = max(ceil(-2*log(T/500)/log(Lambda))+10,20);  % Wilson chain length
        %Nfit = round(-2*log(1e-8)/log(Lambda));         % The Wilson chain site index after which the Wilson chain parameters are obtained via extrapolation
        Nfit = round(-2*log(1e-12)/log(Lambda));         % The Wilson chain site index after which the Wilson chain parameters are obtained via extrapolation
        % nz = 2 : Number of z's to be averaged
        % emin = T : Minimum absolute value of frequency grid
        % emax = 1e3 : Maximum absolute value of frequency grid
        % estep = 250 : # of steps to increase frequency per decade(x10), used in getAdisc and getAcont
                            % Acont: estep, Adisc: 2*estep
        Etrunc = [];
        ETRUNC = [];

        %% System parameters
        phi = pi*phi_div_pi;    % Distortion angle \phi
        t = t_0*cos(phi);       % intra-cell hopping amplitude t
        t_p = t_0*sin(phi);     % intra-cell hopping amplitude t'
        D = 2*V;                % half-bandwidth
        RhoB = @(epsilon) (2 / pi / D)*sqrt(1 - (epsilon/D).^2);     % semi-elliptic DOS for noninteracting Bethe lattice
        
        % (System parameters from parameter files)
        % t_0 : intra-cell hopping amplitude t_0
        % U : on-site Coulomb repulsion
        % D : half bandwidth
        % V : inter-site hopping amplitude

        %% frequency grid
        ocont = getAcont(0, 0, 0, 0, 'emin', emin, 'emax', emax, 'estep', estep);       % logarithmic frequency grid for smoothened spectral functions

        %% define output data: Hybridization functions, Green's functions, SEs, fillings, chemical potentials, Etots, Qtots, Qdiffs, and Wilson chain parameters      
        RhoV2out = zeros(numel(ocont), 2, Ndmft);       % Gamma (= hybridization function) / pi : before Broyden, for whole unit cell (omega, 2, Ndmft)             
        RhoV2out_1e = zeros(numel(ocont), Ndmft);       % Gamma (= hybridization function) / pi : before Broyden, for 1e orbital(1,1 component of RhoV2out), (omega, Ndmft)
        RhoV2out_1o = zeros(numel(ocont), Ndmft);       % Gamma (= hybridization function) / pi : before Broyden, for 1o orbital(4,4 component of RhoV2out), (omega, Ndmft)

        SE_save_e = zeros(numel(ocont),3,3,Ndmft);      % self-energy for the even sector   (omega,3,3,Ndmft)
        SE_save_o = zeros(numel(ocont),2,2,Ndmft);      % self-energy for the odd sector    (omega,2,2,Ndmft)

        RhoV2in = RhoV2out;             % hybridization function after Broyden, for whole unit cell

        Gfloc_e = zeros(3,3,numel(ocont),Ndmft);        % local lattice GF. for even sector
        Gfloc_o = zeros(2,2,numel(ocont),Ndmft);        % local lattice GF. for odd sector
        nures = zeros(5, Ndmft);        % fillings for each iteration [2e,3,2o,1e,1o]
        mures = zeros(1, Ndmft);        % chemical potential for each iteration
        mures(1) = 0;                   % half filling initial condition

        % Eflow data for each DMFT iteration
        Etots = cell(1, Ndmft);
        Qtots = cell(1, Ndmft);
        Qdiffs = cell(1, Ndmft);

        % Wilson chain parameters for each DMFT iteration
        ffs = cell(1, Ndmft);
        ggs = cell(1, Ndmft);
        dffs = cell(1, Ndmft);
        dggs = cell(1, Ndmft);

        %% define local operators

        % local operators in the even/odd basis
        [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge,SU2spin', 'NC', 3);
        [FF_L_eo, ZF_L_eo, EF_L] = setItag('L00', 'op', FF, ZF, IF.E);   % local operators in the L00(impurity) site. [2e, 3, 2o] orbitals
        
        [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge,SU2spin', 'NC', 2);
        [FF_s_eo, ZF_s_eo, EF_s] = setItag('s00', 'op', FF, ZF, IF.E);   % local operators in the s00(first bath) site. [1e, 1o] oribitals

        % number operators in the even/odd basis
        NF_L_eo = quadOp(FF_L_eo, FF_L_eo, []);    % number operators in the L00(impurity) site
        NF_s_eo = quadOp(FF_s_eo, FF_s_eo, []);    % number opators in the s00(first bath) site

        % annihilation ops. in the 'original' basis
        FF_L_orig = QSpace(1,3);        % annihilation ops. for orbital # 2,3,4 (defined at the L00 site)
        FF_s_orig = QSpace(1,2);        % annihilation ops. for oribtal # 1,5 (defined at the s00 site)
        
        FF_L_orig(1) = (FF_L_eo(1) + FF_L_eo(3)) / sqrt(2);
        FF_L_orig(2) = FF_L_eo(2);
        FF_L_orig(3) = (FF_L_eo(1) - FF_L_eo(3)) / sqrt(2);
        FF_s_orig(1) = (FF_s_eo(1) + FF_s_eo(2)) / sqrt(2);
        FF_s_orig(2) = (FF_s_eo(1) - FF_s_eo(2)) / sqrt(2);

        % number operators in the 'original' basis
        NF_L_orig = quadOp(FF_L_orig, FF_L_orig, []);       % number operators in the L00(impurity) site
        NF_s_orig = quadOp(FF_s_orig, FF_s_orig, []);       % number opators in the s00(first bath) site

        %% define system Hamiltonian

        % isometry connecting impurity and first bath site
        A0 = getIdentity(EF_L, 2, EF_s, 2, 'K00*', [1 3 2]);

        % Hamiltonian - interacting part(HU)
        HU_L = QSpace();    % part of HU that is defined on the L00 site
        for itL = 1:3
            HU_L = HU_L + (U/2)*(contract(NF_L_orig(itL),'!1',NF_L_orig(itL)) - NF_L_orig(itL));
        end
        HU_L = HU_L - (U/2)*sum(NF_L_orig);

        HU_s = QSpace();    % part of HU that is defined on the s00 site
        for its = 1:2
            HU_s = HU_s + (U/2)*(contract(NF_s_orig(its),'!1',NF_s_orig(its)) - NF_s_orig(its));
        end
        HU_s = HU_s - (U/2)*sum(NF_s_orig);

        HU = contract(A0,'!2*',{HU_s,A0}) + contract(A0,'!2*',{HU_L,A0});   % full interacting part Hamiltonian
        HU = HU + 1e-40*contract(A0,'!2*',A0);      % add infinitesimal*identity to avoid empty qspace when U=0

        % Hamiltonian - intra-cell hopping(H_hop)
        H_hop = sqrt(2)*t_0*contract(FF_L_eo(1),'!2*',FF_L_eo(2),'!2');
        H_hop = contract(A0,'!2*',{H_hop,A0});
        H_hop = H_hop + (t + t_p)*contract(A0,'!2*',{FF_s_eo(1),{FF_L_eo(1),'!2*',A0}});
        H_hop = H_hop + (t - t_p)*contract(A0,'!2*',{FF_s_eo(2),{FF_L_eo(3),'!2*',A0}});
        H_hop = H_hop + H_hop';

        %% Path to save calculation results
        STG = ['/data/',getenv('USER'),'/DMFT_QI_2CK_lat/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep)];

        %% DOS integration and DMFT self consistency to obtain initial hybridization function

        % initial self-energy
        SE_e = initSE{1};
        SE_o = initSE{2};

        % DOS integration to impose DMFT self-consistency condition
        Gfloc_init_e = zeros(3, 3, numel(ocont));     % initial local lattice GF., even sector
        Gfloc_init_o = zeros(2, 2, numel(ocont));     % initial local lattice GF., odd sector

        disptime('Performing DOS integration to obtain initial hybridization function from input self-energy');

        parfor ito = 1:numel(ocont)

            for it1 = 1:3
                for it2 = 1:3
                    Gfloc_init_e(it1, it2, ito) = integral(@(epsilon) RhoB(epsilon).*Gfk_e(epsilon, ocont(ito), reshape(SE_e(ito,:,:),[3,3]), U, V, t_0, t, t_p, mures(1), ieta, it1, it2), -D, +D, 'AbsTol', abstol);
                end
            end

            for it1 = 1:2
                for it2 = 1:2
                    Gfloc_init_o(it1, it2, ito) = integral(@(epsilon) RhoB(epsilon).*Gfk_o(epsilon, ocont(ito), reshape(SE_o(ito,:,:),[2,2]), U, V, t_0, t, t_p, mures(1), ieta, it1, it2), -D, +D, 'AbsTol', abstol);
                end
            end 
        end % parfor
        disptime('DOS integration complete');

        % impose DMFT self-consistency condition(Gfloc --> Hyb)
        initRVres_1e = zeros(numel(ocont),1);
        initRVres_1o = zeros(numel(ocont),1);
        for ito = 1:numel(ocont)
            Gfloc_e_inv = inv(Gfloc_init_e(:,:,ito));     % inverse of local lattice GF., even sector
            Gfloc_o_inv = inv(Gfloc_init_o(:,:,ito));     % inverse of local lattice GF., odd sector

            initRVres_1e(ito) = -imag( ocont(ito) + mures(1) + U/2 - SE_e(ito,1,1) - Gfloc_e_inv(1,1) ) / pi;
            initRVres_1o(ito) = -imag( ocont(ito) + mures(1) + U/2 - SE_o(ito,1,1) - Gfloc_o_inv(1,1) ) / pi;
        end

        RhoV2in(:,:,1) = cat(2, initRVres_1e, initRVres_1o);
        save([STG, '/RhoV2in.mat'], 'RhoV2in');

        %% DMFT loop
        for itD = 1:Ndmft

            num_threads_SL(PE);
            RVres = RhoV2in(:,:,itD);
            muc = mures(itD);           % chemical potential in current DMFT iteration
            epsd_e = -muc*ones(3,3);    % The matrix elements of the quadratic part of the Hamiltonian, even sector. Option for 'SEtrick'
            epsd_o = -muc*ones(2,2);    % The matrix elements of the quadratic part of the Hamiltonian, odd sector. Option for 'SEtrick'

            disp2('-----------------------------------------------------------------');
            disptime(['DMFT iteration No.', sprintf('%d',itD)]);

            %% NRG calculation for effective impurity to obtain self-energy (using impSE)

            % effective impurity Hamiltonian: quadratic part (interacting part: HU)
            Ntot = contract(A0,'!2*',{sum(NF_L_eo),A0}) + contract(A0,'!2*',{sum(NF_s_eo),A0});     % total number operator
            Hmu = H_hop - muc*Ntot;

            % path to save nrgdata
            nrgdata = ['/tmp/',getenv('USER'),'/nrgdata'];

            [ff, gg, dff, dgg] = doZLD(ocont, RVres, Lambda, N, nz, 'Nfit', Nfit);   % Wilson chain parameters(for saving)
            ffs{itD} = ff;
            ggs{itD} = gg;
            dffs{itD} = dff;
            dggs{itD} = dgg;
            save([STG,'/ffs.mat'], 'ffs');
            save([STG,'/ggs.mat'], 'ggs');
            save([STG,'/dffs.mat'], 'dffs');
            save([STG,'/dggs.mat'], 'dggs');

            Gdep = [1,2,0,3,0;
                    -2,4,0,5,0;
                    0,0,6,0,7;
                    -3,-5,0,8,0;
                    0,0,-7,0,9];
            
            Fdep = [1,2,0,3,0;
                    4,5,0,6,0;
                    0,0,7,0,8;
                    9,10,0,11,0;
                    0,0,12,0,13];

            [ocont,SE_tot,nrgdataz] = impSE(Hmu,HU,A0,FF_s_eo,ZF_s_eo,RVres,T,Lambda,nz,Nkeep,nrgdata,'F_L00',FF_L_eo,'Z_L00',ZF_L_eo,'Gdep',Gdep,'Fdep',Fdep, ...
                                            'doZLD',{'Nfit',Nfit},'NRG_SL',{'Etrunc', Etrunc, 'ETRUNC', ETRUNC},'getRhoFDM',{'Rdiag', true}, ...
                                                'getAdisc',{'emin', emin, 'emax', emax, 'estep', 2*estep}, ...
                                                    'getAcont',{'alphaz', sigmarat / nz, 'smin', log(Lambda) / 15, 'emin', emin, 'emax', emax, 'estep', estep}, ...
                                                        'SEtrick',{});

            disp2(reshape(max(abs(SE_tot),[],1), [5,5]));
            SE_e = SE_tot(:,[1,4,5],[1,4,5]);
            SE_o = SE_tot(:,2:3,2:3);
            SE_save_e(:,:,:,itD) = SE_e;
            SE_save_o(:,:,:,itD) = SE_o;
            save([STG,'/ocont.mat'], 'ocont');          % save the logarithmic frequency grid
            save([STG,'/SE_even.mat'], 'SE_save_e');    % save self-energy of the even sector(3x3)
            save([STG,'/SE_odd.mat'], 'SE_save_o');     % save self-energy of the odd sector(2x2)

            nures_L_eo = getEpVal(nrgdataz, NF_L_eo);        % filling in even/odd orbitals[2e, 3, 2o], L00(impurity) site
            nures_s_eo = getEpVal(nrgdataz, NF_s_eo);        % filling in even/odd orbitals[1e, 1o], s00(first bath) site 
            nures_eo = [nures_L_eo, nures_s_eo];            % filling in all even/odd orbitals [2e,3,2o,1e,1o]
            nures(:,itD) = nures_eo(:);                     % update fillings for current DMFT iteration
            nuressum = reshape(sum(nures,1), [1,Ndmft]);    % total filling for each DMFT iteration
            save([STG,'/nures.mat'], 'nures');              % save fillings

            [Etot,Qtot,Qdiff] = plotE(nrgdataz{1},'Emax',10,'legmax',25);       % data for Eflow diagram
            Etots{itD} = Etot;                                              % update Eflow data for current DMFT iteration
            Qtots{itD} = Qtot;
            Qdiffs{itD} = Qdiff;   
            save([STG,'/Etots.mat'], 'Etots');          % save Eflow data
            save([STG,'/Qtots.mat'], 'Qtots');
            save([STG,'/Qdiffs.mat'], 'Qdiffs');


            %% DOS integration to impose DMFT self-consistency condition
            Gftmp_e = zeros(3, 3, numel(ocont));     % temporary local lattice GF., even sector
            Gftmp_o = zeros(2, 2, numel(ocont));     % temporary local lattice GF., odd sector

            disptime('Starting DOS integration');
            parfor ito = 1:numel(ocont)

                for it1 = 1:3
                    for it2 = 1:3
                        Gftmp_e(it1, it2, ito) = integral(@(epsilon) RhoB(epsilon).*Gfk_e(epsilon, ocont(ito), reshape(SE_e(ito,:,:),[3,3]), U, V, t_0, t, t_p, muc, ieta, it1, it2), -D, +D, 'AbsTol', abstol);
                    end
                end

                for it1 = 1:2
                    for it2 = 1:2
                        Gftmp_o(it1, it2, ito) = integral(@(epsilon) RhoB(epsilon).*Gfk_o(epsilon, ocont(ito), reshape(SE_o(ito,:,:),[2,2]), U, V, t_0, t, t_p, muc, ieta, it1, it2), -D, +D, 'AbsTol', abstol);
                    end
                end 
            end % parfor
            disptime('DOS integration complete');
            Gfloc_e(:,:,:,itD) = Gftmp_e;
            Gfloc_o(:,:,:,itD) = Gftmp_o;

            %% impose DMFT self-consistency condition(Gfloc --> Hyb)
            for ito = 1:numel(ocont)
                Gfloc_e_inv = inv(Gftmp_e(:,:,ito));     % inverse of local lattice GF., even sector
                Gfloc_o_inv = inv(Gftmp_o(:,:,ito));     % inverse of local lattice GF., odd sector

                RhoV2out_1e(ito,itD) = -imag( ocont(ito) + muc + U/2 - SE_e(ito,1,1) - Gfloc_e_inv(1,1) ) / pi;
                RhoV2out_1o(ito,itD) = -imag( ocont(ito) + muc + U/2 - SE_o(ito,1,1) - Gfloc_o_inv(1,1) ) / pi;
            end

            RhoV2out = cat(2, reshape(RhoV2out_1e, [numel(ocont),1,Ndmft]), reshape(RhoV2out_1o, [numel(ocont),1,Ndmft]));
            save([STG,'/RhoV2out.mat'], 'RhoV2out');

            [iscvg, RVres, mutmp] = updateHyb(RhoV2in(:,:,(1:itD)), RhoV2out(:,:,(1:itD)), 3, mures(1:itD), nuressum(1:itD), ndfix, 'cvgth', [cvgth, cvgth], '-v', 'amix', amix, 'Broyden', Bro);
            % RVres : The hybridization function to be used as the input of the next DMFT iteration
            % numel(ocont)x2 matrix, first column representing the hybridization function for the 1e channel, 
            % and the second column representing the hybridization function for the 1o channel

            if ~iscvg && (itD < Ndmft)
                mures(itD + 1) = mutmp;          % Store chemical potential of itD+1 iteration
                RhoV2in(:, :, itD + 1) = RVres;  % Store hyb ftn. for itD+1 iteration
            elseif iscvg
                dispbox('DMFT converged.');
                break;
            end

            save([STG, '/RhoV2in.mat'], 'RhoV2in');
        end

    catch Err
        disp2(getReport(Err));
        rethrow(Err);
    end
end