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

        %% define output data: Hybridization functions, Green's functions, SEs, (improved) impurity spectral functions, fillings, chemical potentials, Etots, Qtots, Qdiffs, and Wilson chain parameters      
        RhoV2out = zeros(numel(ocont), 2, Ndmft);       % Gamma (= hybridization function) / pi : before Broyden, for whole unit cell (omega, 2, Ndmft)             
        RhoV2out_1e = zeros(numel(ocont), Ndmft);       % Gamma (= hybridization function) / pi : before Broyden, for 1e orbital(1,1 component of RhoV2out), (omega, Ndmft)
        RhoV2out_1o = zeros(numel(ocont), Ndmft);       % Gamma (= hybridization function) / pi : before Broyden, for 1o orbital(4,4 component of RhoV2out), (omega, Ndmft)

        SE_save_e = zeros(numel(ocont),3,3,Ndmft);      % self-energy for the even sector   (omega,3,3,Ndmft)
        SE_save_o = zeros(numel(ocont),2,2,Ndmft);      % self-energy for the odd sector    (omega,2,2,Ndmft)

        Aimp_save_e = SE_save_e;
        Aimp_save_o = SE_save_o;

        RhoV2in = RhoV2out;             % hybridization function after Broyden, for whole unit cell
        RhoV2in_1e = RhoV2out_1e;       % hybridization function after Broyden, for even sector
        RhoV2in_1o = RhoV2out_1o;       % hybridization function after Broyden, for odd sector

        Gfloc_e = zeros(3,3,numel(ocont),Ndmft);        % local lattice GF. for even sector
        Gfloc_o = zeros(2,2,numel(ocont),Ndmft);        % local lattice GF. for odd sector
        nures = zeros(5, Ndmft);        % fillings for each iteration
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
        [FF_eo, ZF_eo, SF_eo, IF_eo] = getLocalSpace('FermionS', 'Acharge,SU2spin', 'NC', 5);
        [FF_eo, ZF_eo, EF_eo] = setItag('s00', 'op', FF_eo, ZF_eo, IF_eo.E);
        % local operators in the even/odd basis
        % FF_eo = [FF_1e, FF_2e, FF_3, FF_1o, FF_2o]
        NF_eo = quadOp(FF_eo, FF_eo, []);    % number operators in the even/odd basis

        % annihilation ops. in the 'original' basis
        FF(1) = (FF_eo(1) + FF_eo(4)) / sqrt(2);
        FF(2) = (FF_eo(2) + FF_eo(5)) / sqrt(2);
        FF(3) = FF_eo(3);
        FF(4) = (FF_eo(2) - FF_eo(5)) / sqrt(2);
        FF(5) = (FF_eo(1) - FF_eo(4)) / sqrt(2);
        NF = quadOp(FF, FF, []);      % number operators in the 'original' basis

        %% define system Hamiltonian
        % Hamiltonian - interacting part(HU)
        HU = QSpace();
        for itm = (1:5)
            HU = HU + (U/2)*(contract(NF(itm),'!1',NF(itm)) - NF(itm));
        end
        HU = HU - (U/2)*sum(NF) + 1e-40*EF_eo;

        % Hamiltonian - intra-cell hopping(H_hop)
        H_hop = sqrt(2)*t_0*contract(FF_eo(2),'!2*',FF_eo(3),'!2');
        H_hop = H_hop + (t + t_p)*contract(FF_eo(1),'!2*',FF_eo(2),'!2');
        H_hop = H_hop + (t - t_p)*contract(FF_eo(4),'!2*',FF_eo(5),'!2');
        H_hop = H_hop + H_hop';
        
        % isometry connecting impurity and first bath site
        A0 = getIdentity(setItag('L00', getvac(EF_eo)), 2, EF_eo, 2, 'K00*', [1 3 2]); % isometry

        % commutator [FK,HU], to be used for the self-energy trick
        for ito = (1:numel(FF_eo))
            FHU_eo(ito) = contract(FF_eo(ito), '!1', HU, [1 3 2]) - contract(HU, '!1', FF_eo(ito)); % commutator [FF_eo,HU], to be used for the self-energy trick
        end

        %% Path to save calculation results
        STG = ['/data/',getenv('USER'),'/DMFT_QI_2CK_lat/',JobName,'_Nkeep=',sprintf('%.15g',Nkeep)];
        save([STG,'/ocont.mat'], 'ocont');

        %% DOS integration and DMFT self consistency to obtain initial hybridization function

        % initial self-energy
        SE_e = initSE{1};
        SE_o = initSE{2};

        % DOS integration to impose DMFT self-consistency condition
        Gfloc_init_e = zeros(3, 3, numel(ocont));     % temporary local lattice GF., even sector
        Gfloc_init_o = zeros(2, 2, numel(ocont));     % temporary local lattice GF., odd sector

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

            num_threads_SL(8);
            RVres = RhoV2in(:,:,itD);
            muc = mures(itD);       % chemical potential in current DMFT iteration
            epsd_e = -muc*ones(3,3);  % The matrix elements of the quadratic part of the Hamiltonian, even sector. Option for 'SEtrick'
            epsd_o = -muc*ones(2,2);  % The matrix elements of the quadratic part of the Hamiltonian, odd sector. Option for 'SEtrick'

            disp2('-----------------------------------------');
            disptime(['DMFT iteration No.', sprintf('%d',itD)]);

            %% NRG calculation for effective impurity to obtain self-energy

            % effective impurity Hamiltonian
            H_imp = HU + H_hop - muc*sum(NF_eo);
            H0 = contract(A0, '!2*', {H_imp, '!1', A0}) + 1e-40*getIdentity(A0, 2);    % Add infinitesimal*identity to keep all symmetry sectors

            RVres_exp = [RVres(:,1), zeros(size(RVres,1),2), RVres(:,2), zeros(size(RVres,1), 1)];
            % RVres_exp: expanded RVres. Zeros filled in for the orbitals that do not hybridize with bath
            % numel(ocont)x5 matrix, first column representing the hybridization function for the 1e channel, 
            % and the fourth column representing the hybridization function for the 1o channel
            % other columns represent hybridization function for the 2e, 3, 2o, which are zero(nonexistent)

            [ff, gg, dff, dgg] = doZLD(ocont, RVres_exp, Lambda, N, nz, 'Nfit', Nfit);   % Wilson chain parameters
            ffs{itD} = ff;
            ggs{itD} = gg;
            dffs{itD} = dff;
            dggs{itD} = dgg;
            save([STG,'/ffs.mat'], 'ffs');
            save([STG,'/ggs.mat'], 'ggs');
            save([STG,'/dffs.mat'], 'dffs');
            save([STG,'/dggs.mat'], 'dggs');

            nrgdata = cell(1,nz);

            for itz = (1:nz)
                nrgdata{itz} = NRG_SL([], H0, A0, Lambda, ff{itz}, FF_eo, ZF_eo, ...
                                    gg{itz}, NF_eo, 'Nkeep', Nkeep, 'Etrunc', Etrunc, 'ETRUNC', ETRUNC, ...
                                        'dff', dff{itz}, 'dgg', dgg{itz});
                nrgdata{itz} = getRhoFDM(nrgdata{itz}, T, '-v', 'Rdiag', true);   % calculating the full density matrix(FDM)
            end

            nures_tmp = getEpVal(nrgdata{1}, NF_eo);        % filling in even/odd orbitals
            nures(:,itD) = nures_tmp(:);                    % update fillings for current DMFT iteration
            nuressum = reshape(sum(nures,1), [1,Ndmft]);    % total filling for each DMFT iteration
            
            save([STG,'/nures.mat'], 'nures');

            [Etot,Qtot,Qdiff] = plotE(nrgdata{1},'Emax',10,'legmax',25);       % Data for Eflow diagram
            Etots{itD} = Etot;          % update Eflow data for current DMFT iteration
            Qtots{itD} = Qtot;
            Qdiffs{itD} = Qdiff;   
            save([STG,'/Etots.mat'], 'Etots');
            save([STG,'/Qtots.mat'], 'Qtots');
            save([STG,'/Qdiffs.mat'], 'Qdiffs');

            %% calculate susceptibilities for SEtrick

            % define input variables for SEtrick to be calculated using getAdisc and getAcont
            Acont1_e = zeros(numel(ocont), 3, 3);   % Acont1 for SEtrick, even sector
            Acont2_e = Acont1_e;                    % Acont2 for SEtrick, even sector
            Acont3_e = Acont1_e;                    % Acont3 for SEtrick, even sector
            Adsic2sum_e = zeros(3,3);                 % Adisc2sum for SEtrick, even sector

            Acont1_o = zeros(numel(ocont), 2, 2);   % Acont1 for SEtrick, odd sector
            Acont2_o = Acont1_o;                    % Acont2 for SEtrick, odd sector
            Acont3_o = Acont1_o;                    % Acont3 for SEtrick, odd sector
            Adsic2sum_o = zeros(2,2);                 % Adisc2sum for SEtrick, odd sector

            % define operators for susceptibility calculations
            FF_e = FF_eo(1:3);
            FF_o = FF_eo(4:5);

            FHU_e = FHU_eo(1:3);
            FHU_o = FHU_eo(4:5);

            % calculate Acont1 for SEtrick
            % only calculate on- and upper-diagonal elements of the spectral function and use Hermicity to fill in the lower-diagonal elements
            disptime('Calculating Acont1 for SEtrick');
            Adiscs = cell(9, nz);
            Aconts = cell(1, size(Adiscs,1));

            % operators for getAdisc
            Ops1_e = [repmat(FF_e(1), [1,3]), repmat(FF_e(2), [1,2]), FF_e(3)].';
            Ops2_e = [FF_e(1:3), FF_e(2:3), FF_e(3)].';
            Ops1_o = [repmat(FF_o(1), [1,2]), FF_o(2)].';
            Ops2_o = [FF_o(1:2), FF_o(2)].';

            for itz = 1:nz
                [odisc, Adiscs(:, itz), sigmak] = getAdisc(nrgdata{itz}, [Ops1_e; Ops1_o], [Ops2_e; Ops2_o], ZF_eo, 'emin', emin, 'emax', emax, 'estep', 2*estep);
            end

            for ita = 1:size(Adiscs,1)
                Adisc = mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3);
                [~, Aconts{ita}] = getAcont(odisc, Adisc, sigmak, T / 5, 'alphaz', sigmarat / nz, ...
                                                'smin', log(Lambda) / 15, 'emin', emin, 'emax', emax, 'estep', estep);
            end
            Aconts_e = Aconts(1:6);
            Aconts_o = Aconts(7:9);

            idx = [1,2,3;
                   0,4,5;
                   0,0,6];
            
            for it1 = 1:3
                for it2 = it1:3
                    Acont1_e(:, it1, it2) = Aconts_e{idx(it1, it2)};
                    Acont1_e(:, it2, it1) = conj(Aconts_e{idx(it1, it2)});
                end
            end

            idx = [1,2;
                   0,3];

            for it1 = 1:2
                for it2 = it1:2
                    Acont1_o(:, it1, it2) = Aconts_o{idx(it1, it2)};
                    Acont1_o(:, it2, it1) = conj(Aconts_o{idx(it1, it2)});
                end
            end

            % calculate Acont2 and Adisc2sum for SEtrick
            disptime('Calculate Acont2 and Adisc2sum for SEtrick');
            Adiscs = cell(13, nz);
            Aconts = cell(1, size(Adiscs,1));

            % operators for getAdisc
            Ops1_e = [repmat(FHU_e(1), [3,1]); repmat(FHU_e(2), [3,1]); repmat(FHU_e(3), [3,1])];
            Ops2_e = [FF_e(:); FF_e(:); FF_e(:)];
            Ops1_o = [repmat(FHU_o(1), [2,1]); repmat(FHU_o(2), [2,1])];
            Ops2_o = [FF_o(:); FF_o(:)];

            for itz = 1:nz
                [odisc, Adiscs(:, itz), sigmak] = getAdisc(nrgdata{itz}, [Ops1_e; Ops1_o], [Ops2_e; Ops2_o], ZF_eo, 'emin', emin, 'emax', emax, 'estep', 2*estep);
            end
            Adiscs_e = Adiscs(1:9, :);
            Adiscs_o = Adiscs(10:13, :);

            for ita = 1:size(Adiscs,1)
                Adisc = mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3);
                [~, Aconts{ita}] = getAcont(odisc, Adisc, sigmak, T / 5, 'alphaz', sigmarat / nz, ...
                                                'smin', log(Lambda) / 15, 'emin', emin, 'emax', emax, 'estep', estep);
            end
            Aconts_e = Aconts(1:9);
            Aconts_o = Aconts(10:13);

            for it1 = 1:3
                for it2 = 1:3
                    Acont2_e(:, it1, it2) = Aconts_e{3*(it1-1)+it2};
                end
            end

            for it1 = 1:2
                for it2 = 1:2
                    Acont2_o(:, it1, it2) = Aconts_o{2*(it1-1)+it2};
                end
            end

            Adisc2_e = cell(3, 3, nz);
            for it1 = 1:3
                for it2 = 1:3
                    Adisc2_e(it1, it2, :) = Adiscs_e(3*(it1-1)+it2, :);
                end
            end
            Adisc2sum_e = mean(cellfun(@(x) sum(x(:)), Adisc2_e), 3);

            Adisc2_o = cell(2, 2, nz);
            for it1 = 1:2
                for it2 = 1:2
                    Adisc2_o(it1, it2, :) = Adiscs_o(2*(it1-1)+it2, :);
                end
            end
            Adisc2sum_o = mean(cellfun(@(x) sum(x(:)), Adisc2_o), 3);

            % calculate Acont3 for SEtrick
            disptime('Calculate Acont3 for SEtrick');
            Adiscs = cell(13, nz);
            Aconts = cell(1, size(Adiscs,1));

            % operators for getAdisc
            Ops1_e = [repmat(FHU_e(1), [3,1]); repmat(FHU_e(2), [3,1]); repmat(FHU_e(3), [3,1])];
            Ops2_e = [FHU_e(:); FHU_e(:); FHU_e(:)];
            Ops1_o = [repmat(FHU_o(1), [2,1]); repmat(FHU_o(2), [2,1])];
            Ops2_o = [FHU_o(:); FHU_o(:)];

            for itz = 1:nz
                [odisc, Adiscs(:, itz), sigmak] = getAdisc(nrgdata{itz}, [Ops1_e; Ops1_o], [Ops2_e; Ops2_o], ZF_eo, 'emin', emin, 'emax', emax, 'estep', 2*estep);
            end

            for ita = 1:size(Adiscs,1)
                Adisc = mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3);
                [~, Aconts{ita}] = getAcont(odisc, Adisc, sigmak, T / 5, 'alphaz', sigmarat / nz, ...
                                                'smin', log(Lambda) / 15, 'emin', emin, 'emax', emax, 'estep', estep);
            end
            Aconts_e = Aconts(1:9);
            Aconts_o = Aconts(10:13);

            for it1 = 1:3
                for it2 = 1:3
                    Acont2_e(:, it1, it2) = Aconts_e{3*(it1-1)+it2};
                end
            end

            for it1 = 1:2
                for it2 = 1:2
                    Acont2_o(:, it1, it2) = Aconts_o{2*(it1-1)+it2};
                end
            end


            %% calculate self-energy

            RVres_full_e = zeros(numel(ocont), 3, 3);
            RVres_full_o = zeros(numel(ocont), 2, 2);
            RVres_full_e(:,1,1) = RVres(:,1);
            RVres_full_o(:,1,1) = RVres(:,2);

            [SE_e, Aimp_e, ~] = SEtrick(ocont, Acont1_e, Acont2_e, Adisc2sum_e, Acont3_e, 'ozin', ocont, 'RhoV2in', RVres_full_e, 'epsd', epsd_e);
            [SE_o, Aimp_o, ~] = SEtrick(ocont, Acont1_o, Acont2_o, Adisc2sum_o, Acont3_o, 'ozin', ocont, 'RhoV2in', RVres_full_o, 'epsd', epsd_o);
            SE_save_e(:,:,:,itD) = SE_e;
            SE_save_o(:,:,:,itD) = SE_o;
            Aimp_save_e(:,:,:,itD) = Aimp_e;
            Aimp_save_o(:,:,:,itD) = Aimp_o;
            save([STG,'/SE_even.mat'], 'SE_save_e');
            save([STG,'/SE_odd.mat'], 'SE_save_o');
            save([STG,'/Aimp_even.mat'], 'Aimp_save_e');
            save([STG,'/Aimp_odd.mat'], 'Aimp_save_o');
            disptime('Self-energy and (improved) impurity spectral function saved');


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