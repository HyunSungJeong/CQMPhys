function DMFT_QI_2CK_lat_fullHyb(parfn,varargin)

    setenv('RC_STORE',go('rcs'));           % '/project/hyunsung/RCStore'
    partot = job_func_preamble(parfn, varargin{:});

    [PE, ndfix, T, Lambda, nz, sigmarat, Nkeep, D, U, V, t_0, phi_div_pi, initS] = loadvar(partot, ...
            {'PE', 'ndfix', 'T', 'Lambda', 'nz', 'sigmarat', 'Nkeep', 'D', 'U', 'V', 't_0', 'phi_div_pi', 'initS'}, ...
                {[], [], [], [], [], [], [], [], [], [], [], [], []});

    parpool('IdleTimeout',3500);

    %% DMFT parameters
    cvgth = 1e-3;   % convergence threshold, option for updateHyb (Default: 1e-3)
    amix = 1;       % linear mixing parameter for updateHyb (0<=amix<=1, Default: 1)
    Bro = 0;        % Broyden parameter for updateHyb (# of last iterations to be used for Broyden scheme, Default: 0)
    Ndmft = 50;     % maximum number of DMFT iterations before breaking the loop
    ieta = 1i * 1e-4;   % convergence generating factor
    abstol = 1e-8;      % absolute error tolerance for DOS integration. Option for MATLAB function "integral"
    fullHyb = true;     % if true, use full 5x5 hybridization function. if false, assume that the hybridization function is diagonal(as it should be)
    
    %% Variables used in DMFT loop
    RVres = cell(2,1);      % hybridization function to be used as the input of the next DMFT iteration

    %% NRG parameters
    Lambda = 4;
    N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);  % Wilson chain length
    nz = 2;                                         % Number of z's to be averaged
    Nfit = round(-2*log(1e-8)/log(Lambda));         % The Wilson chain site index after which the Wilson chain parameters are obtained via extrapolation
    emin = T;           % Minimum absolute value of frequency grid
    emax = 1e3;         % Maximum absolute value of frequency grid
    estep = 250;        % # of steps to increase frequency per decade(x10), used in getAdisc and getAcont
                        % Acont: estep, Adisc: 2*estep
    Etrunc = [];
    ETRUNC = [];

    %% System parameters
    t = t_0*cos(phi);       % intra-cell hopping amplitude t
    t_p = t_0*sin(phi);     % intra-cell hopping amplitude t'
    RhoB = @(epsilon) heaviside(D - abs(epsilon))*(2 / pi / D)*sqrt(1 - (epsilon/D)^2);     % semi-elliptic DOS for noninteracting Bethe lattice
    phi = pi*phi_div_pi;    % Distortion angle \phi

    % (System parameters from parameter files)
    % t_0 : intra-cell hopping amplitude t_0
    % U : on-site Coulomb repulsion
    % D : half bandwidth
    % V : inter-site hopping amplitude

    %% frequency grid
    ocont = getAcont(0, 0, 0, 0, 'emin', emin, 'emax', emax, 'estep', estep);       % logarithmic frequency grid for smoothened spectral functions
    docon = zeros(numel(ocont), 1);             % Delta-omega's on discrete frequency grid ocont, for omega integration
    docon(1) = 0.5 * (ocont(2) - ocont(1));
    docon(2:end - 1) = 0.5 * (ocont(3:end) - ocont(1:end - 2));
    docon(end) = 0.5 * (ocont(end) - ocont(end - 1));

    %% TODO: define SEs, RhoV2s, and system Hamiltonian       
    RhoV2out = zeros(numel(ocont), 5, 5, Ndmft);    % Gamma (= hybridization function) / pi : before Broyden, for whole unit cell (omega, 5, 5, Ndmft)             
    RhoV2out_1e = zeros(numel(ocont), Ndmft);       % Gamma (= hybridization function) / pi : before Broyden, for 1e orbital(1,1 component of RhoV2out), (omega, Ndmft)
    RhoV2out_1o = zeros(numel(ocont), Ndmft);       % Gamma (= hybridization function) / pi : before Broyden, for 1o orbital(4,4 component of RhoV2out), (omega, Ndmft)

    RhoV2in = RhoV2out;             % hybridization function after Broyden, for whole unit cell
    RhoV2in_1e = RhoV2out_1e;       % hybridization function after Broyden, for even sector
    RhoV2in_1o = RhoV2out_1o;       % hybridization function after Broyden, for odd sector

    Gfloc_e = zeros(3,3,numel(ocont),Ndmft);        % local lattice GF. for even sector
    Gfloc_o = zeros(2,2,numel(ocont),Ndmft);        % local lattice GF. for odd sector
    SEs = RhoV2out; % self-energy
    nures = zeros(2, 2, Ndmft); %sublattice, spin
    mures = zeros(1, Ndmft);        % Chemical Potential for each iteration
    mures(1) = 0;                   % half filling initial condition

    %% define local operators
    [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge,SU2spin', 'NC', 5);
    [FF, ZF, EF] = setItag('s00', 'op', FF, ZF, IF.E);
    NF = quadOp(FF, FF, []);

    % Hamiltonian - interacting part(HU)
    HU = QSpace();
    for itm = (1:5)
        HU = HU + (U/2)*(contract(NF(itm),'!1',NF(itm)) - NF(itm));
    end
    HU = HU - (U/2)*sum(NF);
    
    A0 = getIdentity(setItag('L00', getvac(EF)), 2, EF, 2, 'K00*', [1 3 2]); % isometry
    for ito = (1:numel(FF))
        FHU(ito) = contract(FF(ito), '!1', HU, [1 3 2]) - contract(HU, '!1', FF(ito)); % commutator [FK,HU], to be used for the self-energy trick
    end

    %% TODO: Initial Self energy
    if isempty(initS)
        SEs{1}(:, :, 1) = -0.1i * 1; % Hilbert transform
        SEs{2}(:, :, 1) = SEs{1}(:, :, 1);
    else
        SEs{1}(:, :, 1) = interp1(initS.ocont, initS.SE{1}, ocont);
        SEs{2}(:, :, 1) = interp1(initS.ocont, initS.SE{2}, ocont);
    end

    %% DMFT loop
    for itD = 1:Ndmft

        num_threads_SL(PE);
        muc = mures(it);

        disp2('-----------------------------------------');
        disptime(['DMFT iteration No.', sprintf('%02d',it)]);

        Gftmp_e = zeros(3, 3, numel(ocont));     % temporary local lattice GF., even sector
        Gftmp_o = zeros(2, 2, numel(ocont));     % temporary local lattice GF., odd sector

        disptime('Starting DOS integration');
        parfor ito = 1:numel(ocont)

            for it1 = 1:3
                for it2 = 1:3
                    Gftmp_e(it1, it2, ito) = integral(@(epsilon) RhoB(epsilon)*Gfk_e(epsilon, ocont(ito), SE_e, U, V, t_0, t, t_p, muc, ieta, it1, it2), -D, +D, 'AbsTol', abstol);
                end
            end

            for it1 = 1:2
                for it2 = 1:2
                    Gftmp_o(it1, it2, ito) = integral(@(epsilon) RhoB(epsilon)*Gfk_o(epsilon, ocont(ito), SE_e, U, V, t_0, t, t_p, muc, ieta, it1, it2), -D, +D, 'AbsTol', abstol);
                end
            end 
        end % parfor
        disptime('DOS integration complete');
        Gfloc_e(:,:,:,it) = Gftmp_e;
        Gfloc_o(:,:,:,it) = Gftmp_o;

        % impose DMFT self-consistency condition(Gfloc --> Hyb)
        if fullHyb          % use full(5x5) hybridization function
            SE_full = blkdiag(SE_e, SE_o);      % full 5x5 self-energy

            G_SEl_inv = @(omega) [omega + ieta + muc + U/2, -(t + t_p), 0, 0, 0;
                                -(t + t_p), omega + ieta + muc + U/2, -sqrt(2)*t_0, 0, 0;
                                0, -sqrt(2)*t_0, omega + ieta + muc + U/2, 0, 0;
                                0, 0, 0, omega + ieta + muc + U/2, -(t - t_p);
                                0, 0, 0, -(t - t_p), omega + ieta + muc + U/2] - SE_full;
            % (inverse of local lattice GF. without SE.)

            Gfloc_e_inv = inv(Gftmp_e);     % inverse of local lattice GF., even sector
            Gfloc_o_inv = inv(Gftmp_o);     % inverse of local lattice GF., odd sector
            RhoV2out(:,:,:,it) = -imag( (ocont(:) + zeros(1,5,5)) + muc + U/2 - SE_e(1,1) - Gfloc_e_inv(1,1) ) / pi;

        else        % Assume that the hybridization function is diagonal(as it should be)
            Gfloc_e_inv = inv(Gftmp_e);     % inverse of local lattice GF., even sector
            Gfloc_o_inv = inv(Gftmp_o);     % inverse of local lattice GF., odd sector

            RhoV2out_1e(:,it) = -imag( ocont + muc + U/2 - SE_e(1,1) - Gfloc_e_inv(1,1) ) / pi;
            RhoV2out_1o(:,it) = -imag( ocont + muc + U/2 - SE_o(1,1) - Gfloc_o_inv(1,1) ) / pi;
        end

        

    end

end