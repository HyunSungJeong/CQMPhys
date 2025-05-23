function [rho, taus, lnZs] = run_XTRG(MPO, tau0, Nkeep, Nsweep, tauMax, varargin)
    % <Description>
    % Performs imaginary-time evolution using XTRG on a 1D system
    % (contour: -1i * Tau )
    %
    % <Input>
    % MPO : [cell vector] Tensors constituting the MPO representation of the system's Hamiltonian
    %
    % tau0 : [numeric] Initial value of the imaginary time
    %
    % tauMax : [numeric] Maximum (imaginary) time to be reached
    %
    % Nkeep : [numeric] Number of states to keep in variational squaring of density matrices (2-site update)
    %
    % Nsweep : [numeric] Number of sweeps (pairs of right --> left & left --> right sweeps) in the varational squaring of density matrices
    %
    % <Output>
    % rho : [cell vector] Tensors constituting the density matrix at final temperature
    %
    % taus : [numeric vector] The imaginary time(inverse temperatures) sequences that the XTRG passed
    %
    % lnZs : [numeric] Logarithms of the partition function measured at time instances given by 'taus'.

    %% Parse inputs

    if ~iscell(MPO)
        error('ERR: ''MPO'' must be a cell vector of tensors forming the MPO Hamiltonian');
    end

    if ~isnumeric(tau0) || ~isscalar(tau0)
        error('ERR: ''tau0'' must be a positive real number');
    elseif tau0 <= 0
        error('ERR: ''tau0'' must be positive');
    end

    if ~isnumeric(tauMax) || ~isscalar(tauMax)
        error('ERR: ''tauMax'' must be a positive real number');
    elseif tauMax <= 0
        error('ERR: ''tauMax'' must be positive');
    end

    if ~isnumeric(Nkeep) || ~isscalar(Nkeep)
        error('ERR: ''Nkeep'' must be a positive integer');
    elseif mod(Nkeep,1) ~= 0 || Nkeep < 1
        error('ERR: ''Nkeep'' must be a positive integer');
    end

    if ~isnumeric(Nsweep) || ~isscalar(Nsweep)
        error('ERR: ''Nsweep'' must be a positive integer');
    elseif mod(Nsweep,1) ~= 0 || Nsweep < 1
        error('ERR: ''Nsweep'' must be a positive integer');
    end

    tobj = tic2;

    LocDim = size(MPO{1}, 1);       % local physical dimension
    ChainLen = numel(MPO);          % chain length

    %% Initialize density matrix at tau0 as I-tau0*H
    rho = cell(ChainLen, 1);

    BondDim = size(MPO{1}, 4);  % Bond dimension of the Hamiltonian MPO. They are assumed to be uniform

    % leftmost MPO for initial density matrix
    rho{1} = zeros(LocDim, LocDim, 1, 2*BondDim);
    rho{1}(:,:,1,1) = eye(LocDim);
    rho{1}(:,:,:, BondDim+1:2*BondDim) = - MPO{1}*tau0;

    % rightmost MPO for initial density matrix
    rho{end} = zeros(LocDim, LocDim, 2*BondDim, 1);
    rho{end}(:,:,1,1) = eye(LocDim);
    rho{end}(:,:, BondDim+1:2*BondDim, 1) = MPO{end};

    % other MPOs for initial density matrix
    for itN = 2:ChainLen-1
        rho{itN} = zeros(LocDim, LocDim, 2*BondDim, 2*BondDim);
        for it = 1:BondDim
            rho{itN}(:,:,it,it) = eye(LocDim);
        end
        rho{itN}(:, :, BondDim+1:2*BondDim, BondDim+1:2*BondDim) = MPO{itN};
    end

    %% Iteratively square density matrix for imaginary-time evolution

    Nstep = round(log2(tauMax)-log2(tau0));
    taus = tau0 * (2.^(1:Nstep));
    lnZs = zeros(size(taus));

    for itS = 1:Nstep

        rho_dag = cell(ChainLen, 1);
        for itN = 1:ChainLen
            rho_dag{itN} = permute(conj(rho{itN}), [2,1,3,4]);
        end
        
        rho = varMul_MPO(rho, rho_dag, Nkeep, Nsweep);

        % compute partition function
        Z = 1;
        for itN = 1:ChainLen
            Z = contract(Z, 2, 2, rho{itN} ,4 ,3);
            Z = contract(Z, 4, [2 3], getIdentity(Z,3), 2, [2 1]);
        end
        lnZs(itS) = log(Z);

        % add the result from the last iteration, since the MPO gets normalized at every iteration to avoid divergence    
        if itS > 1
            lnZs(itS) = lnZs(itS) + 2*lnZs(itS-1);
        end

        % normalize the density matrix
        rho{1} = rho{1} / Z;

        % display information
        disptime(['Sweep #', sprintf('%d',itS), '/', sprintf('%d',Nstep), ': tau = ', sprintf('%.4g',taus(itS)), '/', ...
                                            sprintf('%.4g',taus(end)), ', lnZ = ', sprintf('%.4g', lnZs(itS))]);

    end % itS

    toc2(tobj,'-v');
    chkmem;

end
