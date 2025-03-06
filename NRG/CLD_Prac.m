function [ff,gg] = CLD(ozin,RhoV2in,Lambda,N,varargin)


    % default parameters
    estep = 10;
    emax = max(abs(ozin));
    emin = 1e-20;
    minrho = 0.01;
    fftol = 0.01;

    while ~isempty(varargin)
        switch varargin{1}
            case 'estep'
                estep = varargin{2};
                varargin(1:2) = [];
            case 'emax'
                emax = varargin{2};
                varargin(1:2) = [];
            case 'emin'
                emin = varargin{2};
                varargin(1:2) = [];
            case 'minrho'
                minrho = varargin{2};
                varargin(1:2) = [];
            case 'fftol'
                fftol = varargin{2};
                varargin(1:2) = [];
            otherwise
                error('ERR: check input!');
        end
    end

    % parsing input
    ozin = ozin(:);
    RhoV2in = RhoV2in(:);

    if isempty(ozin)
        error('ERR: Empty frequency input (1st input)');
    elseif isempty(RhoV2in)
        error('ERR: Empty hybridization input (2nd input)');
    elseif numel(ozin) ~= numel(RhoV2in)
        error('ERR: Different # of elements between frequency and hybridization inputs (1st & 2nd)');
    end

    % logarithmic frequency grid on which the hybridization function to be discretized
    xs = (ceil(log(emax)/log(Lambda)*estep):-1:floor(log(emin)/log(Lambda)*estep))/estep;
    xs = flipud(xs(:));
    oz = Lambda.^xs;

    rho1 = interp1(ozin,RhoV2in,+oz,'linear','extrap');
    rho1(rho1<0) = 0;
    [repE1,repT1] = CLD_1side(oz,rho1,estep,minrho);

    rho2 = interp1(ozin,RhoV2in,-oz,'linear','extrap');
    rho2(rho2<0) = 0;
    [repE2,repT2] = CLD_1side(oz,rho2,estep,minrho);

    if (numel(repE1)+numel(repE2)) < N
        fprintf(['WRN: Number of discretization intervals (= ', ...
            sprintf('%i',numel(repE1)+numel(repE2)),') is smaller than the chain length N (= ', ...
            sprintf('%i',N),'\n']);
        N2 = numel(repE1) + numel(repE2);
    else
        N2 = N;
    end

    ff = zeros(N2,1);
    gg = zeros(N2,1);

    % Lanczos tridiagonalization
    Xis = [repE1; -repE2];
    Gammas = [sqrt(repT1);sqrt(repT2)];
    H = [0,Gammas'; Gammas,diag(Xis)];

    U = zeros(size(H,1)*[1,1]);
    U(1,1) = 1;

    for itN = (1:N2)
        v = H*U(:,itN);
        v = v - U(:,1:itN)*(U(:,1:itN)'*v);
        v = v - U(:,1:itN)*(U(:,1:itN)'*v);

        ff(itN) = norm(v);

        if ff(itN) < (Lambda^(-itN/2)*fftol)
            itN = itN-1;
            break;
        end

        U(:,itN+1) = v/ff(itN);
        gg(itN) = U(:,itN)'*H*U(:,itN);     % In the solution, U(:,itN+1) instead of U(:,itN) ... WHY??
    end

    ff(itN+1:end) = [];
    gg(itN+1:end) = [];
end

function [repE,repT] = CLD_1side (oz,rho,estep,minrho)
    % Obtain the representative energies (repE, \mathcal{E} in Campo2005) and
    % the integral of the hybridization function (repT) for each discretization
    % interval, for either positive or negative energy side.
    
    ids0 = find(rho >= minrho*max(rho),1,'last')-estep;
    
    % oz(ids) are the discretization grid points at which the input
    % hybridization function is splot
    ids = [numel(oz),(ids0:-estep:1)]; % The tail outside the effective band edge is integrated into the first 
    
    repT = zeros(numel(ids)-1,1);
    repE = zeros(size(repT));
    
    for itx = (1:numel(repT))
        ozp = oz(ids(itx+1):ids(itx));
        rhop = rho(ids(itx+1):ids(itx));
    
        % determine repE and repT by using numerical integration
        repT(itx) = sum((rhop(2:end)+rhop(1:end-1)).*(ozp(2:end)-ozp(1:end-1)))/2;
        %repE(itx) = sum((rhop(2:end)./ozp(2:end)+rhop(1:end-1)./ozp(1:end-1)).*(ozp(2:end)-ozp(1:end-1)))/2;
        repE(itx) = (rhop(end)-rhop(1)) + ...
        sum( (ozp(2:end).*rhop(1:end-1) - ozp(1:end-1).*rhop(2:end)) ./ ...
            (ozp(2:end) - ozp(1:end-1)) .* log(abs(ozp(2:end)./ozp(1:end-1))) );    % WHY??
    end
    
    repE = repT./repE;
    
    % % % % TODO (end) % % % %
    
end