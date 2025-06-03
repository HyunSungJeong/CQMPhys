function Corr = getSpinCorr(MaxDist, Delta, ChainLen, Nkeep, varargin)
    % <Description>
    % Computes spin-spin correlation function from the ground state MPS 
    % of a single-orbital spin-1/2 quantum chain chain obtained by DMRG
    %
    % <Input>
    % MaxDist : [numeric] Maximum distance to compute spin-spin correlation
    %
    % Delta : [numeric] The spin-z interaction of the spin-1/2 XXZ Heisenberg chain
    %
    % ChainLen : [numeric] Chain length of the spin-1/2 XXZ Heisenberg chain
    %
    % Nkeep : [numeric] Number of states kept in the DMRG calculation
    %
    % <Option>
    % 'range', ... : [numeric vector] two-component vector containing the
    %                       leftmost and rightmost site number to compute correlators.
    %                   (Default: [1, ChainLen])
    %
    % 'NumSamples', ... : [numeric] Number of sample extracted in the DMRG calculation
    %                       (Default: 10000)
    %
    % 'Nsweep', ... : [numeric] Number of DMRG sweeps
    %                       (Default: 10)
    %
    % <Output>
    % Corr : [cell] Cell vector containing spin-spin correlators
    %               Corr{n} is a vector containing the distance n spin-spin correlators,
    %               starting from left to right

    %% Parse inputs

    if ~isnumeric(MaxDist)
        error('ERR: ''MaxDist'' must be a number');
    elseif ~isscalar(MaxDist)
        error('ERR: ''MaxDist'' must be a single number');
    elseif mod(MaxDist, 1) ~= 0 || MaxDist <= 0
        error('ERR: ''MaxDist'' must be a positive integer');
    end

    if ~isnumeric(Delta)
        error('ERR: ''Delta'' must be a number');
    elseif ~isscalar(Delta)
        error('ERR: ''Delta'' must be a single real number');
    end

    if ~isnumeric(ChainLen)
        error('ERR: ''ChainLen'' must be a number');
    elseif ~isscalar(ChainLen)
        error('ERR: ''ChainLen'' must be a single number');
    elseif mod(ChainLen, 1) ~= 0 || ChainLen < 2
        error('ERR: ''ChainLen'' must be an integer larger than one');
    end

    if ~isnumeric(Nkeep) || ~isscalar(Nkeep)
        error('ERR: ''Nkeep'' must be a positive integer');
    elseif mod(Nkeep,1) ~= 0 || Nkeep < 1
        error('ERR: ''Nkeep'' must be a positive integer');
    end

    %% Parse options

    % default options
    range = [1, ChainLen];
    NumSamples = 10000;
    Nsweep = 10;

    while ~isempty(varargin)
        switch varargin{1}
            case 'range'
                if ~isnumeric(varargin{2})
                    error('ERR: ''range'' must be a numeric vector');
                elseif numel(varargin{2}) ~= 2
                    error('ERR: ''range'' must be a length 2 vector');
                elseif mod(varargin{2}(1),1) ~= 0 || mod(varargin{2}(2),1) ~= 0
                    error('ERR: The components of ''range'' must be integers');
                elseif varargin{2}(1) >= varargin{2}(2)
                    error('ERR: The first component of ''range'' must be smaller than the second component');
                elseif varargin{2}(1) <= 0 || varargin{2}(2) <= 0
                    error('ERR: The components of ''range'' must be positive');
                else
                    range = varargin{2};
                    varargin(1:2) = [];
                end

            case 'NumSamples'
                if ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
                    error('ERR: ''NumSamples'' must be a positive integer');
                elseif mod(varargin{2},1) ~= 0 || varargin{2} < 1
                    error('ERR: ''NumSamples'' must be a positive integer');
                else
                    NumSamples = varargin{2};
                    varargin(1:2) = [];
                end

            case 'Nsweep'
                if ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
                    error('ERR: ''Nsweep'' must be a positive integer');
                elseif mod(varargin{2},1) ~= 0 || varargin{2} < 1
                    error('ERR: ''Nsweep'' must be a positive integer');
                else
                    Nsweep = varargin{2};
                    varargin(1:2) = [];
                end

            otherwise
                if ischar(varargin{1})
                    error(['ERR: unknown input ''',varargin{1},'''']);
                else
                    error('ERR: unknown input');
                end

        end % switch-case
    end % while

    %% Compute correlators

    path = '/data/hyunsung/DMRG_XXZ/MPSdata';

    FilePath = [path, filesep, 'MPS_GS_Delta=', sprintf('%.15g',Delta) ,'_ChainLen=', sprintf('%d',ChainLen), '_NumSamples=', ...
                                    sprintf('%d',NumSamples), '_Nkeep=', sprintf('%d',Nkeep), '_Nsweep=', sprintf('%d',Nsweep), '.mat'];

    % load MPS
    MPS = load(FilePath);
    MPS = MPS.MPS_GS;

    % define local operators
    Sz = [1,0;0,-1]/2;      % spin-z operator
    S_r = [0,1;0,0];        % spin raising operator
    S_l = [0,0;1,0];        % spin lowering operator
    I = eye(2);             % identity operator

    % cell variable to store spin-spin correlation functions
    Corr = cell(1, MaxDist);

    for itL = range(1):range(2)-1

        % show status
        disptime(['Computing correlation function starting from site # ', sprintf('%d', itL), '/', sprintf('%d', ChainLen)]);

        % convert to site-canonical form
        MPS = getCanonForm(MPS, itL);

        % Norm of the MPS
        Norm = updateLeft(getIdentity(MPS{itL}, 1), 2, MPS{itL}, [], [], MPS{itL});
        Norm = contract(Norm, 2, [1,2], getIdentity(Norm,2),  2, [1,2]);

        % left part of the MPS-operator-MPS overlap
        Cleft = getIdentity(MPS{itL}, 1);           
        Cleft = updateLeft(Cleft, 2, MPS{itL}, Sz, 3, MPS{itL});

        for itR = itL+1 : min(range(2), itL+MaxDist)

            OverLap = updateLeft(Cleft, 3, MPS{itR}, Sz, 3, MPS{itR});
            OverLap = contract(OverLap, 2, [1,2], getIdentity(MPS{itR}, 2), 2, [1,2]);

            Corr{itR-itL} = [Corr{itR-itL}, OverLap/Norm];

            Cleft = updateLeft(Cleft, 3, MPS{itR}, [], [], MPS{itR});

        end % itR
    end % itL

end