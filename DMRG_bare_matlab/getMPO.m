function MPO = getMPO(ChainLen, varargin)
    % <Description>
    % Generates an MPO form of the given 1D quantum system
    %
    % <Input>
    % ChainLen : [integer] chain length of the 1D system
    %
    % <Options>
    % 'TBchain', ... : [numeric vector] length L-1 vector containing the
    %                       hopping parameters of a length L tight-binding chain
    % '1DHubbard', ... : [numeric] Hubbard U of the 1D Hubbard model. 
    %                       The hopping parameters are taken to be all 1
    % 'XXZ', ... [numeric] Anisotropy \Delta of the XXZ Heisenberg chain
    % 'PBC' : If used, periodic boundary condition is applied.
    %           (Default: not used)
    %
    % <Output>
    % MPO : [cell array of tensors] length L cell array containing the MPOs constituting the Hamiltonian

    %% Parse input

    if isnumeric(ChainLen)
        if mod(ChainLen, 1) == 0
            if ChainLen <= 1
                error('ERR: chain length must be larger than 1');
            end 
        else
            error('ERR: chain length must be an integer larger than 1');
        end
    else
        error('ERR: chain length must be an integer larger than 1');
    end

    %% Parse options

    PBC = false;    % Default: do not apply periodic boundary condition
    System = NaN;

    while ~isempty(varargin)
        switch varargin{1}
            case 'TBchain'
                if isnumeric(varargin{2})
                    if numel(varargin{2}) == ChainLen-1 || numel(varargin{2}) == ChainLen
                        System = 'TBchain';
                        t_hop = varargin{2};
                        t_hop = t_hop(:);
                        varargin(1:2) = [];

                    elseif isscalar(varargin{2})
                        System = 'TBchain';
                        t_hop = varargin{2};
                        varargin(1:2) = [];

                    else
                        error('ERR: the number of hopping parameters must be consistent with the chain length');
                    end
                else
                    error('ERR: The hopping parameters for ''TBchain'' must be a numeric vector');
                end

            case '1DHubbard'
                if isnumeric(varargin{2})
                    System = '1DHubbard';
                    U = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: Hubbard U for ''1DHubbard'' must be a number');
                end

            case 'XXZ'
                if isnumeric(varargin{2})
                    System = 'XXZ';
                    Delta = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: Delta for ''1DHubbard'' must be a number');
                end

            case 'PBC'
                PBC = true;

            otherwise
                if ischar(varargin{1})
                    error(['ERR: unknown input ''',varargin{1},'''']);
                else
                    error('ERR: unknown input');
                end

        end % switch-case
    end % while

    if isnan(System)
        error('ERR: Choose the physical system: either ''TBchain'' or ''1DHubbard''');
    end

    if isequal(System, 'TBchain')
        if PBC
            if isscalar(t_hop)
                t_hop = t_hop*ones(ChainLen, 1);
            elseif numel(t_hop) ~= ChainLen
                error('ERR: the number of hopping parameters must be consistent with the chain length');
            end
    
        else
            if isscalar(t_hop)
                t_hop = t_hop*ones(ChainLen-1, 1);
            elseif numel(t_hop) ~= ChainLen-1
                error('ERR: the number of hopping parameters must be consistent with the chain length');
            end
    
        end
    end

    %% Construct MPOs

    switch System
        case 'TBchain'      % tight-binding chain

            [F,Z,I] = getLocalSpace('Fermion');
            MPO = cell(1, ChainLen);

            if ~PBC
                t_hop = [t_hop; 0];     
                % put dummy hopping parameter at the end for convenience 
                % this parameter won't be present in the final MPO
            end

            for it = 1:ChainLen
                MPO{it} = zeros(2,2,4,4);
                MPO{it}(:,:,1,1) = I;
                MPO{it}(:,:,2,1) = Z*F;
                MPO{it}(:,:,3,1) = F'*Z;
                MPO{it}(:,:,4,2) = t_hop(it) * F';
                MPO{it}(:,:,4,3) = conj(t_hop(it)) * F;
                MPO{it}(:,:,4,4) = I;
            end

            if ~PBC
                MPO{1} = MPO{1}(:,:,end,:);
                MPO{end} = MPO{end}(:,:,:,1);
            end

        case '1DHubbard'    % 1D Hubbard model

            [F,Z,~,I] = getLocalSpace('FermionS');
            MPO = cell(1, ChainLen);
            N = zeros(4,4,2);   % number operator for spin up/down
            N(:,:,1) = F(:,:,1)' * F(:,:,1);
            N(:,:,2) = F(:,:,2)' * F(:,:,2);

            for it = 1:ChainLen
                MPO{it} = zeros(4,4,6,6);
                MPO{it}(:,:,1,1) = I;
                MPO{it}(:,:,2,1) = Z*F(:,:,1);
                MPO{it}(:,:,3,1) = F(:,:,1)'*Z;
                MPO{it}(:,:,4,1) = Z*F(:,:,2);
                MPO{it}(:,:,5,1) = F(:,:,2)'*Z;
                MPO{it}(:,:,6,1) = (U/2)*(sum(N,3) - I)^2;

                MPO{it}(:,:,6,2) = F(:,:,1)';
                MPO{it}(:,:,6,3) = F(:,:,1);
                MPO{it}(:,:,6,4) = F(:,:,2)';
                MPO{it}(:,:,6,5) = F(:,:,2);
                MPO{it}(:,:,6,6) = I;
            end

            if ~PBC
                MPO{1} = MPO{1}(:,:,end,:);
                MPO{end} = MPO{end}(:,:,:,1);
            end

        case 'XXZ'    % 1D Hubbard model

            Sz = [1,0;0,-1]/2;      % spin-z operator
            S_r = [0,1;0,0];        % spin raising operator
            S_l = [0,0;1,0];        % spin lowering operator
            I = eye(2);             % identity operator
            MPO = cell(1, ChainLen);

            for it = 1:ChainLen
                MPO{it} = zeros(2,2,5,5);
                MPO{it}(:,:,1,1) = I;
                MPO{it}(:,:,2,1) = S_l / sqrt(2);
                MPO{it}(:,:,3,1) = Sz;
                MPO{it}(:,:,4,1) = S_r / sqrt(2);

                MPO{it}(:,:,5,2) = S_r / sqrt(2);
                MPO{it}(:,:,5,3) = Delta*Sz;
                MPO{it}(:,:,5,4) = S_l / sqrt(2);
                MPO{it}(:,:,5,5) = I;
            end

            if ~PBC
                MPO{1} = MPO{1}(:,:,end,:);
                MPO{end} = MPO{end}(:,:,:,1);
            end

    end
end