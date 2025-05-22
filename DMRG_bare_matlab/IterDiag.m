function [MPS, E0] = IterDiag(Nkeep, ChainLen, varargin)
    % <Description>
    % Iteratively diagonalizes the Hamiltonian of a given 1-D quantum system
    % and yields its lowest energy states in the form of left-normalized MPS
    % 
    % <Input>
    % Nkeep : [numeric] number of states to keep in the iterative diagonalization
    % ChainLen : [numeric] chain length of the 1D system
    %
    % <Option>
    % 'TBchain', ... : [numeric vector] length L-1 vector containing the
    %                       hopping parameters of a length L tight-binding chain
    %
    % '1DHubbard', ... : [numeric] Hubbard U of the 1D Hubbard model. 
    %                       The hopping parameters are taken to be all 1
    %
    % 'deps', ... : [numeric] The numerical tolerance for restoring degeneracy.
    %                       If the spacing between two energy levels are less than this tolerance, they are regarded as degenerate
    %                           (Default: 1e-13)
    % 
    % 'chooseGS', : If used, the rightmost MPS of the output MPS is dummy, so the MPS represents the ground state 
    %                           (Default: not used)
    %
    % '-v' : If used, information of the iterative diagonalization is displayed in every iteration
    %                           (Defualt: not used)
    %
    % <Output>
    % MPS : [cell vector] cell vector containing the left-normalized isometries 
    %                       forming MPS representing the lowest energy states
    % E0 : [numeric vector] vector of length = Nkeep, containing the energies 
    %                       of the lowest-lying states in increasing order

    %% Parse inputs

    if ~isnumeric(Nkeep)
        error('ERR: ''Nkeep'' must be a positive integer');
    elseif mod(Nkeep,1) ~= 0 || Nkeep < 1
        error('ERR: ''Nkeep'' must be a positive integer');
    end

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

    % default options
    System = NaN;
    deps = 1e-13;
    showInfo = false;
    chooseGS = false;

    while ~isempty(varargin)
        switch varargin{1}
            case 'TBchain'
                if isnumeric(varargin{2})
                    if numel(varargin{2}) == ChainLen-1
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

            case 'deps'
                if isnumeric(varargin{2})
                    deps = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''deps'' must be a number');
                end

            case 'chooseGS'
                chooseGS = true;
                varargin(1) = [];

            case '-v'
                showInfo = true;
                varargin(1) = [];

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
        if isscalar(t_hop)
            t_hop = t_hop*ones(ChainLen-1, 1);
        elseif numel(t_hop) ~= ChainLen-1
            error('ERR: the number of hopping parameters must be consistent with the chain length');
        end
    end


    MPS = cell(ChainLen, 1);
    
    if isequal(System, 'TBchain')   % tight-binding chain
        %% Perform iterative diagonalization for the tight-binding chain

        % define local operators
        [F,Z,I] = getLocalSpace('Fermion');

        % initial Hamiltonian and isometry
        Hprev = 0;
        Aprev = 1;

        for it = 1:ChainLen
            % Add local site to the Hamiltonian and define the identity
            % tensor connecting the previous Hilbert space and the new local space
            Anow = getIdentity(Aprev, 2, I, 2, [1 3 2]);
            Hnow = updateLeft(Hprev, 2, Anow, [], [], Anow);

            if it > 1   
                % add intersite hopping
                Hhop = conj(t_hop(it-1))*updateLeft(Fprev, 3, Anow, F'*Z, 3, Anow);
                Hnow = Hnow + Hhop + Hhop';
            end

            % diagonalize Hnow
            [V,E0] = eig((Hnow+Hnow')/2);

            [E0,idx] = sort(diag(E0), 'ascend');
            V = V(:,idx);

            E_tr = E0(min(numel(E0),Nkeep)) + deps;    % truncation energy
            % truncate and define isometry
            V = V(:, E0<E_tr);
            E0 = E0(E0<E_tr);

            % update isometry and Hamiltonian
            Aprev = contract(Anow, 3, 2, V, 2, 1, [1,3,2]);
            if it == ChainLen && chooseGS
                MPS{it} = contract(Anow, 3, 2, V(:,1), 2, 1, [1,3,2]);
            else
                MPS{it} = Aprev;
            end
            Hprev = diag(E0);

            % update 'Hprev' for the next iteration
            Fprev = updateLeft([], [], Aprev, F, 3, Aprev);

            if showInfo
                disptime(['#', sprintf('%02d/%02d',[it, ChainLen]), ' : ', 'Nkeep=', sprintf('%d/%d', [size(Hprev,2), size(Hnow,2)])]);
            end

        end

    else
        %% Perform iterative diagonalization for the 1D Hubbard chain

        % define local operators
        [F,Z,~,I] = getLocalSpace('FermionS');

        N = zeros(4,4,2);   % number operator for spin up/down
        N(:,:,1) = F(:,:,1)' * F(:,:,1);    % spin up 
        N(:,:,2) = F(:,:,2)' * F(:,:,2);    % spin down
        Ntot = sum(N, 3);                   % total number op. (spin up+down)

        % initial Hamiltonian and isometry
        Hprev = 0;
        Aprev = 1;

        for it = 1:ChainLen
            % Add local site to the Hamiltonian and define the identity
            % tensor connecting the previous Hilbert space and the new local space
            Anow = getIdentity(Aprev, 2, I, 2, [1 3 2]);
            Hnow = updateLeft(Hprev, 2, Anow, [], [], Anow);

            if it > 1   
                % add intersite hopping
                Hhop = updateLeft(Fprev(:,:,1), 3, Anow, F(:,:,1)'*Z, 3, Anow);
                Hhop = Hhop + updateLeft(Fprev(:,:,2), 3, Anow, F(:,:,2)'*Z, 3, Anow);
                Hnow = Hnow + Hhop + Hhop';

                % add local Hubbard repulsion
                HU = updateLeft([], [], Anow, (Ntot-I)^2, 2, Anow);
                Hnow = Hnow + (HU + HU')/2;
            end

            % diagonalize Hnow
            [V,E0] = eig((Hnow+Hnow')/2);

            [E0,idx] = sort(diag(E0), 'ascend');
            V = V(:,idx);

            E_tr = E0(min(numel(E0),Nkeep)) + deps;    % truncation energy
            % truncate and define isometry
            V = V(:, E0<E_tr);
            E0 = E0(E0<E_tr);

            % update isometry and Hamiltonian
            Aprev = contract(Anow, 3, 2, V, 2, 1, [1,3,2]);
            if it == ChainLen && chooseGS
                MPS{it} = contract(Anow, 3, 2, V(:,1), 2, 1, [1,3,2]);
            else
                MPS{it} = Aprev;
            end
            Hprev = diag(E0);

            % update 'Hprev' for the next iteration
            Fprev = updateLeft([], [], Aprev, F, 3, Aprev);

            if showInfo
                disptime(['#', sprintf('%02d/%02d',[it, ChainLen]), ' : ', 'Nkeep=', sprintf('%d/%d', [size(Hprev,2), size(Hnow,2)])]);
            end

        end

    end

end