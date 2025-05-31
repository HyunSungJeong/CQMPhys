function [E0, E_iter, MPS, Sv] = DMRG_GS(MPO, Nkeep, Nsweep, varargin)
    % <Description>
    % Performs Lanczos tridiagonalization
    %
    % <Input>
    % MPS_init : [cell vector of tensors] Tensors forming the initial MPS, from left to right
    %
    % MPO : [cell vector of tensors] Tensors forming the MPO Hamiltonian, form left to right
    %
    % Nkeep : [numeric] Maximum bond dimension of the MPS to be kept
    %
    % Nsweep : [numeric] Number of round trips of sweeps to iterate in the ground state search
    %                   The right --> left and left --> right sweeps will be done Nsweep times, respectively
    %
    % <Options>
    % 'MPSinit', ... : [char or cell vector] (Default: 'IterDiag')
    %                       1. character input;
    %                           If 'IterDiag', use iterative diagonalization result to initialize MPS
    %                           If 'Rand', use random MPS as initial MPS
    %                       2. cell vector input;
    %                           The cell vector must contain the rank-4 tensors constituting the MPO
    %
    % 'KrylovDim', ... : [numeric] The number of Krylov vectors to be used in the Lanczos method
    %                           (Default: 5)
    % 'update', ... : [char] '1site' or '2site'. 
    %                       If '1site', 1-site update DMRG is performed 
    %                       If '2site', 2-site update DMRG is performed
    %                           (Default: 1-site update)
    % 'KrylovTol', ... : [numeric] The tolerance for the norm in the Lanczos method
    %                           (Default: 1e-8)
    %
    % 'margin', ... : [numeric] In each iteration of iterative diagonalization, up to Nkeep*(1+margin) states are kept
    %                       This option is used in the iterative diagonalization, 
    %                       so it is meaningless to apply this option when 'MPSinit' is 'IterDiag'.
    %                       If the spacing between two energy levels are less than this tolerance, they are regarded as degenerate
    %                           (Default: 0.1)
    % 
    % 'Boundary', ... [char] 'open' or 'periodic'.
    %                       If 'open', open boundary condition is imposed
    %                       If 'periodic', periodic boundary condition is imposed
    %                           (Default: 'open')
    % 
    % '-v' : is used, the DMRG sweep information are displayed
    %                       (Default: not used)
    %
    % <Output>
    % E0 : [numeric] The energy of the final ground state MPS
    %
    % E_iter : [numeric array] E_iter(N, M) is the energy of the MPS in
    %                       the N-th iteration in the M-th sweep
    %
    % MPS : [cell vector of tensors] Tensors forming the ground state MPS
    %
    % Sv : [cell vector] length ChainLen+1 cell vector. Sv{n} is a vector containing 
    %                   the singular values on the bond between site n-1 and site n. 
    %                   Sv{1} and Sv{ChainLen+1} are the norm of the MPS in 1-site update.
    %                   In 2-site update, Sv{1} and Sv{ChainLen+1} are empty
    
    %% Parse inputs

    if ~iscell(MPO)
        error('ERR: ''MPO'' must be a cell vector of tensors forming the MPO Hamiltonian');
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

    %% Parse options

    % Default options
    MPSinit_Opt = 'IterDiag';
    KrylovDim = 5;
    updateOpt = '1site';
    KrylovTol = 1e-8;
    margin = 0.1;
    margin_specified = false;
    Boundary = 'open';
    display = false;

    while ~isempty(varargin)
        switch varargin{1}
            case 'MPSinit'
                if ischar(varargin{2})
                    if ~ismember(varargin{2}, {'IterDiag', 'Rand'})
                        if ischar(varargin{2})
                            error(['ERR: Unknown MPS initialization option ''',varargin{2},'''']);
                        else
                            error('ERR: Unknown MPS initialization option');
                        end
                    else
                        MPSinit_Opt = varargin{2};
                        varargin(1:2) = [];
                    end

                elseif iscell(varargin{2})
                    if numel(MPO) ~= numel(varargin{2})
                        error('ERR: Length of MPO and initial MPS is inconsistent');
                    else
                        MPS_init = varargin{2};
                        varargin(1:2) = [];
                        MPSinit_Opt = [];
                    end

                else
                    error('ERR: ''MPS_init'' must be either ''IterDiag'', ''Rand'', or a cell vector of MPS tensors');
                end

            case 'KrylovDim'
                if ~isnumeric(KrylovDim) || ~~isscalar(KrylovDim)
                    error('ERR: ''Krylov'' must be a numeric scalar');
                else
                    KrylovDim = varargin{2};
                    varargin(1:2) = [];
                end

            case 'update'
                if ~ismember(varargin{2}, {'1site', '2site'})
                    if ischar(varargin{2})
                        error(['ERR: Unknown update option ''',varargin{2},'''']);
                    else
                        error('ERR: Unknown update option');
                    end
                else
                    updateOpt = varargin{2};
                    varargin(1:2) = [];
                end

            case 'KrylovTol'
                if ~isnumeric(varargin{2})
                    error('ERR: ''KrylovTol'' must be a positive number');
                elseif ~isscalar(varargin{2})
                    error('ERR: ''KrylovTol'' must be a single positive number');
                else
                    KrylovTol = varargin{2};
                    varargin(1:2) = [];
                end

            case 'margin'
                if isnumeric(varargin{2})
                    margin = varargin{2};
                    varargin(1:2) = [];
                    margin_specified = true;
                else
                    error('ERR: ''margin'' must be a number');
                end

            case 'Boundary'
                if ischar(varargin{2})
                    if ~ismember(varargin{2}, {'open', 'periodic'})
                        if ischar(varargin{2})
                            error(['ERR: Unknown boundary condition ''',varargin{2},'''']);
                        else
                            error('ERR: Unknown boundary condition');
                        end
                    else
                        Boundary = varargin{2};
                        varargin(1:2) = [];
                    end

                else
                    error('ERR: Unknown option for ''Boundary''');
                end

            case '-v'
                display = true;
                varargin(1) = [];

            otherwise 
                if ischar(varargin{1})
                    error(['ERR: Unknown input ''',varargin{1},'''']);
                else
                    error('ERR: Unknown input');
                end
        end % switch-case
    end % while

    ChainLen = numel(MPO);

    if isempty(MPSinit_Opt) || isequal(MPSinit_Opt, 'Rand')
        if margin_specified
            disp('WRN: In the current option of initializing MPS, ''margin'' option is meaningless');
        end
    end

    if isequal(Boundary, 'periodic') && isequal(MPSinit_Opt, 'IterDiag')
        error('ERR: when periodic boundary condition is imposed, the MPS cannot be initialized by iterative diagonalization');
    end

    tobj = tic2;


    %% Initialize MPS

    if ~isempty(MPSinit_Opt)       % if initial MPS is not given as an input

        MPS_init = cell(ChainLen, 1);

        if isequal(MPSinit_Opt, 'IterDiag')     % initialize MPS by iterative diagonalization

            % initial Hamiltonian and isometry
            Aprev = 1;
            Oprev = 1;
    
            for itN = 1:ChainLen
                % Add local site to the Hamiltonian and define the identity
                % tensor connecting the previous Hilbert space and the new local space
                Anow = getIdentity(Aprev, 2, MPO{itN}, 1, [1,3,2]);
    
                Hnow = updateLeft(Oprev, 3, Anow, MPO{itN}(:,:,:,1), 4, Anow);
    
                % diagonalize Hnow
                [U,E0] = eig((Hnow+Hnow')/2);
    
                [E0,id] = sort(diag(E0), 'ascend');
                U = U(:,id);
    
                if numel(E0) < Nkeep
                    Ntr = numel(E0);
                else
                    keepMax = min(round(Nkeep*(1+margin)), numel(E0));
                    Ediff = E0(Nkeep+1 : keepMax) - E0(Nkeep : keepMax-1);
                    [~, id] = max(Ediff);
                    Ntr = Nkeep + id - 1;
                end
                % truncate and define isometry
                U = U(:, 1:Ntr);
                E0 = E0(1:Ntr);
    
                % update isometry and Hamiltonian
                if itN  == ChainLen
                    MPS_init{itN} = contract(Anow, 3, 2, U(:,1), 2, 1, [1,3,2]);
                else
                    Aprev = contract(Anow, 3, 2, U, 2, 1, [1,3,2]);
                    MPS_init{itN} = Aprev;
                    Oprev = updateLeft(Oprev, 3, MPS_init{itN}, MPO{itN}, 4, MPS_init{itN});
                end
    
            end
    
        else    % initialize MPS by a random MPS

            MPS_init = cell(ChainLen,1);

            for itL = 1:ChainLen

                if itL == 1 && isequal(Boundary, 'open')
                    MPS_init{itL} = randn(1, Nkeep, size(MPO{itL}, 1));
                elseif itL == ChainLen && isequal(Boundary, 'open')
                    MPS_init{itL} = randn(Nkeep, 1, size(MPO{itL}, 1));
                else
                    MPS_init{itL} = randn(Nkeep, Nkeep, size(MPO{itL}, 1));
                end
            end % itL

        end
        
    end

    % Convert the initial MPS to left-canonical form
    [MPS, ~, ~] = getCanonForm(MPS_init, ChainLen, 'Nkeep', Nkeep);  


    if isequal(updateOpt, '1site')
        %% Perform 1-site DMRG for ground state search

        H_LR = cell(ChainLen+2, 1);     % left and right parts of the effective Hamiltonian
        H_LR{1} = 1;                    % left 'dummy' space
        H_LR{end} = 1;                  % right 'dummy' space
        E_iter = zeros(ChainLen, 2*Nsweep);     % array to store energy of the MPS over iterations
        Sv = cell(ChainLen+1, 1);               % cell array to store singular values on MPS bonds
    
        for itN = 1:ChainLen-1      % compute left part of the Hamiltonian from the initial left-canonical MPS
            H_LR{itN+1} = updateLeft(H_LR{itN}, 3, MPS{itN}, MPO{itN}, 4, MPS{itN});
        end

        for itS = 1:Nsweep
            
            % right to left sweep
            for itN = ChainLen:-1:1
    
                % Update the current site using the Lanczos method
                [MPS{itN}, E_iter(ChainLen+1-itN, 2*itS-1)] = update_1site(H_LR{itN}, MPO{itN}, H_LR{itN+2}, MPS{itN}, KrylovDim, KrylovTol);
    
                % SVD the updated tensor to move one site left
                [Ul, Sv{itN}, MPS{itN}] = svdTr(MPS{itN}, 3, 1, Nkeep, 0);
    
                if itN > 1
                    % Move one site left
                    MPS{itN-1} = contract(MPS{itN-1}, 3, 2, Ul*diag(Sv{itN}), 2, 1, [1,3,2]);
                else
                    % Sv{itN} must be 1, since it is equal to the MPS norm
                    % absorb Ul into the leftmost tensor
                    MPS{itN} = contract(Ul, 2, 2, MPS{itN}, 3, 1);
                end
    
                % update the right part of the effective Hamiltonian
                H_LR{itN+1} = updateLeft(H_LR{itN+2}, 3, permute(MPS{itN}, [2,1,3]), ...
                                            permute(MPO{itN}, [1,2,4,3]), 4, permute(MPS{itN}, [2,1,3]));
    
            end % itN
    
            % display information of the right to left sweep
            if display
                disptime(['Sweep #', sprintf('%d/%d',[2*itS-1, 2*Nsweep]), ' (right -> left) : Energy = ', sprintf('%.7g',E_iter(ChainLen, 2*itS-1)), ', Norm = ', sprintf('%.5g',Sv{1})]);
            end

            % left to right sweep
            for itN = 1:ChainLen
    
                % Update the current site using the Lanczos method
                [MPS{itN}, E_iter(itN, 2*itS)] = update_1site(H_LR{itN}, MPO{itN}, H_LR{itN+2}, MPS{itN}, KrylovDim, KrylovTol);
    
                % SVD the updated tensor to move one site right
                [MPS{itN}, Sv{itN+1}, Vd] = svdTr(MPS{itN}, 3, [1,3], Nkeep, 0);
                MPS{itN} = permute(MPS{itN}, [1,3,2]);
    
                if itN < ChainLen
                    % Move one site right
                    MPS{itN+1} = contract(diag(Sv{itN+1})*Vd, 2, 2, MPS{itN+1}, 3, 1);
                else
                    % Sv{itN+1} must be 1, since it is equal to the MPS norm
                    % absorb Vd into the rightmost tensor
                    MPS{itN} = contract(MPS{itN}, 3, 2, Vd, 2, 1, [1,3,2]);
                end
                
                % update the left part of the effective Hamiltonian
                H_LR{itN+1} = updateLeft(H_LR{itN}, 3, MPS{itN}, MPO{itN}, 4, MPS{itN});
    
            end % itN
    
            % display information of the left to right sweep
            if display
                disptime(['Sweep #', sprintf('%d/%d',[2*itS, 2*Nsweep]), ' (left -> right) : Energy = ', sprintf('%.9g',E_iter(ChainLen, 2*itS)), ', Norm = ', sprintf('%.5g', Sv{ChainLen+1})]);
            end
            
        end % itS
    
        % The energy of the final ground state MPS
        E0 = E_iter(ChainLen, 2*itS);


    else
        %% Perform 2-site DMRG for ground state search

        H_LR = cell(ChainLen+2, 1);     % left and right parts of the effective Hamiltonian
        H_LR{1} = 1;                    % left 'dummy' space
        H_LR{end} = 1;                  % right 'dummy' space
        E_iter = zeros(ChainLen-1, 2*Nsweep);       % array to store energy of the MPS over iterations
        Sv = cell(ChainLen+1, 1);                   % cell array to store singular values on MPS bonds

        for itN = 1:ChainLen-2      % compute left part of the Hamiltonian from the initial left-canonical MPS
            H_LR{itN+1} = updateLeft(H_LR{itN}, 3, MPS{itN}, MPO{itN}, 4, MPS{itN});
        end

        for itS = 1:Nsweep
            
            % right to left sweep
            for itN = ChainLen-1:-1:1
    
                % Update the current site using the Lanczos method
                C_old = contract(MPS{itN}, 3, 2, MPS{itN+1}, 3, 1, [1,3,2,4]);
                [C_new, E_iter(ChainLen-itN, 2*itS-1)] = update_2site(H_LR{itN}, MPO{itN}, MPO{itN+1}, H_LR{itN+3}, C_old, KrylovDim, KrylovTol);
    
                % SVD the updated tensor to move one site left
                [Ul, Sv{itN+1}, MPS{itN+1}] = svdTr(C_new, 4, [1,3], Nkeep, 0);
    
                % update the orthogonality center
                MPS{itN} = contract(Ul, 3, 3, diag(Sv{itN+1}), 2, 1, [1,3,2]);
                
                % update the right part of the effective Hamiltonian
                H_LR{itN+2} = updateLeft(H_LR{itN+3}, 3, permute(MPS{itN+1}, [2,1,3]), ...
                                            permute(MPO{itN+1}, [1,2,4,3]), 4, permute(MPS{itN+1}, [2,1,3]));
    
            end % itN
    
            % display information of the right to left sweep
            if display
                disptime(['Sweep #', sprintf('%d/%d',[2*itS-1, 2*Nsweep]), ' (right -> left) : Energy = ', sprintf('%.7g',E_iter(ChainLen-1, 2*itS-1))]);
            end
    
            % left to right sweep
            for itN = 1:ChainLen-1
    
                % Update the current site using the Lanczos method
                C_old = contract(MPS{itN}, 3, 2, MPS{itN+1}, 3, 1, [1,3,2,4]);
                [C_new, E_iter(itN, 2*itS)] = update_2site(H_LR{itN}, MPO{itN}, MPO{itN+1}, H_LR{itN+3}, C_old, KrylovDim, KrylovTol);
    
                % SVD the updated tensor to move one site right
                [MPS{itN}, Sv{itN+1}, Vd] = svdTr(C_new, 4, [1,3], Nkeep, 0);
                MPS{itN} = permute(MPS{itN}, [1,3,2]);
    
                % update the orthogonality center
                MPS{itN+1} = contract(diag(Sv{itN+1}), 2, 2, Vd, 3, 1);
                
                % update the left part of the effective Hamiltonian
                H_LR{itN+1} = updateLeft(H_LR{itN}, 3, MPS{itN}, MPO{itN}, 4, MPS{itN});
    
            end % itN
    
            % display information of the left to right sweep
            if display
                disptime(['Sweep #', sprintf('%d/%d',[2*itS, 2*Nsweep]), ' (left -> right) : Energy = ', sprintf('%.9g',E_iter(ChainLen-1, 2*itS))]);
            end
            
        end % itS

        % The energy of the final ground state MPS
        E0 = E_iter(ChainLen-1, 2*itS);
        
    end

    toc2(tobj,'-v');


    %% Define 1-site update function

    function [C_new, E_new] = update_1site(Hleft, Hcen, Hright, C_old, KrylovMax, tol)
        % <Description>
        % Update the orthogonality centor by diagonalizing the effective Hamiltonian using the Lanczos method
        %
        % <Input>
        % Hleft : [tensor] left part of the effective Hamiltonian, rank 3 tensor
        % Hcen : [tensor] center part of the effective Hamiltonian, rank 4 tensor
        %                   (i.e., local MPO)
        % Hright : [tensor] right part of the effective Hamiltonian, rank 3 tensor
        % C_old : [tensor] the local MPS tensor to be updated, rank 3 tensor
        % KrylovDim : [tensor] the dimension of the Krylov subspace to be used in the Lanczos method
        % tol : [numeric] the tolerance for beta
        %
        % <Output>
        % C_new : [tensor] the updated local MPS tensor
        % E_new : [numeric] the eigenvalue of the updated tensor

        V = zeros(size(C_old,1), size(C_old,2), size(C_old,3), KrylovMax);
        V(:,:,:,1) = C_old;
    
        alpha = zeros(KrylovMax,1);
        beta = zeros(KrylovMax-1,1);
        KrylovNum = 0;  % counter for number of Krylov vectors
    
        for itf = 1:KrylovMax

            % Calculate norm & normalize
            normV = sqrt(abs(contract(conj(V(:,:,:,itf)), 3, [1,2,3], V(:,:,:,itf), 3, [1,2,3])));

            if itf > 1
                beta(itf-1) = normV;
            end

            if normV < tol
                break;
            end
            V(:,:,:,itf) = V(:,:,:,itf) / normV;
    
            % matrix-vector multiplication
            MatVec = contract(Hleft, 3, 2, V(:,:,:,itf), 3, 1);
            MatVec = contract(MatVec, 4, [2,4], Hcen, 4, [3,2]);
            MatVec = contract(MatVec, 4, [2,4], Hright, 3, [2,3]);
            MatVec = permute(MatVec, [1,3,2]);

            % compute alpha
            alpha(itf) = real(contract(conj(V(:,:,:,itf)), 3, [1,2,3], MatVec, 3, [1,2,3]));

            % orthogonalize
            for itO = 1:2   % orthogonalize twice to reduce numerical noise
                Proj = contract(conj(V(:,:,:,1:itf)), 4, [1,2,3], MatVec, 3, [1,2,3]);
                Proj = contract(V(:,:,:,1:itf), 4, 4, Proj, 2, 1);
                MatVec = MatVec - Proj;
            end

            V(:,:,:,itf+1) = MatVec;

            KrylovNum = KrylovNum + 1;

        end % itf
        V = V(:,:,:,1:KrylovNum);
        alpha = alpha(1:KrylovNum);
        beta = beta(1:KrylovNum-1);

        % construct and diagonalize the tridiagonal matrix
        Tri = diag(alpha) + diag(beta, 1) + diag(beta, -1);
        [C, E] = eig(Tri);
        [~, idx] = min(diag(E));

        GS_Krylov = C(:,idx);
        C_new = contract(V, 4, 4, GS_Krylov, 2, 1);

        % compute the expectation value of the effective Hamiltonian with respect to C_new
        tmp = updateLeft(Hleft, 3, C_new, Hcen, 4, C_new);
        E_new = real(contract(tmp, 3, [1,2,3], Hright, 3, [1,2,3]));
    end


    %% Define 2-site update function

    function [C_new, E_new] = update_2site(Hleft, Hcen1, Hcen2, Hright, C_old, KrylovMax, tol)
        % <Description>
        % Update the orthogonality centor by diagonalizing the effective Hamiltonian using the Lanczos method
        %
        % <Input>
        % Hleft : [tensor] left part of the effective Hamiltonian, rank 3 tensor
        % Hcen1 : [tensor] left of the two center parts of the effective Hamiltonian, rank 4 tensor
        %                   (i.e., local MPO)
        % Hcen2 : [tensor] right of the two center parts of the effective Hamiltonian, rank 4 tensor
        %                   (i.e., local MPO)
        % Hright : [tensor] right part of the effective Hamiltonian, rank 3 tensor
        % C_old : [tensor] the contraction of the two local MPS tensors to be updated, rank 4 tensor
        % KrylovDim : [tensor] the dimension of the Krylov subspace to be used in the Lanczos method
        % tol : [numeric] the tolerance for beta
        %
        % <Output>
        % C_new : [tensor] the updated local MPS tensor
        % E_new : [numeric] the eigenvalue of the updated tensor

        V = zeros(size(C_old,1), size(C_old,2), size(C_old,3), size(C_old,4), KrylovMax);
        V(:,:,:,:,1) = C_old;
    
        alpha = zeros(KrylovMax,1);
        beta = zeros(KrylovMax-1,1);
        KrylovNum = 0;  % counter for number of Krylov vectors
    
        for itf = 1:KrylovMax

            % Calculate norm & normalize
            normV = sqrt(abs(contract(conj(V(:,:,:,:,itf)), 4, [1,2,3,4], V(:,:,:,:,itf), 4, [1,2,3,4])));

            if itf > 1
                beta(itf-1) = normV;
            end

            if normV < tol
                break;
            end
            V(:,:,:,:,itf) = V(:,:,:,:,itf) / normV;
    
            % matrix-vector multiplication
            MatVec = contract(Hleft, 3, 2, V(:,:,:,:,itf), 4, 1);
            MatVec = contract(MatVec, 5, [2,4], Hcen1, 4, [3,2]);
            MatVec = contract(MatVec, 5, [3,5], Hcen2, 4, [2,3]);
            MatVec = contract(MatVec, 5, [2,5], Hright, 3, [2,3]);
            MatVec = permute(MatVec, [1,4,2,3]);

            % compute alpha
            alpha(itf) = real(contract(conj(V(:,:,:,:,itf)), 4, [1,2,3,4], MatVec, 4, [1,2,3,4]));

            % orthogonalize
            for itO = 1:2   % orthogonalize twice to reduce numerical noise
                Proj = contract(conj(V(:,:,:,:,1:itf)), 5, [1,2,3,4], MatVec, 4, [1,2,3,4]);
                Proj = contract(V(:,:,:,:,1:itf), 5, 5, Proj, 2, 1);
                MatVec = MatVec - Proj;
            end

            V(:,:,:,:,itf+1) = MatVec;

            KrylovNum = KrylovNum + 1;

        end % itf
        V = V(:,:,:,:,1:KrylovNum);
        alpha = alpha(1:KrylovNum);
        beta = beta(1:KrylovNum-1);

        % construct and diagonalize the tridiagonal matrix
        Tri = diag(alpha) + diag(beta, 1) + diag(beta, -1);
        [C, E] = eig(Tri);
        [~, idx] = min(diag(E));

        GS_Krylov = C(:,idx);
        C_new = contract(V, 5, 5, GS_Krylov, 2, 1);

        % compute the expectation value of the effective Hamiltonian with respect to C_new
        tmp = contract(Hleft, 3, 2, C_new, 4, 1);
        tmp = contract(tmp, 5, [2,4], Hcen1, 4, [3,2]);
        tmp = contract(tmp, 5, [3,5], Hcen2, 4, [2,3]);
        tmp = contract(tmp, 5, [1,3,4], conj(C_new), 4, [1,3,4]);
        E_new = real(contract(tmp, 3, [1,2,3], Hright, 3, [2,3,1]));
    end
    
end