function MPO_prod = varMul_MPO(MPO1, MPO2, Nkeep, Nsweep, varargin)
    % <Description>
    % Variationally multiplies two MPOs of equal sizes. MPO_prod =  MPO2*MPO1
    %
    % <Input>
    % MPO1 : [cell vector] Tensors constituting the 'top' MPO in the product
    %
    % MPO2 : [cell vector] Tensors constituting the 'bottom' MPO in the product
    % 
    % Nkeep : [numeric] The bond dimension of the MPO to be kept in the variational multiplication
    %
    % Nsweep : [numeric] The number of sweeps(left --> right & right --> left pairs) to perform in the variational multiplication
    % 
    % <Option>
    %
    % '-v' : If used, the information of sweep is displayed
    %           (Default: not used)
    %
    % <Output>
    % MPO_prod : [cell vector] Tensors constituting the variational product of two MPOs

    %% Parse inputs

    if ~iscell(MPO1)
        error('ERR: ''MPO1'' must be a cell vector of tensors forming the MPO Hamiltonian');
    end

    if ~iscell(MPO2)
        error('ERR: ''MPO2'' must be a cell vector of tensors forming the MPO Hamiltonian');
    end

    if numel(MPO1) ~= numel(MPO2)
        error('ERR: ''MPO1'' and ''MPO2'' must have same lengths');
    else
        ChainLen = numel(MPO1);
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

    % default options
    display = false;

    while ~isempty(varargin)
        switch varargin{1}
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

    %% Perform variational product of two MPOs

    % initialize MPO_prod as MPO1 in right-canonical form
    MPO_prod = getCanonForm_MPO(MPO1, 1);

    Op_LR = cell(ChainLen+2, 1);     % left and right parts of the MPO_prod * MPO2 * MPO1 ('triple') contraction
    Op_LR{1} = 1;                    % left 'dummy' space
    Op_LR{end} = 1;                  % right 'dummy' space

    % calculate the right part of the 'triple' contraction
    for itN = ChainLen:-1:3
        % 'zipper' contraction
        tmp = contract(Op_LR{itN+2}, 3, 2, MPO1{itN}, 4, 4);
        tmp = contract(tmp, 5, [2,3], MPO2{itN}, 4, [4,2]);
        Op_LR{itN+1} = contract(tmp, 5, [1,2,4], MPO_prod{itN}, 4, [4,1,2], [3,1,2]);
    end

    for itS = 1:Nsweep

        %% left --> right sweep

        % display information of the left --> right sweep
        if display
            disptime(['Performing sweep #', sprintf('%d/%d',[2*itS-1, 2*Nsweep]), ' (left -> right)']);
        end

        for itN = 1:ChainLen-1
            
            % matrix-vector multiplication
            tmp = contract(Op_LR{itN}, 3, 2, MPO1{itN}, 4, 3);
            tmp = contract(tmp, 5, [2,3], MPO2{itN}, 4, [3,2]);
            tmp = contract(tmp, 5, 3, MPO1{itN+1}, 4, 3);
            tmp = contract(tmp, 7, [4,5], MPO2{itN+1}, 4, [3,2]);
            C_new = contract(tmp, 7, [5,7], Op_LR{itN+3}, 3, [2,3]);

            % SVD and move one site right
            [MPO_prod{itN}, S, Vd] = svdTr(C_new, 6, [1,2,3], Nkeep, []);
            MPO_prod{itN} = permute(MPO_prod{itN}, [3,2,1,4]);
            MPO_prod{itN+1} = contract(diag(S), 2, 2, Vd, 4, 1, [3,2,1,4]);

            % update the left part of the 'triple' contraction
            tmp = contract(Op_LR{itN}, 3, 2, MPO1{itN}, 4, 3);
            tmp = contract(tmp, 5, [2,3], MPO2{itN}, 4, [3,2]);
            Op_LR{itN+1} = contract(tmp, 5, [1,2,4], MPO_prod{itN}, 4, [3,1,2], [3,1,2]);
        end

        %% right --> left sweep
        
        % display information of the right --> left sweep
        if display
            disptime(['Performing sweep #', sprintf('%d/%d',[2*itS, 2*Nsweep]), ' (right -> left)']);
        end

        for itN = ChainLen-1:-1:1

            % matrix-vector multiplication
            tmp = contract(Op_LR{itN}, 3, 2, MPO1{itN}, 4, 3);
            tmp = contract(tmp, 5, [2,3], MPO2{itN}, 4, [3,2]);
            tmp = contract(tmp, 5, 3, MPO1{itN+1}, 4, 3);
            tmp = contract(tmp, 7, [4,5], MPO2{itN+1}, 4, [3,2]);
            C_new = contract(tmp, 7, [5,7], Op_LR{itN+3}, 3, [2,3]);

            % SVD and move one site left
            [U, S, MPO_prod{itN+1}] = svdTr(C_new, 6, [1,2,3], Nkeep, []);
            MPO_prod{itN+1} = permute(MPO_prod{itN+1}, [3,2,1,4]);
            MPO_prod{itN} = contract(U, 4, 4, diag(S), 2, 1, [3,2,1,4]);

            % update the right part of the 'triple' contraction
            tmp = contract(Op_LR{itN+3}, 3, 2, MPO1{itN+1}, 4, 4);
            tmp = contract(tmp, 5, [2,3], MPO2{itN+1}, 4, [4,2]);
            Op_LR{itN+2} = contract(tmp, 5, [1,2,4], MPO_prod{itN+1}, 4, [4,1,2], [3,1,2]);
        end

    end

end