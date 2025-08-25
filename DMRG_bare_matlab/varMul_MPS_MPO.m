function MPS_prod = varMul_MPS_MPO(MPS, MPO, Nkeep, Nsweep, varargin)
    % <Description>
    % variationally multiplies an MPO to an MPS to obtain an MPS
    %
    % <Input>
    % MPS : [cell] tensors constituting the initial MPS
    %
    % MPO : [cell] tensors consitituting the MPO to be multiplied on the MPS
    %
    % Nkeep : [numeric] the maximum bond dimension to be kept
    %
    % Nsweep : [numeric] the number of variational sweeps
    %
    % <Options>
    % '1site' : If used, 1-site update is used
    %               (Default: used)
    %
    % '2site' : If used, 2-site update is used
    %               (Default: not used)
    %
    % '-v' : if used, the information of variational multiplication sweeps are shown as output
    %               (Default: not used)
    %
    % <Output>
    % MPS_prod : [cell] tensors consituting the MPS obtained by the MPS-MPO product

    %% Parse inputs

    if ~iscell(MPS)
        error('ERR: ''MPS'' must be a cell vector of tensors forming the initial MPS');
    end

    if ~iscell(MPO)
        error('ERR: ''MPO'' must be a cell vector of tensors forming the MPO');
    end

    if numel(MPS) ~= numel(MPO)
        error('ERR: ''MPS'' and ''MPO'' must have same lengths');
    else
        ChainLen = numel(MPS);
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
    update = '1site';

    while ~isempty(varargin)
        switch varargin{1}
            case '-v'
                display = true;
                varargin(1) = [];

            case '1site'
                update = '1site';
                varargin(1) = [];

            case '2site'
                update = '2site';
                varargin(1) = [];

            otherwise 
                if ischar(varargin{1})
                    error(['ERR: Unknown input ''',varargin{1},'''']);
                else
                    error('ERR: Unknown input');
                end
        end % switch-case
    end % while

    %% Perform variational MPS-MPO multiplication

    % initialize MPS_prod as initial MPS in the right-canonical form
    MPS_prod = getCanonForm(MPS, 0);

    Op_LR = cell(ChainLen+2,1);
    Op_LR{1} = 1;
    Op_LR{end} = 1;

    % calculate the right part of the 'triple' contraction
    switch update
        case '1site'
            END1 = 2;
            END2 = ChainLen;
        case '2site'
            END1 = 3;
            END2 = CHainLen-1;
    end

    for itN = ChainLen:-1:END1
        Op_LR{itN+1} = updateLeft(Op_LR{itN+2}, 3, permute(MPS_prod{itN}, [2,1,3]), ...
                                    permute(MPO{itN}, [1,2,4,3]), 4, permute(MPS{itN}, [2,1,3]));
    end

    for itS = 1:Nsweep

        %% Left --> right sweep

        % display information of the left --> right sweep
        if display
            disptime(['Performing sweep #', sprintf('%d/%d',[2*itS-1, 2*Nsweep]), ' (left -> right)']);
        end

        for itN = 1:END2

            switch update

                case '1site'

                    % matrix-vector multiplication
                    tmp = contract(Op_LR{itN}, 3, 2, MPS{itN}, 3, 1);
                    tmp = contract(tmp, 4, [2,4], MPO{itN}, 4, [3,2]);
                    C_new = contract(tmp, 4, [2,4], Op_LR{itN+2}, 3, [2,3]);

                    if itN < ChainLen
                        % SVD and move one site right
                        [MPS_prod{itN}, S, Vd] = svdTr(C_new, 3, [1,2], Nkeep, []);
                        MPS_prod{itN} = permute(MPS_prod{itN}, [1,3,2]);
                        tmp = contract(diag(S), 2, 2, Vd, 2, 1);
                        MPS_prod{itN+1} = contract(tmp, 2, 2, MPS_prod{itN+1}, 3, 1);
                    else
                        MPS_prod{itN} = permute(C_new, [1,3,2]);
                    end

                    % update the left part of the 'triple' contraction
                    % Op_LR{2} is empty after itN == 1 iteration
                    Op_LR{itN+1} = updateLeft(Op_LR{itN}, 3, MPS_prod{itN}, MPO{itN}, 4, MPS{itN});

                case '2site'
            
                    % matrix-vector multiplication
                    tmp = contract(Op_LR{itN}, 3, 2, MPS{itN}, 3, 1);
                    tmp = contract(tmp, 4, [2,4], MPO{itN}, 4, [3,2]);
                    tmp = contract(tmp, 4, 2, MPS{itN+1}, 3, 1);
                    tmp = contract(tmp, 5, [3,5], MPO{itN+1}, 4, [3,2]);
                    C_new = contract(tmp, 5, [2,5], Op_LR{itN+3}, 3, [2,3]);
        
                    % SVD and move one site right
                    [MPS_prod{itN}, S, Vd] = svdTr(C_new, 4, [1,2], Nkeep, []);
                    MPS_prod{itN} = permute(MPS_prod{itN}, [1,3,2]);
                    MPS_prod{itN+1} = contract(diag(S), 2, 2, Vd, 3, 1, [1,3,2]);
        
                    % update the left part of the 'triple' contraction
                    Op_LR{itN+1} = updateLeft(Op_LR{itN}, 3, MPS_prod{itN+1}, MPO{itN+1}, 4, MPS{itN+1});

            end % switch-case
        end % itN

        %% Right --> left sweep

        % display information of the right --> left sweep
        if display
            disptime(['Performing sweep #', sprintf('%d/%d',[2*itS, 2*Nsweep]), ' (right -> left)']);
        end

        for itN = END2:-1:1

            switch update

                case '1site'

                    % matrix-vector multiplication
                    tmp = contract(Op_LR{itN}, 3, 2, MPS{itN}, 3, 1);
                    tmp = contract(tmp, 4, [2,4], MPO{itN}, 4, [3,2]);
                    C_new = contract(tmp, 4, [2,4], Op_LR{itN+2}, 3, [2,3]);

                    if itN > 1
                        % SVD and move one site right
                        [U, S, MPS_prod{itN}] = svdTr(C_new, 3, 1, Nkeep, []);
                        MPS_prod{itN} = permute(MPS_prod{itN}, [1,3,2]);
                        tmp = contract(U, 2, 2, diag(S), 2, 1);
                        MPS_prod{itN-1} = contract(MPS_prod{itN-1}, 3, 2, tmp, 2, 1, [1,3,2]);
                    else
                        MPS_prod{itN} = permute(C_new, [1,3,2]);
                    end

                    % update the right part of the 'triple' contraction
                    Op_LR{itN+1} = updateLeft(Op_LR{itN+2}, 3, permute(MPS_prod{itN},[2,1,3]), ...
                                                    permute(MPO{itN},[1,2,4,3]), 4, permute(MPS{itN},[2,1,3]));

                case '2site'
            
                    % matrix-vector multiplication
                    tmp = contract(Op_LR{itN}, 3, 2, MPS{itN}, 3, 1);
                    tmp = contract(tmp, 4, [2,4], MPO{itN}, 4, [3,2]);
                    tmp = contract(tmp, 4, 2, MPS{itN+1}, 3, 1);
                    tmp = contract(tmp, 5, [3,5], MPO{itN+1}, 4, [3,2]);
                    C_new = contract(tmp, 5, [2,5], Op_LR{itN+3}, 3, [2,3]);
        
                    % SVD and move one site left
                    [U, S, MPS_prod{itN+1}] = svdTr(C_new, 4, [1,2], Nkeep, []);
                    MPS_prod{itN+1} = permute(MPS_prod{itN+1}, [1,3,2]);
                    MPS_prod{itN} = contract(U, 3, 3, diag(S), 2, 1, [1,3,2]);
        
                    % update the right part of the 'triple' contraction
                    Op_LR{itN+2} = updateLeft(Op_LR{itN+3}, 3, permute(MPS_prod{itN+1},[2,1,3]), ...
                                                permute(MPO{itN+1},[1,2,4,3]), 4, permute(MPS{itN+1},[2,1,3]));

            end % switch-case
        end % itN
    end % itS

end