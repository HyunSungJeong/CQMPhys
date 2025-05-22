function [MPS_canon, Sv, dw] = getCanonForm(MPS, idx, varargin)
    % <Description>
    % converts the input MPS into a canonical form
    % 
    % <Input>
    % MPS : [cell vector] cell vector containing the tensors forming the MPS to be converted to canonical form
    % idx : [numeric] Index specifying the position of the diamond tensor in the bond canonical form.
    %               The diamond tensor containing the singular values is located between site # idx and site # idx+1
    %           1. idx == 0
    %               Convert to right canonical form
    %           2. idx == numel(MPS)
    %               Convert to left canonical form
    %           3. otherwise
    %               Convert to bond canonical form
    % 
    % <Option>
    % 'Nkeep', ... : [numeric] Maximum number of singluar values to keep at each SVD
    %                   (Default: Inf. i.e., no truncation)
    % 'Skeep', ... : [numeric] Minimum absolute value of the singular value to keep
    %                   (Default: 10*eps(Sv(1)), where Sv(1) is the maximum singular value
    %                               i.e., same as the default option of Skeep in the function 'svdTr')
    % 'Normalize', ... : [logical] Whether to normalize or not normalize the given MPS
    %                           This option cannot be turned off if idx == 0 or idx == numel(MPS) (i.e., right- or left- canonical)
    %                           (Default: true (normalize))
    %
    % <Output>
    % MPS_canon : [cell vector] cell vector containing the left- and right- normalized tensors forming the MPS in the canonical form
    % Sv : [numeric vector] Vector containing the diagonal elements(singular values) of the diamond tensor.
    %                       Sorted in descending order
    % dw : [numeric vector] dw(it) : discarded weights in the SVD between MPS{n} and MPS{n+1}

    %% Parse inputs

    if ~iscell(MPS)
        error('ERR: ''MPS'' must be a cell array');
    end

    if ~isnumeric(idx)
        error('ERR: ''idx'' must be a nonnegative integer');
    elseif mod(idx, 1) ~= 0
        error('ERR: ''idx'' must be a nonnegative integer');
    elseif idx < 0 || idx > numel(MPS)
        error('ERR: ''idx'' must be an integer between 0 and chain length');
    end

    %% Parse options
    
    Nkeep = [];
    Skeep = [];
    Normalize = true;

    while ~isempty(varargin)
        switch varargin{1}
            case 'Nkeep'
                if ~isnumeric(varargin{2})
                    error('ERR: ''Nkeep'' must be an integer');
                else
                    Nkeep = varargin{2};
                    varargin(1:2) = [];
                end

            case 'Skeep'
                if ~isnumeric(varargin{2})
                    error('ERR: ''Skeep'' must be a number');
                elseif varargin{2} <= 0
                    error('ERR: ''Skeep'' must be positive');
                else
                    Skeep = varargin{2};
                    varargin(1:2) = [];
                end

            case 'Normalize'
                if islogical(varargin{2})
                    Normalize = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''Normalize'' must be either true or false');
                end

            otherwise
                if ischar(varargin{1})
                    error(['ERR: unknown input ''',varargin{1},'''']);
                else
                    error('ERR: unknown input');
                end
        end % switch-case
    end % while

    %% Convert the MPS to canonical form
    
    MPS_canon = cell(numel(MPS), 1);
    dw = zeros(numel(MPS)-1, 1);

    % Convert the left part of the MPS to the left-canonical form
    for it = 1:idx-1
        [MPS_canon{it}, S, Vd, dw(it)] = svdTr(MPS{it}, 3, [1,3], Nkeep, Skeep);
        MPS_canon{it} = permute(MPS_canon{it}, [1,3,2]);

        MPS{it+1} = contract(diag(S)*Vd, 2, 2, MPS{it+1}, 3, 1);
    end

    % Convert the right part of the MPS to the right-canonical form
    for it = numel(MPS):-1:idx+2
        [U, S, MPS_canon{it}, dw(it)] = svdTr(MPS{it}, 3, 1, Nkeep, Skeep);

        MPS{it-1} = contract(MPS{it-1}, 3, 2, U*diag(S), 2, 1, [1,3,2]);
    end

    % Treat the tensors near the orthogonality center
    if idx == 0     % right-canonical form
        [U, Sv, MPS_canon{1}, dw(1)] = svdTr(MPS{1}, 3, 1, Nkeep, Skeep);
        MPS_canon{1} = contract(U, 2, 2, MPS_canon{1}, 3, 1);   % absorb the overall phase factor to MPS{1}

    elseif idx == numel(MPS)    % left-canonical form
        [MPS_canon{idx}, Sv, Vd, dw(idx)] = svdTr(MPS{idx}, 3, [1,3], Nkeep, Skeep);
        MPS_canon{idx} = contract(MPS_canon{idx}, 3, 3, Vd, 2, 1, [1,3,2]);

    else    % bond-canonical form
        tmp = contract(MPS{idx}, 3, 2, MPS{idx+1}, 3, 1);   % rank-4 tensor to be svd'ed
        [MPS_canon{idx}, Sv, MPS_canon{idx+1}] = svdTr(tmp, 4, [1,2], Nkeep, Skeep);
        if Normalize
            Sv = Sv / norm(Sv);
        end
        MPS_canon{idx} = permute(MPS_canon{idx}, [1,3,2]);
    end


end