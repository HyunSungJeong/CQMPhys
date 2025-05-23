function [MPO_canon, Sv, dw] = getCanonForm_MPO(MPO, idx, varargin)
    % <Description>
    % converts the input MPO into a canonical form
    % 
    % <Input>
    % MPO : [cell vector] cell vector containing the tensors forming the MPO to be converted to canonical form
    %
    % idx : [numeric] Index specifying the center of the site-canonical form of the MPO
    % 
    % <Option>
    % 'Nkeep', ... : [numeric] Maximum number of singluar values to keep at each SVD
    %                          If empty([]), Nkeep is regarded as infinite
    %                               (Default: Inf. i.e., no truncation)
    %
    % 'Skeep', ... : [numeric] Minimum absolute value of the singular value to keep
    %                          If empty([]), Skeep is set as its default value
    %                               (Default: 10*eps(Sv(1)), where Sv(1) is the maximum singular value
    %                                         i.e., same as the default option of Skeep in the function 'svdTr')
    %
    % <Output>
    % MPO_canon : [cell vector] cell vector containing the left- and right- normalized tensors forming the MPO in the canonical form
    %
    % Sv : [numeric vector] Vector containing the diagonal elements(singular values) of the diamond tensor.
    %                       Sorted in descending order
    %
    % dw : [numeric vector] dw(it) : discarded weights in the SVD between MPO{n} and MPO{n+1}

    %% Parse inputs

    if ~iscell(MPO)
        error('ERR: ''MPO'' must be a cell array');
    end

    if ~isnumeric(idx)
        error('ERR: ''idx'' must be a nonnegative integer');
    elseif mod(idx, 1) ~= 0
        error('ERR: ''idx'' must be a nonnegative integer');
    elseif idx < 0 || idx > numel(MPO)
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

            otherwise
                if ischar(varargin{1})
                    error(['ERR: unknown input ''',varargin{1},'''']);
                else
                    error('ERR: unknown input');
                end
        end % switch-case
    end % while

    %% Convert the MPO to canonical form
    
    MPO_canon = cell(numel(MPO), 1);
    dw = zeros(numel(MPO)-1, 1);

    % Convert the left part of the MPS to the left-canonical form
    for it = 1:idx-1
        [MPO_canon{it}, S, Vd, dw(it)] = svdTr(MPO{it}, 4, [1,2,3], Nkeep, Skeep);

        MPO{it+1} = contract(diag(S)*Vd, 2, 2, MPO{it+1}, 4, 3, [2,3,1,4]);
    end

    % Convert the right part of the MPS to the right-canonical form
    for it = numel(MPO):-1:idx+2
        [U, S, MPO_canon{it}, dw(it)] = svdTr(MPO{it}, 4, 3, Nkeep, Skeep);
        MPO_canon{it} = permute(MPO_canon{it}, [2,3,1,4]);

        MPO{it-1} = contract(MPO{it-1}, 4, 4, U*diag(S), 2, 1);
    end

    % Treat the tensors near the orthogonality center
    if idx < numel(MPO)
        tmp = contract(MPO{idx}, 4, 4, MPO{idx+1}, 4, 3);   % rank-4 tensor to be svd'ed
        [MPO_canon{idx}, Sv, MPO_canon{idx+1}] = svdTr(tmp, 6, [1,2,3], Nkeep, Skeep);

        MPO_canon{idx} = contract(MPO_canon{idx}, 4, 4, diag(Sv), 2, 2);
        MPO_canon{idx+1} = permute(MPO_canon{idx+1}, [2,3,1,4]);
    end


end