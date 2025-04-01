function plotColorgrid(Cmesh, colormap,varargin)
    % <Description>
    % divides plot into rectangular meshes specified by 'Cmesh'
    % and colors them according to the 'colormap'
    %
    % <Input>
    % Cmesh : [cell array] Cell array with two elements. 
    %                   First cell element is a vector defining the verical grid,
    %                   and the second cell element is a vector defining the horizontal grid
    %
    % colormap : [numeric array] 3-D array specifying the color of each patches in the mesh. 
    %                       The (n,m,:) components define the rgb color of each patch defined by 'Cmesh'.
    %                       If Cmesh{1} is nx1 and Cmesh{2} is mx1, 'colormap' must be (n-1)x(m-1)x3.
    %
    % <Option>
    % 'Opacity', .. : [numeric] numeric value in range [0,1]. 
    %                       'FaceAlpha' of patches will be set to 1-Opacity
    %                           (Default: 0)
    %
    % <Output>
    % A plot colored according to 'Cmesh' and 'colormap'

    %% Check inputs

    if ~iscell(Cmesh)
        error('ERR: ''Cmesh'' must be a cell array');
    elseif numel(Cmesh) ~= 2
        error('ERR: ''Cmesh'' must have two elements');
    elseif  ~ismatrix(Cmesh{1}) || ~ismatrix(Cmesh{2})
        error('ERR: elements of ''Cmesh'' must be numeric vectors');
    elseif ~issorted(Cmesh{1})
        error('ERR: ''Cmesh{1}'' must be in ascending order');
    elseif ~issorted(Cmesh{2})
        error('ERR: ''Cmesh{2}'' must be in ascending order');
    end

    if ~isnumeric(colormap)
        error('ERR: ''colormap'' must be a numeric array');
    elseif ndims(colormap) ~= 3
        error('ERR: ''colormap'' must be a 3-dimensional array');
    elseif ~isequal(size(colormap), [numel(Cmesh{1})-1, numel(Cmesh{2})-1, 3])
        error('ERR: ''colormap'' and ''Cmesh'' must have compatible sizes');
    elseif ~(all(colormap>=0, 'all') && all(colormap<=1, 'all'))
        error('ERR: ''colormap'' must specify RGB colors: its array elements must be between 0 and 1');
    end

    %% Parse options

    % default options
    Opacity = 0;

    % input options
    while ~isempty(varargin)
        switch varargin{1}
            case 'Opacity'
                if ~isnumeric(varargin{2})
                    error('ERR: ''Opacity'' must be a real number');
                elseif ~(varargin{2} >= 0 && varargin{2} <= 1)
                    error('ERR: ''Opacity'' must be in range [0,1]');
                else
                    Opacity = varargin{2};
                    varargin(1:2) = [];
                end
            case 'opacity'
                if ~isnumeric(varargin{2})
                    error('ERR: ''Opacity'' must be a real number');
                elseif ~(varargin{2} >= 0 && varargin{2} <= 1)
                    error('ERR: ''Opacity'' must be in range [0,1]');
                else
                    Opacity = varargin{2};
                    varargin(1:2) = [];
                end
            otherwise
                if ischar(varargin{1})
                    error(['ERR: unknown option ''', varargin{1}, '''']);
                else
                    error('ERR: unknown input');
                end
        end
    end

    %% plot

    for it1 = 1:size(colormap,1)
        for it2 = 1:size(colormap,2)
            patchX = [Cmesh{1}(it1), Cmesh{1}(it1+1), Cmesh{1}(it1+1), Cmesh{1}(it1)];
            patchY = [Cmesh{2}(it2), Cmesh{2}(it2), Cmesh{2}(it2+1), Cmesh{2}(it2+1)];

            patch(patchX, patchY, colormap(it1,it2,:), 'FaceAlpha', 1-Opacity, 'linestyle', 'none');
        end
    end

end