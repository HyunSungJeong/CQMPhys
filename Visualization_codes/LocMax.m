function [Maxima, MaxIdx, MaxPos] = LocMax(x,y,varargin)
    %
    % <Description>
    %
    % Find 'reasonable' local minima in noisy data
    % local minima search tuned by 
    % 'sparseness', 'minpeak', 'maxsharpness','gloMinRange'
    % 
    % The function finds "reasonable" local maxima according to the following steps :
    %
    % 1. The function 'sparsifies' the input data by extracting data points in
    %   interval of 'sparseness'(Default: 10).
    % 2. The function finds local maxima from this 'sparse' data.
    % 3. The function weeds out local maxima with peak height smaller than 'minpeak' 
    %   (Default: 0.1)
    % 4. The function weeds out local maxima that are too 'sharp'
    %   A local maximum is regarded as too 'sharp' if the difference
    %   between the peak and its adjacent (sparse) data point exceeds
    %   a certain proportion(defined by 'maxsharpness') of its peak height
    %   (Default: 0.1 (10%))
    % 5. The function takes the orginal data points in the neighborhood of
    %   each local maximum of the "sparse" data and finds their maximum.
    %   If a local maximum of the 'sparse' data is the N th data of the
    %   original data, its "neighborhood" is defined by data in index range
    %   [ N - gloMinRange , N + gloMinRange ] 
    %
    % <Input>
    % X : [numeric vector] x-coordinates of data
    % Y : [numeric vector] y-coordinates of data
    % 
    % <Output>
    % Maxima : [numeric vector] vector of local maximum values
    % MaxIdx : [numeric vector] indices corresponding to local maxima 
    %                           of y-coordinates
    % MaxPos : [numeric vector] x-coordinates corresponding to local maxima
    %                           of y-coordinates
    %
    % <Option>
    % 'sparseness',... : sparseness of the "sparsified" data (Default: 10)
    % 'minpeak',... : minimum peak height that are regarded as relevant (Default: 0.1)
    % 'maxsharpness',... : maximum sharpness of peaks that are regarded as relevant
    %                   read description for the definition of "sharpness"
    %                   (Default: 0.1)
    % 'gloMinRange',... : The half-width of the range of original data to
    %                   be searched when finding maxima in original data
    %                   (Default: 20)
    
    %% Parsing inputs
    if numel(x) ~= numel(y)
        error('ERR: the number of x and y data must be the same');
    end

    %% Parsing options

    % Default values of options
    sparseness = 10;
    minpeak = 0.1;
    maxsharpness = 0.1;
    gloMinRange = 20;

    while ~isempty(varargin)
        switch varargin{1}
            case 'sparseness'
                if isnumeric(varargin{2})
                    sparseness = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''sparseness'' must be a positive integer');
                end

            case 'minpeak'
                if isnumeric(varargin{2})
                    minpeak = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''minpeak'' must be a positive real number');
                end

            case 'maxsharpness'
                if isnumeric(varargin{2})
                    maxsharpness = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''maxsharpness'' must be a positive real number');
                end

            case 'gloMinRange'
                if isnumeric(varargin{2})
                    gloMinRange = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''gloMinRange'' must be a positive integer');
                end

        end
    end

    x_sparse = x(1:sparseness:end);
    y_sparse = y(1:sparseness:end);
    sparsemaxIdx = islocalmax(y_sparse);
    near_locmax = [];
    
    for it = (1:numel(sparsemaxIdx))
	    if sparsemaxIdx(it)
		    if abs(y_sparse(it)) < minpeak
			    sparsemaxIdx(it) = false;
		    elseif abs(y_sparse(it) - y_sparse(it+1)) > abs(y_sparse(it))*maxsharpness
			    sparsemaxIdx(it) = false;
		    else
			    near_locmax = cat(2, near_locmax, sparseness*(it-1)+1);
		    end
	    end
    end
    
    Maxima = zeros(numel(near_locmax), 1);
    MaxIdx = zeros(numel(near_locmax), 1);
    MaxPos = zeros(numel(near_locmax), 1);
    
    for it = (1:numel(near_locmax))
	    [tmp1, tmp2] = max(y( max(near_locmax(it)-gloMinRange,1) : near_locmax(it)+gloMinRange));
	    Maxima(it) = tmp1;
        MaxIdx(it) = near_locmax(it) - gloMinRange + tmp2 - 1;
	    MaxPos(it) = x(MaxIdx(it));
    end

end