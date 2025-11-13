function [axes, varargout] = SLplot(X, Y, varargin)
    
    % <Description>
    %
    % Symmetric Log Plot
    %
    % <Input>
    % X : [numeric vector / cell vector of numeric vectors]
    %       if X is a numeric vector, 
    %           X is a vector of x-coordinates of data points
    %       if X is a cell vector of numeric vectors,
    %           each cell element is a vector of x-coordinates of data points
    %
    % Y : [numeric vector / cell vector of numeric vectors]
    %       if Y is a numeric vector, 
    %           Y is a vector of y-coordinates of data points
    %       if Y is a cell vector of numeric vectors,
    %           each cell element is a vector of y-coordinates of data points
    %
    % (The sizes of cell arrays X, Y and the lengths of vectors color, legend
    % must all be identical)
    %
    % <Option>
    % 'XScale',... : [char] x-axis scale. 'linear', 'log', or 'symlog' (Default: 'symlog')
    % 'YScale',... : [char] y-axis scale. 'linear', 'log', or 'symlog' (Default: 'linear')
    % 'XexpLim',... : [length 2 numeric vector] minimum and maximum natural logarithmic shown on x-axis in 'symlog' scale, respectively
    %                   (default: caculated from input data)
    % 'YexpLim',... : [length 2 numeric vector] minimum and maximum natural logarithmic shown on y-axis in 'symlog' scale, respectively
    %                   (default: caculated from input data)
    % 'XcontAtzero' : if used, the data is shown continuously near x = 0 (Default: not used)
    %                   (Only works when XScale is symlog and when the X data contains both positive and negative values)
    % 'YcontAtzero' : if used, the data is shown continuously near y = 0 (Default: not used)
    %                   (Only works when YScale is symlog and when the Y data contains both positive and negative values)
    % 'xzerowidth',... : [numeric] the ratio of distance between minimum absolute
    %                       value Xticks with respect to total x-axis length 
    %                       (Default: 0.08) (0<Xzerowidth<0.9)
    % 'yzerowidth',... : [numeric] the ratio of distance between minimum absolute
    %                       value Yticks with respect to total x-axis length
    %                       (Default: 0.08) (0<Yzerowidth<0.9)
    %
    % <Output>
    % 'axes' : x and y axis for symlog plot
    %  
    % if 'XScale' == 'symlog'
    % 'lin2sym_X' : function handle to convert linear X data to symmetric log X data
    %                   input X data must be within the interval specified by 'XexpLim'
    % 'sym2lin_X' : function handle to convert symmetric log X data to linear X data
    %
    % if 'YScale' == 'symlog'
    % 'lin2sym_Y' : function handle to convert linear Y data to symmetric log Y data
    %                   input Y data must be within the interval specified by 'YexpLim'
    % 'sym2lin_Y' : function handle to convert symmetric log Y data to linear Y data

    %% checking inputs
    if isequal(class(X),'double')
        X = {X};
        lenX = 1;
    elseif isequal(class(X),'cell')
        lenX = numel(X);
    else
        error('ERR: X must be either a numeric vector or a cell vector of numeric vectors');
    end

    if isequal(class(Y),'double')
        Y = {Y};
        lenY = 1;
    elseif isequal(class(Y),'cell')
        lenY = numel(Y);
    else
        error('ERR: Y must be either a numeric vector or a cell vector of numeric vectors');
    end

    if ~isequal(lenX, lenY)
        error('ERR: The number of data types, colors, and legends must be consistent');
    end


    %% Parsing options

    % Default values of options
    XScale = 'symlog';
    YScale = 'linear';
    Xzerowidth = 0.08;
    Yzerowidth = 0.08;
    XcontAtzero = false;
    YcontAtzero = false;
    XexpLim_spec = false;   % is XexpLim specified?
    YexpLim_spec = false;   % is YexpLim specified?

    while ~isempty(varargin)
        switch varargin{1}
            case 'XScale'
                if ismember(varargin{2}, {'linear', 'log', 'symlog'})
                    XScale = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''XScale'' must be ''linear'', ''log'', or ''symlog''');
                end

            case 'YScale'
                if ismember(varargin{2}, {'linear', 'log', 'symlog'})
                    YScale = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''XScale'' must be ''linear'', ''log'', or ''symlog''');
                end
            
            case 'XexpLim'
                
                if isnumeric(varargin{2})
                    XexpLim_spec = true;
                    XexpLim = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''XexpLim'' must be a numeric vector of length 2');
                end
                if ~isequal(XScale, 'symlog')
                    disp('Alert: ''XexpLim'' is not needed in current XScale');
                end

            case 'YexpLim'
                
                if isnumeric(varargin{2})
                    YexpLim_spec = true;
                    YexpLim = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''YexpLim'' must be a numeric vector of length 2');
                end
                if ~isequal(YScale, 'symlog')
                    disp('Alert: ''YexpLim'' is not needed in current YScale');
                end

            case 'xzerowidth'
                if isnumeric(varargin{2})
                    if varargin{2} >= 0 && varargin{2} <0.9
                        Xzerowidth = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: xzerowidth must be a between 0 and 0.9');
                    end
                else
                    error('ERR: xzerowidth must be a number');
                end

            case 'yzerowidth'
                if isnumeric(varargin{2})
                    if varargin{2} >= 0 && varargin{2} <0.9
                        Yzerowidth = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: yzerowidth must be a nonnegative real number smaller than 0.9');
                    end
                else
                    error('ERR: yzerowidth must be a number');
                end

            case 'XcontAtzero'
                XcontAtzero = true;
                varargin(1) = [];
                Xzerowidth = 0;

            case 'YcontAtzero'
                YcontAtzero = true;
                varargin(1) = [];
                Yzerowidth = 0;

            otherwise
                if ischar(varargin{1})
                    error(['ERR: Unknown option ',varargin{1}]);
                else
                    error('ERR: Unknown option');
                end
        end
    end

    if XcontAtzero && ~isequal(XScale,'symlog')
        error('''XcontAtzero'' can only be used when XScale is symlog')
    end

    if YcontAtzero && ~isequal(YScale,'symlog')
        error('''YcontAtzero'' can only be used when XScale is symlog')
    end

    %% Set X and Y ticks and convert data for plotting
    varargout = cell(0,0);
    axes = gca;
    axes.XAxis.TickLabelInterpreter = 'latex';
    axes.YAxis.TickLabelInterpreter = 'latex';

    if isequal(XScale,'symlog')

        sign_X = '';
        
        minExp_arr = zeros(numel(X), 1);
        maxExp_arr = zeros(numel(X), 1);

        for it = (1:numel(X))

            X_temp = X{it};
            
            % Update the sign for the dataset:
            % "positive" if the dataset X contains only positive values
            % "negative" if the dataset X contains only negative values
            % "undefined" if the dataset X contains both positive and negative values
            if all(X_temp>0)
                minExp_arr(it) = floor(log(min(X_temp))/log(10));   
                % largest integer whose power by 10 is smaller than the absolute value of any element of X_temp
                % i.e. all(abs(X_temp) > 10^minExp_arr(it))
                maxExp_arr(it) = ceil(log(max(X_temp))/log(10));
                % smallest integer whose power by 10 is larger than the absolute value of any element of X_temp
                % i.e. all(abs(X_temp) < 10^maxExp_arr(it))

                if ~isequal(sign_X,'undefined')
                    if isequal(sign_X,'negative')
                        sign_X = 'undefined';
                    else
                        sign_X = 'positive';
                    end
                end

            elseif all(X_temp<0)
                minExp_arr(it) = floor(log(-max(X_temp))/log(10));
                % largest integer whose power by 10 is smaller than the absolute value of any element of X_temp
                % i.e. all(abs(X_temp) > 10^minExp_arr(it))
                maxExp_arr(it) = ceil(log(-min(X_temp))/log(10));
                % smallest integer whose power by 10 is larger than the absolute value of any element of X_temp
                % i.e. all(abs(X_temp) < 10^maxExp_arr(it))
                
                if ~isequal(sign_X,'undefined')
                    if isequal(sign_X,'positive')
                        sign_X = 'undefined';
                    else
                        sign_X = 'negative';
                    end
                end

            elseif ~isequal(X_temp,0)
                pos_X = X_temp(X_temp>0);
                neg_X = X_temp(X_temp<0);
                minExp_arr(it) = min( floor(log(min(pos_X))/log(10)), floor(log(-max(neg_X))/log(10)) );
                % largest integer whose power by 10 is smaller than the absolute value of any element of X_temp
                % i.e. all(abs(X_temp) > 10^minExp_arr(it))
                maxExp_arr(it) = max( ceil(log(max(pos_X))/log(10)), ceil(log(-min(neg_X))/log(10)) );
                % smallest integer whose power by 10 is larger than the absolute value of any element of X_temp
                % i.e. all(abs(X_temp) < 10^maxExp_arr(it))
                sign_X = 'undefined';

            else
                sign_X = 'undefined';
            end
        end

        XminExp = min(minExp_arr);
        % largest integer whose power by 10 is smaller than the absolute value of any input X
        % i.e. all(abs(X_temp) > 10^minExp_arr(it))
        XmaxExp = max(maxExp_arr);
        % smallest integer whose power by 10 is larger than the absolute value of any input X
        % i.e. all(abs(X) < 10^maxExp_arr(it))

        if XexpLim_spec
            XminExp = XexpLim(1);
            XmaxExp = XexpLim(2);
        end

        if XmaxExp - XminExp < 5
            XminTick = XminExp;
            XmaxTick = XmaxExp;
            XtickInterval = 1;
        elseif ceil(XmaxExp/2) - ceil(XminExp/2) < 5
            XminTick = 2*ceil(XminExp/2);
            XmaxTick = 2*ceil(XmaxExp/2);
            XtickInterval = 2;
        else
            XminTick = 5*ceil(XminExp/5);
            XmaxTick = 5*floor(XmaxExp/5);
            XtickInterval = 5*ceil((XmaxTick - XminTick)/25);
        end

        % Function handles for linear <--> symlog transformation of x-coordinates
        lin2sym_X = @(x_lin) (heaviside(log(abs(x_lin)+1e-50)/log(10)-XminExp+1e-10).*sign(x_lin).*(Xzerowidth + (1-Xzerowidth-0.05)*(log(abs(x_lin)+1e-50)/log(10) - floor(XminExp))/(ceil(XmaxExp) - floor(XminExp))));

        sym2lin_X = @(x_sl) ((sign(x_sl).*power(10, (abs(x_sl)-Xzerowidth)*(ceil(XmaxExp) - floor(XminExp))/(1-Xzerowidth-0.05) + floor(XminExp)) ) );

        % generate xtick positions and their labels
        temp = (XminTick : XtickInterval : XmaxTick);
        XtickPositions = [];
        XminorTickPositions = [];

        for it = (1:numel(temp))
            XtickPositions = [XtickPositions, ...
                                Xzerowidth + (1-Xzerowidth-0.05)*(temp(it) - floor(XminExp))/(ceil(XmaxExp) - floor(XminExp))];
        end
        
        for it = (floor(XminExp):ceil(XmaxExp))
            if ~ismember(it, temp)
                XminorTickPositions = [XminorTickPositions, ...
                                        Xzerowidth + (1-Xzerowidth-0.05)*(it - floor(XminExp))/(ceil(XmaxExp) - floor(XminExp))];
            end
        end
        
        XtickLabels = cell(0,0);

        sign_X = 'undefined';   % To be fixed!!!

        switch sign_X
            case 'undefined'    % sign of dataset X is "undefined"
                if XcontAtzero
                    XtickPositions = [-flip(XtickPositions(2:end)), XtickPositions];
                    XminorTickPositions = [-flip(XminorTickPositions), XminorTickPositions];
                    
                    temp = flip(temp);
                    for it = (1:numel(temp)-1)
                        XtickLabels = cat(2,XtickLabels,'$-10^{'+string(temp(it))+'}$');
                    end
                    XtickLabels = cat(2,XtickLabels,'$\pm10^{'+string(temp(end))+'}$');
                    temp = flip(temp);
                    for it = (2:numel(temp))
                        XtickLabels = cat(2,XtickLabels,'$10^{'+string(temp(it))+'}$');
                    end
                else
                    XtickPositions = [-flip(XtickPositions), 0, XtickPositions];
                    XminorTickPositions = [-flip(XminorTickPositions), XminorTickPositions];

                    temp = flip(temp);
                    for it = (1:numel(temp))
                        XtickLabels = cat(2,XtickLabels,'$-10^{'+string(temp(it))+'}$');
                    end
                    XtickLabels = cat(2,XtickLabels,'0');
                    temp = flip(temp);
                    for it = (1:numel(temp))
                        XtickLabels = cat(2,XtickLabels,'$10^{'+string(temp(it))+'}$');
                    end
                end

                

            case 'positive'     % sign of dataset X is "positive"

                if XcontAtzero
                    error('''XcontAtzero'' can only be used when the X data contains both positive and negative values');
                end
                XtickPositions = [0, XtickPositions];

                for it = (1:numel(temp))
                    XtickLabels = cat(2,XtickLabels,'$10^{'+string(temp(it))+'}$');
                end
                XtickLabels = cat(2,'0',XtickLabels);

            case 'negative'     % sign of dataset X is "negative"

                if XcontAtzero
                    error('''XcontAtzero'' can only be used when the X data contains both positive and negative values');
                end
                XtickPositions = [-flip(XtickPositions), 0];

                temp = flip(temp);
                for it = (1:numel(temp))
                    XtickLabels = cat(2,XtickLabels,'$-10^{'+string(temp(it))+'}$');
                end
                XtickLabels = cat(2,XtickLabels,'0');
        end
        
        axes.XAxis.Scale = 'linear';
        if isequal(sign_X,'negative')
            axes.XAxis.Limits = [-1,0.1];
        elseif isequal(sign_X,'positive')
            axes.XAxis.Limits = [-0.1,1];
        else
            axes.XAxis.Limits = [-1,1];
        end
        axes.XAxis.TickValues = XtickPositions;
        axes.XAxis.TickLabels = XtickLabels;
        axes.XAxis.MinorTickValues = XminorTickPositions;

        varargout = cat(2,varargout,{lin2sym_X, sym2lin_X});

    else
        axes.XAxis.Scale = XScale;
    end



    if isequal(YScale,'symlog')

        sign_Y = '';

        minExp_arr = zeros(numel(Y), 1);
        maxExp_arr = zeros(numel(Y), 1);

        for it = (1:numel(Y))

            Y_temp = Y{it};

            
            % Update the sign for the dataset Y:
            % "positive" if the dataset Y contains only positive values
            % "negative" if the dataset Y contains only negative values
            % "undefined" if the dataset Y contains both positive and negative values
            if all(Y_temp>0)

                minExp_arr(it) = floor(log(min(Y_temp))/log(10));
                % largest integer whose power by 10 is smaller than the absolute value of any element of Y_temp
                % i.e. all(abs(Y_temp) > 10^minExp_arr(it))
                maxExp_arr(it) = ceil(log(max(Y_temp))/log(10));
                % smallest integer whose power by 10 is larger than the absolute value of any element of Y_temp
                % i.e. all(abs(Y_temp) < 10^minExp_arr(it))

                if ~isequal(sign_Y,'undefined')
                    if isequal(sign_Y,'negative')
                        sign_Y = 'undefined';
                    else
                        sign_Y = 'positive';
                    end
                end


            elseif all(Y_temp<0)

                minExp_arr(it) = floor(log(-max(Y_temp))/log(10));
                % largest integer whose power by 10 is smaller than the absolute value of any element of Y_temp
                % i.e. all(abs(Y_temp) > 10^minExp_arr(it))
                maxExp_arr(it) = ceil(log(-min(Y_temp))/log(10));
                % smallest integer whose power by 10 is larger than the absolute value of any element of Y_temp
                % i.e. all(abs(Y_temp) < 10^minExp_arr(it))

                if ~isequal(sign_Y,'undefined')
                    if isequal(sign_Y,'positive')
                        sign_Y = 'undefined';
                    else
                        sign_Y = 'negative';
                    end
                end


            else
                pos_Y = Y_temp(Y_temp>0);
                neg_Y = Y_temp(Y_temp<0);
                minExp_arr(it) = min( floor(log(min(pos_Y))/log(10)), floor(log(-max(neg_Y))/log(10)) );
                % largest integer whose power by 10 is smaller than the absolute value of any element of Y_temp
                % i.e. all(abs(Y_temp) > 10^minExp_arr(it))
                maxExp_arr(it) = max( ceil(log(max(pos_Y))/log(10)), ceil(log(-min(neg_Y))/log(10)) );
                % smallest integer whose power by 10 is larger than the absolute value of any element of Y_temp
                % i.e. all(abs(Y_temp) < 10^minExp_arr(it))

                sign_Y = 'undefined';
            end
        end

        YminExp = min(minExp_arr);
        YmaxExp = max(maxExp_arr);

        if YexpLim_spec
            YminExp = YexpLim(1);
            YmaxExp = YexpLim(2);
        end

        if YmaxExp - YminExp < 5
            YminTick = YminExp;
            YmaxTick = YmaxExp;
            YtickInterval = 1;
        elseif ceil(YmaxExp/2) - ceil(YminExp/2) < 5
            YminTick = 2*ceil(YminExp/2);
            YmaxTick = 2*ceil(YmaxExp/2);
            YtickInterval = 2;
        else
            YminTick = 5*ceil(YminExp/5);
            YmaxTick = 5*floor(YmaxExp/5);
            YtickInterval = 5*ceil((YmaxTick-YminTick)/25);
        end

        % Function handles for linear <--> symlog transformation of y-coordinates
        lin2sym_Y = @(y_lin) (heaviside(log(abs(y_lin)+1e-50)/log(10)-YminExp+1e-10).*sign(y_lin).*(Yzerowidth + (1-Yzerowidth-0.05)*(log(abs(y_lin)+1e-50)/log(10) - floor(YminExp))/(ceil(YmaxExp) - floor(YminExp))));

        sym2lin_Y = @(y_sl) ((sign(y_sl).*power(10, (abs(y_sl)-Yzerowidth)*(ceil(YmaxExp) - floor(YminExp))/(1-Yzerowidth-0.05) + floor(YminExp)) ) );

        % generate ytick positions and their labels
        temp = (YminTick : YtickInterval : YmaxTick);
        YtickPositions = [];
        YminorTickPositions = [];
        for it = (1:numel(temp))
            YtickPositions = [YtickPositions, ...
                                Yzerowidth + (1-Yzerowidth-0.05)*(temp(it) - floor(YminExp))/(ceil(YmaxExp) - floor(YminExp))];
        end

        for it = (YminExp:YmaxExp)
            YminorTickPositions = [YminorTickPositions, ...
                                    Yzerowidth + (1-Yzerowidth-0.05)*(it - floor(YminExp))/(ceil(YmaxExp) - floor(YminExp))];
        end

        YtickLabels = cell(0,0);
        switch sign_Y
            case 'undefined'    % sign of dataset Y is "undefined"
                if YcontAtzero
                    YtickPositions = [-flip(YtickPositions(2:end)), YtickPositions];
                    YminorTickPositions = [-flip(YminorTickPositions), YminorTickPositions];

                    temp = flip(temp);
                    for it = (1:numel(temp)-1)
                        YtickLabels = cat(2,YtickLabels,'$-10^{'+string(temp(it))+'}$');
                    end
                    YtickLabels = cat(2,YtickLabels,'$\pm10^{'+string(temp(end))+'}$');
                    temp = flip(temp);
                    for it = (2:numel(temp))
                        YtickLabels = cat(2,YtickLabels,'$10^{'+string(temp(it))+'}$');
                    end
                else
                    YtickPositions = [-flip(YtickPositions), 0, YtickPositions];
                    YminorTickPositions = [-flip(YminorTickPositions), YminorTickPositions];

                    temp = flip(temp);
                    for it = (1:numel(temp))
                        YtickLabels = cat(2,YtickLabels,'$-10^{'+string(temp(it))+'}$');
                    end
                    YtickLabels = cat(2,YtickLabels,'0');
                    temp = flip(temp);
                    for it = (1:numel(temp))
                        YtickLabels = cat(2,YtickLabels,'$10^{'+string(temp(it))+'}$');
                    end
                end

            case 'positive'     % sign of dataset Y is "positive"

                if YcontAtzero
                    error('''YcontAtzero'' can only be used when the Y data contains both positive and negative values');
                end
                YtickPositions = [0, YtickPositions];

                for it = (1:numel(temp))
                    YtickLabels = cat(2,YtickLabels,'$10^{'+string(temp(it))+'}$');
                end
                YtickLabels = cat(2,'0',YtickLabels);

            case 'negative'     % sign of dataset Y is "negative"

                if YcontAtzero
                    error('''YcontAtzero'' can only be used when the Y data contains both positive and negative values');
                end
                YtickPositions = [-flip(YtickPositions), 0];

                temp = flip(temp);
                for it = (1:numel(temp))
                    YtickLabels = cat(2,YtickLabels,'$-10^{'+string(temp(it))+'}$');
                end
                YtickLabels = cat(2,YtickLabels,'0');
        end        
        
        [YtickPositions, ids] = sort(YtickPositions);
        YtickLabels = YtickLabels(ids);

        axes.YAxis.Scale = 'linear';
        if isequal(sign_Y,'negative')
            axes.YAxis.Limits = [-1,0.1];
        elseif isequal(sign_Y,'positive')
            axes.YAxis.Limits = [-0.1,1];
        else
            axes.YAxis.Limits = [-1,1];
        end
        axes.YAxis.TickValues = YtickPositions;
        axes.YAxis.TickLabels = YtickLabels;
        axes.YAxis.MinorTickValues = YminorTickPositions;

        varargout = cat(2,varargout,{lin2sym_Y, sym2lin_Y});

    else
        axes.YAxis.Scale = YScale;
    end
end
