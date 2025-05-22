function [I, Phase_range, Phase_name] = PhaseRange_TI(J0,K0,varargin)
    % <Description>
    % Returns the range and corresponding temperature range for data of given (J0, K0)
    %
    % <Input>
    % J0 : [numeric] spin-spin coupling for the isotropic 2-channel spin-orbital Kondo(2soK) model 
    %           in which the phase range will be evaluated
    % K0 : [numeric] orbital-orbital coupling for the isotropic 2-channel spin-orbital Kondo(2soK) model
    %           in which the phase range will be evaluated
    %
    % <Options>
    % 'flatThres', .. : [numeric] maximum value of second derivative to be considered as zero
    %                           (Default: 0.05)
    % 'platLenThres', .. : [numeric] minimum length of plateau(in log_{10} scale) to be considered as a well-defined phase
    %                           (Default: 1)
    % 'NFLboundary' : if used, the boundary between NFL phase and FL phase is defined as the lowerest temperature showing NFL power law.
    %                 if not used, the boundary between NFL phase and FL phase is defined as the highest temperature showingFL power law.
    %                           (Default: not used)
    % 'I0min', ... : [numeric] minimum absolute value of I0 to be used in the plot
    %                           (Default: 10^{-14})
    % '-v' : when used, impurity spin susceptibilities are plotted with phase ranges for the used parameter sets
    %       (Default: not used)
    %
    % <Output>
    % I : [numeric vector ] row vector containing the spin-orbital coupling
    %           constants I0 for given (J0, K0)
    %
    % Phase_range : [cell array of numeric arrays] temperature ranges of phases for each I.
    %                   1 x numel(I) cell array. Each cell element is an array with two columns.
    %                   Each row represents the minmum and maximum value of temperature 
    %                   in which a phase is well-defined
    %
    % Phase_name : [cell array of cell arrays] names of phases for each I.
    %                   Each cell element is a cell array containing names of phases specfied in Phase_range.
    %                   Phase{itD}{itP} is the name of the phase defined in temperature range specified by
    %                   Phase_range{itD}(:,itP)
    %

    %% Parse options

    flatThres = 0.05;   % default value of 'flatThres'
    platLenThres = 0.5;   % default value of 'platLenThres'
    plotSusc = false;
    NFLboundary = false;
    I0min = 1e-14;

    while ~isempty(varargin)
        switch varargin{1}
            case 'flatThres'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        flatThres = varargin{2};
                    else
                        error('''flatThres'' must be positive')
                    end
                else
                    error('''flatThres'' must be a number');
                end
                varargin(1:2) = [];

            case 'platLenThres'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        platLenThres = varargin{2};
                    else
                        error('''platLenThres'' must be positive')
                    end
                else
                    error('''platLenThres'' must be a number');
                end
                varargin(1:2) = [];

            case 'NFLboundary'
                NFLboundary = true;
                varargin(1) = [];

            case 'I0min'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        I0min = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: ''I0min'' must be positive');
                    end
                else
                    error('ERR: ''I0min'' must be a number');
                end

            case '-v'
                plotSusc = true;
                varargin(1) = [];

            otherwise
                if isequal(class(varargin{1}), 'char')
                    error(['Unkown option ''',varargin{1},'''']);
                else
                    error('Unkown input');
                end

        end
    end

    %% Load the selected data folder
    path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK_selected';
    FileInfo = dir(path);

    %% Make a list of the names of data folders
    DataFolders = cell(1,numel(FileInfo));      % names of folders that hold data
    cnt = 1;
    for it = 1:numel(FileInfo)      
        DirName = FileInfo(it).name;
        if DirName(1) == 'J'        % add new data folder to cell 'DataFolders'
            DataFolders{cnt} = DirName;
            cnt = cnt + 1;
        end
    end
    DataFolders = DataFolders(1:cnt-1);

    %% Extract available parameter sets from the data folder names
    J = zeros(1,numel(DataFolders));
    K = J;
    I = J;
    T = J;
    for it = 1:numel(DataFolders)
        tmp = sscanf(DataFolders{it}, 'J0=%f_K0=%f_I0=%f_T=%f');
        J(it) = tmp(1);
        K(it) = tmp(2);
        I(it) = tmp(3);
        T(it) = tmp(4);
    end
    [I, idx] = sort(I,'ascend');
    J = J(idx);
    K = K(idx);
    T = T(idx);
    DataFolders = DataFolders(idx);

    %% Load data to be used for plotting
    RelevIdx = (J == J0) & (K == K0) & (abs(I) < 0.2) & (abs(I) >= I0min);       % Indices of folders that are relevant for plotting
    RelevData = DataFolders(RelevIdx);      % folder names of data that are relevant for plotting: ones that are (J, K) = (J0, K0)
    I = I(RelevIdx);                        % parameters 'I' that are available for plotting
    N_Relev = numel(RelevData);             % number of relevant data
    %disp(RelevIdx);
    %disp(J);
    %disp(K);
    %disp(I);

    ocont = cell(1,numel(RelevData));
    ImpSp = ocont;
    ImpOrb = ocont;
    ImpSpOrb = ocont;

    for itSub = 1:numel(RelevData)        % for all data with relevant parameters
        DataPath = [path, '\', RelevData{itSub}];
        FileInfo = dir(DataPath);

        for itD = 1:numel(FileInfo)     % load impurity susceptibilities for each relevant parameter

                switch FileInfo(itD).name
                    case 'ocont.mat'
                        tmp = load([DataPath,'\ocont.mat']);
                        field = fieldnames(tmp);
                        ocont{itSub} = getfield(tmp, field{1});

                    case 'NRG_Op=ImpSp.mat'
                        tmp = load([DataPath,'\NRG_Op=ImpSp.mat']);
                        field = fieldnames(tmp);
                        ImpSp{itSub} = getfield(tmp, field{1});

                    case 'NRG_Op=ImpOrb.mat'
                        tmp = load([DataPath,'\NRG_Op=ImpOrb.mat']);
                        field = fieldnames(tmp);
                        ImpOrb{itSub} = getfield(tmp, field{1});

                    case 'NRG_Op=ImpSpOrb.mat'
                        tmp = load([DataPath,'\NRG_Op=ImpSpOrb.mat']);
                        field = fieldnames(tmp);
                        ImpSpOrb{itSub} = getfield(tmp, field{1});

                end % switch - case
        end % itD

    end % itF

    %% Parse data

    Phase_range = cell(1,N_Relev);      % temperature ranges of phases for each dataset(I)
    Phase_name = cell(1,N_Relev);       % names of phases for each dataset(I).


    for itSub = 1:N_Relev

        %% Find plateaus in second derivatives of impurity dynamic susceptibilities

        ImpDynSusc = {ImpSp{itSub},ImpOrb{itSub}};  % impurity dynamic susceptibilities (spin, orbital)
        Range_Plateau = cell(1,2);              % temperature ranges of plateaus in second derivatives of impurity dynamic susceptibilities (spin, orbital)

        for itS = 1:2   % for spin and orbital dynamic susceptibilities (itS == 1: spin, itS == 2: orbital)
            
            [logT, Susc_2ndDer] = log_Susc_2ndDer(ocont{itSub}, ImpDynSusc{itS});    % second derivative of impurity dynamic susceptibility
    
            Susc_2ndDer_mov = movmean(Susc_2ndDer,[10,10]);       % movmean of dynamic susceptibility: to "flatten" numerical fluctuations
    
            % find plateaus in 2nd derivatives of impurity spin susceptibility
            streak = false;
            for itT = 1:numel(logT)
                if abs(Susc_2ndDer_mov(itT)) < flatThres   % if second derivative is zero
                    if ~streak
                        streak = true;
                        plateauMin = itT;
                    end
    
                else    % if second derivative is nonzero
                    if streak
                        streak = false;
                        plateauMax = itT - 1;
                        
                        if logT(plateauMax) - logT(plateauMin) > platLenThres
                            Range_Plateau{itS} = cat(1, Range_Plateau{itS}, [logT(plateauMin), logT(plateauMax)]);
                        end
                    end
    
                end
            end % itT

        end % itS

        %% Determine phases of the plateau regions found above

        Phase_range{itSub} = [];
        Phase_name{itSub} = {};
        logT = log(ocont{itSub}(ocont{itSub}>0))./log(10);
        logSp = log(ImpSp{itSub}(ocont{itSub}>0))./log(10);
        logOrb = log(ImpOrb{itSub}(ocont{itSub}>0))./log(10);

        for it1 = 1:size(Range_Plateau{1},1)
            range_Sp = Range_Plateau{1}(it1,:);

            for it2 = 1:size(Range_Plateau{2},1)
                range_Orb = Range_Plateau{2}(it2,:);

                maxMin = max(range_Sp(1), range_Orb(1));
                minMax = min(range_Sp(2), range_Orb(2));
                if maxMin < minMax      % if plateau range is well-defined    

                    [ExpSp, ~] = logFit(logT,logSp,[maxMin,minMax]);
                    [ExpOrb, ~] = logFit(logT,logOrb,[maxMin,minMax]);

                    eps = 0.1;  % a small number
                    if abs(ExpSp-1) < eps && abs(ExpOrb-1) < eps
                        Phase_range{itSub} = cat(1,Phase_range{itSub},[-24, minMax]);
                        Phase_name{itSub} = cat(1,Phase_name{itSub},{'Fermi Liquid'});

                    elseif abs(ExpSp) < eps && abs(ExpOrb) < eps
                        Phase_range{itSub} = cat(1,Phase_range{itSub},[maxMin, minMax]);
                        Phase_name{itSub} = cat(1,Phase_name{itSub},{'Fully Overscreened'});
                        
                    elseif ExpSp < -1 && ExpOrb < -1
                        Phase_range{itSub} = cat(1,Phase_range{itSub},[maxMin, max(3,minMax)]);
                        Phase_name{itSub} = cat(1,Phase_name{itSub},{'Unscreened'});
                    end
                    
                        
                end
            end % it2
        end % it1

        for it = 1:size(Range_Plateau{1},1)
            [ExpSp, ~] = logFit(logT,logSp,Range_Plateau{1}(it,:));
            [ExpOrb, ~] = logFit(logT,logOrb,Range_Plateau{1}(it,:));
            
            if abs(ExpSp) < eps && abs(ExpOrb) > eps
                Phase_range{itSub} = cat(1, Phase_range{itSub}, Range_Plateau{1}(it,:));
                Phase_name{itSub} = cat(1, Phase_name{itSub}, {'Spin Overscreened'});
            end
        end

        for it = 1:size(Range_Plateau{2},1)
            [ExpSp, ~] = logFit(logT,logSp,Range_Plateau{2}(it,:));
            [ExpOrb, ~] = logFit(logT,logOrb,Range_Plateau{2}(it,:));

            if abs(ExpOrb) < eps && abs(ExpSp) > eps
                Phase_range{itSub} = cat(1, Phase_range{itSub}, Range_Plateau{2}(it,:));
                Phase_name{itSub} = cat(1, Phase_name{itSub}, {'Orbital Overscreened'});
            end
        end 

        if ismember('Fully Overscreened', Phase_name{itSub})
            Idx1 = find(strcmp(Phase_name{itSub}, 'Fully Overscreened'));

            if ismember('Spin Overscreened', Phase_name{itSub})
                Idx2 = find(strcmp(Phase_name{itSub}, 'Spin Overscreened'));
                Phase_range{itSub}(Idx2,1) = Phase_range{itSub}(Idx1,2);

            elseif ismember('Orbital Overscreened', Phase_name{itSub})
                Idx2 = find(strcmp(Phase_name{itSub}, 'Orbital Overscreened'));
                Phase_range{itSub}(Idx2,1) = Phase_range{itSub}(Idx1,2);
            end
        end

        %% Clarify phase boundaries

        if ~isempty(Phase_range{itSub})
            if Phase_range{itSub}(1,1) > -22
                if NFLboundary
                    Phase_range{itSub} = cat(1, [-24, Phase_range{itSub}(1,1)], Phase_range{itSub});
                    Phase_name{itSub} = cat(1, {'Fermi Liquid'}, Phase_name{itSub});
                else
                end

            elseif isequal(Phase_name{itSub}{1}, 'Fermi Liquid') && numel(Phase_name{itSub}) > 1
                if NFLboundary
                    Phase_range{itSub}(1,2) = Phase_range{itSub}(2,1);
                else
                    Phase_range{itSub}(2,1) = Phase_range{itSub}(1,2);
                end
            end
        end

        
    end % itSub

    %% Plot deailed figures, if '-v' option is used

    if plotSusc

        Phase_handles = cell(1,4);
        ExistPhase = false(1,4);
        Nplot = ceil(N_Relev/23);

        for itP = 1:Nplot
            
            if itP < Nplot
                Nsub = 23;
            else
                Nsub = N_Relev - 23*(Nplot-1);
            end

            %% Plot impurity dynamic susceptibilities with phase ranges

            figure;
            hold on;
            sgtitle('$\mathrm{Impurity \ Dynamic \ Susceptibilities}$','Interpreter','latex','FontSize',20);

            for itSub = 1:Nsub
               
                DataNum = itSub + 23*(itP-1);
                ax = subplot(4,6,itSub);
                hold on;
                PosIdx = ImpSp{DataNum} > 0 & ImpOrb{DataNum};
                if ~isempty(ocont{DataNum})
                    Sp = plot(ocont{DataNum}(PosIdx),ImpSp{DataNum}(PosIdx),'color','blue','Linewidth',1,'LineStyle','-');
                    Orb = plot(ocont{DataNum}(PosIdx),ImpOrb{DataNum}(PosIdx),'color','red','Linewidth',1,'LineStyle','--');
                end
                set(ax,'XScale','log','YScale','log');
                title(['$\left(J_{0},K_{0},I_{0}\right) = \left(',sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',SciNot(I(DataNum),'Signif',1),' \right)$'], ...
                        'Interpreter','latex','FontSize',7);
        
                xlim([1e-24,1e3]);
                ylim([1e-15, 1e10]);
                Yrange = ylim;
                for it = 1:numel(Phase_name{DataNum})
                    patch_X = power(10, [Phase_range{DataNum}(it,:),Phase_range{DataNum}(it,2:-1:1)]);
                    patch_Y = [Yrange(1),Yrange(1),Yrange(2),Yrange(2)];
        
                    switch Phase_name{DataNum}{it}
                        case 'Fermi Liquid'
                            Phase_handles{1} = patch(patch_X,patch_Y,'green','FaceAlpha', 0.3, 'linestyle', 'none');
                            ExistPhase(1) = true;
                        case 'Orbital Overscreened'
                            Phase_handles{2} = patch(patch_X,patch_Y,'magenta','FaceAlpha', 0.3, 'linestyle', 'none');
                            ExistPhase(2) = true;
                        case 'Spin Overscreened'
                            Phase_handles{3} = patch(patch_X,patch_Y,'blue','FaceAlpha', 0.3, 'linestyle', 'none');
                            ExistPhase(3) = true;
                        case 'Fully Overscreened'
                            Phase_handles{4} = patch(patch_X,patch_Y,'cyan','FaceAlpha', 0.3, 'linestyle', 'none');
                            ExistPhase(4) = true;
                    end
                    
                
                end % it
                
                ax.XAxis.FontSize = 5;
                ax.YAxis.FontSize = 5;
                hold off;
            end % itSub

            handles = [Sp, Orb, Phase_handles{ExistPhase}];
            legends = {'$\chi_{\mathrm{sp}}$', '$\chi_{\mathrm{orb}}$', '$\mathrm{Fermi \ Liquid}$', ... 
                        '$\mathrm{Orbital \ Overscreened}$', '$\mathrm{Spin \ Overscreened}$', '$\mathrm{Fully \ Overscreened}$'};
            for it = 4:-1:1
                if ~ExistPhase(it)
                    legends(it+2) = [];
                end
            end
            legend(handles,legends,'Interpreter','latex','Location','bestoutside','FontSize',15);
            hold off;


            %% Plot second derivatives of impurity dynamic susceptibilities with phase ranges

            figure;
            hold on;
            sgtitle('$\mathrm{Second \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities}$','Interpreter','latex','FontSize',20);

            for itSub = 1:Nsub

                DataNum = itSub + 23*(itP-1);
                PosIdx = ImpSp{DataNum} > 0 & ImpOrb{DataNum};
                [logT, Sp_2ndDer] = log_Susc_2ndDer(ocont{DataNum}(PosIdx), ImpSp{DataNum}(PosIdx));    % second derivative of impurity spin dynamic susceptibility
                [~, Orb_2ndDer] = log_Susc_2ndDer(ocont{DataNum}(PosIdx), ImpOrb{DataNum}(PosIdx));    % second derivative of impurity spin dynamic susceptibility
                Sp_2ndDer_mov = movmean(Sp_2ndDer,[10,10]);
                Orb_2ndDer_mov = movmean(Orb_2ndDer,[10,10]);

                ax = subplot(4,6,itSub);
                hold on;
                
                if ~isempty(ocont{DataNum})
                    Sp = plot(power(10,logT),Sp_2ndDer_mov,'color','blue','Linewidth',1,'LineStyle','-');
                    Orb = plot(power(10,logT),Orb_2ndDer_mov,'color','red','Linewidth',1,'LineStyle','--');
                end
                set(ax,'XScale','log','YScale','linear');
                title(['$\left(J_{0},K_{0},I_{0}\right) = \left(',sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',SciNot(I(DataNum),'Signif',1),' \right)$'], ...
                        'Interpreter','latex','FontSize',7);
        
                xlim([1e-24,1e3]);
                ylim([-1,1]);
                Yrange = ylim;
                for it = 1:numel(Phase_name{DataNum})
                    patch_X = power(10, [Phase_range{DataNum}(it,:),Phase_range{DataNum}(it,2:-1:1)]);
                    patch_Y = [Yrange(1),Yrange(1),Yrange(2),Yrange(2)];
        
                    switch Phase_name{DataNum}{it}
                        case 'Fermi Liquid'
                            Phase_handles{1} = patch(patch_X,patch_Y,'green','FaceAlpha', 0.3, 'linestyle', 'none');
                        case 'Orbital Overscreened'
                            Phase_handles{2} = patch(patch_X,patch_Y,'magenta','FaceAlpha', 0.3, 'linestyle', 'none');
                        case 'Spin Overscreened'
                            Phase_handles{3} = patch(patch_X,patch_Y,'blue','FaceAlpha', 0.3, 'linestyle', 'none');
                        case 'Fully Overscreened'
                            Phase_handles{4} = patch(patch_X,patch_Y,'cyan','FaceAlpha', 0.3, 'linestyle', 'none');
                    end
                    
                
                end % it
                
                ax.XAxis.FontSize = 5;
                ax.YAxis.FontSize = 5;
                hold off;
            end % itSub

            handles = [Sp, Orb, Phase_handles{ExistPhase}];
            legends = {'$\chi_{\mathrm{sp}}$', '$\chi_{\mathrm{orb}}$', '$\mathrm{Fermi \ Liquid}$', ... 
                        '$\mathrm{Orbital \ Overscreened}$', '$\mathrm{Spin \ Overscreened}$', '$\mathrm{Fully \ Overscreened}$'};
            for it = 4:-1:1
                if ~ExistPhase(it)
                    legends(it+2) = [];
                end
            end
            legend(handles,legends,'Interpreter','latex','Location','bestoutside','FontSize',15);
            hold off;


        end % itP
    end




    %% function for calculating the log-2nd-derivative of dynamic susceptibilities

    function [logT_grid_2ndDer, log_Susc_2ndDer] = log_Susc_2ndDer(T_grid, Susc)    

        logT_grid = log(T_grid(T_grid>0))./log(10);      % log temperature grid
        log_Susc = log(Susc(T_grid>0))./log(10);

        log_Susc_1stDer = diff(log_Susc,1)./diff(logT_grid,1);      % first derivative
        logT_grid_1stDer = movmean(logT_grid, [0,1]);           % log temperature grid for first derivative
        logT_grid_1stDer = logT_grid_1stDer(1:end-1);

        log_Susc_2ndDer = diff(log_Susc_1stDer,1)./diff(logT_grid_1stDer,1);    % Second derivative
        logT_grid_2ndDer = movmean(logT_grid_1stDer, [0,1]);                    % log temperature grid for second derivative
        logT_grid_2ndDer = logT_grid_2ndDer(1:end-1);   
    end

    %% function for fitting impurity dynamic susceptibilities

    function [Exp, Rsq] = logFit(logT, ImpDyn, fit_range)

        x = logT(logT > fit_range(1));
        y = ImpDyn(logT > fit_range(1));
        y = y(x < fit_range(2));
        x = x(x < fit_range(2));
        coeff = polyfit(x,y,1);
        Exp = coeff(1);                     % fitted critical exponent
        yfit = polyval(coeff, x);          % Estimated  Regression Line
        SStot = sum((y-mean(y)).^2);        % Total Sum-Of-Squares
        SSres = sum((y-yfit).^2);           % Residual Sum-Of-Squares
        Rsq = 1-SSres/SStot;                % R^2
    end

end