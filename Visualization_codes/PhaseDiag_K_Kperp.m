function PhaseDiag_K_Kperp(J0, I0)
    
    useEnt = true;

    %% Load the selected data folder
    path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK_Aniso_selected';
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
    K_perp = J;
    K_z = J;
    I = J;
    T = J;
    Nkeep = J;
    for it = 1:numel(DataFolders)
        tmp = sscanf(DataFolders{it}, 'J0=%f_K_perp=%f_K_z=%f_I0=%f_T=%f_Nkeep=%f');
        J(it) = tmp(1);
        K_perp(it) = tmp(2);
        K_z(it) = tmp(3);
        I(it) = tmp(4);
        T(it) = tmp(5);
        Nkeep(it) = tmp(6);
    end

    %% Load data to be used for plotting
    RelevIdx = (J == J0) & (I == I0) & (Nkeep == 6000);       % Indices of folders that are relevant for plotting
    RelevData = DataFolders(RelevIdx);      % folder names of data that are relevant for plotting: ones that are (J, K) = (J0, K0)
    J = J(RelevIdx);                        % parameters 'J' that are available for plotting
    K_perp = K_perp(RelevIdx);
    K_z = K_z(RelevIdx);
    N_Relev = numel(RelevData);             % number of relevant data

    ocont = cell(1,numel(RelevData));
    ImpSp = ocont;
    ImpOrb_plus = ocont;
    ImpOrb_z = ocont;
    Sent = ocont;
    DiscData = ocont;

    for itSub = 1:numel(RelevData)        % for all data with relevant parameters
        DataPath = [path, filesep, RelevData{itSub}];
        FileInfo = dir(DataPath);

        for itD = 1:numel(FileInfo)     % load impurity susceptibilities and entropy for each relevant parameter

                switch FileInfo(itD).name
                    case 'ocont.mat'
                        tmp = load([DataPath,'\ocont.mat']);
                        field = fieldnames(tmp);
                        ocont{itSub} = getfield(tmp, field{1});

                    case 'NRG_Op=ImpOrb_plus.mat'
                        tmp = load([DataPath,'\NRG_Op=ImpOrb_plus.mat']);
                        field = fieldnames(tmp);
                        ImpOrb_plus{itSub} = getfield(tmp, field{1});

                    case 'NRG_Op=ImpOrb_z.mat'
                        tmp = load([DataPath,'\NRG_Op=ImpOrb_z.mat']);
                        field = fieldnames(tmp);
                        ImpOrb_z{itSub} = getfield(tmp, field{1});

                    case 'NRG_Op=ImpSp.mat'
                        tmp = load([DataPath,'\NRG_Op=ImpSp.mat']);
                        field = fieldnames(tmp);
                        ImpSp{itSub} = getfield(tmp, field{1});

                    case 'DiscData.mat'
                        tmp = load([DataPath,'\DiscData.mat']);
                        field = fieldnames(tmp);
                        DiscData{itSub} = getfield(tmp, field{1});

                    case 'Temps.mat'
                        tmp = load([DataPath,'\Temps.mat']);
                        field = fieldnames(tmp);
                        Temps = getfield(tmp, field{1});

                    case 'Sent_imp_beta=1.5.mat'
                        tmp = load([DataPath,'\Sent_imp_beta=1.5.mat']);
                        field = fieldnames(tmp);
                        Sent{itSub} = getfield(tmp, field{1});

                end % switch - case
        end % itD

    end % itF

    %% extract impurity orbital_z dynamic susceptibility power law

    if ~useEnt
        fit_range = [-6.5, -3];
        Orbz_exp = zeros(1, N_Relev);
        OrbP_exp = zeros(1, N_Relev);
    
        for itD = 1:N_Relev
            
            log_T = log10(ocont{itD}(ocont{itD}>0));          % log temperatures
            log_Susc = log10(ImpOrb_z{itD}(ocont{itD}>0));
    
            log_Susc = log_Susc(log_T>fit_range(1) & log_T<fit_range(2));
            log_T = log_T(log_T>fit_range(1) & log_T<fit_range(2));
    
            % fit 
            f = polyfit(log_T, log_Susc, 1);
            Orbz_exp(itD) = f(1);
    
            log_T = log10(ocont{itD}(ocont{itD}>0));          % log temperatures
            log_Susc = log10(ImpOrb_plus{itD}(ocont{itD}>0));
    
            log_Susc = log_Susc(log_T>fit_range(1) & log_T<fit_range(2));
            log_T = log_T(log_T>fit_range(1) & log_T<fit_range(2));
    
            % fit 
            f = polyfit(log_T, log_Susc, 1);
            OrbP_exp(itD) = f(1);
    
        end

    end

    %% Plot phase diagram with entropy colormap
    figure('Position', [100, 200, 550, 400]);
    hold on;
    ax = gca;
    pos = [0.12, 0.15, 0.8, 0.82];
    set(ax, 'Position', pos);
    set(ax, 'XScale', 'linear', 'YScale', 'linear');
    set(ax, 'XTick', -0.5:0.2:0.5);
    set(ax, 'YTick', 0:0.1:0.5);
    set(ax, 'FontSize', 13);
    set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(ax, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);
    ax.TickLength = [0.015, 0.002];      % ticksize : [major, minor]

    Kz_min = -0.5;
    Kz_max = 0.5;

    K_perp_min = 0;
    K_perp_max = 0.5;

    xlim([Kz_min, Kz_max]);
    ylim([K_perp_min, K_perp_max]);


    % define x- and y-labels
    xlabel('$K_{z}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$K_{\perp}$', 'Interpreter', 'latex', 'FontSize', 20);
    hx = get(ax, 'XLabel');
    hy = get(ax, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';
    
    % Modify x- and y-label positions manually
    hx.Position = hx.Position + [0, 0, 0];
    hy.Position = hy.Position + [-0.01, 0.02, 0];
    
    % redefine x-tick labels using annotation()
    tickLabelFont = 18;
        
    ax.XTickLabel = [];
    xticks = Kz_min : 0.2 : Kz_max;
    xticks = 0.1*round(10*xticks+1e-10);
        
    xtickWidth = 0.09;
    xtickHeight = 0.02;
    for itx = 1:numel(xticks)
        X_pos = pos(1) + pos(3) * (xticks(itx) - Kz_min) / 1.15;
        Y_pos = pos(2);
        X_pos = X_pos - 0.05;
        Y_pos = Y_pos - 0.015;
        annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', ['$',sprintf('%.15g',xticks(itx)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end
    
    % redefine y-tick labels using annotation()
    ax.YTickLabel = [];
    yticks = K_perp_min : 0.1 : K_perp_max;
    yticks = 0.1*round(10*yticks+1e-10);
        
    ytickWidth = 0.09;
    ytickHeight = 0.03;
    for ity = 1:numel(yticks)
        X_pos = pos(1);
        Y_pos = pos(2) + pos(4) * (yticks(ity) - K_perp_min) / (K_perp_max - K_perp_min);
        X_pos = X_pos - 0.09;
        Y_pos = Y_pos + 0.01;
        annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',sprintf('%.15g', yticks(ity)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end

    % define and plot colormap
    eps = 1e-5;
    Kp_intvl = 0.025;
    N_Kp = 0.5/Kp_intvl + 1;
    Cmap = nan(N_Kp, 21);
    PB = zeros(3,21);   % K_perp at phase boundary

    for itz = 1:21
        PB_found = false(1,4);

        for itp = 1:N_Kp

            Kp = Kp_intvl * itp - Kp_intvl;
            Kz = 0.05*itz - 0.55;
            idx = find(J == J0 & abs(K_perp-Kp) < eps & abs(K_z-Kz) < eps);

            if useEnt
                if ~isempty(idx)
                    if ~isempty(Sent{idx})
                        Cmap(itp,itz) = exp(Sent{idx}(end));
                        
                        for itT = 1:3
                            if exp(Sent{idx}(50*itT)) < 2*sqrt(2) - 0.05 && ~PB_found(itT)
                                PB(itT,itz) = Kp;
                                PB_found(itT) = true;
                            end
                        end

                    else
                        Cmap(itp,itz) = 0;    
                    end
                end
            else
                %Cmap(itp, itz) = Orbz_exp(idx) - OrbP_exp(idx);
                Cmap(itp, itz) = Orbz_exp(idx);
                %Cmap(itp, itz) = OrbP_exp(idx);
            end 
        end
    end

    if useEnt
        for itz = 1:21
            for itp = 1:N_Kp

                Kp = Kp_intvl*itp - Kp_intvl;
                
                if Cmap(itp,itz) == 0 || isnan(Cmap(itp,itz))
                    step = 0.05;
                    min_val = 0;
                    max_val = 0.5;
                    
                    % Round to nearest multiple of 0.05
                    if Kp > PB(end,itz)
                        rounded = ceil(Kp/step) * step;
                    else
                        rounded = floor(Kp/step) * step;
                    end
                    
                    % Clamp within bounds
                    clamped = min(max(rounded, min_val), max_val);
                    itp_closest = round(clamped/Kp_intvl) + 1;

                    Cmap(itp,itz) = Cmap(itp_closest,itz);
                end

            end % itp
        end % itz
    end

    patch_purple = [.937, .910, .973];
    line_purple = [.455, .165, .765];
    patch_teal = [.627, .843, .843];
    line_teal = [.100, .600, .600];
    
    % Number of steps
    n = 256;
    
    % Interpolation
    r = linspace(patch_purple(1), patch_teal(1), n)';
    g = linspace(patch_purple(2), patch_teal(2), n)';
    b = linspace(patch_purple(3), patch_teal(3), n)';
    
    custom_cmap = [r, g, b];

    [arrow_X, arrow_Y] = meshgrid(-0.525:0.05:0.525, -Kp_intvl/2 : Kp_intvl : 0.5+Kp_intvl/2);
    Z = ones(size(arrow_X));  % some positive number for height
    surf(arrow_X, arrow_Y, Z, Cmap, 'EdgeColor', 'none');
    shading flat;      
    colormap(custom_cmap);
    cb = colorbar;
    yl = ylabel(cb, '$\exp(S_{\mathrm{imp}})$', 'Interpreter', 'latex');
    yl.Units = 'Normalized';
    yl.Position = yl.Position + [-0.01, 0, 0];

    cb.Ticks = [2:0.2:2.6, 2*sqrt(2)];
    cb.TickLabels = {'$2$', '$2.2$', '$2.4$', '$2.6$', '$2 \sqrt{2}$'};
    set(cb, 'TickLabelInterpreter', 'latex', 'FontSize', tickLabelFont-4);
    cb.Label.String = '$\exp(S_{\mathrm{imp}})$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = 18;


    Kz = -0.5:0.05:0.5;
    Z = 5*ones(1,21);
    legends = cell(1,4);
    handles = [];
    colors = [.466, .674, .188;
              .929, .694, .125;
              .850, .325, .098;
              0, .447, .741];
    linestyles = {'-.', '--', '-.', '-'};

    for itT = 1:3
        PB_now = PB(itT,:);
        hdl = plot3(Kz(PB_now>0), PB_now(PB_now>0), Z(PB_now>0), linestyles{itT}, 'Color', colors(itT,:), 'LineWidth', 2);
        handles = [handles,hdl];
        legends{itT} = ['$ T = ', SciNot(sqrt(2.5)^-(50*itT)), '$'];
    end

    %plot3(Kz(PB>0), PB(PB>0), Z(PB>0), '-', 'Color', 'red', 'LineWidth', 2);
    hdl = plot3(Kz, -Kz, Z, '--', 'Color', [.2,.2,.2], 'LineWidth', 1.5);
    handles = [handles, hdl];
    legends{4} = '$\mathrm{ poor \ man''s \ scaling}$';

    legend(handles, legends, 'Interpreter', 'latex', 'FontSize', 15);

    FS = 15;
    text_X = 0.02;
    text_Y = 0.25;
    text(text_X, text_Y, 5, '$\mathrm{Fully \ Overscreened}$', 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);

    text_X = -0.38;
    text_Y = 0.1;
    text(text_X, text_Y, 5, '$\mathrm{Orbital}$', 'Color', line_teal, 'Interpreter', 'latex', 'FontSize', FS);
    text_X = -0.45;
    text_Y = 0.06;
    text(text_X, text_Y, 5, '$\mathrm{Ferromagnetic}$', 'Color', line_teal, 'Interpreter', 'latex', 'FontSize', FS);

    if useEnt
        clim([2,2*sqrt(2)]);            % Set color value limits
    else
        %clim([0,0.15]);
        clim([-1.2,-0.8]);
        %clim([-1.2,-0.8]);
    end

    view(2);
    set(gca, 'Layer', 'top');

    fig = gcf;

    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'Units', 'inches');
    fig_pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperSize', fig_pos(3:4));
    
    % Specify full path
    output_path = 'C:\Users\hsjun\OneDrive\Physics\Research\TsoK_publication\Figures\Phase_Diagram_Kperp_Kz.pdf';
    set(fig, 'Renderer', 'painters');
    print(fig, output_path, '-dpdf');

    hold off;


end