clear;

figureHandle = figure('Position', [50, -30, 550, 850]);

% RGflow with I0 == 0
J_0 = [0.1, 1, 0.1, 0.15];
K_0 = [1, 0.1, 0.15, 0.1];

J_0 = [J_0, -0.9, -0.9, -0.9, -0.9];
K_0 = [K_0, 1, 0.7, 0.3, 0.15];

J_0 = [J_0, 1, 0.7, 0.3, 0.15];
K_0 = [K_0, -0.9, -0.9, -0.9, -0.9];

J_0 = [J_0, -0.9, -0.03, -0.9];
K_0 = [K_0, -0.9, -0.9, -0.03];

I_0 = zeros(1,15);

PM_RGflow_zeroI(J_0, K_0, zeros(1,15));

% RGflow with I0 ~= 0
J_0 = [0.03, 0.2, 0.2];
K_0 = [0.2, 0.2, 0.03];

J_0 = [J_0, -0.5, -0.5, -0.5, -0.5];
K_0 = [K_0, 1, 0.7, 0.3, 0.03];

J_0 = [J_0, 1, 0.7, 0.3, 0.03];
K_0 = [K_0, -0.5, -0.5, -0.5, -0.5];

J_0 = [J_0, -0.5, -0.03, -0.5];
K_0 = [K_0, -0.5, -0.5, -0.03];

I_0 = 0.01*ones(1,15);

PM_RGflow_finiteI(J_0, K_0, I_0)

% subplot labels
%{}
annotation(figureHandle, 'textbox', [0.06, 0.975, 0, 0], ...
    'String', '$\mathrm{(a)}$', 'Interpreter', 'latex', 'FontSize', 18, ...
    'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'LineStyle', 'none');

annotation(figureHandle, 'textbox', [0.06, 0.515, 0, 0], ...
    'String', '$\mathrm{(b)}$', 'Interpreter', 'latex', 'FontSize', 18, ...
    'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'LineStyle', 'none');
%}

fig = gcf;

set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'inches');
fig_pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', fig_pos(3:4));

% Specify full path
output_path = 'C:\Users\hsjun\OneDrive\Physics\Research\TsoK_publication\Figures\Poorman_RG_JK_subplots.pdf';
set(fig, 'Renderer', 'painters');
print(fig, output_path, '-dpdf');

%% Function to plot RG flow for I0 == 0

function PM_RGflow_zeroI(J0, K0, I0)

    Nflow = numel(J0);

    Kmin = -0.5;
    Kmax = 1.05;
    Jmin = -0.5;
    Jmax = 1.05;

    subplot1 = subplot(2, 1, 1); 
    set(subplot1, 'Position', [0.15, 0.51, 0.75, 0.44]); % [left bottom width height]
    hold on;
    set(gca, 'Layer', 'top');

    % define plot positions
    ax = gca;
    pos = ax.Position;  % [left bottom width height]
    pos(2) = pos(2) + 0.02;
    ax.Position = pos;

    set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(ax, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

    %xlabel('$K_{0}$', 'Interpreter', 'latex', 'FontSize', 24);
    ylabel('$J_{0}$', 'Interpreter', 'latex', 'FontSize', 20);

    hx = get(ax, 'XLabel');
    hy = get(ax, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';

    % Modify x- and y-label positions manually
    hx.Position = hx.Position + [0, -0.005, 0];
    hy.Position = hy.Position + [-0.02, 0, 0];

    % redefine x-tick labels using annotation()
    tickLabelFont = 18;
    ax.XTickLabel = [];
    xticks = [-0.4, -0.2, 0, 0.3, 0.6, 1];
    ax.XTick = xticks;

    % redefine y-tick labels using annotation()
    ax.YTickLabel = [];
    yticks = [-0.4, -0.2, 0, 0.3, 0.6, 1];
    ax.YTick = yticks;
    
    for ity = 1:numel(yticks)
        text(Kmin-0.025, yticks(ity)+0.04, ...
        sprintf('$%.15g$', yticks(ity)), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickLabelFont);
    end

    xlim([Kmin, Kmax]);
    ylim([Jmin, Jmax]);

    %% Plot phase patches

    % define custom colors
    patch_blue = [.910, .906, .969];
    patch_purple = [.937, .910, .973];
    patch_orange = [.984, .957, .910];
    patch_rose = [0.957, 0.835, 0.875];
    patch_gray = 0.85 * [1, 1, 1];
    

    % plot phase patches
    X_FO = [0, Kmax, Kmax, 0];  Y_FO = [0, 0, Jmax, Jmax];
    X_SO = [Kmin, 0, 0, Kmin];  Y_SO = [0, 0, Jmax, Jmax];
    X_OO = [0, Kmax, Kmax, 0];  Y_OO = [0, 0, Jmin, Jmin];
    X_FI = [Kmin, 0, 0, Kmin];  Y_FI = [0, 0, Jmin, Jmin];

    FaceAlpha = 1;
    patch(X_FO, Y_FO, patch_purple, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);
    patch(X_SO, Y_SO, patch_orange, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);
    patch(X_OO, Y_OO, patch_rose, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);
    patch(X_FI, Y_FI, patch_gray, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);

    %% Calculate and plot RG flows from poor man's scaling results
    cmap = turbo(Nflow);

    for itN = 1:Nflow
        [J,K,I,log_D] = RG_PoorMan(J0(itN), K0(itN), I0(itN), 'order', 3);

        Nsteps = 75;
        Steps = round(linspace(10, numel(log_D)-1, Nsteps));

        for itS = 1:Nsteps
            NumD = Steps(itS);
            dK = K(NumD+1) - K(NumD-1);
            dJ = J(NumD+1) - J(NumD-1);
            dS = sqrt(dK^2 + dJ^2);

            %height = 0.1 / abs(log(dS));
            height = min(0.1, 45*dS);
            angle = atan(dJ/dK) - pi/2;

            [X, Y] = getArrow(K(NumD), J(NumD), 0.045, 0.045, 0.015, height, angle);

            patch(X, Y, cmap(itN,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
            %patch(X, Y, 'black', 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
        
    end

    %% Mark RG fixed points

    FI = [0, 0];

    plot(FI(1), FI(2), 's', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.39, 0.62, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FI}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    SO = [0, 1];

    plot(SO(1), SO(2), 'X', ...
    'MarkerSize', 15, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.39, 0.91, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{SO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    J_OO = [1, 0];

    plot(J_OO(1), J_OO(2), '^', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.77, 0.62, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{OO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    FO = [1, 1];

    plot(FO(1), FO(2), 'o', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.88, 0.9, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

end

%% Function to plot RG flow for finite I0

function PM_RGflow_finiteI(J0, K0, I0)

    Nflow = numel(J0);

    Kmin = -0.5;
    Kmax = 1.05;
    Jmin = -0.5;
    Jmax = 1.05;

    subplot2 = subplot(2, 1, 2);
    set(subplot2, 'Position', [0.15, 0.05, 0.75, 0.44]); % [left bottom width height]
    hold on;
    set(gca, 'Layer', 'top');

    % define plot positions
    ax = gca;
    pos = ax.Position;  % [left bottom width height]
    pos(2) = pos(2) + 0.02;
    ax.Position = pos;

    set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(ax, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

    xlabel('$K_{0}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$J_{0}$', 'Interpreter', 'latex', 'FontSize', 20);

    hx = get(ax, 'XLabel');
    hy = get(ax, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';

    % Modify x- and y-label positions manually
    hx.Position = hx.Position + [0, -0.005, 0];
    hy.Position = hy.Position + [-0.02, 0, 0];

    % redefine x-tick labels using annotation()
    tickLabelFont = 18;
    ax.XTickLabel = [];
    xticks = [-0.4, -0.2, 0, 0.3, 0.6, 1];
    ax.XTick = xticks;
    
    for itx = 1:numel(xticks)
        text(xticks(itx), Jmin - 0.005, ...
        sprintf('$%.15g$', xticks(itx)), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickLabelFont);
    end

    % redefine y-tick labels using annotation()
    ax.YTickLabel = [];
    yticks = [-0.4, -0.2, 0, 0.3, 0.6, 1];
    ax.YTick = yticks;
    
    for ity = 1:numel(yticks)
        text(Kmin-0.025, yticks(ity)+0.04, ...
        sprintf('$%.15g$', yticks(ity)), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickLabelFont);
    end

    xlim([Kmin, Kmax]);
    ylim([Jmin, Jmax]);

    %% Plot phase patches

    %{}
    % define custom colors
    patch_blue = [.910, .906, .969];
    patch_purple = [.937, .910, .973];
    patch_orange = [.984, .957, .910];
    patch_rose = [0.957, 0.835, 0.875];
    patch_gray = 0.85 * [1, 1, 1];
    patch_green = [.925, .973, .910];

    % plot phase patches
    X_FL = [Kmin, Kmax, Kmax, Kmin];   Y_FL = [Jmin, Jmin, Jmax, Jmax];
    X_FI = [Kmin, 0, 0, Kmin];  Y_FI = [0, 0, Jmin, Jmin];

    FaceAlpha = 1;
    patch(X_FL, Y_FL, patch_green, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);
    patch(X_FI, Y_FI, patch_gray, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);

    %% Calculate and plot RG flows from poor man's scaling results
    cmap = turbo(Nflow);

    for itN = 1:Nflow
        [J,K,~,log_D] = RG_PoorMan(J0(itN), K0(itN), I0(itN), 'order', 3);

        Nsteps = 100;
        Steps = round(linspace(10, numel(log_D)-1, Nsteps));

        for itS = 1:Nsteps
            NumD = Steps(itS);
            dK = K(NumD+1) - K(NumD-1);
            dJ = J(NumD+1) - J(NumD-1);
            dS = sqrt(dK^2 + dJ^2);

            %height = 0.1 / abs(log(dS));
            height = min(0.1, 45*dS);
            angle = atan(dJ/dK) - pi/2;

            [X, Y] = getArrow(K(NumD), J(NumD), 0.045, 0.045, 0.015, height, angle);

            patch(X, Y, cmap(itN,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
            %patch(X, Y, 'black', 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
        
    end

    %% Mark RG fixed points

    FI = [0, 0];

    plot(FI(1), FI(2), 's', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.39, 0.16, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FI}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    FL = [1, 1];

    plot(FL(1), FL(2), 'diamond', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.89, 0.43, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FL}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

end

%% Function to plot RG flow on K0-I0 section

function PM_RGflow_KI(J0, K0, I0)

    Nflow = numel(J0);

    Kmin = -0.5;
    Kmax = 1.05;
    Imin = 0;
    Imax = 1.05;

    subplot1 = subplot(2, 1, 1);
    set(subplot1, 'Position', [0.15, 0.51, 0.75, 0.44]);
    hold on;
    set(gca, 'Layer', 'top');

    % define plot positions
    ax = gca;
    pos = ax.Position;  % [left bottom width height]
    pos(2) = pos(2) + 0.02;
    ax.Position = pos;

    set(ax, 'XScale', 'linear', 'YScale', 'linear');

    set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(ax, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

    xlabel('$I_{0}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$K_{0}$', 'Interpreter', 'latex', 'FontSize', 20);

    hx = get(ax, 'XLabel');
    hy = get(ax, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';

    % Modify x- and y-label positions manually
    hx.Position = hx.Position + [0, -0.005, 0];
    hy.Position = hy.Position + [-0.02, 0, 0];

    xlim([Imin, Imax]);
    ylim([Kmin, Kmax]);

    %% Plot phase patches

    % define custom colors
    patch_blue = [.910, .906, .969];
    patch_purple = [.937, .910, .973];
    patch_orange = [.984, .957, .910];
    patch_rose = [0.957, 0.835, 0.875];
    patch_gray = 0.85 * [1, 1, 1];
    

    % plot phase patches
    %X_FO = [0, Kmax, Kmax, 0];  Y_FO = [0, 0, Jmax, Jmax];
    %X_SO = [Kmin, 0, 0, Kmin];  Y_SO = [0, 0, Jmax, Jmax];
    %X_OO = [0, Kmax, Kmax, 0];  Y_OO = [0, 0, Jmin, Jmin];
    %X_FI = [Kmin, 0, 0, Kmin];  Y_FI = [0, 0, Jmin, Jmin];

    FaceAlpha = 1;
    %patch(X_FO, Y_FO, patch_purple, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);
    %patch(X_SO, Y_SO, patch_rose, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);
    %patch(X_OO, Y_OO, patch_orange, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);
    %patch(X_FI, Y_FI, patch_gray, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);

    %% Calculate and plot RG flows from poor man's scaling results
    cmap = turbo(Nflow);

    for itN = 1:Nflow
        [J,K,I,log_D] = RG_PoorMan(J0(itN), K0(itN), I0(itN), 'order', 3);

        Nsteps = 75;
        Steps = round(linspace(10, numel(log_D)-1, Nsteps));

        for itS = 1:Nsteps
            NumD = Steps(itS);
            dI = I(NumD+1) - I(NumD-1);
            dK = K(NumD+1) - K(NumD-1);
            dS = sqrt(dI^2 + dK^2);

            %height = 0.1 / abs(log(dS));
            height = min(0.1, 45*dS);
            angle = atan(dK/dI) - pi/2;

            [X, Y] = getArrow(I(NumD), K(NumD), 0.045, 0.045, 0.015, height, angle);

            patch(X, Y, cmap(itN,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
            %patch(X, Y, 'black', 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
        
    end

    %% Mark RG fixed points

    FI = [0, 0];

    plot(FI(1), FI(2), 's', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.39, 0.62, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FI}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    SO = [0, 1];

    plot(SO(1), SO(2), 'X', ...
    'MarkerSize', 15, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.39, 0.91, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{SO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    J_OO = [1, 0];

    plot(J_OO(1), J_OO(2), '^', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.77, 0.62, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{OO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    FO = [1, 1];

    plot(FO(1), FO(2), 'o', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.88, 0.9, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

end



function PM_RGflow_legend(J0, K0, I0)
    
    Nflow = numel(J0);
    legends = cell(1,Nflow);

    cmap = turbo(Nflow);

    for itN = 1:Nflow
        legends{itN} = ['$(K_{0}, J_{0}, I_{0})=(',sprintf('%.15g',K0(itN)),',', ...
                            sprintf('%.15g',J0(itN)),',',sprintf('%.15g',I0(itN)),')$'];
    end

    %% Plot custom legends using patch and annotation
    figure;
    hold on;
    set(gcf, 'Position', [100, 100, 800, 400]);
    xlim([0,2]);
    ylim([0,1]);

    X_arrow = [0.13, 0.77, 1.43];
    Y_arrow = linspace(0.5, 0.77, 5);

    X_label = X_arrow + 0.01;
    Y_label = 0.635*ones(1,3);

    for it1 = 1:3

        lgdStr = '$\begin{array}{l@{\hskip 1pt}l@{\hskip 1pt}l}';
        for it2 = 1:5

            itN = 5*(it1-1) + it2;

            [X, Y] = getArrow(X_arrow(it1), Y_arrow(it2), 0.04, 0.04, 0.015, 0.07, -pi/2);
            patch(X, Y, cmap(itN,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
    
            if it2 < 5
                lgdStr = [lgdStr, '\mathbf{c}_{0} = (',sprintf('%.15g',J0(itN)),', &',sprintf('%.15g',K0(itN)), ...
                            ', &',sprintf('%.15g',I0(itN)),') \\ '];
            else
                lgdStr = [lgdStr, '\mathbf{c}_{0} = (',sprintf('%.15g',J0(itN)),', &',sprintf('%.15g',K0(itN)), ...
                            ', &',sprintf('%.15g',I0(itN)),')', ...
                                '\end{array}$'];
            end
        end

        text(X_label(it1), Y_label(it1), lgdStr, 'Interpreter', 'latex', 'FontSize', 13.9);
    end

    hold off;
end



function [X_all, Y_all] = getArrow(cx, cy, tri_width, tri_height, rect_width, rect_height, angle)
    % Returns coordinates for an arrow composed of:
    % - An isosceles triangle with base width `tri_width` and height `tri_height`
    % - A rectangle with width rect_width` and height `rect_height`
    % - Rotated by 'angle' around the circumcenter, placed so the rectangle
    %   sits flush below the triangle's base
    % - The overall shape is centered at the triangle's circumcenter at (cx, cy)
    
    % === STEP 1: Build triangle with circumcenter at origin ===
    x1 = -tri_width/2; y1 = 0;
    x2 =  tri_width/2; y2 = 0;
    x3 = 0;            y3 = tri_height;
    
    % Compute circumcenter of triangle
    D = 2*(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
    Ux = ((x1^2 + y1^2)*(y2 - y3) + ...
          (x2^2 + y2^2)*(y3 - y1) + ...
          (x3^2 + y3^2)*(y1 - y2)) / D;
    Uy = ((x1^2 + y1^2)*(x3 - x2) + ...
          (x2^2 + y2^2)*(x1 - x3) + ...
          (x3^2 + y3^2)*(x2 - x1)) / D;
    
    % Triangle coordinates centered at origin
    X_tri = [x1, x2, x3] - Ux;
    Y_tri = [y1, y2, y3] - Uy;
    
    % === STEP 2: Create rectangle just below triangle base ===
    % Rectangle is centered horizontally, and directly below base of triangle
    % Its top edge lies along y = 0 (same as triangle base)
    rect_x1 = -rect_width/2;
    rect_x2 =  rect_width/2;
    rect_y1 = -rect_height;
    rect_y2 = 0;
    
    X_rect = [rect_x1, rect_x2, rect_x2, rect_x1] - Ux;
    Y_rect = [rect_y2, rect_y2, rect_y1, rect_y1] - Uy;
    
    % === STEP 3: Rotate both triangle and rectangle ===
    R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    
    XY_tri_rot = R * [X_tri; Y_tri];
    XY_rect_rot = R * [X_rect; Y_rect];
    
    % === STEP 4: Translate everything to target circumcenter (cx, cy) ===
    X_tri_final = XY_tri_rot(1, :) + cx;
    Y_tri_final = XY_tri_rot(2, :) + cy;
    
    X_rect_final = XY_rect_rot(1, :) + cx;
    Y_rect_final = XY_rect_rot(2, :) + cy;
    
    % Combine shapes into one patchable polygon
    X_all = [X_rect_final(4:-1:2), X_tri_final([2,3,1]), X_rect_final(1)];
    Y_all = [Y_rect_final(4:-1:2), Y_tri_final([2,3,1]), Y_rect_final(1)];
end
