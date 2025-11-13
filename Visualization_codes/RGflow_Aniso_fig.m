clear;

%{
J_0 = [0.1*ones(1,5), -0.9*ones(1,4), 0.03, -0.03, -0.9];
K_0 = [-0.9, -0.1, 0.15, 0.4, 1, -0.9, 0.15, 0.3, 1, 0, -0.9, 0];
I_0 = zeros(1,12);
%}

Kp0 = [0.3, 0.5, 0.8];
Kz0 = [-0.9, -0.9, -0.9];

Kp0 = [Kp0, 0.02, 0.02, 0.02];
Kz0 = [Kz0, 0.3, 0.5, 0.8];

Kp0 = [Kp0, 0.95, 0.95, 0.95];
Kz0 = [Kz0, -0.85, -0.7, -0.55];

J0 = 0.3*ones(1,numel(Kp0));
Ip0 = zeros(1,numel(Kp0));
Iz0 = zeros(1,numel(Kp0));

I_0 = zeros(1,15);

PM_RGflow(J0, Kp0, Kz0, Ip0, Iz0);
%PM_RGflow_legend(J_0, K_0, I_0);

function PM_RGflow(J0, Kp0, Kz0, Ip0, Iz0)

    Nflow = numel(J0);

    Kz_min = -1.05;
    Kz_max = 1.05;
    Kp_min = 0;
    Kp_max = 1.05;

    figure;
    hold on;
    set(gca, 'Layer', 'top');
    set(gcf, 'Position', [100, 100, 550, 400]);

    % define plot positions
    ax = gca;
    pos = ax.Position;  % [left bottom width height]
    pos(2) = pos(2) + 0.04;
    ax.Position = pos;

    set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(ax, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

    xlabel('$K_{z}$', 'Interpreter', 'latex', 'FontSize', 24);
    ylabel('$K_{\perp}$', 'Interpreter', 'latex', 'FontSize', 24);

    hx = get(ax, 'XLabel');
    hy = get(ax, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';

    % Modify x- and y-label positions manually
    hx.Position = hx.Position + [0, 0, 0];
    hy.Position = hy.Position + [-0.02, 0, 0];

    % redefine x-tick labels using annotation()
    tickLabelFont = 18;
    ax.XTickLabel = [];
    xticks = -1:0.4:1;
    ax.XTick = xticks;
    
    for itx = 1:numel(xticks)
        text(xticks(itx), Kp_min - 0.005, ...
        sprintf('$%.15g$', xticks(itx)), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickLabelFont);
    end

    % redefine y-tick labels using annotation()
    ax.YTickLabel = [];
    yticks = 0:0.2:1;
    ax.YTick = yticks;
    
    for ity = 1:numel(yticks)
        text(Kz_min-0.025, yticks(ity)+0.04, ...
        sprintf('$%.15g$', yticks(ity)), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickLabelFont);
    end

    xlim([Kz_min, Kz_max]);
    ylim([Kp_min, Kp_max]);

    %% Plot phase patches

    % define custom colors
    patch_purple = [.937, .910, .973];
    line_purple = [.455, .165, .765];
    patch_teal = [.627, .843, .843];
    line_teal = [.100, .600, .600];

    % plot phase patches
    X_OF = [Kz_min, 0, Kz_min];  Y_OF = [0, 0, Kp_max];
    X_FO = [0, Kz_max, Kz_max, Kz_min];  Y_FO = [0, 0, Kp_max, Kp_max];

    FaceAlpha = 1;
    patch(X_OF, Y_OF, patch_teal, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);
    patch(X_FO, Y_FO, patch_purple, 'EdgeColor', 'none', 'FaceAlpha', FaceAlpha);

    plot([0,Kz_max], [0,Kp_max], '--', 'Color', 'black', 'LineWidth', 2);


    %% Calculate and plot RG flows from poor man's scaling results
    cmap = turbo(Nflow);

    for itN = 1:Nflow

        if itN < 4
            [J,Kp,Kz,Ip,Iz,log_D] = RG_PoorMan_Aniso(J0(itN), Kp0(itN), Kz0(itN), Ip0(itN), Iz0(itN), 'order', 2, 'Temp', 1e-18);
            Nsteps = 75;
        else
            [J,Kp,Kz,Ip,Iz,log_D] = RG_PoorMan_Aniso(J0(itN), Kp0(itN), Kz0(itN), Ip0(itN), Iz0(itN), 'order', 2, 'Temp', 1e-12);
            
            idx = Kp < 1 & Kz < 1;
            Kp = Kp(idx);
            Kz = Kz(idx);

            Nsteps = 15;
        end

        Steps = round(linspace(10, numel(Kp)-1, Nsteps));

        for itS = 1:Nsteps
            NumD = Steps(itS);
            dKp = Kp(NumD+1) - Kp(NumD-1);
            dKz = Kz(NumD+1) - Kz(NumD-1);
            dS = sqrt(dKz^2 + dKp^2);

            %height = 0.1 / abs(log(dS));
            height = min(0.07, 30*dS);
            angle = atan(dKp/dKz) - pi/2;

            [X, Y] = getArrow(Kz(NumD), Kp(NumD), 0.045, 0.045, 0.015, height, angle);

            patch(X, Y, cmap(itN,:), 'EdgeColor', 'none', 'FaceAlpha', 1);
            %patch(X, Y, 'black', 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
        
    end

    %% Mark RG fixed points

    %{
    FI = [0, 0];

    plot(FI(1), FI(2), 's', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.37, 0.36, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FI}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    SO = [0, 1];

    plot(SO(1), SO(2), 'X', ...
    'MarkerSize', 15, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.36, 0.9, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{SO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    J_OO = [1, 0];

    plot(J_OO(1), J_OO(2), '^', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.78, 0.36, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{OO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    FO = [1, 1];

    plot(FO(1), FO(2), 'o', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.86, 0.9, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);
    %}
    
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
