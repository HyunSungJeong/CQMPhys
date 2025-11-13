clear;

%% Extract crossover scale T_c from NRG data

path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK_TI_selected';

J0_NRG = 0.3*ones(1,8);
K0_NRG = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2];
I0_NRG = 1e-6*ones(1,8);
Tc_NRG = zeros(1,8);

figure;
hold on;
ylim([-2,2]);
handles = zeros(1,8);
leg = cell(1,8);

for itN = 1:numel(J0_NRG)
    folder = ['\J0=',sprintf('%.15g',J0_NRG(itN)),'_K0=',sprintf('%.15g',K0_NRG(itN)),'_I0=',sprintf('%.15g',I0_NRG(itN)),'_T=1e-24_Nkeep=3000_Lambda=2.5'];
    
    DiscData = load([path, folder, '\DiscData.mat']);
    field = fieldnames(DiscData);
    DiscData = getfield(DiscData, field{1});

    odisc = DiscData.odisc;
    sigmak = DiscData.sigmak;
    Adiscs = DiscData.Adiscs;
    nz = DiscData.nz;
    emin = DiscData.emin;
    T = 1e-24;

    Adisc = mean(cell2mat(reshape(Adiscs(2,:),[1 1 nz])),3);
            
    [ocont, ImpOrb] = getAcont(odisc,Adisc,sigmak,T/5,'alphaz',0.9,'emin',emin);
    
    T_grid = ocont(ocont > 0);
    OrbSusc = ImpOrb(ocont > 0);

    [logT, Orb_2Der] = log_Susc_2ndDer(T_grid, OrbSusc);
    [~, ~, MaxPos] = LocMax(logT, Orb_2Der);
    
    Idx = find(MaxPos > -10 & MaxPos < -2, 1, 'first');
    
    h = plot(logT, Orb_2Der);
    handles(itN) = h;
    plot(MaxPos(Idx), Orb_2Der(find(abs(logT - MaxPos(Idx)) < 1e-2 & logT < -2 & logT > -10, 1, 'first')), 'o', 'MarkerSize', 7);

    Tc_NRG(itN) = MaxPos(Idx);

    leg{itN} = ['$T=',SciNot(K0_NRG(itN)),'$'];
end

legend(handles, leg, 'Interpreter', 'latex');
hold off;

Tc_NRG = power(10, Tc_NRG);

%% compute the crossover temperature T_c^PM from poor man's scaling and plot figures

K_0 = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1];

J_0 = 0.3*ones(1,numel(K_0));

I_0 = 1e-6*ones(1,numel(K_0));

range = 3;

figureHandle = figure('Position', [50, 50, 550, 900]);

Tc = PM_RGflow_KI(J_0, K_0, I_0, range);
Tc = Tc(:, K_0 < 1e-1);
K_0 = K_0(K_0 < 1e-1);

subplot2 = subplot(2,1,2);
hold on;
pos2 = [0.15, 0.058, 0.8, 0.435];
set(subplot2, 'Position', pos2);
set(gca, 'Layer', 'top');
ax = gca;

Imin = 5e-6;
Imax = 1e-1;
Tmin = 1e-12;
Tmax = 1e-1;

set(ax, 'Layer', 'top');
set(ax, 'XScale', 'log', 'YScale', 'log');
set(ax, 'FontSize', 18);

set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
set(ax, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

xlabel('$K_{0}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\mathrm{temperature / energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 20);

hx = get(ax, 'XLabel');
hy = get(ax, 'YLabel');
hx.Units = 'normalized';
hy.Units = 'normalized';

% Modify x- and y-label positions manually
hx.Position = hx.Position + [0, 0.04, 0];
hy.Position = hy.Position + [0.02, 0, 0];

% redefine x-tick labels using annotation()
tickLabelFont = 18;
ax.XTickLabel = [];
xticks = logspace(-5,-1,3);
xtick_labels = {'$10^{-5}$', '$10^{-3}$', '$10^{-1}$'};
ax.XTick = xticks;
    
xtickWidth = 0.09;
xtickHeight = 0.06;
xtickShift = [-0.04, -0.04, -0.045];
for itx = 1:numel(xticks)
    X_pos = pos2(1) + pos2(3) * (log10(xticks(itx)) - log10(Imin)) / (log10(Imax) - log10(Imin)) + xtickShift(itx);
    Y_pos = pos2(2);
    Y_pos = Y_pos - 0.058;
    annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', xtick_labels{itx}, 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'left', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end

% redefine y-tick labels using annotation()
ax.YTickLabel = [];
yticks = logspace(-12, -3, 4);
ax.YTick = yticks;
    
ytickWidth = 0.09;
ytickHeight = 0.06;
for ity = 1:numel(yticks)
    X_pos = pos2(1);
    Y_pos = pos2(2) + pos2(4) * (log10(yticks(ity)) - log10(Tmin)) / (log10(Tmax) - log10(Tmin));
    X_pos = X_pos - 0.08;
    Y_pos = Y_pos - 0.04;
    annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',SciNot(yticks(ity)),'$'], 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end

xlim([Imin, Imax]);
ylim([Tmin, Tmax]);

h2 = plot(K_0, Tc(1,:), 'o-', 'Color', [0, .447, .741], 'LineWidth', 2);

idx = ~isnan(Tc(2,:)) & ~isnan(Tc(3,:));

X_patch = [K_0(idx), fliplr(K_0(idx))];
Y_patch = [Tc(2,idx), fliplr(Tc(3,idx))];

h3 = fill(X_patch, Y_patch, [0, .447, .741], 'LineStyle', 'none', 'FaceAlpha',0.3);

h1 = plot(K0_NRG, Tc_NRG, 'x-', 'LineWidth', 2, 'MarkerSize', 10);

set(gca, 'FontSize', 18);

legends = {'$T^{\mathrm{NRG}}_{\mathrm{c}}$', '$T^{\mathrm{PM}}_{\mathrm{c}} = T(K_{0} = I_{0})$', '$ T(1/3 \leq K_{0}/I_{0} \leq 3)$'};
legend([h1, h2, h3], legends, 'Interpreter', 'latex', 'Location', 'northeast');

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'inches');
fig_pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', fig_pos(3:4));

% subplot labels
annotation(figureHandle, 'textbox', [0.035, 1, 0, 0], ...
    'String', '$\mathrm{(a)}$', 'Interpreter', 'latex', 'FontSize', 18, ...
    'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'LineStyle', 'none');
    
annotation(figureHandle, 'textbox', [0.035, 0.5, 0, 0], ...
    'String', '$\mathrm{(b)}$', 'Interpreter', 'latex', 'FontSize', 18, ...
    'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'LineStyle', 'none');

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'inches');
fig_pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', fig_pos(3:4));

output_path = 'C:\Users\hsjun\OneDrive\Physics\Research\TsoK_publication\Figures\Tc_comparison.pdf';
set(fig, 'Renderer', 'painters');
print(fig, output_path, '-dpdf');

hold off;

function T_c = PM_RGflow_KI(J0, K0, I0, range)

    Nflow = numel(J0);

    Kmin = 5e-7;
    Kmax = 1.05;
    Imin = 5e-7;
    Imax = 4.05;

    subplot1 = subplot(2, 1, 1);
    hold on;
    pos1 = [0.15, 0.56, 0.8, 0.435]; % [left bottom width height]
    set(subplot1, 'Position', pos1); 
    set(gca, 'Layer', 'top');
    ax = gca;

    set(ax, 'Layer', 'top');
    set(ax, 'XScale', 'log', 'YScale', 'log');
    set(ax, 'FontSize', 18);
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
    hx.Position = hx.Position + [0, 0.055, 0];
    hy.Position = hy.Position + [0.025, 0, 0];

    % redefine x-tick labels using annotation()
    tickLabelFont = 18;
    ax.XTickLabel = [];
    xticks = [logspace(-6,0,4), 4];
    xtick_labels = {'$10^{-6}$', '$10^{-4}$', '$10^{-2}$', '$10^{0}$', '4'};
    ax.XTick = xticks;
    
    xtickWidth = 0.09;
    xtickHeight = 0.06;
    for itx = 1:numel(xticks)
        X_pos = pos1(1) + pos1(3) * (log10(xticks(itx)) - log10(Imin)) / (log10(Imax) - log10(Imin));
        Y_pos = pos1(2);
        if itx == 5
            X_pos = X_pos - 0.025;
        else
            X_pos = X_pos - 0.04;
        end
        Y_pos = Y_pos - 0.058;
        annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', xtick_labels{itx}, 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'left', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end

    % redefine y-tick labels using annotation()
    ax.YTickLabel = [];
    yticks = logspace(-6,0,4);
    ax.YTick = yticks;
    
    ytickWidth = 0.09;
    ytickHeight = 0.06;
    for ity = 1:numel(yticks)
        X_pos = pos1(1);
        Y_pos = pos1(2) + pos1(4) * (log10(yticks(ity)) - log10(Kmin)) / (log10(Kmax) - log10(Kmin));
        X_pos = X_pos - 0.08;
        Y_pos = Y_pos - 0.045;
        annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',SciNot(yticks(ity)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end

    xlim([Imin, Imax]);
    ylim([Kmin, Kmax]);

    %% Plot phase patches

    % define custom colors
    patch_blue = [.910, .906, .969];
    patch_purple = [.937, .910, .973];
    patch_orange = [.984, .957, .910];
    patch_rose = [0.957, 0.835, 0.875];
    patch_gray = 0.85 * [1, 1, 1];

    %% Calculate and plot RG flows from poor man's scaling results
    cmap = turbo(6);
    T_c = nan(3,Nflow);

    for itN = 1:Nflow
        [J,K,I,log_D] = RG_PoorMan(J0(itN), K0(itN), I0(itN), 'order', 3);

        Nsteps = 75;
        Steps = round(linspace(10, numel(log_D)-1, Nsteps));

        if mod(itN, 2) == 1
            for itS = 1:Nsteps
                NumD = Steps(itS);
                dI = log10(I(NumD+1)) - log10(I(NumD-1));
                dK = log10(K(NumD+1)) - log10(K(NumD-1));
                dS = sqrt(dI^2 + dK^2);
    
                %height = 0.1 / abs(log(dS));
                height = min(0.2, 50*dS);
                angle = atan(dK/dI) - pi/2;
    
                [X, Y] = getArrow(log10(I(NumD)), log10(K(NumD)), 0.15, 0.15, 0.05, height, angle);
                X = power(10, X);
                Y = power(10, Y);
    
                patch(X, Y, cmap(ceil(itN/2),:), 'EdgeColor', 'none', 'FaceAlpha', 1);
                %patch(X, Y, 'black', 'EdgeColor', 'none', 'FaceAlpha', 1);
            end
        end

        % plot reference lines
        I_ref = logspace(log10(Imin), log10(Imax), 10);

        X_patch = [I_ref, fliplr(I_ref)];
        Y_patch = [range*I_ref, fliplr(I_ref/range)];
        
        patch(X_patch, Y_patch, [.6 .6 .6], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        plot(I_ref, I_ref, 'Color', 'black', 'LineWidth', 1.5);

        % find the temperatures at which K0 == I0

        eps = 1e-2; % a small number

        ratio = [1,1/range,range];

        for itR = 1:3

            for it = 1:numel(log_D)
                if abs(K(it)/I(it) - ratio(itR)) < eps && isnan(T_c(itR,itN))

                    T_c(itR,itN) = power(10, log_D(it));
                    I_c = I(it);
                    K_c = K(it);
                end
            end % it
            %plot(I_c, K_c, 'o', 'MarkerSize', 7, 'LineWidth', 2);
        end % itR
        
    end % itN

    %% Mark RG fixed points
    SO = [Imin, Kmin];

    plot(SO(1), SO(2), 'X', ...
    'MarkerSize', 15, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.06, 0.47, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{SO}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

    FL = [4, 1];

    plot(FL(1), FL(2), 'diamond', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2);

    annotation('textbox', [0.91, 0.885, 0.1, 0.1], 'String', '$\mathbf{c}_{\mathrm{FL}}^{*}$', 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', 24);

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