
figureHandle = figure('Position', [50, 10, 550, 900]);

%% Subplot 1: Temperature-K0 phase diagram with J0 = 0.3, I0 = 0;

J0 = 0.3;
I0 = 0;

T_thres = -22;
Tmin = 1e-18;
Tmax = 1e-2;   

Tscale_labels = false;
phase_labels = true;

% Custom colors

patch_blue = [.910, .906, .969];
patch_purple = [.937, .910, .973];
patch_orange = [.984, .957, .910];
patch_green = [.925, .973, .910];
patch_gray = 0.85 * [1, 1, 1];

line_blue = [.188, .137, .761];
line_purple = [.455, .165, .765];
line_orange = [.902, .612, .247];
line_green = [.431, .780, .259];
line_gray = 0.5 * [1, 1, 1];

% Define temperature scales for each phase

[K0, Phase_range, Phase_name, Temps, Sent] = PhaseRange_TK(J0, I0, 'getEnt', 'showAllBound', 'InflecBound');

T_FL_1 = nan(1, numel(K0));
T_FO_1 = nan(1, numel(K0));
T_OO_1 = nan(1, numel(K0));
T_SO_1 = nan(1, numel(K0));
T_F_1 = nan(1, numel(K0));

T_FL_2 = nan(1, numel(K0));
T_FO_2 = nan(1, numel(K0));
T_OO_2 = nan(1, numel(K0));
T_SO_2 = nan(1, numel(K0));
T_F_2 = nan(1, numel(K0));

for itD = 1:numel(K0)

    for itP = 1:numel(Phase_name{itD})

        switch Phase_name{itD}{itP}

            case 'Fermi Liquid'
                T_FL_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_FL_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            case 'Fully Overscreened'
                T_FO_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_FO_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            case 'Orbital Overscreened'
                T_OO_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_OO_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            case 'Spin Overscreened'
                T_SO_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_SO_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            case 'Unscreened'
                T_F_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_F_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            otherwise

        end % switch-case
            
    end % it2
end % it1

K0min = min(K0);
K0max = max(K0);

subplot1 = subplot(2, 1, 1);
hold on;
pos1 = [0.15, 0.58, 0.8, 0.42]; % [left bottom width height]
set(subplot1, 'Position', pos1); 
set(gca, 'Layer', 'top');
ax = gca; 

set(ax, 'XScale', 'linear', 'YScale', 'log', 'FontSize', 18);
set(ax, 'XTick', -0.4:0.1:0.3);
set(ax, 'YTick', 10.^(log10(Tmin):3:log10(Tmax)));
set(ax, 'FontSize', 13);
set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
set(ax, 'XMinorTick', 'off', 'YMinorTick', 'on');
set(ax, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);
ax.TickLength = [0.015, 0.002];      % ticksize : [major, minor]

% define x- and y-labels
xlabel('$K_{0}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\mathrm{temperature / energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 20);

hx = get(ax, 'XLabel');
hy = get(ax, 'YLabel');
hx.Units = 'normalized';
hy.Units = 'normalized';

% Modify x- and y-label positions manually
hx.Position = hx.Position + [0, 0, 0];
hy.Position = hy.Position + [-0.09, -0.03, 0];

% redefine x-tick labels using annotation()
tickLabelFont = 15;
    
ax.XTickLabel = [];
xticks = K0min : 0.1 : K0max;
xticks = 0.1*round(10*xticks+1e-10);
    
xtickWidth = 0.09;
xtickHeight = 0.02;
for itx = 1:numel(xticks)
    X_pos = pos1(1) + pos1(3) * (xticks(itx) - K0min) / (K0max - K0min);
    Y_pos = pos1(2);
    X_pos = X_pos - 0.047;
    Y_pos = Y_pos - 0.02;
    annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', ['$',sprintf('%.15g',xticks(itx)),'$'], 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'center', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end

% redefine y-tick labels using annotation()
ax.YTickLabel = [];
yticks = 3*ceil(log10(Tmin)/3) : 3 : log10(Tmax);
yticks = power(10, yticks);
    
ytickWidth = 0.09;
ytickHeight = 0.03;
for ity = 1:numel(yticks)
    X_pos = pos1(1);
    Y_pos = pos1(2) + pos1(4) * (log10(yticks(ity)) - log10(Tmin)) / (log10(Tmax) - log10(Tmin));
    X_pos = X_pos - 0.08;
    Y_pos = Y_pos - 0.015;
    annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',SciNot(yticks(ity)),'$'], 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end

% define fully overscreened phase patch
idx1 = find(~isnan(T_FO_2), 1);

FO_x = [K0(idx1), K0max, K0(end:-1:idx1)];
FO_y = [Tmin, Tmin, T_FO_2(end:-1:idx1)];

% define spin overscreened phase patch
idx2 = find(~isnan(T_SO_2), 1, 'last');
SO_x = [K0min, K0(idx1:end), K0max, K0(idx2:-1:1)];
SO_y = [Tmin, T_FO_2(idx1:end), T_SO_2(idx2), T_SO_2(idx2:-1:1)];

% define unscreened phase patch
F_x = [K0(1:idx2), K0max, K0max, K0min];
F_y = [T_SO_2(1:idx2), T_SO_2(idx2), Tmax, Tmax];

% plot phase patches
FO_patch = patch(FO_x, FO_y, patch_purple, 'EdgeColor', 'none');
SO_patch = patch(SO_x, SO_y, patch_orange, 'EdgeColor', 'none');
F_patch = patch(F_x, F_y, patch_gray, 'EdgeColor', 'none');

% plot phase boundaries
MS = 10;
LW = 1;
LS = '-x';
plot(K0, T_FO_2, LS, 'color', line_purple, 'Markersize', MS, 'LineWidth', LW);
plot(K0, T_SO_2, LS, 'color', line_orange, 'Markersize', MS, 'LineWidth', LW);

% plot phase boundary markers
FO_text = '$T_{\mathrm{FO}}$';
SO_text = '$T_{\mathrm{SO}}$';
FS = 15;
LW = 1;
HL = 7;
HW = 7;

% Spin Overscreened scale
arrow_X = [0.5, 0.47];
arrow_Y = [0.93, 0.96];
text_X = -0.09;
text_Y = 1.5e-5;

annotation('arrow', arrow_X, arrow_Y, 'Color', line_orange, 'LineWidth', LW, 'HeadLength', HL, 'HeadWidth', HW);
text(text_X, text_Y, SO_text, 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);

% Fully Overscreened scale
arrow_X = [0.825, 0.79];
arrow_Y = [0.845, 0.87];
text_X = 0.195;
text_Y = 9e-9;

annotation('arrow', arrow_X, arrow_Y, 'Color', line_purple, 'LineWidth', LW, 'HeadLength', HL, 'HeadWidth', HW);
text(text_X, text_Y, FO_text, 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);

% insert phase labels
FS = 18;

% fully overscreened
text_X = 0.14;
text_Y = 1e-14;
text(text_X, text_Y, '$\mathrm{Fully}$', 'Interpreter', 'latex', 'Color', line_purple, 'FontName', 'Calibri', 'FontSize', FS);
text_X = 0.08;
text_Y = 1e-15;
text(text_X, text_Y, '$\mathrm{Overscreened}$', 'Interpreter', 'latex', 'Color', line_purple, 'FontName', 'Calibri', 'FontSize', FS);
    
% orbital overscreened
text_X = -0.25;
text_Y = 1e-7;
text(text_X, text_Y, '$\mathrm{Spin \ Overscreened}$', 'Interpreter', 'latex', 'Color', line_orange, 'FontName', 'Calibri', 'FontSize', FS);

% fit fully overscreened scale
set(gca, 'Layer', 'top');
xlim([K0min, K0max]);
ylim([Tmin, Tmax]);

% fit T - 1/K0
idx = ~isnan(T_FO_2);
X = 1./K0(K0>0 & idx);
Y = log(T_FO_2(K0>0 & idx));
f = polyfit(X, Y, 1);

yfit = polyval(f, Y);           % Estimated  Regression Line
SStot = sum((Y-mean(Y)).^2);    % Total Sum-Of-Squares
SSres = sum((Y-yfit).^2);       % Residual Sum-Of-Squares
Rsq = 1- SSres/SStot;           % R^2

% plot inset
inset_pos = [0.27 0.65 0.37 0.18];  % Adjust these values as needed
inset_axes = axes('Position', inset_pos);
hold on;
box on;     % Draw box around inset
set(inset_axes, 'XScale', 'linear', 'YScale', 'log', 'FontSize', 12);
inset_axes.TickLength = [0.015, 0.002];      % ticksize : [major, minor]
set(inset_axes, 'XMinorTick', 'off', 'YMinorTick', 'off');

plot(X, T_FO_2(K0>0 & idx), '-x', 'Color', line_purple, 'LineWidth', 1.5, 'MarkerSize', 5);
text2 = ['$\log (T) = - \frac{',sprintf('%.2g',-f(1)),'}{K_{0}}',sprintf('%.2g',f(2)),'$'];
text(4.5, 1e-3, text2,'Interpreter','latex','FontSize',13);

% define x- and y-labels
xlabel(inset_axes, '$1/K_{0}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(inset_axes, '$T_{\mathrm{FO}}$', 'Interpreter', 'latex', 'FontSize', 14);
    
hx = get(inset_axes, 'XLabel');
hy = get(inset_axes, 'YLabel');
hx.Units = 'normalized';
hy.Units = 'normalized';
    
% Modify x- and y-label positions manually
hx.Position = hx.Position + [0, 0.07, 0];
hy.Position = hy.Position + [0.05, 0.02, 0];
    
% redefine x-tick labels using annotation()
tickLabelFont = 12;
       
Xmin = 0;
Xmax = 20;
xlim([Xmin,Xmax]);
xticks = Xmin : 5 : Xmax;
inset_axes.XTick = xticks;
inset_axes.XTickLabel = [];
        
xtickWidth = 0.09;
xtickHeight = 0.06;
for itx = 1:numel(xticks)
    X_pos = inset_pos(1) + inset_pos(3) * (xticks(itx) - Xmin) / (Xmax - Xmin);
    Y_pos = inset_pos(2);
    X_pos = X_pos - 0.024;
    Y_pos = Y_pos - 0.055;
    annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', ['$',sprintf('%.15g',xticks(itx)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'left', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end
    
% redefine y-tick labels using annotation()
Ymin = 1e-20;
Ymax = 1;
ylim([Ymin,Ymax]);
yticks = log10(Ymin) : 5 : log10(Ymax);
yticks = power(10, yticks);
inset_axes.YTick = yticks;
inset_axes.YTickLabel = [];
        
ytickWidth = 0.09;
ytickHeight = 0.15;
for ity = 1:numel(yticks)
    X_pos = inset_pos(1);
    Y_pos = inset_pos(2) + inset_pos(4) * (log10(yticks(ity)) - log10(Ymin)) / (log10(Ymax) - log10(Ymin));
    X_pos = X_pos - 0.08;
    Y_pos = Y_pos - 0.135;
    annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',SciNot(yticks(ity)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end



%% Subplot 2: Temperature-I0 phase diagram with J0 = 0.3, K0 = 0.16;

J0 = 0.16;
K0 = 0.3;
SpOrbFlip = true;
InflecBound = true;

I0min = 1e-9;
I0max = 1e-2;
T_thres = -22;
Tmin = 1e-12;
Tmax = 1e-2;

Tscale_labels = true;
phase_labels = true;
plotEnt = false;

% Define temperatures scales for each phase
[I0, Phase_range, Phase_name, Temps, Sent] = PhaseRange_TI(J0, K0, 'I0min', I0min, 'showAllBound', 'flatThres', 0.2, 'alphaz', 1, 'getEnt');

idx = abs(2*log10(I0) - round(2*log10(I0))) < 1e-3 & I0 >= I0min & I0 <= I0max;
%idx = I0 >= I0min & I0 <=I0max;
I0 = I0(idx);
Phase_range = Phase_range(idx);
Phase_name = Phase_name(idx);
Sent = Sent(idx);

I0min = min(I0);

T_FL_1 = nan(1, numel(I0));
T_FO_1 = nan(1, numel(I0));
T_OO_1 = nan(1, numel(I0));
T_SO_1 = nan(1, numel(I0));
T_F_1 = nan(1, numel(I0));

T_FL_2 = nan(1, numel(I0));
T_FO_2 = nan(1, numel(I0));
T_OO_2 = nan(1, numel(I0));
T_SO_2 = nan(1, numel(I0));
T_F_2 = nan(1, numel(I0));

for itD = 1:numel(I0)

    for itP = 1:numel(Phase_name{itD})

        switch Phase_name{itD}{itP}

            case 'Fermi Liquid'
                T_FL_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_FL_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            case 'Fully Overscreened'
                T_FO_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_FO_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            case 'Orbital Overscreened'
                T_OO_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_OO_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            case 'Spin Overscreened'
                T_SO_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_SO_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            case 'Unscreened'
                T_F_2(itD) = power(10, Phase_range{itD}(itP,2) );

                if Phase_range{itD}(itP,1) > T_thres
                    T_F_1(itD) = power(10, Phase_range{itD}(itP,1) );
                end

            otherwise

        end % switch-case
            
    end % itP

end % itD

NFL_bound = nan(1, numel(I0));

% Find phase boundaries using inflection points
[~, Phase_range_inflec, Phase_name_inflec] = PhaseRange_TI(J0, K0, 'I0min', I0min, 'flatThres', 0.2, 'InflecBound', 'alphaz', 1);
for itD = 1:numel(I0)

    if InflecBound
        for itP = 1:numel(Phase_name_inflec{itD})
            
            switch Phase_name_inflec{itD}{itP}
    
                case 'Orbital Overscreened'
                    T_OO_2(itD) = power(10, Phase_range_inflec{itD}(itP,2) );

                case 'Fermi Liquid'
                    T_FL_2(itD) = power(10, Phase_range_inflec{itD}(itP,2) );
    
                otherwise
    
            end % switch-case
        end % itP
    end

    for itP = 1:numel(Phase_name_inflec{itD})

        switch Phase_name_inflec{itD}{itP}

            case 'Fully Overscreened'
                NFL_bound(itD) = power(10, Phase_range_inflec{itD}(itP,2));

            otherwise

        end % switch-case
            
    end % itP
end % itD

% Find phase boundaries using inflection points - using small alphaz for small \omega
[~, Phase_range_inflec, Phase_name_inflec] = PhaseRange_TI(J0, K0, 'I0min', I0min, 'flatThres', 0.2, 'InflecBound', 'alphaz', 0.8);

Thres = 1e-8;
for itD = 1:numel(I0)

    if InflecBound
        for itP = 1:numel(Phase_name_inflec{itD})
            
            switch Phase_name_inflec{itD}{itP}
    
                case 'Orbital Overscreened'
                    if T_OO_2(itD) < Thres
                        T_OO_2(itD) = power(10, Phase_range_inflec{itD}(itP,2) );
                    end

                case 'Fermi Liquid'
                    if T_FL_2(itD) < Thres
                        %T_FL_2(itD) = power(10, Phase_range_inflec{itD}(itP,2) );
                    end

                case 'Fully Overscreened'
                    if NFL_bound(itD) < Thres
                        NFL_bound(itD) = power(10, Phase_range_inflec{itD}(itP,2));
                    end
    
                otherwise
    
            end % switch-case
        end % itP
    end

end % itD

NFL_bound = NFL_bound(~isnan(NFL_bound));

subplot2 = subplot(2,1,2);
hold on;
pos2 = [0.15, 0.07, 0.8, 0.42];
set(subplot2, 'Position', pos2);
set(gca, 'Layer', 'top');
ax = gca;
set(ax, 'Position', pos2);
set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 20);
set(ax, 'XTick', 10.^(-8:2:-2));
set(ax, 'YTick', 10.^(-12:2:-2));
set(ax, 'FontSize', 20);
set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on');
ax.TickLength = [0.015, 0.002];      % ticksize : [major, minor]
set(gca, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

xlabel('$I_{0}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\mathrm{temperature / energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 20);

hx = get(ax, 'XLabel');
hy = get(ax, 'YLabel');
hx.Units = 'normalized';
hy.Units = 'normalized';

% Modify x- and y-label positions manually
hx.Position = hx.Position + [0, -0.04, 0];
hy.Position = hy.Position + [-0.09, -0.03, 0];

% redefine x-tick labels using annotation()
tickLabelFont = 15;
ax.XTickLabel = [];
xticks = 10.^(-2:-2:-8);
    
xtickWidth = 0.09;
xtickHeight = 0.02;
for itx = 1:numel(xticks)
    X_pos = pos2(1) + pos2(3) * (log10(xticks(itx)) - log10(I0min)) / (log10(I0max) - log10(I0min));
    Y_pos = pos2(2);
    X_pos = X_pos - 0.04;
    Y_pos = Y_pos - 0.02;
    annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', ['$',SciNot(xticks(itx)),'$'], 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'left', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end

% redefine y-tick labels using annotation()
ax.YTickLabel = [];
yticks = 10.^(-2:-2:-12);
    
ytickWidth = 0.09;
ytickHeight = 0.03;
for ity = 1:numel(yticks)
    X_pos = pos2(1);
    Y_pos = pos2(2) + pos2(4) * (log10(yticks(ity)) - log10(Tmin)) / (log10(Tmax) - log10(Tmin));
    X_pos = X_pos - 0.08;
    Y_pos = Y_pos - 0.015;
    annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',SciNot(yticks(ity)),'$'], 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end

% extrapolate T_FO
FO_Tmax = T_FO_2(find(I0 == I0min));
idx1 = find(T_FO_1 > Tmin, 1);
idx2 = find(~isnan(T_FO_1), 1, 'last');
idx3 = find(~isnan(T_OO_1), 1, 'last');

X = log10(I0(idx2+1:idx3));
Y = log10(T_OO_1(idx2+1:idx3));
f = polyfit(X, Y, 1);

I0_FO_extrap = ( log10(FO_Tmax) - f(2) ) / f(1);
I0_FO_extrap = power(10, I0_FO_extrap);
 
% extrapolate FO - OO phase boundary (extrapolate inflection point)
NFL_bound_mean = mean(log10(NFL_bound));
I0_NFL_bound_extrap = ( NFL_bound_mean - f(2) ) / f(1);
I0_NFL_bound_extrap = power(10, I0_NFL_bound_extrap);
NFL_bound_mean = power(10,NFL_bound_mean);

% extrapolate FL - OO - free phase boundaries
OO_T_extrap = T_OO_2(find(~isnan(T_OO_2), 1));
I0_FL_extrap = ( log10(OO_T_extrap) - f(2) ) / f(1);
I0_FL_extrap = power(10, I0_FL_extrap);

%%
%{}

% define orbital overscreened phase patch
idx4 = find(I0 > I0_FO_extrap, 1);
idx5 = find(~isnan(T_OO_2), 1);
idx6 = find(~isnan(T_OO_1), 1, 'last');
idx7 = find(~isnan(T_OO_2), 1, 'last');


if InflecBound
    % define fully overscreened phase patch
    FO_x = [I0min, I0(idx1), I0(idx1:idx2+1), I0_NFL_bound_extrap, I0(numel(NFL_bound):-1:1)];
    FO_y = [Tmin, Tmin, T_FO_1(idx1:idx2), T_OO_1(idx2+1), NFL_bound_mean, NFL_bound(end:-1:1)];

    % define orbital overscreened phase patch
    OO_x = [I0(1:numel(NFL_bound)), I0_NFL_bound_extrap, I0(idx4:idx6), I0_FL_extrap, I0(idx7:-1:idx5), I0min];
    OO_y = [NFL_bound, NFL_bound_mean, T_OO_1(idx4:idx6), OO_T_extrap, T_OO_2(idx7:-1:idx5), OO_T_extrap];
else
     % define fully overscreened phase patch
    FO_x = [I0min, I0(idx1), I0(idx1:idx2+1), I0_FO_extrap, I0min];
    FO_y = [Tmin, Tmin, T_FO_1(idx1:idx2), T_OO_1(idx2+1), FO_Tmax, FO_Tmax];

    % define orbital overscreened phase patch
    OO_x = [I0min, I0_FO_extrap, I0(idx4:idx6), I0_FL_extrap, I0(idx7:-1:idx5), I0min];
    OO_y = [FO_Tmax, FO_Tmax, T_OO_1(idx4:idx6), OO_T_extrap, T_OO_2(idx7:-1:idx5), OO_T_extrap];
end


% define Fermi liquid phase patch
FL_x = [I0(idx1), I0max, I0max, I0_FL_extrap, I0(idx6:-1:idx2+1), I0(idx2:-1:idx1)];
FL_y = [Tmin, Tmin, OO_T_extrap, OO_T_extrap, T_OO_1(idx6:-1:idx2+1), T_FO_1(idx2:-1:idx1)];

% define unscreened phase patch
F_x = [I0min, I0max, I0max, I0min];
F_y = [OO_T_extrap, OO_T_extrap, Tmax, Tmax];

% plot phase patches
FO_patch = patch(FO_x, FO_y, patch_purple, 'EdgeColor', 'none');
OO_patch = patch(OO_x, OO_y, patch_orange, 'EdgeColor', 'none');
FL_patch = patch(FL_x, FL_y, patch_green, 'EdgeColor', 'none');
F_patch = patch(F_x, F_y, patch_gray, 'EdgeColor', 'none');

% define extrapolated parts of phase boundaries
I0_UO_ext_bound = [I0min, I0(idx5), I0(idx7), I0_FL_extrap, I0max];
T_UO_ext_bound = repmat(OO_T_extrap, [1,5]);

I0_FL_ext_bound = [I0(idx6), I0_FL_extrap];
T_FL_ext_bound = [T_OO_1(idx6), OO_T_extrap];

% plot phase boundaries
MS = 5;
LW = 1;
LS = '-';

plot(I0, T_OO_2, LS, 'color', line_orange, 'Markersize', MS, 'LineWidth', LW);
plot(I0, T_FL_2, LS, 'color', line_green, 'Markersize', MS, 'LineWidth', LW);

plot(I0, T_OO_2, 'x', 'color', line_orange, 'Markersize', MS, 'LineWidth', LW);
plot(I0, T_FL_2, 'x', 'color', line_green, 'Markersize', MS, 'LineWidth', LW);

if InflecBound
    I0_FO_ext_bound = [I0(idx2), I0_NFL_bound_extrap];
    T_FO_ext_bound = [NFL_bound_mean, NFL_bound_mean];
    plot(I0(1:numel(NFL_bound)), NFL_bound, LS, 'color', line_purple, 'Markersize', MS, 'LineWidth', LW);
    plot(I0(1:numel(NFL_bound)), NFL_bound, 'x', 'color', line_purple, 'Markersize', MS, 'LineWidth', LW);
else
    I0_FO_ext_bound = [I0(idx2), I0_FO_extrap];
    T_FO_ext_bound = [FO_Tmax, FO_Tmax];
    plot(I0, T_FO_2, LS, 'color', line_purple, 'Markersize', MS, 'LineWidth', LW);
    plot(I0, T_FO_2, 'x', 'color', line_purple, 'Markersize', MS, 'LineWidth', LW);
end
    
% plot extrapolated parts of phase boundaries
plot(I0_UO_ext_bound(1:2), T_UO_ext_bound(1:2), '--', 'color', line_orange, 'Markersize', MS, 'LineWidth', 2);
plot(I0_UO_ext_bound(3:4), T_UO_ext_bound(3:4), '--', 'color', line_orange, 'Markersize', MS, 'LineWidth', 2);
plot(I0_FO_ext_bound, T_FO_ext_bound, '--', 'color', line_purple, 'Markersize', MS, 'LineWidth', 2);

% plot phase boundary markers
FO_text = '$T_{\mathrm{FO}}$';
SO_text = '$T_{\mathrm{SO}}$';
OO_text = '$T_{\mathrm{OO}}$';
FL_text = '$T_{\mathrm{FL}}$';
PHB_text = '$T_{\mathrm{PHB}}$';
LW = 1;
HL = 7;
HW = 7;

FS = 15;
% Fermi liquid scale
if J0 == 0.2 && ~InflecBound
    arrow_X = [0.74, 0.69];
    arrow_Y = [0.26, 0.27];
    text_X = 1.5e-4;
    text_Y = 3e-8;
elseif J0 == 0.16 && InflecBound
    arrow_X = [0.72, 0.67];
    arrow_Y = [0.28, 0.29];
    text_X = 1e-4;
    text_Y = 9e-8;
end
annotation('arrow', arrow_X, arrow_Y, 'Color', line_green, 'LineWidth', LW, 'HeadLength', HL, 'HeadWidth', HW);
text(text_X, text_Y, FL_text, 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);

% fully overscreened phase scale
if J0 == 0.2 && ~InflecBound
    arrow_X = [0.44, 0.39];
    arrow_Y = [0.34, 0.32];
    text_X = 3.5e-7;
    text_Y = 2.5e-6;
elseif J0 == 0.16 && InflecBound
    arrow_X = [0.45, 0.4];
    arrow_Y = [0.335, 0.315];
    text_X = 4.5e-7;
    text_Y = 2e-6;
end
annotation('arrow', arrow_X, arrow_Y, 'Color', line_purple, 'LineWidth', LW, 'HeadLength', HL, 'HeadWidth', HW);
text(text_X, text_Y, FO_text, 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
    
% orbital overscreened phase scale
if J0 == 0.2 && ~InlfecBound
    arrow_X = [0.45, 0.4];
    arrow_Y = [0.425, 0.405];
    text_X = 4.5e-7;
    text_Y = 2.5e-4;
elseif J0 == 0.16 && InflecBound
    arrow_X = [0.45, 0.4];
    arrow_Y = [0.465, 0.445];
    text_X = 4.5e-7;
    text_Y = 2.5e-3;
end
annotation('arrow', arrow_X, arrow_Y, 'Color', line_orange, 'LineWidth', LW, 'HeadLength', HL, 'HeadWidth', HW);
if SpOrbFlip
    text(text_X, text_Y, SO_text, 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
else
    text(text_X, text_Y, OO_text, 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
end

FS = 18;
            
% Fermi liquid
text_X = 1e-5;
text_Y = 1e-10;
text(text_X, text_Y, '$\mathrm{Fermi \ Liquid}$', 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);
        
% fully overscreened
text_X = 6e-9;
text_Y = 8e-8;
text(text_X, text_Y, '$\mathrm{Fully}$', 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
text_X = 1.5e-9;
text_Y = 2e-8;
text(text_X, text_Y, '$\mathrm{Overscreened}$', 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
    
% orbital overscreened
text_X = 1.5e-9;
text_Y = 1e-5;
if SpOrbFlip
    text(text_X, text_Y, '$\mathrm{Spin \ Overscreened}$', 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
else
    text(text_X, text_Y, '$\mathrm{Orbital \ Overscreened}$', 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
end

% linear fit to find power laws

X = I0(idx1:idx3);

Y = T_FL_2(idx1:idx3);
idx = X >= 3e-9 & X < 1e-6;
f = polyfit(log10(X(idx)), log10(Y(idx)), 1);
Xfit = linspace(-9+log10(3), -6, 10);
Yfit = polyval(f, Xfit);
Xfit = power(10, Xfit);     Yfit = power(10, Yfit)/2;
plot(Xfit, Yfit, '-', 'Color', line_green, 'LineWidth', 1);
textX = 3e-7;
textY = 3e-10;
text(textX, textY, ['$|I_{0}|^{',sprintf('%d',round(f(1))),'}$'], 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);

set(gca, 'Layer', 'top');
xlim([I0min, I0max]);
ylim([Tmin, Tmax]);

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

% save as pdf

fig = gcf;

set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'inches');
fig_pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', fig_pos(3:4));

% Specify full path
output_path = 'C:\Users\hsjun\OneDrive\Physics\Research\TsoK_publication\Figures\Phase_Diagrams_Iso.pdf';
set(fig, 'Renderer', 'painters');
print(fig, output_path, '-dpdf');

hold off;