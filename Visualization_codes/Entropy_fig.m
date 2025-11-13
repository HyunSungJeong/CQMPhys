%% Load data
EntData = cell(1,3);
tmp = load('C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK\J0=0.3_K0=0.08_I0=0_T=1e-24_Nkeep=3000_Lambda=2.5\EntData.mat');
EntData{1} = tmp.EntData;
tmp = load('C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK\J0=0.3_K0=-0.08_I0=0_T=1e-24_Nkeep=3000_Lambda=2.5\EntData.mat');
EntData{2} = tmp.EntData;
tmp = load('C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK\J0=0.3_K0=0.08_I0=1e-06_T=1e-24_Nkeep=3000_Lambda=2.5\EntData.mat');
EntData{3} = tmp.EntData;

Tmin = 1e-20;
Tmax = 1;
Smin = -0.03;
Smax = log(4)+0.03;

%% Plot entropy figure
figure;
hold on;
set(gca, 'Layer', 'top');
set(gca,'XScale','log','YScale','linear');
set(gcf, 'Position', [100, 100, 550, 400]);

% define plot positions
ax = gca;
pos = ax.Position;  % [left bottom width height]
pos = pos + [0.06, 0.04, 0, 0];
ax.Position = pos;

set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
set(ax, 'XMinorTick', 'off', 'YMinorTick', 'off');
set(ax, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

xlabel('$T$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$S_{\mathrm{imp}}$', 'Interpreter', 'latex', 'FontSize', 20);

hx = get(ax, 'XLabel');
hy = get(ax, 'YLabel');
hx.Units = 'normalized';
hy.Units = 'normalized';

% Modify x- and y-label positions manually
hx.Position = hx.Position + [0, 0, 0];
hy.Position = hy.Position + [-0.08, 0.02, 0];

% redefine x-tick labels using annotation()
tickLabelFont = 15;
ax.XTickLabel = [];
xtick_pos = 10.^(-20:5:0);

xticks(xtick_pos);
    
xtickWidth = 0.09;
xtickHeight = 0.02;
for itx = 1:numel(xtick_pos)
    X_pos = pos(1) + pos(3) * (log10(xtick_pos(itx)) - log10(Tmin)) / (log10(Tmax) - log10(Tmin));
    Y_pos = pos(2);
    X_pos = X_pos - 0.04;
    Y_pos = Y_pos - 0.02;
    annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', ['$',SciNot(xtick_pos(itx)),'$'], 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'left', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
end

% redefine y-tick labels using annotation()
ax.YTickLabel = [];
ytick_pos = [0, log(2), 3/2*log(2), 2*log(2)];
ytick_labels = {'$0$', '$\ln 2$', '$\ln(2 \sqrt{2})$', '$\ln 4$'};

yticks(ytick_pos);
    
ytickWidth = 0.09;
ytickHeight = 0.03;
for ity = 1:numel(ytick_pos)
    X_pos = pos(1);
    Y_pos = pos(2) + pos(4) * (ytick_pos(ity) - Smin) / (Smax - Smin);
    X_pos = X_pos - 0.085;
    Y_pos = Y_pos + 0.01;
    annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ytick_labels{ity}, 'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    plot([Tmin,Tmax], [ytick_pos(ity), ytick_pos(ity)], '--', 'Color', [.7 .7 .7], 'LineWidth', 1);
end

xlim([Tmin, Tmax]);
ylim([Smin, Smax]);

beta = 1.3;
colors = [0, .447, .741;
          .85, .325, .098;
          .929, .694, .125];
handles = zeros(1,numel(EntData));

for itN = 1:numel(EntData)

    idx = find(abs(EntData{itN}.beta - beta) < 1e-3);
    Temps = EntData{itN}.Temps{idx};
    S_imp = EntData{itN}.S_imp{idx};

    S_imp = S_imp(Temps >= Tmin & Temps <= Tmax);
    Temps = Temps(Temps >= Tmin & Temps <= Tmax);

     h = plot(Temps, S_imp, 'Color', colors(itN,:), 'LineWidth', 2);
     handles(itN) = h;
end

text(1e-17, 0.05, '$\mathrm{FL}$', 'Interpreter', 'latex', 'FontSize', 15);
text(1e-18, 0.75, '$\mathrm{FO}$', 'Interpreter', 'latex', 'FontSize', 15);
text(1e-15, 1.1, '$\mathrm{SO}$', 'Interpreter', 'latex', 'FontSize', 15);

X_J0 = 5e-8;
X_K0 = 1e-5;
X_I0 = 2e-3;

Y_FO = 0.4;
Y_SO = 0.26;
Y_FL = 0.12;

text(X_J0*1.5, 0.52, '$J_{0}$', 'Interpreter', 'latex', 'FontSize', 15);
text(X_K0*1.5, 0.52, '$K_{0}$', 'Interpreter', 'latex', 'FontSize', 15);
text(X_I0*2, 0.52, '$I_{0}$', 'Interpreter', 'latex', 'FontSize', 15);

text(X_J0, Y_FO, '$0.3$', 'Interpreter', 'latex', 'FontSize', 15);
text(X_K0, Y_FO, '$0.08$', 'Interpreter', 'latex', 'FontSize', 15);
text(X_I0*2, Y_FO, '$0$', 'Interpreter', 'latex', 'FontSize', 15);

text(X_J0, Y_SO, '$0.3$', 'Interpreter', 'latex', 'FontSize', 15);
text(X_K0/4, Y_SO, '$-0.08$', 'Interpreter', 'latex', 'FontSize', 15);
text(X_I0*2, Y_SO, '$0$', 'Interpreter', 'latex', 'FontSize', 15);

text(X_J0, Y_FL, '$0.3$', 'Interpreter', 'latex', 'FontSize', 15);
text(X_K0, Y_FL, '$0.08$', 'Interpreter', 'latex', 'FontSize', 15);
text(X_I0, Y_FL, '$10^{-6}$', 'Interpreter', 'latex', 'FontSize', 15);

lgd = legend(handles, {'','',''}, 'Location', 'southeast', 'FontSize', 18, 'Box', 'off');
lgd.Position = [0.55, 0.22, 0.2, 0.2];


hold off;