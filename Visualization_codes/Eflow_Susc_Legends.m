%% Isotropic case

QnumKeys = [
        0, 1, 1;
        -1, 0, 0;
        1, 0, 0
        -1, 2, 0;
        1, 2, 0;
        1, 0, 2
        -1, 0, 2;
        2, 1, 1;
        -2, 1, 1
        ];

ColorDict = [
    .486, .282, .125;
    .466, .674, .188;
    1, 0, 0;
    0, .447, .741;
    .929, .694, .125;
    1, 0, 1;
    0, .5, .3;
    .475, .435, .655;
    .745, .792, .259
    ];

LinestyleDict = {'-', '-', '--', '-', '--', '--', '-', '--', '-'};

% define 'HalfInt' function
function Str = HalfInt(X)

    if rem(X*2,1) ~= 0
        error('ERR: Input of ''HalfInt'' must be a half integers');
    end

    if rem(X*2, 2) == 0
        Str = sprintf('%d',X);
    else
        Str = ['\frac{', sprintf('%d',2*X), '}{2}'];
    end
end

% Define signed_int_str()

function s = signed_int_str(x)
    if x == 0
        s = '0';
    else
        s = sprintf('%+d', x);
    end
end


% Plot the legend for Eflow plots
figure;
hold on;

legends = cell(1,size(QnumKeys,1));

for itN = 1:size(QnumKeys,1)
    plot(nan, nan, LinestyleDict{itN}, 'Color', ColorDict(itN,:), 'LineWidth', 2);
    Q_Sp = HalfInt(QnumKeys(itN,2)/2);
    Q_orb = HalfInt(QnumKeys(itN,2)/2);
    legends{itN} = ['$\left(', signed_int_str(QnumKeys(itN,1)), ', ' Q_Sp, ', ', Q_orb, '\right)$'];
end

legend(legends, 'Interpreter', 'latex', 'NumColumns', 3, 'Location', 'best', 'FontSize', 15, 'Box', 'off');


% Plot the legend for dynamic susceptibilitiy plots
figure;
hold on;

plot(nan, nan, 'Color', [.466, .674, .188], 'LineWidth', 2);
plot(nan, nan, 'Color', [.850, .325, .098], 'LineWidth', 2);
plot(nan, nan, 'Color', [0, .447, .741], 'LineWidth', 2);

legend({'$\chi_{\mathrm{sp}}$', '$\chi_{\mathrm{orb}}$', '$\chi_{\mathrm{sp-orb}}$'}, ...
                'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 16, 'Box', 'on', 'Color', 'none');



%% Anisotropic case

QnumKeys = [
        0, -1, 1;
        -1, 0, 1;
        -1, -1, 0;
        0, -2, 0;
        1, -1, 0;
        -2, 0, 2;
        0, 0, 2;
        1, -1, 2;
        -1, -1, 2
        ];

    ColorDict = [
        .466, .674, .188;
        .635, .078, .184;
        0, .447, .741;
        1, 0, 0;
        .929, .694, .125;
        1, 0, 1;
        0, .5, .3;
        .475, .435, .655;
        .745, .792, .259
        ];

LinestyleDict = {'-', '-', '--', '-', '--', '--', '-', '--', '-'};

% Plot the legend for Eflow plots
figure;
hold on;

legends = cell(1,size(QnumKeys,1));

for itN = 1:size(QnumKeys,1)
    plot(nan, nan, LinestyleDict{itN}, 'Color', ColorDict(itN,:), 'LineWidth', 2);
    Q_Sp = HalfInt(QnumKeys(itN,3)/2);
    legends{itN} = ['$\left(', signed_int_str(QnumKeys(itN,1)), ', ' signed_int_str(QnumKeys(itN,2)), ', ' Q_Sp, '\right)$'];
end

legend(legends, 'Interpreter', 'latex', 'NumColumns', 3, 'Location', 'best', 'FontSize', 15, 'Box', 'off');


% Plot the legend for dynamic susceptibilitiy plots
figure;
hold on;

plot(nan, nan, 'Color', [.466, .674, .188], 'LineWidth', 2);
plot(nan, nan, 'Color', [.850, .325, .098], 'LineWidth', 2);
plot(nan, nan, 'Color', [0, .447, .741], 'LineWidth', 2);

legend({'$\chi_{\mathrm{sp}}$', '$\chi^{+}_{\mathrm{orb}}$', '$\chi^{z}_{\mathrm{orb}}$'}, ...
                'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 16, 'Box', 'on', 'Color', 'none');