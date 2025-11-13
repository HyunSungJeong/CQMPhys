clear;

Lor = @(x,width) (width ./ (x.^2+width^2) ./ pi);

figure;
hold on;
set(gca,'XScale','log','YScale','log');

X_lim = [1e-10,1];

xlim(X_lim);

N_bins = 100;
bins = logspace(log10(X_lim(1)), log10(X_lim(2)), N_bins);

N_points = 1000;

X_max = zeros(1,N_bins);
Y_max = zeros(1,N_bins);

colors = [0, .447, .741;
          .85, .325, .098;
          .929, .694, .125];
Exp = [0,-1,1];

handles = zeros(1,3);

for itX = 1:numel(Exp)

    for itN = 1:N_bins
        
        X = logspace(log10(bins(itN))-1, log10(bins(itN))+1, N_points);
        Y = Lor(X, bins(itN)/5) * power(bins(itN), Exp(itX));
    
        if itN == 1
            h = plot(X, Y, 'Color', colors(itX,:));
            handles(itX) = h;
        else
            plot(X, Y, 'Color', colors(itX,:));
        end

        [Y_max(itN), idx] = max(Y);
        X_max(itN) = X(idx);
    end % itN
end % Exp

legend(handles, {'$\mathrm{A(\omega) \propto \omega^{0}}$', '$\mathrm{A(\omega) \propto \omega^{-1}}$', '$\mathrm{A(\omega) \propto \omega^{1}}$'}, ...
                'Interpreter', 'latex', 'FontSize', 20);

%plot(X_max, Y_max, 'Color', [0, .447, .741]);

hold off;