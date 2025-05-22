
% (J0, K_perp, I0) = (0.1, 0.1, 0)
% Data : (K_z, Exponent)
% Fit: 10^{-9} - 10^{-14}

clear;
K_perp = 0.1;

Data = [-0.5, -0.46;
        -0.4, -0.58;
        -0.35, -0.64;
        -0.3, -0.70;
        -0.25, -0.76;
        -0.2, -0.82;
        -0.15, -0.88;
        -0.1, -0.96;
        -0.05, -1.06
        ];

Data = Data(abs(Data(:,1)) > K_perp, :);

K_z_star = -sqrt(Data(:,1).^2 - K_perp^2);
SuscExp = Data(:,2);

coeff = polyfit(K_z_star, SuscExp, 1);
Xfit = linspace(-0.6, 0, 100);
Yfit = polyval(coeff, Xfit);
FitEq = ['$y = ', sprintf('%.2g', coeff(1)), 'x ', sprintf('%.2g', coeff(2)), '$'];

figure;
hold on;
set(gca,'FontSize',25);
xlabel('$K_z^{*}$','Interpreter','latex','FontSize',30);
ylabel('$\mathrm{Exponent}$','Interpreter','latex','FontSize',30);
title('$T^{z} \ \mathrm{exponent} - K_{z}$', 'Interpreter', 'latex', 'FontSize', 20);
plot(Xfit, Yfit, '--', 'Color', 'red');
plot(K_z_star, SuscExp, '.', 'Color', 'blue', 'MarkerSize', 15);
text(-0.2, -0.5, FitEq, 'Interpreter', 'latex', 'FontSize', 20);
hold off;