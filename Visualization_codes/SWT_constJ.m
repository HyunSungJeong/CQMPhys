function SWT_constJ(U, J, Hyb)

    JL = linspace(-U,U,100);
    
    % Compute candidate bounds for mu
    expr1 = U - 2*J + JL/4;
    expr2 = U - (3/4)*JL;
    expr3 = U + JL/4;
        
    minExpr = min(min(expr1, expr2), expr3);
        
    % Define mu from condition in the problem
    mu = 0.5 * (minExpr - JL/4);
    gap = 0.5 * (minExpr + JL/4);

    %eps = sqrt(Hyb);
    %idx = find(abs(U-2*J+JL-mu) < eps);
    idx = 50;
        
    % Check inequality alpha:  -JL/4 < mu < minExpr
    cond1 = (-JL/4 < mu);
    cond2 = (mu < minExpr);
    valid = cond1 & cond2;
    JL = JL(valid);
    mu = mu(valid);
    gap = gap(valid);
        
    % Define I
    I0 = 1 ./ (U - 2*J + JL/4 - mu) + ...
    1 ./ (U - 3*JL/4 - mu) + ...
    2 ./ (JL/4 + mu);

    J0 = -1 ./ (U-2*J + JL/4 - mu) ...
        +2 ./ (U - 3*JL/4 - mu) ...
        +1 ./ (U + JL/4 - mu) ...
        +2 ./ (JL/4 + mu);
    J0 = J0 ./ 2;
    
    K_perp = 3 ./ (U-2*J + JL/4 - mu) ...
             -1 ./ (U + JL/4 - mu) ...
             +2 ./ (JL/4 + mu);
    K_perp = K_perp ./ 2;
    
    K_z = 3 ./ (U-2*J + JL/4 - mu) ...
          -2 ./ (U - 3*JL/4 - mu) ...
          +1 ./ (U + JL/4 - mu) ...
          +2 ./ (JL/4 + mu);
    K_z = K_z ./ 2;
        
    I0 = Hyb*I0;  J0 = Hyb*J0;    K_perp = Hyb*K_perp;   K_z = Hyb*K_z;

    linestyle = '-';
    colors = {[0 .447 .741];
              [.850 .325 .098];
              [.929 .694 .125];
              [.494 .184 .556]};
    range1 = 1:idx(1)-1;
    range2 = idx(end)+1:numel(JL);
    LW = 2;
    MS = 4;

    figure;
    hold on;
    set(gca, 'XScale', 'linear', 'YScale', 'linear');

    xlim([2*JL(1)/3-1, JL(end)/3+0.1]);
    ylim([1.1*min(K_perp), 1.1*max(K_perp)]);

    plot(repmat(JL(idx(ceil(end/2))), [1,2]), [-10,10], '--', 'color', [.7 .7 .7], 'LineWidth', 1.5); 

    h1 = plot(JL(range1), J0(range1), linestyle, 'Color', colors{1}, 'MarkerFaceColor', colors{1}, 'LineWidth', LW, 'MarkerSize', MS);
    plot(JL(range2), J0(range2), linestyle, 'Color', colors{1}, 'MarkerFaceColor', colors{1}, 'LineWidth', LW, 'MarkerSize', MS);

    h2 = plot(JL(range1), K_perp(range1), linestyle, 'Color', colors{2}, 'MarkerFaceColor', colors{2}, 'LineWidth', LW, 'MarkerSize', MS);
    plot(JL(range2), K_perp(range2), linestyle, 'Color', colors{2}, 'MarkerFaceColor', colors{2}, 'LineWidth', LW, 'MarkerSize', MS);
    
    h3 = plot(JL(range1), K_z(range1), linestyle, 'Color', colors{3}, 'MarkerFaceColor', colors{3}, 'LineWidth', LW, 'MarkerSize', MS);
    plot(JL(range2), K_z(range2), linestyle, 'Color', colors{3}, 'MarkerFaceColor', colors{3}, 'LineWidth', LW, 'MarkerSize', MS);

    h4 = plot(JL, I0, linestyle, 'Color', colors{4}, 'MarkerFaceColor', colors{4}, 'LineWidth', LW, 'MarkerSize', MS);

    legends = {'$J_{0}$', '$K_{\perp}$', '$K_{z}$', '$I_{0}$'};
    legend([h1, h2, h3, h4], legends, 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 18);

    hold off;
end
    