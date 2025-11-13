function RG_plot(J0, K0, I0, order)

    Dmin = 5e-6; 
    D0 = 1;
    
    [J,K,I,log_D] = RG_PoorMan(J0,K0,I0,'order',order,'Rho',1);

    RGB = orderedcolors('gem');
    blue = RGB(1,:);
    orange = RGB(2,:);
    green = RGB(5,:);

    figure('Position', [100, 100, 800, 350]);   % left, bottom, width, height
    hold on;

    if I0 == 0
        yticks = 0:0.5:1;
        ymax = 1;
    else
        yticks = 0:4;
        ymax = 4;
    end

    set(gca, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 18);
    set(gca, 'XTick', 10.^(-5:1:0));
    set(gca, 'YTick', yticks);
    set(gca, 'FontSize', 24);
    set(gca, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(gca, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);
    xlim([Dmin, D0]);
    ylim([min([J,K,I])-0.1,ymax+0.4]);

    LW = 3;
    J_handle = plot(exp(log_D), J, 'Color', blue, 'LineWidth', LW);
    K_handle = plot(exp(log_D), K, 'Color', orange, 'LineWidth', LW);
    I_handle = plot(exp(log_D), I, 'Color', green, 'LineWidth', LW);

    legends = {'$J$', '$K$', '$I$'};
    handles = [J_handle, K_handle, I_handle];
    
    lgd = legend(handles, legends, 'Interpreter', 'latex', 'location', 'northeast', 'Box', 'off', 'FontSize', 37);
    lgd.Position = [0.5, 0.52, 0.1, 0.1];

    xlabel('$D/D_{0}$', 'Interpreter', 'latex', 'FontSize', 25);
    title(['$ J_{0}=', sprintf('%.15g',J0),',K_{0}=',sprintf('%.15g',K0),',I_{0}=',sprintf('%.15g',I0),'$'],'Interpreter','latex','FontSize',25);

    hold off;
end