function RG_plot_Aniso(J0, Kp0, Kz0, Ip0, Iz0)

    Dmin = 1e-9; 
    D0 = 1;
    
    [J,Kp,Kz,Ip,Iz,log_D] = RG_PoorMan_Aniso(J0,Kp0,Kz0,Ip0,Iz0);

    RGB = orderedcolors('gem');
    blue = RGB(1,:);
    orange = RGB(2,:);
    green = RGB(5,:);

    figure('Position', [100, 100, 800, 350]);   % left, bottom, width, height
    hold on;

    ymax = 4;
    yticks = 0:4;

    set(gca, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 18);
    set(gca, 'XTick', 10.^(-9:3:0));
    set(gca, 'YTick', yticks);
    set(gca, 'FontSize', 24);
    set(gca, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(gca, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);
    xlim([Dmin, D0]);
    ylim([min([Kz,Ip])-0.1,ymax+0.4]);

    LW = 3;
    J_handle = plot(exp(log_D), J, '-', 'Color', RGB(1,:), 'LineWidth', LW);
    Kp_handle = plot(exp(log_D), Kp, '-', 'Color', RGB(2,:), 'LineWidth', LW);
    Kz_handle = plot(exp(log_D), Kz, ':', 'Color', RGB(3,:), 'LineWidth', LW);
    Ip_handle = plot(exp(log_D), Ip, '-', 'Color', RGB(4,:), 'LineWidth', LW);
    Iz_handle = plot(exp(log_D), Iz, '--', 'Color', RGB(5,:), 'LineWidth', LW);

    legends = {'$J$', '$K_{\perp}$', '$K_{z}$', '$I_{\perp}$', '$I_{z}$'};
    handles = [J_handle, Kp_handle, Kz_handle, Ip_handle, Iz_handle];
    
    lgd = legend(handles, legends, 'Interpreter', 'latex', 'location', 'northeast', 'Box', 'off', 'FontSize', 24);
    lgd.Position = [0.5, 0.52, 0.1, 0.1];

    xlabel('$D/D_{0}$', 'Interpreter', 'latex', 'FontSize', 25);
    title(['$ J_{0}=', sprintf('%.15g',J0),',K_{\perp}=',sprintf('%.15g',Kp0),',K_{z}=',sprintf('%.15g',Kz0), ...
            ',I_{\perp}=',sprintf('%.15g',Ip0),',I_{z}=',sprintf('%.15g',Iz0),'$'],'Interpreter','latex','FontSize',25);

    hold off;
end