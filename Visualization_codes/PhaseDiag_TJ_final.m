function PhaseDiag_TJ_final(K0, I0, varargin)

    T_thres = -22;
    Tmin = 1e-18;
    Tmax = 1e-2;

    Tscale_labels = false;
    phase_labels = true;

    %% Custom colors

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

    %% Define temperatures scales for each phase

    [J0, Phase_range, Phase_name, Temps, Sent] = PhaseRange_TJ(K0, I0, 'getEnt', 'showAllBound', 'InflecBound');
    %[J0, Phase_range, Phase_name, Temps, Sent] = PhaseRange_TJ(K0, I0, 'getEnt', 'showAllBound');

    T_FL_1 = nan(1, numel(J0));
    T_FO_1 = nan(1, numel(J0));
    T_OO_1 = nan(1, numel(J0));
    T_SO_1 = nan(1, numel(J0));
    T_F_1 = nan(1, numel(J0));

    T_FL_2 = nan(1, numel(J0));
    T_FO_2 = nan(1, numel(J0));
    T_OO_2 = nan(1, numel(J0));
    T_SO_2 = nan(1, numel(J0));
    T_F_2 = nan(1, numel(J0));

    for itD = 1:numel(J0)

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

    %% Plot phase diagram with entropy colormap
    figure;
    hold on;
    set(gca, 'XScale', 'linear', 'YScale', 'log', 'FontSize', 20);
    set(gca, 'XTick', -0.4:0.1:0.3);
    set(gca, 'YTick', 10.^(log10(Tmin):3:log10(Tmax)));
    set(gca, 'FontSize', 13);
    set(gca, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(gca, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

    xlabel('$\mathrm{J_{0}}$', 'Interpreter', 'latex', 'FontSize', 25);
    ylabel('$\mathrm{energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 25);
            

    % define and plot colormap for entropy
    Cmap = zeros(numel(Temps), numel(J0));

    for it = 1:numel(J0)
        Cmap(:,it) = exp(Sent{it});
    end

    J0min = min(J0);
    J0max = max(J0);
    xlim([J0min, J0max]);
    ylim([Tmin, Tmax]);

    [arrow_X, arrow_Y] = meshgrid(J0, Temps);
    surf(arrow_X, arrow_Y, Cmap, 'EdgeColor', 'none');
    view(2);      
    colormap('jet');
    cb = colorbar;
    ylabel(cb, '$\exp(S_{\mathrm{imp}})$', 'Interpreter', 'latex'); 
    clim([0,4]);            % Set color value limits

    Z = 5*ones(1,numel(J0));
    % plot phase temperature scales
    plot3(J0, T_FO_1, Z, 'o-', 'color', 'red', 'Markersize', 5, 'LineWidth', 2);
    plot3(J0, T_SO_1, Z, 'o-', 'color', 'green', 'Markersize', 5, 'LineWidth', 2);
    plot3(J0, T_OO_1, Z, 'o-', 'color', 'blue', 'Markersize', 5, 'LineWidth', 2);
    plot3(J0, T_FL_1, Z, 'o-', 'color', patch_blue, 'Markersize', 5, 'LineWidth', 2);
    plot3(J0, T_F_1, Z, 'o-', 'color', [.7,.7,.7], 'Markersize', 5, 'LineWidth', 2);

    plot3(J0, T_FO_2, Z, 'o-', 'color', 'red', 'Markersize', 5, 'LineWidth', 2);
    plot3(J0, T_SO_2, Z, 'o-', 'color', 'green', 'Markersize', 5, 'LineWidth', 2);
    plot3(J0, T_OO_2, Z, 'o-', 'color', 'blue', 'Markersize', 5, 'LineWidth', 2);
    plot3(J0, T_FL_2, Z, 'o-', 'color', patch_blue, 'Markersize', 5, 'LineWidth', 2);
    plot3(J0, T_F_2, Z, 'o-', 'color', [.7,.7,.7], 'Markersize', 5, 'LineWidth', 2);

    plot3(J0, sqrt(T_FL_2 .* T_OO_1), Z, 'o-', 'color', 'black', 'Markersize', 5, 'LineWidth', 2);

    set(gca, 'Layer', 'top');
    hold off;


    %% Plot phase diagram without entropy colormap

    figure;
    hold on;
    set(gca, 'XScale', 'linear', 'YScale', 'log', 'FontSize', 18);
    set(gca, 'XTick', -0.4:0.1:0.3);
    set(gca, 'YTick', 10.^(log10(Tmin):3:log10(Tmax)));
    set(gca, 'FontSize', 13);
    set(gca, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on');
    set(gca, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

    xlabel('$J_{0}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$\mathrm{Temperature / energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 20);

    % define fully overscreened phase patch
    idx1 = find(~isnan(T_FO_2), 1);

    FO_x = [J0(idx1), J0max, J0(end:-1:idx1)];
    FO_y = [Tmin, Tmin, T_FO_2(end:-1:idx1)];

    % define orbital overscreened phase patch
    idx2 = find(~isnan(T_OO_2), 1, 'last');
    OO_x = [J0min, J0(idx1:end), J0max, J0(idx2:-1:1)];
    OO_y = [Tmin, T_FO_2(idx1:end), T_OO_2(idx2), T_OO_2(idx2:-1:1)];

    % define unscreened phase patch
    F_x = [J0(1:idx2), J0max, J0max, J0min];
    F_y = [T_OO_2(1:idx2), T_OO_2(idx2), Tmax, Tmax];

    % plot phase patches
    FO_patch = patch(FO_x, FO_y, patch_purple, 'EdgeColor', 'none');
    OO_patch = patch(OO_x, OO_y, patch_orange, 'EdgeColor', 'none');
    F_patch = patch(F_x, F_y, patch_gray, 'EdgeColor', 'none');

    % plot phase boundaries
    MS = 10;
    LW = 1;
    LS = '-x';
    plot(J0, T_FO_2, LS, 'color', line_purple, 'Markersize', MS, 'LineWidth', LW);
    plot(J0, T_OO_2, LS, 'color', line_orange, 'Markersize', MS, 'LineWidth', LW);

    % plot phase boundary markers
    FO_text = '$\mathrm{T_{FO}}$';
    OO_text = '$\mathrm{T_{OO}}$';
    FS = 13;
    LW = 1.3;

    if K0 == 0.3

        if Tscale_labels
            % fully overscreened phase scale
            arrow_X = [0.38, 0.32];
            arrow_Y = [0.69, 0.71];
            text_X = 6e-10;
            text_Y = 1e-8;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_purple, 'LineWidth', LW);
            text(text_X, text_Y, FO_text, 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
    
            % orbital overscreened phase scale
            arrow_X = [0.41, 0.35];
            arrow_Y = [0.865, 0.84];
            text_X = 1e-9;
            text_Y = 5e-5;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_orange, 'LineWidth', LW);
            text(text_X, text_Y, OO_text, 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
        end

        if phase_labels

            FS = 15;
            % fully overscreened
            text_X = 0.16;
            text_Y = 1e-14;
            text(text_X, text_Y, '$\mathrm{Fully}$', 'Interpreter', 'latex', 'Color', line_purple, 'FontName', 'Calibri', 'FontSize', FS);
            text_X = 0.1;
            text_Y = 1e-15;
            text(text_X, text_Y, '$\mathrm{Overscreened}$', 'Interpreter', 'latex', 'Color', line_purple, 'FontName', 'Calibri', 'FontSize', FS);
    
            % orbital overscreened
            text_X = -0.3;
            text_Y = 1e-7;
            text(text_X, text_Y, '$\mathrm{Orbital \ Overscreened}$', 'Interpreter', 'latex', 'Color', line_orange, 'FontName', 'Calibri', 'FontSize', FS);
        end

        set(gca, 'Layer', 'top');
        xlim([J0min, J0max]);
        ylim([Tmin, Tmax]);

        % fit T - 1/J0
        idx = ~isnan(T_FO_2);
        X = 1./J0(J0>0 & idx);
        Y = log(T_FO_2(J0>0 & idx));
        f = polyfit(X, Y, 1);

        yfit = polyval(f, Y);           % Estimated  Regression Line
        SStot = sum((Y-mean(Y)).^2);    % Total Sum-Of-Squares
        SSres = sum((Y-yfit).^2);       % Residual Sum-Of-Squares
        Rsq = 1- SSres/SStot;           % R^2

        % plot inset
        inset_pos = [0.27 0.25 0.32 0.32];  % Adjust these values as needed
        inset_axes = axes('Position', inset_pos);
        box on  % Draw box around inset
        
        plot(X, T_FO_2(J0>0 & idx), '-x', 'Color', line_purple, 'LineWidth', 1.5, 'MarkerSize', 5);
        text2 = ['$\log (T) = - \frac{',sprintf('%.2g',-f(1)),'}{J_{0}}',sprintf('%.2g',f(2)),'$'];
        text(5.5, 5e-4, text2,'Interpreter','latex','FontSize',10.5);
        xlabel(inset_axes, '$\mathbf{1/J_{0}}$', 'Interpreter', 'latex', 'FontSize', 10);
        ylabel(inset_axes, '$\mathbf{T_{FO}}$', 'Interpreter', 'latex', 'FontSize', 10);
        set(inset_axes, 'XScale', 'linear');
        set(inset_axes, 'YScale', 'log');
        

    elseif K0 == -0.3
    end

    

    hold off;

end