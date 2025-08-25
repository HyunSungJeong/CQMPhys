function PhaseDiag_TK_final(J0, I0, varargin)

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

    %% Define temperature scales for each phase

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

    xlabel('$\mathrm{K_{0}}$', 'Interpreter', 'latex', 'FontSize', 25);
    ylabel('$\mathrm{energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 25);

    % define and plot colormap for entropy
    Cmap = zeros(numel(Temps), numel(K0));

    for it = 1:numel(K0)
        Cmap(:,it) = exp(Sent{it});
    end

    K0min = min(K0);
    K0max = max(K0);
    xlim([K0min, K0max]);
    ylim([Tmin, Tmax]);

    [arrow_X, arrow_Y] = meshgrid(K0, Temps);
    surf(arrow_X, arrow_Y, Cmap, 'EdgeColor', 'none');
    view(2);      
    colormap('jet');
    cb = colorbar;
    ylabel(cb, '$\exp(S_{\mathrm{imp}})$', 'Interpreter', 'latex'); 
    clim([0,4]);            % Set color value limits

    Z = 5*ones(1,numel(K0));
    % plot phase temperature scales
    plot3(K0, T_FO_1, Z, 'o-', 'color', 'red', 'Markersize', 5, 'LineWidth', 2);
    plot3(K0, T_SO_1, Z, 'o-', 'color', 'green', 'Markersize', 5, 'LineWidth', 2);
    plot3(K0, T_OO_1, Z, 'o-', 'color', 'blue', 'Markersize', 5, 'LineWidth', 2);
    plot3(K0, T_FL_1, Z, 'o-', 'color', patch_blue, 'Markersize', 5, 'LineWidth', 2);
    plot3(K0, T_F_1, Z, 'o-', 'color', [.7,.7,.7], 'Markersize', 5, 'LineWidth', 2);

    plot3(K0, T_FO_2, Z, 'o-', 'color', 'red', 'Markersize', 5, 'LineWidth', 2);
    plot3(K0, T_SO_2, Z, 'o-', 'color', 'green', 'Markersize', 5, 'LineWidth', 2);
    plot3(K0, T_OO_2, Z, 'o-', 'color', 'blue', 'Markersize', 5, 'LineWidth', 2);
    plot3(K0, T_FL_2, Z, 'o-', 'color', patch_blue, 'Markersize', 5, 'LineWidth', 2);
    plot3(K0, T_F_2, Z, 'o-', 'color', [.7,.7,.7], 'Markersize', 5, 'LineWidth', 2);

    plot3(K0, sqrt(T_FL_2 .* T_OO_1), Z, 'o-', 'color', 'black', 'Markersize', 5, 'LineWidth', 2);

    set(gca, 'Layer', 'top');
    hold off;


    %% Plot phase diagram without entropy colormap

    figure('Position', [50, 200, 550, 400]);
    hold on;
    ax = gca;
    pos = [0.15, 0.15, 0.775, 0.815];
    set(ax, 'Position', pos);
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
    ylabel('$\mathrm{Temperature / energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 20);

    hx = get(ax, 'XLabel');
    hy = get(ax, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';

    % Modify x- and y-label positions manually
    hx.Position = hx.Position + [0, 0.02, 0];
    hy.Position = hy.Position + [-0.09, 0, 0];

    % redefine x-tick labels using annotation()
    tickLabelFont = 15;
    
    ax.XTickLabel = [];
    xticks = K0min : 0.1 : K0max;
    xticks = 0.1*round(10*xticks);
    
    xtickWidth = 0.09;
    xtickHeight = 0.06;
    for itx = 1:numel(xticks)
        X_pos = pos(1) + pos(3) * (xticks(itx) - K0min) / (K0max - K0min);
        Y_pos = pos(2);
        X_pos = X_pos - 0.05;
        Y_pos = Y_pos - 0.06;
        annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', ['$',sprintf('%.15g',xticks(itx)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'left', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end

    % redefine y-tick labels using annotation()
    ax.YTickLabel = [];
    yticks = 3*ceil(log10(Tmin)/3) : 3 : log10(Tmax);
    yticks = power(10, yticks);
    
    ytickWidth = 0.09;
    ytickHeight = 0.15;
    for ity = 1:numel(yticks)
        X_pos = pos(1);
        Y_pos = pos(2) + pos(4) * (log10(yticks(ity)) - log10(Tmin)) / (log10(Tmax) - log10(Tmin));
        X_pos = X_pos - 0.08;
        Y_pos = Y_pos - 0.11;
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
    FO_text = '$\mathrm{T_{FO}}$';
    SO_text = '$\mathrm{T_{SO}}$';
    FS = 13;
    LW = 1.3;

    if J0 == 0.3

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
            text(text_X, text_Y, SO_text, 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
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
            text_Y = 1e-6;
            text(text_X, text_Y, '$\mathrm{Spin \ Overscreened}$', 'Interpreter', 'latex', 'Color', line_orange, 'FontName', 'Calibri', 'FontSize', FS);
        end

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

        %% plot inset
        inset_pos = [0.27 0.32 0.37 0.37];  % Adjust these values as needed
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
        hy.Position = hy.Position + [0.05, 0, 0];
    
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
            X_pos = X_pos - 0.022;
            Y_pos = Y_pos - 0.05;
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
            Y_pos = Y_pos - 0.12;
            annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',SciNot(yticks(ity)),'$'], 'Interpreter', 'latex', ...
                            'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
        end
        

    elseif J0 == -0.3
    end

    

    hold off;

end