function PhaseDiag_TI_final(J0, K0, varargin)
    % <Description>
    % Draws T-I phase diagram from NRG data
    %
    % <Input>
    % J0 : [numeric] spin-spin coupling strength in the isotropic 2soK model
    % K0 : [numeric] pseudospin-pseudospin coupling strength in the isotropic 2soK model
    %
    % <Option>
    % 'SpOrbFlip' : If used, the spin and orbital sectors are flipped
    %
    % <Output>
    % Temperature - I0 phase diagram

    SpOrbFlip = false;

    while ~isempty(varargin)

        switch varargin{1}
            case 'SpOrbFlip'
                SpOrbFlip = true;
                varargin(1) = [];
            otherwise
        end
    end % while

    I0min = 1e-9;
    %I0min = 1e-13;
    I0max = 2e-2;
    T_thres = -22;
    Tmin = 1e-12;
    %Tmin = 1e-15;
    Tmax = 1e-2;

    Tscale_labels = true;
    phase_labels = true;
    plotEnt = true;

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

    [I0, Phase_range, Phase_name, Temps, Sent] = PhaseRange_TI(J0, K0, 'I0min', I0min, 'getEnt', 'showAllBound', 'flatThres', 0.2, 'alphaz', 0.7);
    %[I0, Phase_range, Phase_name] = PhaseRange_TI(J0, K0, 'I0min', I0min, 'showAllBound', 'flatThres', 0.2);

    idx = I0 >= I0min & I0 <= I0max;
    I0 = I0(idx);
    Phase_range = Phase_range(idx);
    Phase_name = Phase_name(idx);

    I0min = min(I0);
    %I0max = max(I0);

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
            
        end % it2
    end % it1


    %% Plot phase diagram with entropy colormap

    if plotEnt

        figure;
        hold on;
        set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 20);
        ax = gca;
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;
        xlabel('$\mathrm{I_{0}}$', 'Interpreter', 'latex', 'FontSize', 25);
        ylabel('$\mathrm{energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 25);
                
    
        % define and plot colormap for entropy
        Cmap = zeros(numel(Temps), numel(I0));
    
        Tsqrt2 = nan(1, numel(I0));
    
        for it = 1:numel(I0)
            Cmap(:,it) = exp(Sent{it});
            
            T1 = []; T2 = [];
            Idx1 = find(Sent{it} < log(2)/2, 1, 'first');
            Idx2 = find(Sent{it} > log(2)/2, 1, 'last');
            T1 = Temps(Idx1);
            T2 = Temps(Idx2);
    
            if ~isempty(T1)
    
                d1 = abs(sqrt(2) - exp(Sent{it}(Idx1)));
                d2 = abs(sqrt(2) - exp(Sent{it}(Idx2)));
    
                % Internal division point weighted by distances
                Tsqrt2(it) = exp( (d2 * log(T1) + d1 * log(T2)) / (d1 + d2) );
            end
        end
    
        xlim([I0min, I0max]);
        ylim([Tmin, Tmax]);
    
        [arrow_X, arrow_Y] = meshgrid(I0, Temps);
        surf(arrow_X, arrow_Y, Cmap, 'EdgeColor', 'none');
        view(2);      
        colormap('jet');
        cb = colorbar;
        ylabel(cb, '$\exp(S_{\mathrm{imp}})$', 'Interpreter', 'latex'); 
        clim([0,4]);            % Set color value limits
    
        Z = 5*ones(1,numel(I0));
        % plot phase temperature scales
        plot3(I0, T_FO_1, Z, 'o-', 'color', 'red', 'Markersize', 5, 'LineWidth', 2);
        plot3(I0, T_SO_1, Z, 'o-', 'color', 'green', 'Markersize', 5, 'LineWidth', 2);
        plot3(I0, T_OO_1, Z, 'o-', 'color', 'blue', 'Markersize', 5, 'LineWidth', 2);
        plot3(I0, T_FL_1, Z, 'o-', 'color', patch_blue, 'Markersize', 5, 'LineWidth', 2);
        plot3(I0, T_F_1, Z, 'o-', 'color', [.7,.7,.7], 'Markersize', 5, 'LineWidth', 2);
    
        plot3(I0, T_FO_2, Z, 'o-', 'color', 'red', 'Markersize', 5, 'LineWidth', 2);
        plot3(I0, T_SO_2, Z, 'o-', 'color', 'green', 'Markersize', 5, 'LineWidth', 2);
        plot3(I0, T_OO_2, Z, 'o-', 'color', 'blue', 'Markersize', 5, 'LineWidth', 2);
        plot3(I0, T_FL_2, Z, 'o-', 'color', patch_blue, 'Markersize', 5, 'LineWidth', 2);
        plot3(I0, T_F_2, Z, 'o-', 'color', [.7,.7,.7], 'Markersize', 5, 'LineWidth', 2);
    
        Idx = ~isnan(Tsqrt2);
        plot3(I0(Idx), Tsqrt2(Idx), Z(Idx), '--', 'color', [.7,.7,.7], 'LineWidth', 1);
    
        %plot3(I0, sqrt(T_FL_2 .* T_OO_1), Z, 'o-', 'color', 'black', 'Markersize', 5, 'LineWidth', 2);
    
        set(gca, 'Layer', 'top');
        hold off;
    end


    %% Plot phase diagram without entropy colormap

    figure('Position', [400, 300, 550, 400]);
    hold on;
    ax = gca;
    pos = [0.15, 0.15, 0.775, 0.815];
    set(ax, 'Position', pos);
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 20);
    set(ax, 'XTick', 10.^(-8:2:-2));
    set(ax, 'YTick', 10.^(-12:2:-2));
    set(ax, 'FontSize', 20);
    set(ax, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on');
    ax.TickLength = [0.015, 0.002];      % ticksize : [major, minor]
    set(gca, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);
    %set(gcf, 'Position', [400, 300, 550, 400]);  % Width = 567 px, Height = 420 px

    ax = gca;

    xlabel('$I_{0}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$\mathrm{Temperature / energy \ scale}$', 'Interpreter', 'latex', 'FontSize', 20);

    hx = get(ax, 'XLabel');
    hy = get(ax, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';

    % Modify x- and y-label positions manually
    hx.Position = hx.Position + [0, -0.05, 0];
    hy.Position = hy.Position + [-0.09, 0, 0];

    pos = ax.Position;  % [left bottom width height]
    %{
    pos(1) = pos(1) + 0.08;
    pos(3) = pos(3) - 0.08;     % shrink width by 0.1 to make space on the left
    pos(2) = pos(2) + 0.05;
    pos(4) = pos(4) - 0.02;     % shrink height by 0.02 to make space on the top
    ax.Position = pos;
    %}

    % redefine x-tick labels using annotation()
    tickLabelFont = 15;
    ax.XTickLabel = [];
    xticks = 10.^(-2:-2:-8);
    
    xtickWidth = 0.09;
    xtickHeight = 0.06;
    for itx = 1:numel(xticks)
        X_pos = pos(1) + pos(3) * (log10(xticks(itx)) - log10(I0min)) / (log10(I0max) - log10(I0min));
        Y_pos = pos(2);
        X_pos = X_pos - 0.05;
        Y_pos = Y_pos - 0.06;
        annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', ['$',SciNot(xticks(itx)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'left', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end

    % redefine y-tick labels using annotation()
    ax.YTickLabel = [];
    yticks = 10.^(-2:-2:-12);
    
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

    % extrapolate fully overscreened -  orbital overscreened phase boundary
    FO_Tmax = T_FO_2(find(I0 == I0min));
    idx1 = find(T_FO_1 > Tmin, 1);
    idx2 = find(~isnan(T_FO_1), 1, 'last');
    idx3 = find(~isnan(T_OO_1), 1, 'last');

    X = log10(I0(idx2+1:idx3));
    Y = log10(T_OO_1(idx2+1:idx3));
    f = polyfit(X, Y, 1);
    I0_FO_extrap = ( log10(FO_Tmax) - f(2) ) / f(1);
    I0_FO_extrap = power(10, I0_FO_extrap);

    % extrapolate FL - OO - free phase boundaries
    OO_T_extrap = T_OO_2(find(~isnan(T_OO_2), 1));
    I0_FL_extrap = ( log10(OO_T_extrap) - f(2) ) / f(1);
    I0_FL_extrap = power(10, I0_FL_extrap);

    % define fully overscreened phase patch
    FO_x = [I0min, I0(idx1), I0(idx1:idx2+1), I0_FO_extrap, I0min];
    FO_y = [Tmin, Tmin, T_FO_1(idx1:idx2), T_OO_1(idx2+1), FO_Tmax, FO_Tmax];
    
    % define orbital overscreened phase patch
    idx4 = find(I0 > I0_FO_extrap, 1);
    idx5 = find(~isnan(T_OO_2), 1);
    idx6 = find(~isnan(T_OO_1), 1, 'last');
    idx7 = find(~isnan(T_OO_2), 1, 'last');

    OO_x = [I0min, I0_FO_extrap, I0(idx4:idx6), I0_FL_extrap, I0(idx7:-1:idx5), I0min];
    OO_y = [FO_Tmax, FO_Tmax, T_OO_1(idx4:idx6), OO_T_extrap, T_OO_2(idx7:-1:idx5), OO_T_extrap];

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
    I0_UO_ext_bound = [I0min, I0(idx5), I0(idx6), I0_FL_extrap, I0max];
    T_UO_ext_bound = repmat(OO_T_extrap, [1,5]);

    I0_FO_ext_bound = [I0(idx2+1), I0_FO_extrap];
    T_FO_ext_bound = [FO_Tmax, FO_Tmax];

    I0_FL_ext_bound = [I0(idx6), I0_FL_extrap];
    T_FL_ext_bound = [T_OO_1(idx6), OO_T_extrap];

    % plot phase boundaries
    MS = 5;
    LW = 1;
    LS = '-';
    %plot(I0, T_FO_1, LS, 'color', line_gray, 'Markersize', MS, 'LineWidth', LW);
    %plot(I0(idx2+1:end), T_OO_1(idx2+1:end), LS, 'color', line_gray, 'Markersize', MS, 'LineWidth', LW);

    plot(I0, T_FO_2, LS, 'color', line_purple, 'Markersize', MS, 'LineWidth', LW);
    plot(I0, T_OO_2, LS, 'color', line_orange, 'Markersize', MS, 'LineWidth', LW);
    plot(I0, T_FL_2, LS, 'color', line_green, 'Markersize', MS, 'LineWidth', LW);

    showIdx = ismember(I0, 10.^(-14:0.2:-2));
    plot(I0(showIdx), T_FO_2(showIdx), 'x', 'color', line_purple, 'Markersize', MS, 'LineWidth', LW);
    plot(I0(showIdx), T_OO_2(showIdx), 'x', 'color', line_orange, 'Markersize', MS, 'LineWidth', LW);
    plot(I0(showIdx), T_FL_2(showIdx), 'x', 'color', line_green, 'Markersize', MS, 'LineWidth', LW);
    
    % plot extrapolated parts of phase boundaries
    plot(I0_UO_ext_bound(1:2), T_UO_ext_bound(1:2), '--', 'color', line_orange, 'Markersize', MS, 'LineWidth', 2);
    plot(I0_UO_ext_bound(3:4), T_UO_ext_bound(3:4), '--', 'color', line_orange, 'Markersize', MS, 'LineWidth', 2);
    %plot(I0_UO_ext_bound(4:5), T_UO_ext_bound(4:5), '--', 'color', line_gray, 'Markersize', MS, 'LineWidth', 2);
    plot(I0_FO_ext_bound, T_FO_ext_bound, '--', 'color', line_purple, 'Markersize', MS, 'LineWidth', 2);
    %plot(I0_FL_ext_bound, T_FL_ext_bound, '--', 'color', line_gray, 'Markersize', MS, 'LineWidth', 2);

    % plot phase boundary markers
    FO_text = '$T_{\mathrm{FO}}$';
    SO_text = '$T_{\mathrm{SO}}$';
    OO_text = '$T_{\mathrm{OO}}$';
    FL_text = '$T_{\mathrm{FL}}$';
    PHB_text = '$T_{\mathrm{PHB}}$';
    LW = 1.3;

    if J0 == 0.2

        if Tscale_labels
            FS = 15;
            % Fermi liquid scale
            arrow_X = [0.73, 0.68];
            arrow_Y = [0.54, 0.55];
            text_X = 2.3e-4;
            text_Y = 4.5e-8;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_green, 'LineWidth', LW);
            text(text_X, text_Y, FL_text, 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);

            %{
            % particle-hole symmetry breaking scale
            arrow_X = [0.71, 0.66];
            arrow_Y = [0.665, 0.65];
            text_X = 1.8e-4;
            text_Y = 4.5e-6;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_gray, 'LineWidth', LW);
            text(text_X, text_Y, PHB_text, 'Color', line_gray, 'Interpreter', 'latex', 'FontSize', FS);
            %}

            % fully overscreened phase scale
            arrow_X = [0.52, 0.47];
            arrow_Y = [0.66, 0.64];
            text_X = 1.5e-6;
            text_Y = 2e-6;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_purple, 'LineWidth', LW);
            text(text_X, text_Y, FO_text, 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
    
            % orbital overscreened phase scale
            arrow_X = [0.45, 0.4];
            arrow_Y = [0.815, 0.795];
            text_X = 3e-7;
            text_Y = 2.2e-4;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_orange, 'LineWidth', LW);
            if SpOrbFlip
                text(text_X, text_Y, SO_text, 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
            else
                text(text_X, text_Y, OO_text, 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
            end
        end

        if phase_labels
            FS = 15;
            
            % Fermi liquid
            text_X = 1e-5;
            text_Y = 1e-10;
            text(text_X, text_Y, '$\mathrm{Fermi \ Liquid}$', 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);
        
            % fully overscreened
            text_X = 6e-9;
            text_Y = 2.3e-7;
            text(text_X, text_Y, '$\mathrm{Fully}$', 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
            text_X = 1.5e-9;
            text_Y = 6e-8;
            text(text_X, text_Y, '$\mathrm{Overscreened}$', 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
    
            % orbital overscreened
            text_X = 1.5e-9;
            text_Y = 1e-5;
            if SpOrbFlip
                text(text_X, text_Y, '$\mathrm{Spin \ Overscreened}$', 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
            else
                text(text_X, text_Y, '$\mathrm{Orbital \ Overscreened}$', 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
            end
        end

        %% linear fit to find power laws

        X = I0(idx1:idx3);
        %Y = [T_FO_1(idx1:idx2), T_OO_1(idx2+1:idx3)];
    
        %{
        idx = X >= 1e-7 & X < 5e-5;
        f = polyfit(log10(X(idx)), log10(Y(idx)), 1);
        Xfit = linspace(-7, -5+log10(5), 10);
        Yfit = polyval(f, Xfit);
        Xfit = power(10, Xfit);     Yfit = power(10, Yfit)/2;
        plot(Xfit, Yfit, '-', 'Color', line_gray, 'LineWidth', 1);
        textX = 3e-6;
        textY = 8e-8;
        text(textX, textY, ['$\omega^{',sprintf('%.1f',f(1)),'}$'], 'Color', line_gray, 'Interpreter', 'latex', 'FontSize', FS);
        %}

        Y = T_FL_2(idx1:idx3);
        idx = X >= 1e-7 & X < 5e-5;
        f = polyfit(log10(X(idx)), log10(Y(idx)), 1);
        Xfit = linspace(-7, -5+log10(5), 10);
        Yfit = polyval(f, Xfit);
        Xfit = power(10, Xfit);     Yfit = power(10, Yfit)/2;
        plot(Xfit, Yfit, '-', 'Color', line_green, 'LineWidth', 1);
        textX = 3e-6;
        textY = 2e-9;
        text(textX, textY, ['$\omega^{',sprintf('%.1f',f(1)),'}$'], 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);
        
    elseif J0 == 0.1

        if Tscale_labels
            % Fermi liquid scale
            arrow_X = [0.72, 0.67];
            arrow_Y = [0.52, 0.53];
            text_X = 2.3e-4;
            text_Y = 4.5e-8;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_green, 'LineWidth', LW);
            text(text_X, text_Y, FL_text, 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);

            %{
            % particle-hole symmetry breaking scale
            arrow_X = [0.71, 0.66];
            arrow_Y = [0.665, 0.65];
            text_X = 1.8e-4;
            text_Y = 4.5e-6;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_gray, 'LineWidth', LW);
            text(text_X, text_Y, PHB_text, 'Color', line_gray, 'Interpreter', 'latex', 'FontSize', FS);
            %}

            % fully overscreened phase scale
            arrow_X = [0.48, 0.43];
            arrow_Y = [0.64, 0.63];
            text_X = 1.5e-6;
            text_Y = 2e-6;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_purple, 'LineWidth', LW);
            text(text_X, text_Y, FO_text, 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
    
            % orbital overscreened phase scale
            arrow_X = [0.4, 0.35];
            arrow_Y = [0.795, 0.775];
            text_X = 3e-7;
            text_Y = 2.2e-4;
            annotation('arrow', arrow_X, arrow_Y, 'Color', line_orange, 'LineWidth', LW);
            text(text_X, text_Y, OO_text, 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
        end

        if phase_labels
            FS = 15;
            
            % Fermi liquid
            text_X = 1e-5;
            text_Y = 1e-10;
            text(text_X, text_Y, '$\mathrm{Fermi \ Liquid}$', 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);
        
            % fully overscreened
            text_X = 6e-9;
            text_Y = 1e-7;
            text(text_X, text_Y, '$\mathrm{Fully}$', 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
            text_X = 1.5e-9;
            text_Y = 3e-8;
            text(text_X, text_Y, '$\mathrm{Overscreened}$', 'Color', line_purple, 'Interpreter', 'latex', 'FontSize', FS);
    
            % orbital overscreened
            text_X = 1.5e-9;
            text_Y = 8e-6;
            text(text_X, text_Y, '$\mathrm{Orbital \ Overscreened}$', 'Color', line_orange, 'Interpreter', 'latex', 'FontSize', FS);
        end

        %% linear fit to find power laws

        X = I0(idx1:idx3);
        Y = [T_FO_1(idx1:idx2), T_OO_1(idx2+1:idx3)];
    
        %{
        idx = X >= 1e-7 & X < 5e-5;
        f = polyfit(log10(X(idx)), log10(Y(idx)), 1);
        Xfit = linspace(-7, -5+log10(5), 10);
        Yfit = polyval(f, Xfit);
        Xfit = power(10, Xfit);     Yfit = power(10, Yfit)/2;
        plot(Xfit, Yfit, '-', 'Color', line_gray, 'LineWidth', 1);
        textX = 3e-6;
        textY = 8e-8;
        text(textX, textY, ['$\omega^{',sprintf('%.1f',f(1)),'}$'], 'Color', line_gray, 'Interpreter', 'latex', 'FontSize', FS);
        %}

        Y = T_FL_2(idx1:idx3);
        idx = X >= 1e-7 & X < 5e-5;
        f = polyfit(log10(X(idx)), log10(Y(idx)), 1);
        Xfit = linspace(-7, -5+log10(5), 10);
        Yfit = polyval(f, Xfit);
        Xfit = power(10, Xfit);     Yfit = power(10, Yfit)/2;
        plot(Xfit, Yfit, '-', 'Color', line_green, 'LineWidth', 1);
        textX = 3e-6;
        textY = 2e-9;
        text(textX, textY, ['$\omega^{',sprintf('%.1f',f(1)),'}$'], 'Color', line_green, 'Interpreter', 'latex', 'FontSize', FS);

    end

    set(gca, 'Layer', 'top');
    xlim([I0min, I0max]);
    ylim([Tmin, Tmax]);

    hold off;

end