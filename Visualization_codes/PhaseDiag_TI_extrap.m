function PhaseDiag_TI_extrap(J0, K0)
    % <Description>
    % Draws T-I phase diagram from NRG data
    %
    % <Input>
    % J0 : [numeric] spin-spin coupling strength in the isotropic 2soK model
    % K0 : [numeric] pseudospin-pseudospin coupling strength in the isotropic 2soK model
    %

    NFLboundary = true;
    I0min = 1e-11;
    Ymin = 1e-20;
    NeedExtrap = true;
    ShowAll = true;
    ShowMarkers = true;

    %% Define phase boundaries using FL scale

    [I0, Phase_range, Phase_name] = PhaseRange_TI(J0,K0,'I0min',I0min);

    T_K = nan(1,numel(I0));
    T_K_diff = nan(1,numel(I0));
    T_2 = nan(1,numel(I0));
    T_K_min = -25;

    X_FL_L = [];
    X_FL_R = [];
    Y_FL_L = []; 
    Y_FL_R = [];
    X_FO = [];
    Y_FO = [];

    for it1 = 1:numel(I0)
        for it2 = 1:numel(Phase_name{it1})

            if it2 < numel(Phase_name{it1})

                if isequal(Phase_name{it1}{it2}, 'Fully Overscreened') && isequal(Phase_name{it1}{it2+1}, 'Orbital Overscreened')
    
                    T_2(it1) = Phase_range{it1}(it2,2);
                    T_K_min = Phase_range{it1}(it2,2);
    
                    X_FO = [X_FO, I0(it1)];
                    Y_FO = cat(2, Y_FO, [Phase_range{it1}(it2,1); Phase_range{it1}(it2,2)] );

                end
    
            elseif isequal(Phase_name{it1}{it2}, 'Fully Overscreened')
    
                T_2(it1) = Phase_range{it1}(it2,2);
                T_K_min = Phase_range{it1}(it2,2);
    
                X_FO = [X_FO, I0(it1)];
                Y_FO = cat(2, Y_FO, [Phase_range{it1}(it2,1); Phase_range{it1}(it2,2)] );
            end

            if isequal(Phase_name{it1}{it2}, 'Fermi Liquid')
                    
                T_K(it1) = Phase_range{it1}(it2,2);

                if I0(it1) < 0
                    X_FL_L = [X_FL_L, I0(it1)];
                    Y_FL_L = cat(2, Y_FL_L, [Phase_range{it1}(it2,1); Phase_range{it1}(it2,2) ] );
                else
                    X_FL_R = [X_FL_R, I0(it1)];
                    Y_FL_R = cat(2, Y_FL_R, [Phase_range{it1}(it2,1); Phase_range{it1}(it2,2) ] );
                end
            end

        end % it2
    end % it1

    %% Extrapolate FL boundaries

    if NeedExtrap

        logYmin = log10(Ymin);

        EraseIdx = Y_FO(1,:) > logYmin & X_FO < 0;
        X_FO(EraseIdx) = [];
        Y_FO(:, EraseIdx) = [];
        EraseIdx = Y_FO(1,:) > logYmin & X_FO > 0;
        X_FO(EraseIdx) = [];
        Y_FO(:, EraseIdx) = [];
        
        % coordinates of upper boundary of FL region to be linearly fitted
        X_FL_L_fit = log10(abs(X_FL_L(Y_FL_L(2,:) < T_K_min)));
        Y_FL_L_fit = Y_FL_L(2, Y_FL_L(2,:) < T_K_min);

        % linear fit
        coeff_L = polyfit(X_FL_L_fit, Y_FL_L_fit, 1);

        % extrapolation from the end of FL region to FO region
        I0_Extrap_L = linspace(log10(abs(X_FL_L(end)))-0.1, log10(abs(X_FO(1))), 10);
        X_FL_L_Extrap = -power(10, I0_Extrap_L);
        Y_FL_L_Extrap = [logYmin*ones(1,10); polyval(coeff_L, I0_Extrap_L)];

        % discard extrapolated boundary under Ymin
        X_FL_L_Extrap = X_FL_L_Extrap(Y_FL_L_Extrap(2,:) > logYmin);
        Y_FL_L_Extrap = Y_FL_L_Extrap(:, Y_FL_L_Extrap(2,:) > logYmin);


        X_FL_L = [X_FL_L, X_FL_L_Extrap];
        Y_FL_L = cat(2, Y_FL_L, Y_FL_L_Extrap);
    
        X_FL_R_fit = log10(abs(X_FL_R(Y_FL_R(2,:) < T_K_min)));
        Y_FL_R_fit = Y_FL_R(2, Y_FL_R(2,:) < T_K_min);
        coeff_R = polyfit(X_FL_R_fit, Y_FL_R_fit, 1);
        I0_Extrap_R = linspace(log10(X_FO(end)), log10(abs(X_FL_R(1)))-0.1, 10);
        X_FL_R_Extrap = power(10, I0_Extrap_R);
        Y_FL_R_Extrap = [logYmin*ones(1,10); polyval(coeff_R, I0_Extrap_R)];
        X_FL_R_Extrap = X_FL_R_Extrap(Y_FL_R_Extrap(2,:) > logYmin);
        Y_FL_R_Extrap = Y_FL_R_Extrap(:, Y_FL_R_Extrap(2,:) > logYmin);
        X_FL_R = [X_FL_R_Extrap, X_FL_R];
        Y_FL_R = cat(2, Y_FL_R_Extrap, Y_FL_R);
    end


    %% Adjust NFL boundary to the extrapolated FL boundary

    if NeedExtrap
            
        X_FL_L_under = X_FL_L(Y_FL_L(2,:) <= T_K_min);
        Y_FL_L_under = Y_FL_L(2, Y_FL_L(2,:) <= T_K_min);
        X_FL_R_under = X_FL_R(Y_FL_R(2,:) <= T_K_min);
        Y_FL_R_under = Y_FL_R(2, Y_FL_R(2,:) <= T_K_min);
        
        X_FO = [X_FL_L_under, X_FO, X_FL_R_under];
        Y_FO = cat(2, [Y_FL_L_under; T_K_min*ones(1, numel(Y_FL_L_under))], Y_FO);
        Y_FO = cat(2, Y_FO, [Y_FL_R_under; T_K_min*ones(1, numel(Y_FL_R_under))]);
        
        Idx_L = find(Y_FL_L(2,:) > T_K_min, 1, 'last');
        Idx_R = find(Y_FL_R(2,:) > T_K_min, 1, 'first');
        X_FO = [X_FL_L(Idx_L), X_FO, X_FL_R(Idx_R)];
        Y_FO = cat(2, repmat(Y_FL_L(2,Idx_L), [2,1]), Y_FO);
        Y_FO = cat(2, Y_FO, repmat(Y_FL_R(2,Idx_R), [2,1]));
    end

    %% Define phase region patches
    X_FL_L = [X_FL_L, fliplr(X_FL_L)];
    Y_FL_L = [Y_FL_L(1,:), fliplr(Y_FL_L(2,:))];
    X_FL_R = [X_FL_R, fliplr(X_FL_R)];
    Y_FL_R = [Y_FL_R(1,:), fliplr(Y_FL_R(2,:))];
    X_FO = [X_FO, fliplr(X_FO)];
    Y_FO = [Y_FO(1,:), fliplr(Y_FO(2,:))];

    %% Exponentiate
    T_K = power(10,T_K);
    T_K_diff = power(10,T_K_diff);
    T_2 = power(10,T_2);
    T_K_min = power(10,T_K_min);

    %T_K(T_K < T_K_min) = nan;
    T_2(T_K < T_K_min) = T_K_min;
    
    Y_FL_L = power(10,Y_FL_L);
    Y_FL_R = power(10,Y_FL_R);
    Y_FO = power(10,Y_FO);
    
    figure;
    hold on;
    ylim([Ymin, 1e-1]);
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.XAxis.MinorTick = 'on';
    ax.Layer = 'top';
    title(['$J_{0} = ',sprintf('%.15g',J0),', K_{0} = ',sprintf('%.15g',K0),'$'], 'interpreter', 'latex', 'fontsize', 21);
    xlabel('$I_{0}$', 'interpreter', 'latex', 'fontsize', 25);
    ylabel('$\mathrm{T}$', 'interpreter', 'latex', 'fontsize', 25);
    [~,lin2sym_X,~] = SLplot(I0, T_K, 'XScale', 'symlog', 'YScale', 'log');

    X_FL_L = lin2sym_X(X_FL_L);
    X_FL_R = lin2sym_X(X_FL_R);
    X_FO = lin2sym_X(X_FO);
    patch(X_FL_L, Y_FL_L, [.753,.878,.753], 'FaceAlpha', 1, 'linestyle', 'none');
    patch(X_FL_R, Y_FL_R, [.753,.878,.753], 'FaceAlpha', 1, 'linestyle', 'none');
    patch(X_FO, Y_FO, 1.1*[.745 .682 .898], 'FaceAlpha', 1, 'linestyle', 'none');

    if ShowMarkers
        I0_sym = lin2sym_X(I0);
        plot(I0_sym(T_K > T_K_min), T_K(T_K > T_K_min), '.', 'color', 'red', 'MarkerSize', 10);
        plot(I0_sym(T_K < T_K_min), T_K(T_K < T_K_min), '.', 'color', 'black', 'MarkerSize', 10);
        plot(I0_sym, T_K_diff, '.', 'Color', 'black', 'MarkerSize', 10);
        plot(I0_sym, T_2, '.', 'Color', 'blue', 'MarkerSize', 10);
    end

    if NeedExtrap
        X_FL_L_Extrap_sym = lin2sym_X(X_FL_L_Extrap);
        X_FL_R_Extrap_sym = lin2sym_X(X_FL_R_Extrap);
        plot(X_FL_L_Extrap_sym, power(10, Y_FL_L_Extrap(2,:)), '--', 'Color', 'black', 'LineWidth', 1);
        plot(X_FL_R_Extrap_sym, power(10, Y_FL_R_Extrap(2,:)), '--', 'Color', 'black', 'LineWidth', 1);
    end
    



    %% Define phase boundaries using NFL scales

    [I0, Phase_range, Phase_name] = PhaseRange_TI(J0,K0,'NFLboundary','I0min',I0min);

    T_K = nan(1,numel(I0));
    T_K_diff = nan(1,numel(I0));
    T_2 = nan(1,numel(I0));
    T_K_min = -25;

    X_FL_L = [];
    X_FL_R = [];
    Y_FL_L = []; 
    Y_FL_R = [];
    X_FO = [];
    Y_FO = [];

    for it1 = 1:numel(I0)
        for it2 = 1:numel(Phase_name{it1})

            if it2 < numel(Phase_name{it1})
                if isequal(Phase_name{it1}{it2}, 'Fermi Liquid') && ismember(Phase_name{it1}{it2+1}, {'Orbital Overscreened', 'Spin Overscreened', 'Fully Overscreened'})
                        
                    T_K(it1) = Phase_range{it1}(it2+1,1);
    
                    if I0(it1) < 0
                        X_FL_L = [X_FL_L, I0(it1)];
                        Y_FL_L = cat(2, Y_FL_L, [Phase_range{it1}(it2,1); Phase_range{it1}(it2+1,1) ] );
                    else
                        X_FL_R = [X_FL_R, I0(it1)];
                        Y_FL_R = cat(2, Y_FL_R, [Phase_range{it1}(it2,1); Phase_range{it1}(it2+1,1) ] );
                    end
                        
                elseif isequal(Phase_name{it1}{it2}, 'Fully Overscreened') && isequal(Phase_name{it1}{it2+1}, 'Orbital Overscreened')
    
                    T_2(it1) = Phase_range{it1}(it2,2);
                    T_K_min = Phase_range{it1}(it2,2);
    
                    X_FO = [X_FO, I0(it1)];
                    Y_FO = cat(2, Y_FO, [Phase_range{it1}(it2,1); Phase_range{it1}(it2,2)] );

                end
    
            elseif isequal(Phase_name{it1}{it2}, 'Fermi Liquid')
                T_K_diff(it1) = Phase_range{it1}(it2,2);
    
            elseif isequal(Phase_name{it1}{it2}, 'Fully Overscreened')
    
                T_2(it1) = Phase_range{it1}(it2,2);
                T_K_min = Phase_range{it1}(it2,2);
    
                X_FO = [X_FO, I0(it1)];
                Y_FO = cat(2, Y_FO, [Phase_range{it1}(it2,1); Phase_range{it1}(it2,2)] );
            end

        end % it2
    end % it1

    %% Define phase regions

    %{}
    Idx_L = find(Y_FL_L(2,:) > T_K_min, 1, 'last');
    Idx_R = find(Y_FL_R(2,:) > T_K_min, 1, 'first');
    X_FO = [X_FL_L(Idx_L), X_FO, X_FL_R(Idx_R)];
    Y_FO = cat(2, repmat(Y_FL_L(2,Idx_L), [2,1]), Y_FO);
    Y_FO = cat(2, Y_FO, repmat(Y_FL_R(2,Idx_R), [2,1]));
    %}

    %% Define phase region patches
    X_FL_L = [X_FL_L, fliplr(X_FL_L)];
    Y_FL_L = [Y_FL_L(1,:), fliplr(Y_FL_L(2,:))];
    X_FL_R = [X_FL_R, fliplr(X_FL_R)];
    Y_FL_R = [Y_FL_R(1,:), fliplr(Y_FL_R(2,:))];
    X_FO = [X_FO, fliplr(X_FO)];
    Y_FO = [Y_FO(1,:), fliplr(Y_FO(2,:))];

    %% Exponentiate
    T_K = power(10,T_K);
    T_K_diff = power(10,T_K_diff);
    T_2 = power(10,T_2);
    T_K_min = power(10,T_K_min);
    T_2(T_K < T_K_min) = T_K_min;
    
    Y_FL_L = power(10,Y_FL_L);
    Y_FL_R = power(10,Y_FL_R);
    Y_FO = power(10,Y_FO);

    ylim([Ymin, 1e-1]);
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.XAxis.MinorTick = 'on';
    ax.Layer = 'top';
    title(['$J_{0} = ',sprintf('%.15g',J0),', K_{0} = ',sprintf('%.15g',K0),'$'], 'interpreter', 'latex', 'fontsize', 21);
    xlabel('$I_{0}$', 'interpreter', 'latex', 'fontsize', 25);
    ylabel('$\mathrm{T}$', 'interpreter', 'latex', 'fontsize', 25);
    [~,lin2sym_X,~] = SLplot(I0, T_K, 'XScale', 'symlog', 'YScale', 'log');

    if ShowAll
        X_FL_L = lin2sym_X(X_FL_L);
        X_FL_R = lin2sym_X(X_FL_R);
        X_FO = lin2sym_X(X_FO);
        %{
        patch(X_FL_L, Y_FL_L, [.753,.878,.753], 'FaceAlpha', 1, 'linestyle', 'none');
        patch(X_FL_R, Y_FL_R, [.753,.878,.753], 'FaceAlpha', 1, 'linestyle', 'none');
        %}
        patch(X_FO, Y_FO, 0.9*[.745 .682 .898], 'FaceAlpha', 1, 'linestyle', 'none');
    
        if ShowMarkers
            I0_sym = lin2sym_X(I0);
            plot(I0_sym(T_K > T_K_min), T_K(T_K > T_K_min), '.', 'color', 'red', 'MarkerSize', 10);
            plot(I0_sym(T_K > T_K_min & I0 < 0), T_K(T_K > T_K_min & I0 < 0), '--', 'color', 'red', 'MarkerSize', 10);
            
            plot(I0_sym(T_K > T_K_min & I0 > 0), T_K(T_K > T_K_min & I0 > 0), '--', 'color', 'red', 'MarkerSize', 10);
            plot(I0_sym(T_K < T_K_min), T_K(T_K < T_K_min), '.', 'color', 'black', 'MarkerSize', 10);
            %plot(I0_sym, T_2, '.', 'Color', 'blue', 'MarkerSize', 10);
        end
    
        xlim([-0.8709, 0.8709]);
        hold off;
    end

end