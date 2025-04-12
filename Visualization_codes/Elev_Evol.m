function Elev_Evol(N, N_elev, Nshow)
    if ~ismember(N_elev, 1:N)
        error('ERR: N_elev must be a natural number that does not exceed N')
    end

    %{
    n = 30;  % Number of colors
    colors = parula(n);  % 'jet', 'hsv', 'parula', etc., also work
    colors = mat2cell(colors, ones(1, n), 3); % Convert to cell array
    %}

    colors = {...
    [1, 0, 0],        % Red
    [0, 1, 0],        % Green
    [0, 0, 1],        % Blue
    [1, 1, 0],        % Yellow
    [1, 0, 1],        % Magenta
    [0, 1, 1],        % Cyan
    [0.5, 0.5, 0.5],  % Gray
    [0.8, 0.4, 0],    % Orange
    [0.5, 0, 0.5],    % Purple
    [0, 0.5, 0],      % Dark Green
    [0.5, 0.5, 0],    % Olive
    [0, 0.5, 0.5],    % Teal
    [0.5, 0, 0],      % Dark Red
    [0.7, 0.7, 0.7],  % Light Gray
    [0.3, 0.3, 0.3],  % Dark Gray
    [1, 0.5, 0.5],    % Light Red
    [0.5, 1, 0.5],    % Light Green
    [0.5, 0.5, 1],    % Light Blue
    [1, 0.7, 0],      % Golden Yellow
    [0.7, 0, 1],      % Violet
    [0, 1, 0.7],      % Mint Green
    [1, 0, 0.7],      % Pink
    [0, 0.7, 1],      % Sky Blue
    [0.6, 0.3, 0],    % Brown
    [0.9, 0.6, 0],    % Amber
    [0.8, 0.1, 0.6],  % Deep Pink
    [0.2, 0.6, 0.3],  % Forest Green
    [0.3, 0.2, 0.8],  % Indigo
    [0.8, 0.3, 0.2],  % Rust
    [0.9, 0.8, 0],    % Mustard
    [0, 0.8, 0.2],    % Bright Green
    [0.4, 0.2, 0],    % Chocolate
    [0.6, 0.2, 0.8],  % Orchid
    [0.3, 0.7, 0.9],  % Turquoise
    [0.9, 0.5, 0.3],  % Peach
    [0.3, 0.9, 0.5],  % Light Sea Green
    [0.7, 0.2, 0.5],  % Rose
    [0.5, 0.8, 0.2],  % Apple Green
    [0.2, 0.2, 0.2]   % Almost Black
    };

    linestyles = {'--', '-', ':', '-.'};
    
    K_z = 0:0.1:1;
    Q = 0.01*ones(numel(K_z));
    %{
    K_z = [0.1:0.1:0.5, 0.51:0.01:0.59, 0.6:0.1:1];
    K_z = K_z(1:7);
    Q = 0.1*ones(numel(K_z));
    Lbd_z_fix = [0.05236, 0.1029, 0.1539, 0.2057, 0.2642, 0.2725, 0.2849];
    %}
    %{
    K_z = [0.1:0.1:0.4, 0.41:0.01:0.49, 0.5:0.1:1];
    K_z = K_z(1:9);
    Q = 0.2*ones(numel(K_z));
    Lbd_z_fix = [0.05948, 0.1119, 0.1660, 0.2262, 0.2333, 0.2407, 0.2487, 0.2585, 0.2701];
    %}
    %{
    K_z = [0.1:0.1:0.4, 0.41:0.01:0.49, 0.5:0.1:1];
    K_z = K_z(1:3);
    Q = 0.3*ones(numel(K_z));
    Lbd_z_fix = [0.07201, 0.1277, 0.1890];
    %}

    %K_z = 0.1:0.1:1;
    %Q = zeros(1,10);
    Elev_ev = cell(1,numel(K_z));

    for it = 1:numel(K_z)
        path = 'C:\Users\82104\Documents\Physics\Research\data\Quartic\K_z=';
        path = [path, num2str(K_z(it))];
        path = [path, '_Q=', num2str(Q(it)), '_T=1e-24_Nkeep=3000'];

        Etot = load([path, '\Etot.mat']);
        Etot = Etot.Etot;

        Qtot = load([path, '\Qtot.mat']);
        Qtot = Qtot.Qtot;

        if it == 1
            Elev_ev{it} = Elev_label(Etot, Qtot, N, N_elev);
        else
            Elev_ev{it} = Elev_label(Etot, Qtot, N, N_elev*5);
        end
    end

    cnt = 0;
    Qnum_groups = zeros(1,5);
    Qnum_gp_colors = cell(1,1);
    Qnum_gp_data = cell(1,1);
    legends = cell(1,1);
    Esep = [];
    Lbd_z_fit = [];
    for it1 = 1:numel(Elev_ev)
        found = [false, false];

        for it2 = 1:size(Elev_ev{it1},1)
            Qnums = round( Elev_ev{it1}(it2, 2:6) );
            Elev = Elev_ev{it1}(it2,1);
            if ismember(Qnums, Qnum_groups, 'rows')
                [~,idx] = ismember(Qnums, Qnum_groups, 'rows');
                Qnum_gp_data{idx} = [ Qnum_gp_data{idx}; [K_z(it1), Elev] ];
            else
                if it1 == 1
                    cnt = cnt + 1;
                    Qnum_groups = [Qnum_groups; Qnums];
                    Qnum_gp_colors = cat(2, Qnum_gp_colors, colors{rem(cnt,39)+1});
                    Qnum_gp_data = cat(2, Qnum_gp_data, {[K_z(it1), Elev]} );
                    if cnt <= Nshow
                        legends = cat(1, legends, {['[',num2str(Qnums(1)), ',', num2str(Qnums(2)), ',', num2str(Qnums(3)), ',', num2str(Qnums(4)), ']']});
                    end
                end
            end

            if isequal(Qnums, [-1,0,0,1,1])
                E1 = Elev;
                found(1) = true;
            elseif isequal(Qnums, [1,-2,0,1,1])
                E2 = Elev;
                found(2) = true;
            end

        end

        if found(1) && found(2)
            if K_z(it1) < 0.65
                Esep = [Esep, E2-E1];
                Lbd_z_fit = [Lbd_z_fit, K_z(it1)/2];
            end
       end
    end
    Qnum_groups(1,:) = [];
    Qnum_gp_colors(1) = [];
    Qnum_gp_data(1) = [];
    legends(1) = [];

    coeff = polyfit(Lbd_z_fit,Esep,1);
    Efit = polyval(coeff, Lbd_z_fit);
    SStot = sum((Esep-mean(Esep)).^2);
    SSres = sum((Esep-Efit).^2);
    Rsq = 1 - SSres/SStot;
    disp(Esep);
    disp(coeff);
    disp(Rsq);

    figure;
    hold on;
    for it = 1:size(Qnum_groups,1)

        plot(Qnum_gp_data{it}(:,1)/2, Qnum_gp_data{it}(:,2), linestyles{rem(it,4)+1}, 'linewidth', 1.5, 'Color', Qnum_gp_colors{it});
        %plot(Lbd_z_fix, Qnum_gp_data{it}(1:3,2), linestyles{rem(it,4)+1}, 'linewidth', 1.5, 'Color', Qnum_gp_colors{it});
        %{
        if isequal(Qnum_groups(it,:), [-2, 1, 1, 0, 1])
            plot(Qnum_gp_data{it}(:,1)/2, Qnum_gp_data{it}(:,2), linestyles{rem(it,4)+1}, 'linewidth', 1.5, 'Color', Qnum_gp_colors{it});
        end
        %}
    end
    set(gca,'fontsize',20);
    xlabel('$\lambda_{z}$','Interpreter','latex','FontSize',25);
    ylabel('$\mathrm{Rescaled \ energy}$','Interpreter','latex','FontSize',25);
    legend(legends,'Location','eastoutside','FontSize',15, 'Box', 'off');
    title('$\mathrm{Lowest \ Energy \ Levels} \left( \lambda_{x} = 0.01 \right)$', 'Interpreter', 'latex', 'FontSize', 30);
    
    hold off;

end