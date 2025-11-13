function NRG_RGflow_fig(J0, K0, I0, FitInfo, varargin)

    YLims = [1e-6, 1e16;
             1e-17, 5e14];
    SubFig_idx = {'$\mathrm{(a)}$', '$\mathrm{(b)}$', '$\mathrm{(c)}$'};

    alphaz = 0.95;
    while ~isempty(varargin)

        switch varargin{1}
            case 'YLims'
                YLims = varargin{2};
                varargin(1:2) = [];

            case 'SubFig_idx'
                SubFig_idx = varargin{2};
                varargin(1:2) = [];

            case 'alphaz'
                alphaz = varargin{2};
                varargin(1:2) = [];

            otherwise
        end % switch - case
    end % while

    if isscalar(alphaz)
        alphaz = alphaz*[1,1];
    end

    QnumKeys = [
        0, 1, 1;
        -1, 0, 0;
        1, 0, 0
        -1, 2, 0;
        1, 2, 0;
        1, 0, 2
        -1, 0, 2;
        2, 1, 1;
        -2, 1, 1
        ];

    ColorDict = [
        .486, .282, .125;
        .466, .674, .188;
        1, 0, 0;
        0, .447, .741;
        .929, .694, .125;
        1, 0, 1;
        0, .5, .3;
        .475, .435, .655;
        .745, .792, .259
        ];

    LinestyleDict = {'-', '-', '--', '-', '--', '--', '-', '--', '-'};
    
    N_elev = 40;
    N_legends = 20;

    Tmin = 1e-20;
    Tmax = 1;
    Emin = -0.5;
    ImpSuscMin = YLims(1,1);
    ImpSuscMax = YLims(1,2);
    BathSuscMin = YLims(2,1);
    BathSuscMax = YLims(2,2);
    if J0*K0 < 0
        Emax = 2;
    else
        Emax = 1;
    end
    
    %% Load data
    Niter = 120;
    Nkeep = 3000;
    Lambda = 2.5;
    TsoK_path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK_TI_selected';
    Data_path = ['J0=', sprintf('%.15g',J0), '_K0=', sprintf('%.15g', K0), ...
                    '_I0=', sprintf('%.15g', I0),'_T=1e-24_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];

    Full_path = [TsoK_path, filesep, Data_path];
    FileInfo = dir(Full_path);
    
    Etot = load([Full_path, filesep, 'Etot.mat']);      Etot = Etot.Etot;
    Qtot = load([Full_path, filesep, 'Qtot.mat']);      Qtot = Qtot.Qtot;
    if ismember('DiscData.mat', {FileInfo.name})
        DiscData = load([Full_path, filesep, 'DiscData.mat']);      DiscData = DiscData.DiscData;
    else
        DiscData = [];
    end

    ImpSusc = cell(1,3);
    tmp = load([Full_path, filesep, 'NRG_Op=ImpSp.mat']);   ImpSusc{1} = tmp.temp;
    tmp = load([Full_path, filesep, 'NRG_Op=ImpOrb.mat']);  ImpSusc{2} = tmp.temp;
    tmp = load([Full_path, filesep, 'NRG_Op=ImpSpOrb.mat']);  ImpSusc{3} = tmp.temp;

    BathSusc = cell(1,3);
    tmp = load([Full_path, filesep, 'NRG_Op=BathSp.mat']);   BathSusc{1} = tmp.temp;
    tmp = load([Full_path, filesep, 'NRG_Op=BathOrb.mat']);  BathSusc{2} = tmp.temp;
    tmp = load([Full_path, filesep, 'NRG_Op=BathSpOrb.mat']);  BathSusc{3} = tmp.temp;

    ocont = load([Full_path, filesep, 'ocont.mat']);            ocont = ocont.ocont;

    %% Calculate Acont if "DiscData" is available

    if ~isempty(DiscData)

        odisc = DiscData.odisc;
        sigmak = DiscData.sigmak;
        Adiscs = DiscData.Adiscs;
        nz = DiscData.nz;
        emin = DiscData.emin;
        T = 1e-24;
        
        for ita = (1:size(Adiscs,1))
        
            Adisc = mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3);
            
            [ocont, Acont] = getAcont(odisc,Adisc,sigmak,T/5,'alphaz',alphaz(2),'emin',emin,'tol',1e-19);
            
            if ita <=3
                ImpSusc{ita} = Acont;
            else
                BathSusc{ita-3} = Acont;
            end

        end % ita
    end

    for itN = 1:3
        ImpSusc{itN} = ImpSusc{itN}(ocont>0);
        BathSusc{itN} = BathSusc{itN}(ocont>0);
    end
    ocont = ocont(ocont>0);
    

    %% Label and rescale energy spectrum

    Elev = cell(1, N_elev);
    Qnum = cell(1, N_elev);
    
    Elev_plot = cell(1,0);
    Qnum_plot = zeros(0,4);
    
    for itN = 1:floor(Niter/2)
    
        Elev_lab = Elev_label(Etot, Qtot, 2*itN-1, N_elev);
        Elev{itN} = Elev_lab(:,1);
        Qnum{itN} = Elev_lab(:,2:5) + [1,0,0,0];    % shift U(1) charge by +1
        idx = find(ismember(Qnum{itN}, [0,1,1,1], 'rows'));
        Elev{itN} = Elev{itN} - Elev{itN}(idx);
    
        for itL = 1:N_elev
    
            Enow = Elev{itN}(itL);
            Qnow = Qnum{itN}(itL,:);
    
            if ismember(Qnow, Qnum_plot, 'rows')
                idx = find(ismember(Qnum_plot, Qnow, 'rows'));
                Elev_plot{idx}(itN) = Enow;
            else
                Qnum_plot = [Qnum_plot; Qnow];
                Elev_plot = cat(2, Elev_plot, [nan(1,itN-1), Enow, nan(1, floor(Niter/2)-itN)]);
            end
        end
    end

    %% Extract energy scales and phase widths

    [~, Phase_range, Phase_name] = PhaseRange_TI(J0, K0, 'I0', I0, 'showAllBound', 'flatThres', 0.2, 'alphaz', alphaz(1), 'InflecBound');

    E_scale = nan(1,4); % [FO, SO, OO, FL]
    E_width = nan(1,4); % [FO, SO, OO, FL]

    for itP = 1:numel(Phase_name)

        switch Phase_name{itP}
            case 'Fully Overscreened'
                E_scale(1) = power(10, Phase_range(itP,2));
                E_width(1) = Phase_range(itP,2) - Phase_range(itP,1);
            case 'Spin Overscreened'
                E_scale(2) = power(10, Phase_range(itP,2));
                E_width(2) = Phase_range(itP,2) - Phase_range(itP,1);
            case 'Orbital Overscreened'
                E_scale(3) = power(10, Phase_range(itP,2));
                E_width(3) = Phase_range(itP,2) - Phase_range(itP,1);
            case 'Fermi Liquid'
                E_scale(4) = power(10, Phase_range(itP,2));
                E_width(4) = Phase_range(itP,2) - Phase_range(itP,1);
            otherwise
        end % switch - case
    end % itP


    Thres = 1e-8;

    [~, Phase_range, Phase_name] = PhaseRange_TI(J0, K0, 'I0', I0, 'showAllBound', 'flatThres', 0.2, 'alphaz', alphaz(2), 'InflecBound');

    for itP = 1:numel(Phase_name)

        switch Phase_name{itP}
            case 'Fully Overscreened'
                if E_scale(1) < Thres
                    E_scale(1) = power(10, Phase_range(itP,2));
                    E_width(1) = Phase_range(itP,2) - Phase_range(itP,1);
                end

            case 'Spin Overscreened'
                if E_scale(2) < Thres
                    E_scale(2) = power(10, Phase_range(itP,2));
                    E_width(2) = Phase_range(itP,2) - Phase_range(itP,1);
                end

            case 'Orbital Overscreened'
                if E_scale(3) < Thres
                    E_scale(3) = power(10, Phase_range(itP,2));
                    E_width(3) = Phase_range(itP,2) - Phase_range(itP,1);
                end

            case 'Fermi Liquid'
                if E_scale(4) < Thres
                    E_scale(4) = power(10, Phase_range(itP,2));
                    E_width(4) = Phase_range(itP,2) - Phase_range(itP,1);
                end

            otherwise
        end % switch - case
    end % itP



    %% Plot figure

    figureHandle = figure('Position', [100, -50, 550, 1100]);    % Width = 550 px, Height = 1100 px
    hold on;

    % subplot labels
    annotation(figureHandle, 'textbox', [0.04, 0.995, 0, 0], ...
    'String', SubFig_idx{1}, 'Interpreter', 'latex', 'FontSize', 18, ...
    'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'LineStyle', 'none');

    annotation(figureHandle, 'textbox', [0.04, 0.54, 0, 0], ...
        'String', SubFig_idx{2}, 'Interpreter', 'latex', 'FontSize', 18, ...
        'FitBoxToText', 'on', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'LineStyle', 'none');

    annotation(figureHandle, 'textbox', [0.04, 0.31, 0, 0], ...
        'String', SubFig_idx{3}, 'Interpreter', 'latex', 'FontSize', 18, ...
        'FitBoxToText', 'on', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'LineStyle', 'none');

    %% plot subplot 1 ========================================================================
    
    ax1 = subplot(3,1,1);
    pos1 = [0.17, 0.68, 0.8, 0.3]; % [left bottom width height]
    set(ax1, 'Position', pos1); 
    hold on;

    % set axis properties
    set(ax1, 'XScale', 'log', 'YScale', 'linear');
    set(ax1, 'XTickLabel', []);
    set(ax1, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    ax1.XTick = 10.^(-20:4:0);
    Ystep = round(2*(Emax-Emin)/3)/2;
    ax1.YTick = Ystep*ceil(Emin/Ystep) : Ystep : Emax;
    ax1.TickLength = [0.015, 0.002];      % ticksize : [major, minor]
    set(ax1, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(ax1, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);
    ax1.XAxis.FontSize = 24;
    ax1.YAxis.FontSize = 24;
    set(ax1, 'Layer', 'top');
    xlim([Tmin, Tmax]);
    ylim([Emin,Emax]);

    % define x- and y-labels
    ylabel('$\mathrm{rescaled \ energy}$', 'Interpreter', 'latex', 'FontSize', 20);

    hx = get(ax1, 'XLabel');
    hy = get(ax1, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';

    % Modify y-label positions manually
    if Emax == 2
        hy.Position = hy.Position + [-0.045, 0, 0];
    else
        hy.Position = hy.Position + [0.04, 0.002, 0];
    end

    % remove x-tick labels
    tickLabelFont = 18;
    ax1.XTickLabel = [];

    % redefine y-tick labels using annotation()
    ax1.YTickLabel = [];
    yticks = Ystep*ceil(Emin/Ystep) : Ystep : Emax;
    
    ytickWidth = 0.09;
    ytickHeight = 0.15;
    for ity = 1:numel(yticks)
        X_pos = pos1(1) - 0.09;
        Y_pos = pos1(2) + pos1(4) * (yticks(ity) - Emin) / (Emax - Emin);
        X_pos = X_pos - 0;
        Y_pos = Y_pos - 0.13;
        annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',sprintf('%.15g',yticks(ity)+0),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end

    %% plot rescaled energy levels and legends
    legends = cell(1, N_legends);
    handles = zeros(1, N_legends);

    cnt_leg = 0;    % counter for legends
        
    for it = 1:numel(Elev_plot)

        if ismember(Qnum_plot(it,1:3), QnumKeys, 'rows')
            key = find(ismember(QnumKeys, Qnum_plot(it,1:3), 'rows'));

            Color = ColorDict(key,:);
            linestyle = LinestyleDict{key};
            linewidth = 2;
        else
            Color = [.7, .7, .7];
            linestyle = '-';
            linewidth = 0.7;
        end
                
        if cnt_leg <= N_legends-1 && ismember(Qnum_plot(it,1:3), QnumKeys, 'rows') && Qnum_plot(it,4) == 1
    
            cnt_leg = cnt_leg + 1;
            h = plot(sqrt(Lambda) .^ -(1:2:Niter), Elev_plot{it}, linestyle, 'Color', Color, 'LineWidth', linewidth);
            Q_Sp = HalfInt(Qnum_plot(it,2)/2);
            Q_orb = HalfInt(Qnum_plot(it,3)/2);
            legends{cnt_leg} = ['$\left(', signed_int_str(Qnum_plot(it,1)), ', ' Q_Sp, ', ', Q_orb, '\right)$'];
            handles(cnt_leg) = h;
        else
            plot(sqrt(Lambda) .^ -(1:2:Niter), Elev_plot{it}, linestyle, 'Color', Color, 'LineWidth', linewidth);
        end
    end

    legends = legends(1:cnt_leg);
    handles = handles(1:cnt_leg);
    
    legend(handles, legends, 'Interpreter', 'latex', 'NumColumns', 3, 'Location', 'best', 'FontSize', 18, 'Box', 'off');

    hold off;

    %% plot subplot 2 ========================================================================

    ax2 = subplot(3,1,2);
    pos2 = [0.17, 0.33, 0.8, 0.2]; % [left bottom width height]
    set(ax2, 'Position', pos2);
    hold on;

    % set axis properties
    set(ax2, 'XScale', 'log', 'YScale', 'log', 'FontSize', 18);
    set(ax2, 'XTickLabel', []);
    set(ax2, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    ax2.XTick = 10.^( 5*ceil(log10(Tmin)/5) : 5 : 5*floor(log10(Tmax)/5) );
    ax2.YTick = 10.^( 5*ceil(log10(ImpSuscMin)/5) : 5 : 5*floor(log10(ImpSuscMax)/5) );
    ax2.TickLength = [0.015, 0.002];      % ticksize : [major, minor]
    set(ax2, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(ax2, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);
    ax2.XAxis.FontSize = 24;
    ax2.YAxis.FontSize = 24;
    set(ax2, 'Layer', 'top');
    xlim([Tmin, Tmax]);
    ylim([ImpSuscMin, ImpSuscMax]);

    % define y-label
    ylabel('$\chi^{\mathrm{imp}}$', 'Interpreter', 'latex', 'FontSize', 20);

    hx = get(ax2, 'XLabel');
    hy = get(ax2, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';

    % Modify y-label position manually
    hy.Position = hy.Position + [0.06, 0, 0];

    % redefine y-tick labels using annotation()
    ax2.YTickLabel = [];
    yticks = 10.^( 5*ceil(log10(ImpSuscMin)/5) : 5 : 5*floor(log10(ImpSuscMax)/5) );
    
    ytickWidth = 0.09;
    ytickHeight = 0.1;
    for ity = 1:numel(yticks)
        X_pos = pos2(1) - 0.09;
        Y_pos = pos2(2) + pos2(4) * (log10(yticks(ity)) - log10(ImpSuscMin)) / (log10(ImpSuscMax) - log10(ImpSuscMin));
        X_pos = X_pos + 0.005;
        Y_pos = Y_pos - 0.085;
        annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',SciNot(yticks(ity)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end
    
    % plot impurity dynamic susceptibilities
    hdl_sp = plot(ocont, ImpSusc{1}, 'Color', [.466, .674, .188], 'LineWidth', 2);
    hdl_orb = plot(ocont, ImpSusc{2}, 'Color', [.850, .325, .098], 'LineWidth', 2);
    hdl_sporb = plot(ocont, ImpSusc{3}, 'Color', [0, .447, .741], 'LineWidth', 2);

    TempFontSize = 16;
    FS = 18;

    % linear fit and show power laws
    for itS = 1:3
        FitRange = FitInfo.Imp{itS}.Range;
        LineShift = FitInfo.Imp{itS}.LineShift;
        TextShift = FitInfo.Imp{itS}.TextShift;

        for itN = 1:numel(LineShift)

            FitMin = FitRange(itN,1);
            FitMax = FitRange(itN,2);
            FitCenter = mean([FitMin,FitMax]);
            LShift = LineShift(itN);
            TShift = TextShift(itN,:);

            % linear fit
            idx = log10(ocont) > FitMin & log10(ocont) < FitMax;
            f = polyfit(log10(ocont(idx)), log10(ImpSusc{itS}(idx)), 1);
            Xfit = linspace(FitMin, FitMax, 10);
            Yfit = polyval(f, Xfit);
            Xfit = power(10, Xfit);     

            CritExp = round(abs(f(1)));  % critical exponent
            if CritExp > 0
                CritExp = CritExp * sign(f(1));
            end
            Yfit = power(10, Yfit) * LShift;
        
            % plot the linear fit result
            plot(Xfit, Yfit, '-', 'Color', 'black', 'LineWidth', 1);
            textX = power(10, FitCenter + TShift(1) );
            textY = power(10, polyval(f, FitCenter) + TShift(2) );
            if itS == 2 && itN == 2 && false
                text(textX, textY, ['$\omega^{',sprintf('%d',CritExp),'} \log^{-1}(\omega / T_{\mathrm{FO}}) $'], 'Color', 'black', 'Interpreter', 'latex', 'FontSize', FS);
            else
                text(textX, textY, ['$\omega^{',sprintf('%d',CritExp),'}$'], 'Color', 'black', 'Interpreter', 'latex', 'FontSize', FS);
            end

        end % itN
    end % itS

    hold off;



    %% plot subplot 3 ========================================================================

    ax3 = subplot(3,1,3);
    pos3 = [0.17, 0.08, 0.8, 0.22]; % [left bottom width height]
    set(ax3, 'Position', pos3);
    hold on;

    % set axis properties
    set(ax3, 'XScale', 'log', 'YScale', 'log', 'FontSize', 18);
    set(ax3, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
    ax3.XTick = 10.^( 5*ceil(log10(Tmin)/5) : 5 : 5*floor(log10(Tmax)/5) );
    ax3.YTick = 10.^( 5*ceil(log10(BathSuscMin)/5) : 5 : 5*floor(log10(BathSuscMax)/5) );
    ax3.TickLength = [0.015, 0.002];      % ticksize : [major, minor]
    set(ax3, 'XMinorTick', 'off', 'YMinorTick', 'off');
    set(ax3, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);
    ax3.XAxis.FontSize = 24;
    ax3.YAxis.FontSize = 24;
    set(ax3, 'Layer', 'top');
    xlim([Tmin, Tmax]);
    ylim([BathSuscMin, BathSuscMax]);

    % define x- and y-labels
    xlabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$\chi^{\mathrm{bath}}$', 'Interpreter', 'latex', 'FontSize', 20);

    hx = get(ax3, 'XLabel');
    hy = get(ax3, 'YLabel');
    hx.Units = 'normalized';
    hy.Units = 'normalized';
    
    % Modify x- and y-label positions manually
    hx.Position = hx.Position + [0, 0.25, 0];
    hy.Position = hy.Position + [0.06, 0, 0];

    % redefine x-tick labels using annotation()
    tickLabelFont = 18;
    ax3.XTickLabel = [];
    xticks = 10.^( 5*ceil(log10(Tmin)/5) : 5 : 5*floor(log10(Tmax)/5) );
    
    xtickWidth = 0.09;
    xtickHeight = 0.06;
    for itx = 1:numel(xticks)
        X_pos = pos3(1) + pos3(3) * (log10(xticks(itx)) - log10(Tmin)) / (log10(Tmax) - log10(Tmin));
        Y_pos = pos3(2);
        X_pos = X_pos - 0.05;
        Y_pos = Y_pos - 0.06;
        annotation('textbox', [X_pos, Y_pos, xtickWidth, xtickHeight], 'String', ['$',SciNot(xticks(itx)),'$'], 'Interpreter', 'latex', ...
                        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end

    % redefine y-tick labels using annotation()
    ax3.YTickLabel = [];
    yticks = 10.^( 5*ceil(log10(BathSuscMin)/5) : 5 : 5*floor(log10(BathSuscMax)/5) );
    
    ytickWidth = 0.09;
    ytickHeight = 0.1;
    for ity = 1:numel(yticks)
        X_pos = pos3(1) - 0.09;
        Y_pos = pos3(2) + pos3(4) * (log10(yticks(ity)) - log10(BathSuscMin)) / (log10(BathSuscMax) - log10(BathSuscMin));
        X_pos = X_pos + 0.005;
        Y_pos = Y_pos - 0.085;
        annotation('textbox', [X_pos, Y_pos, ytickWidth, ytickHeight], 'String', ['$',SciNot(yticks(ity)),'$'], 'Interpreter', 'latex', ...
                        'HorizontalAlignment', 'right', 'FitBoxToText', 'off', 'Units', 'normalized', 'LineStyle', 'none', 'FontSize', tickLabelFont);
    end
    
    % plot impurity dynamic susceptibilities
    plot(ocont, BathSusc{1}, 'Color', [.466, .674, .188], 'LineWidth', 2);
    plot(ocont, BathSusc{2}, 'Color', [.850, .325, .098], 'LineWidth', 2);
    plot(ocont, BathSusc{3}, 'Color', [0, .447, .741], 'LineWidth', 2);

    % Mark characteristic energy scales
    Scale_names = {'FO','SO','OO','FL'};
    
    for itN = 1:4
        if ~isnan(E_scale(itN))

            % identify energy scales and mark them on the plot
            xline_fig(E_scale(itN), '--');

            x_text = pos2(1) + pos2(3) * (log10(E_scale(itN)) - log10(Tmin)) / (log10(Tmax) - log10(Tmin)) - 0.03;
            y_text = 0.6;

            annotation('textbox', [x_text, y_text, 0.1, 0.1], ...
                    'String', ['$T_{\mathrm{',Scale_names{itN},'}}$'], 'Interpreter', 'latex', ...
                    'EdgeColor', 'none', ...          % remove box
                    'BackgroundColor', [1,1,1], ...
                    'FitBoxToText', 'on', ...
                    'FontSize', TempFontSize, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'middle');
        end
    end

    % linear fit and show power laws
    for itS = 1:3
        FitRange = FitInfo.Bath{itS}.Range;
        LineShift = FitInfo.Bath{itS}.LineShift;
        TextShift = FitInfo.Bath{itS}.TextShift;

        for itN = 1:numel(LineShift)

            FitMin = FitRange(itN,1);
            FitMax = FitRange(itN,2);
            FitCenter = mean([FitMin,FitMax]);
            LShift = LineShift(itN);
            TShift = TextShift(itN,:);

            % linear fit
            idx = log10(ocont) > FitMin & log10(ocont) < FitMax;
            f = polyfit(log10(ocont(idx)), log10(BathSusc{itS}(idx)), 1);
            Xfit = linspace(FitMin, FitMax, 10);
            Yfit = polyval(f, Xfit);
            Xfit = power(10, Xfit);     

            CritExp = round(abs(f(1)));  % critical exponent
            if CritExp > 0
                CritExp = CritExp * sign(f(1));
            end
            Yfit = power(10, Yfit) * LShift;
        
            % plot the linear fit result
            h = plot(Xfit, Yfit, '-', 'Color', 'black', 'LineWidth', 1);
            uistack(h, 'bottom');
            textX = power(10, FitCenter + TShift(1) );
            textY = power(10, polyval(f, FitCenter) + TShift(2) );
            if itS == 2 && itN == 2 && false
                text(textX, textY, ['$\omega^{',sprintf('%d',CritExp),'} \log^{-1}(\omega / T_{\mathrm{FO}}) $'], 'Color', 'black', 'Interpreter', 'latex', 'FontSize', FS);
            else
                text(textX, textY, ['$\omega^{',sprintf('%d',CritExp),'}$'], 'Color', 'black', 'Interpreter', 'latex', 'FontSize', FS);
            end

        end % itN
    end % itS

    fig = gcf;

    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'Units', 'centimeters');
    fig_pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', fig_pos(3:4));
    
    % Specify full path
    output_path = 'C:\Users\hsjun\OneDrive\Physics\Research\TsoK_publication\Figures';
    output_path = [output_path, '\Eflow_Susc_J0=', sprintf('%.15g',J0), '_K0=', sprintf('%.15g',K0), '.svg'];
    print(fig, output_path, '-dsvg');

    hold off;

    %{
    lgd = legend([hdl_sp, hdl_orb, hdl_sporb],{'$\chi_{\mathrm{sp}}$', '$\chi_{\mathrm{orb}}$', '$\chi_{\mathrm{sp-orb}}$'}, ...
                        'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 16, 'Box', 'on', 'Color', 'none');

    uistack(lgd, 'top');  % bring the legend to the top
    %}



    %% define 'HalfInt' function
    function Str = HalfInt(X)

        if rem(X*2,1) ~= 0
            error('ERR: Input of ''HalfInt'' must be a half integers');
        end

        if rem(X*2, 2) == 0
            Str = sprintf('%d',X);
        else
            Str = ['\frac{', sprintf('%d',2*X), '}{2}'];
        end
    end

    %% Define 'xline_fig' function
    function xline_fig(xval, style)
        % Add vertical line at xval that spans all subplots vertically
        % xval: x position in data units (log-scale supported)
        % style: linestyle, e.g., 'k--'
    
        fig = gcf;
        all_axes = findall(fig, 'type', 'axes');
        
        % Filter out legend/colorbar axes, keep only visible main ones
        all_axes = flipud(all_axes(arrayfun(@(ax) strcmp(ax.Tag,''), all_axes)));  % bottom to top order
    
        % Use first axes (bottom subplot) for x-normalization
        ax_ref = all_axes(1);
        xscale = get(ax_ref, 'XScale');
        axpos = get(ax_ref, 'Position');
        xlims = get(ax_ref, 'XLim');
    
        if strcmp(xscale, 'log')
            x_norm = (log10(xval) - log10(xlims(1))) / (log10(xlims(2)) - log10(xlims(1)));
        else
            x_norm = (xval - xlims(1)) / (xlims(2) - xlims(1));
        end
        x_fig = axpos(1) + x_norm * axpos(3);  % normalized figure x
    
        % Get tightest y-span across all axes
        y0 = inf; y1 = -inf;
        for ax3 = all_axes'
            pos3 = get(ax3, 'Position');  % [x y w h] in normalized units
            y0 = min(y0, pos3(2));
            y1 = max(y1, pos3(2) + pos3(4));
        end
    
        % Add vertical line
        annotation('line', [x_fig x_fig], [y0 y1], 'Color', 'k', ...
                   'LineStyle', style, 'LineWidth', 1.5);
    end

    %% Define signed_int_str()

    function s = signed_int_str(x)
        if x == 0
            s = '0';
        else
            s = sprintf('%+d', x);
        end
    end

end

