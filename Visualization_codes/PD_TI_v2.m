function PD_TI_v2
    Unscreened_NFL1 = [
        0, -3.02;
        1e-8, -3.02;
        1e-7, -3.02;
        1e-6, -3.02;
        1e-5, -3.025;
        1e-4, -3.025;
        1e-3, -3.025;
        -1e-8, -3.02;
        -1e-7, -3.02;
        -1e-6, -3.02;
        -1e-5, -3.025;
        -1e-4, -3.025;
        -1e-3, -3.025;
    ];

    % avg(green, magenta)
    %{
    NFL1_FL = [
        1e-8, -12.96;
        1e-7, -11.79;
        1e-6, -10.78;
        1e-5, -9.57;
        1e-4, -8.18;
        1e-3, -6.65;
        -1e-8, -12.96;
        -1e-7, -11.79;
        -1e-6, -10.78;
        -1e-5, -9.57;
        -1e-4, -8.18;
        -1e-3, -6.65;
    ];
    %}

    % avg(red, magenta)
    %{}
    NFL1_FL = [
        1e-8, -12.41;
        1e-7, -11.14;
        1e-6, -9.92;
        1e-5, -8.71;
        1e-4, -7.42;
        1e-3, -5.99;
        -1e-8, -12.41;
        -1e-7, -11.14;
        -1e-6, -9.92;
        -1e-5, -8.71;
        -1e-4, -7.42;
        -1e-3, -5.99;
    ];
    %}

    NFL1_NFL2 = [
        0, -9.98;
    ];

    xl = [1e-7, 1];
    yl = [1e-13, 1e-2];
    

    Unscreened_NFL1(abs(Unscreened_NFL1(:,1))<xl(1),:) = [];
    NFL1_FL(abs(NFL1_FL(:,1))<xl(1),:) = [];

    X_bound = log(xl)/log(10);
    A1 = polyfit(Unscreened_NFL1(:,1), Unscreened_NFL1(:,2), 1);

    pos_Idx = NFL1_FL(:,1) > 0;
    A2 = polyfit(log(NFL1_FL(pos_Idx, 1))/log(10), NFL1_FL(pos_Idx, 2), 1);
    X2 = log(xl)/log(10);
    Y2 = power(10, polyval(A2,X2));
    X2 = power(10,X2);

    neg_Idx = NFL1_FL(:,1) < 0;
    A3 = polyfit(log(-NFL1_FL(neg_Idx, 1))/log(10), NFL1_FL(neg_Idx, 2), 1);
    X3 = log(xl)/log(10);
    Y3 = power(10, polyval(A3,X3));
    X3 = -power(10,X3);

    Intsec = [-power(10, (A1(2)-A3(2))/(A3(1)-A1(1))), power(10, (A1(2)-A2(2))/(A2(1)-A1(1))); ...
                    power(10, (A1(2)*A3(1)-A1(1)*A3(2))/(A3(1)-A1(1))), ...
                        power(10, (A1(2)*A2(1)-A1(1)*A2(2))/(A2(1)-A1(1)))];
    
    X1 = [Intsec(1,1), Intsec(1,2)];
    Y1 = power(10, polyval(A1,X1));
        
    figure;
    hold on;
    set(gcf,'position',[0,0,800,1000])

    ylim(yl);
    xlim([-1,1]);

    Xzerowidth = 0.05;
    XzeroRatio = 0.2;

    NFL1_X = [Intsec(1,1), Intsec(1,2), xl(1), xl(1), -xl(1), -xl(1), Intsec(1,1)]; ...
    NFL1_Y = [Intsec(2,1), Intsec(2,2), power(10, polyval(A2,log(xl(1))/log(10)) ), yl(1), yl(1),...
                power(10, polyval(A3,log(xl(1))/log(10)) ), Intsec(2,1)];
    NFL2_X = XzeroRatio*Xzerowidth*[-1,-1,1,1,-1];
    NFL2_Y = [yl(1), power(10, NFL1_NFL2(2)), power(10, NFL1_NFL2(2)), yl(1), yl(1)];

    FL1_X = [-xl(2), -xl(2), -power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), -xl(1), -xl(1), -xl(2)];
    FL1_Y = [yl(1), yl(2), yl(2), power(10, polyval(A3,log(xl(1))/log(10)) ), yl(1), yl(1)];
    FL2_X = [xl(2), xl(2), power(10, (log(yl(2))/log(10)-A2(2))/A2(1)), xl(1), xl(1), xl(2)];
    FL2_Y = [yl(1), yl(2), yl(2), power(10, polyval(A2,log(xl(1))/log(10)) ), yl(1), yl(1)];
    Unscreened_X = [Intsec(1,1), -power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), ...
                            power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), Intsec(1,2), Intsec(1,1)];
    Unscreened_Y = [Intsec(2,1), yl(2), yl(2), Intsec(2,2), Intsec(2,1)];

    [axes, lin2sym_X, sym2lin_X] = SLplot(X1, Y1, 'XScale', 'symlog', 'YScale', 'log', 'XexpLim', [log(xl(1))/log(10), log(xl(2))/log(10)], 'xzerowidth', Xzerowidth);


    NFL1_X = lin2sym_X(NFL1_X);
    FL1_X = lin2sym_X(FL1_X);
    FL2_X = lin2sym_X(FL2_X);
    Unscreened_X = lin2sym_X(Unscreened_X);

    FL1_X = [-1,-1,FL1_X(3:5),-1];
    FL2_X = [1,1,FL2_X(3:5),1];

    axes.XAxis.FontSize = 15;
    axes.YAxis.FontSize = 15;
    axes.XAxis.MinorTick = 'on';
    axes.Layer = 'top';
    title('$J_{0} = 0.1, K_{0} = 0.3$', 'interpreter', 'latex', 'fontsize', 21);
    xlabel(axes,'$I_{0}$', 'interpreter', 'latex', 'fontsize', 25);
    ylabel(axes,'$\mathrm{T}$', 'interpreter', 'latex', 'fontsize', 25);
        
    NFL1_h = patch(NFL1_X, NFL1_Y, [.737,.890,.996], 'FaceAlpha', 1, 'linestyle', 'none');
    FL_h = patch(FL1_X, FL1_Y, [.753,.878,.753], 'FaceAlpha', 1, 'linestyle', 'none');
    patch(FL2_X, FL2_Y, [.753,.878,.753], 'FaceAlpha', 1, 'linestyle', 'none');
    Unscreened_h = patch(Unscreened_X, Unscreened_Y, [.686,.694,.698], 'FaceAlpha', 1, 'linestyle', 'none');
    NFL2_h = patch(NFL2_X,NFL2_Y, [.745 .682 .898], 'FaceAlpha', 1, 'linestyle', 'none');

    patch([-Xzerowidth*[1,1], -XzeroRatio*Xzerowidth*[1,1], -Xzerowidth], ...
            [yl(1), yl(2)*[1,1], yl(1)*[1,1]], [.5,.5,.5], 'FaceAlpha', 0.5, 'linestyle', 'none');

    patch([XzeroRatio*Xzerowidth*[1,1], Xzerowidth*[1,1], XzeroRatio*Xzerowidth], ...
            [yl(1), yl(2)*[1,1], yl(1)*[1,1]], [.5,.5,.5], 'FaceAlpha', 0.5, 'linestyle', 'none');

    % [.988,.769,.8] : pink

    X1 = lin2sym_X(X1);
    X2 = lin2sym_X(X2);
    X3 = lin2sym_X(X3);
    Unscreened_NFL1_X = lin2sym_X(Unscreened_NFL1(:,1));
    NFL1_FL_X = lin2sym_X(NFL1_FL(:,1));

    Unscreened_NFL1_Y = power(10, Unscreened_NFL1(:,2));
    NFL1_FL_Y = power(10, NFL1_FL(:,2));

    %{
    plot(X1,Y1,'--','color','black');
    plot(X2,Y2,'--','color','black');
    plot(X3,Y3,'--','color','black');
    %}
    %{  
    plot(Unscreened_NFL1_X, Unscreened_NFL1_Y,'.','color',[.85, .325, .098],'MarkerSize',15);
    plot(NFL1_FL_X, NFL1_FL_Y,'.','color','blue','MarkerSize',15);
    %}

    handle = [NFL1_h, NFL2_h, FL_h, Unscreened_h];
    legends = {'$\mathrm{NFL}_{1}$', '$\mathrm{NFL}_{2}$', '$\mathrm{FL}$', '$\mathrm{Unscreened}$'};
    legend(handle, legends,'Interpreter','latex','Location','southeast','FontSize',15, 'AutoUpdate', 'off');

        
    hold off;

    saveas(gcf,fullfile('C:\Users\82104\Documents\Physics\Research\Figures','Red_Mag_Phase_Diagram_TI_J0=0.1_K0=0.3_v2.png'),'png');

end