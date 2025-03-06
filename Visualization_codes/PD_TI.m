function PD_TI(J)
    if J == 0.0625

        Unscreened_NFL1 = [
            0, -3.015;
            1e-8, -3.015;
            1e-7, -3.015;
            1e-6, -3.015;
            1e-5, -3.015;
            1e-4, -3.015;
            1e-3, -3.015;
            1e-2, -3.015;
            -1e-8, -3.015;
            -1e-7, -3.015;
            -1e-6, -3.015;
            -1e-5, -3.015;
            -1e-4, -3.015;
            -1e-3, -3.015;
            -1e-2, -3.015;
        ];

        % avg(green, magenta)
        %{
        NFL1_FL = [
            1e-8, -14.93;
            1e-7, -14.48;
            1e-6, -13.02;
            1e-5, -11.23;
            1e-4, -9.33;
            1e-3, -7.40;
            1e-2, -5.49;
            1e-1, -2.88;
            -1e-8, -14.93;
            -1e-7, -14.48;
            -1e-6, -13.02;
            -1e-5, -11.23;
            -1e-4, -9.33;
            -1e-3, -7.40;
            -1e-2, -5.49;
            -1e-1, -2.88;
        ];
        %}

        % avg(red, magenta)
        %{}
        NFL1_FL = [
            1e-8, -14.10;
            1e-7, -12.80;
            1e-6, -11.41;
            1e-5, -9.91;
            1e-4, -8.33;
            1e-3, -6.66;
            1e-2, -4.81;
            1e-1, -2.87;
            -1e-8, -14.10;
            -1e-7, -12.80;
            -1e-6, -11.41;
            -1e-5, -9.91;
            -1e-4, -8.33;
            -1e-3, -6.66;
            -1e-2, -4.81;
            -1e-1, -2.87;
        ];
        %}

        %{}
        NFL1_NFL2 = [
            0, -15.76;
            1e-14, -15.76;
            1e-13, -15.76;
            1e-12, -15.76;
            1e-11, -15.76;
            -1e-14, -15.76;
            -1e-13, -15.76;
            -1e-12, -15.76;
            -1e-11, -15.76;
        ];
        %}

        xl = [1e-9, 1];
        yl = [1e-16, 1e-2];

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

        ylim(yl);
        xlim([-1,1]);

        Xzerowidth = 0.08;

        NFL1_X = [Intsec(1,1), Intsec(1,2), xl(1), xl(1), -xl(1), -xl(1), Intsec(1,1)]; ...
        NFL1_Y = [Intsec(2,1), Intsec(2,2), power(10, polyval(A2,log(xl(1))/log(10)) ), yl(1), yl(1),...
                    power(10, polyval(A3,log(xl(1))/log(10)) ), Intsec(2,1)];

        FL1_X = [-xl(2), -xl(2), -power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), -xl(1), -xl(1), -xl(2)];
        FL1_Y = [yl(1), yl(2), yl(2), power(10, polyval(A3,log(xl(1))/log(10)) ), yl(1), yl(1)];
        FL2_X = [xl(2), xl(2), power(10, (log(yl(2))/log(10)-A2(2))/A2(1)), xl(1), xl(1), xl(2)];
        FL2_Y = [yl(1), yl(2), yl(2), power(10, polyval(A2,log(xl(1))/log(10)) ), yl(1), yl(1)];
        Unscreened_X = [Intsec(1,1), -power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), ...
                                power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), Intsec(1,2), Intsec(1,1)];
        Unscreened_Y = [Intsec(2,1), yl(2), yl(2), Intsec(2,2), Intsec(2,1)];

        [axes, lin2sym_X, sym2lin_X] = SLplot(X1, Y1, 'YScale', 'log', 'XexpLim', [log(xl(1))/log(10), log(xl(2))/log(10)], 'xzerowidth', Xzerowidth);

        NFL1_X = lin2sym_X(NFL1_X);
        FL1_X = lin2sym_X(FL1_X);
        FL2_X = lin2sym_X(FL2_X);
        Unscreened_X = lin2sym_X(Unscreened_X);

        NFL2_X = 0.1*Xzerowidth*[-1,-1,1,1,-1];
        NFL2_Y = [yl(1), yl(2), yl(2), yl(1), yl(1)];

        FL1_X = [-1,-1,FL1_X(3:5),-1];
        FL2_X = [1,1,FL2_X(3:5),1];

        axes.XAxis.FontSize = 11;
        axes.YAxis.FontSize = 11;
        axes.XAxis.MinorTick = 'on';
        axes.Layer = 'top';
        title('$J_{0} = 0.0625, K_{0} = 0.3$', 'interpreter', 'latex', 'fontsize', 21);
        xlabel(axes,'$I_{0}$', 'interpreter', 'latex', 'fontsize', 16);
        ylabel(axes,'$\mathrm{T}$', 'interpreter', 'latex', 'fontsize', 16);
            
        NFL1_h = patch(NFL1_X, NFL1_Y, [.996,.734,.2], 'FaceAlpha', 0.7, 'linestyle', 'none');
        FL_h = patch(FL1_X, FL1_Y, [0,0.9,0.05], 'FaceAlpha', 0.7, 'linestyle', 'none');
        patch(FL2_X, FL2_Y, [0,0.9,0.05], 'FaceAlpha', 0.7, 'linestyle', 'none');
        Unscreened_h = patch(Unscreened_X, Unscreened_Y, [.7,.7,.7], 'FaceAlpha', 0.7, 'linestyle', 'none');

        patch([-Xzerowidth*[1,1], -0.1*Xzerowidth*[1,1], -Xzerowidth], ...
                [yl(1), yl(2)*[1,1], yl(1)*[1,1]], [.5,.5,.5], 'FaceAlpha', 0.5, 'linestyle', 'none');

        patch([0.1*Xzerowidth*[1,1], Xzerowidth*[1,1], 0.1*Xzerowidth], ...
                [yl(1), yl(2)*[1,1], yl(1)*[1,1]], [.5,.5,.5], 'FaceAlpha', 0.5, 'linestyle', 'none');

        % [.988,.769,.8] : pink

        X1 = lin2sym_X(X1);
        X2 = lin2sym_X(X2);
        X3 = lin2sym_X(X3);
        Unscreened_NFL1_X = lin2sym_X(Unscreened_NFL1(:,1));
        NFL1_FL_X = lin2sym_X(NFL1_FL(:,1));

        Unscreened_NFL1_Y = power(10, Unscreened_NFL1(:,2));
        NFL1_FL_Y = power(10, NFL1_FL(:,2));

        plot(X1,Y1,'--','color','black');
        plot(X2,Y2,'--','color','black');
        plot(X3,Y3,'--','color','black');
        plot(Unscreened_NFL1_X, Unscreened_NFL1_Y,'.','color',[.85, .325, .098],'MarkerSize',15);
        plot(NFL1_FL_X, NFL1_FL_Y,'.','color','blue','MarkerSize',15);

        handle = [NFL1_h, FL_h, Unscreened_h];
        legends = {'$NFL_{1}$', '$FL$', '$Unscreened$'};
        legend(handle, legends,'Interpreter','latex','Location','southeast','FontSize',15, 'AutoUpdate', 'off');

            
        hold off;

        saveas(gcf,fullfile('/home/hyunsung/MyWork/Figures','Red_Mag_Phase_Diagram_TI_J0=0.0625_K0=0.3.png'),'png');
    
    elseif J == 0.1
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
            -1e-10, -9.99
            -1e-9, -9.99;
            -2e-9, -9.99;
            -3e-9, -9.99;
            -4e-9, -9.99;
            -5e-9, -10.01;
            -8e-9, -10.01;
            1e-10, -9.99;
            1e-9, -9.99;
            2e-9, -10.00;
            3e-9, -9.99;
            4e-9, -9.99;
            5e-9, -10.01;
            6e-9, -10.01;
            7e-9, -10.02;
            8e-9, -10.05;
        ];

        NFL2_FL = [
            -1e-9, -13.675;
            -2e-9, -13.30;
            -3e-9, -13.07;
            -4e-9, -12.92;
            -5e-9, -12.79;
            -8e-9, -12.53;
            1e-9, -13.675;
            2e-9, -13.30;
            3e-9, -13.07;
            4e-9, -12.92;
            5e-9, -12.79;
            6e-9, -12.69;
            7e-9, -12.61;
            8e-9, -12.53;
        ];

        xl = [1e-9, 1];
        yl = [1e-16, 1e-2];

        Unscreened_NFL1(abs(Unscreened_NFL1(:,1))<xl(1),:) = [];
        NFL1_FL(abs(NFL1_FL(:,1))<xl(1),:) = [];
        NFL1_NFL2(abs(NFL1_NFL2(:,1))<xl(1),:) = [];

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

        ylim(yl);
        xlim([-1,1]);

        Xzerowidth = 0.08;

        NFL1_X = [Intsec(1,1), Intsec(1,2), xl(1), xl(1), -xl(1), -xl(1), Intsec(1,1)]; ...
        NFL1_Y = [Intsec(2,1), Intsec(2,2), power(10, polyval(A2,log(xl(1))/log(10)) ), yl(1), yl(1),...
                    power(10, polyval(A3,log(xl(1))/log(10)) ), Intsec(2,1)];
        FL1_X = [-xl(2), -xl(2), -power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), -xl(1), -xl(1), -xl(2)];
        FL1_Y = [yl(1), yl(2), yl(2), power(10, polyval(A3,log(xl(1))/log(10)) ), yl(1), yl(1)];
        FL2_X = [xl(2), xl(2), power(10, (log(yl(2))/log(10)-A2(2))/A2(1)), xl(1), xl(1), xl(2)];
        FL2_Y = [yl(1), yl(2), yl(2), power(10, polyval(A2,log(xl(1))/log(10)) ), yl(1), yl(1)];
        Unscreened_X = [Intsec(1,1), -power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), ...
                                power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), Intsec(1,2), Intsec(1,1)];
        Unscreened_Y = [Intsec(2,1), yl(2), yl(2), Intsec(2,2), Intsec(2,1)];

        [NFL2_X, Idx] = sort(NFL1_NFL2(:,1).', 'ascend');
        NFL2_Y = power(10, NFL1_NFL2(:,2).');
        NFL2_Y = NFL2_Y(Idx);
        [tmp_X,Idx] = sort(NFL2_FL(:,1).', 'descend');
        tmp_X_pos = tmp_X(tmp_X>0);
        tmp_X_neg = tmp_X(tmp_X<0);
        
        tmp_Y = power(10, NFL2_FL(:,2).');
        tmp_Y = tmp_Y(Idx);

        NFL2_X = [NFL2_X, tmp_X, NFL2_X(1)];
        NFL2_Y = [NFL2_Y, tmp_Y, NFL2_Y(1)];

        [axes, lin2sym_X, sym2lin_X] = SLplot(X1, Y1, 'YScale', 'log', 'XexpLim', [log(xl(1))/log(10), log(xl(2))/log(10)], 'xzerowidth', Xzerowidth);


        NFL1_X = lin2sym_X(NFL1_X);
        FL1_X = lin2sym_X(FL1_X);
        FL2_X = lin2sym_X(FL2_X);
        Unscreened_X = lin2sym_X(Unscreened_X);
        NFL2_X = lin2sym_X(NFL2_X);
        %disp(NFL2_X);
        %disp(log(NFL2_Y)/log(10));

        FL1_X = [-1,-1,FL1_X(3:5),-1];
        FL2_X = [1,1,FL2_X(3:5),1];

        axes.XAxis.FontSize = 11;
        axes.YAxis.FontSize = 11;
        axes.XAxis.MinorTick = 'on';
        axes.Layer = 'top';
        title('$J_{0} = 0.1, K_{0} = 0.3$', 'interpreter', 'latex', 'fontsize', 21);
        xlabel(axes,'$I_{0}$', 'Interpreter', 'latex', 'fontsize', 16);
        ylabel(axes,'$\mathrm{T}$', 'Interpreter', 'latex', 'fontsize', 16);
            
        NFL1_h = patch(NFL1_X, NFL1_Y, [.996,.734,.2], 'FaceAlpha', 0.7, 'linestyle', 'none');
        FL_h = patch(FL1_X, FL1_Y, [0,0.9,0.05], 'FaceAlpha', 0.7, 'linestyle', 'none');
        patch(FL2_X, FL2_Y, [0,0.9,0.05], 'FaceAlpha', 0.7, 'linestyle', 'none');
        Unscreened_h = patch(Unscreened_X, Unscreened_Y, [.7,.7,.7], 'FaceAlpha', 0.7, 'linestyle', 'none');
        NFL2_h = patch(NFL2_X, NFL2_Y, [.235 .000 .392], 'FaceAlpha', 0.5, 'linestyle', 'none');


        patch([-Xzerowidth*[1,1], -0.1*Xzerowidth*[1,1], -Xzerowidth], ...
                [yl(1), yl(2)*[1,1], yl(1)*[1,1]], [.5,.5,.5], 'FaceAlpha', 0.5, 'linestyle', 'none');

        patch([0.1*Xzerowidth*[1,1], Xzerowidth*[1,1], 0.1*Xzerowidth], ...
                [yl(1), yl(2)*[1,1], yl(1)*[1,1]], [.5,.5,.5], 'FaceAlpha', 0.5, 'linestyle', 'none');

        %patch([-0.1*Xzerowidth*[1,1], 0.1*Xzerowidth*[1,1], -0.1*Xzerowidth], [yl(1), power(10, NFL1_NFL2(1,2))*[1,1], yl(1)*[1,1]], ...
        %        [.235 .000 .392], 'FaceAlpha', 0.5, 'linestyle', 'none');

        % [.988,.769,.8] : pink

        X1 = lin2sym_X(X1);
        X2 = lin2sym_X(X2);
        X3 = lin2sym_X(X3);
        Unscreened_NFL1_X = lin2sym_X(Unscreened_NFL1(:,1));
        NFL1_FL_X = lin2sym_X(NFL1_FL(:,1));

        Unscreened_NFL1_Y = power(10, Unscreened_NFL1(:,2));
        NFL1_FL_Y = power(10, NFL1_FL(:,2));

        NFL1_NFL2_X = lin2sym_X(NFL1_NFL2(:,1));
        NFL1_NFL2_Y = power(10, NFL1_NFL2(:,2));

        NFL2_FL_X = lin2sym_X(NFL2_FL(:,1));
        NFL2_FL_Y = power(10, NFL2_FL(:,2));

        %plot(X1,Y1,'--','color','black');
        %plot(X2,Y2,'--','color','black');
        %plot(X3,Y3,'--','color','black');
        plot(Unscreened_NFL1_X, Unscreened_NFL1_Y,'.','color',[.85, .325, .098],'MarkerSize',15);
        plot(NFL1_FL_X, NFL1_FL_Y,'.','color','blue','MarkerSize',15);
        plot(NFL1_NFL2_X, NFL1_NFL2_Y,'.','color',[.494 .184 .556],'MarkerSize',15);
        plot(NFL2_FL_X, NFL2_FL_Y,'.','color','black','MarkerSize',15);

        handle = [NFL1_h, NFL2_h, FL_h, Unscreened_h];
        legends = {'$\mathrm{NFL}_{1}$', '$\mathrm{NFL}_{2}$', '$\mathrm{FL}$', '$\mathrm{Unscreened}$'};
        legend(handle, legends,'Interpreter','latex','Location','southeast','FontSize',15, 'AutoUpdate', 'off');

            
        hold off;

        saveas(gcf,fullfile('/home/hyunsung/MyWork/Figures','Red_Mag_Phase_Diagram_TI_J0=0.1_K0=0.3.png'),'png');

    elseif J == 0.2
        Unscreened_NFL1 = [
            0, -3.065;
            1e-8, -3.065;
            1e-7, -3.065;
            1e-6, -3.065;
            1e-5, -3.065;
            1e-4, -3.065;
            -1e-8, -3.065;
            -1e-7, -3.065;
            -1e-6, -3.065;
            -1e-5, -3.065;
            -1e-4, -3.065;
        ];

        % avg(green, magenta)
        %{
        NFL1_FL = [
            
        ];
        %}

        % avg(red, magenta)
        %{}
        NFL1_FL = [
            1e-5, -7.18;
            1e-4, -5.97;
            -1e-5, -7.18;
            -1e-4, -5.97;
        ];
        %}

        NFL1_NFL2 = [
            0, -4.98;
            -1e-10, -4.98;
            -1e-9, -4.98;
            -1e-8, -4.985;
            -1e-7, -4.985;
            -1e-6, -4.985;
            1e-10, -4.98;
            1e-9, -4.98;
            1e-8, -4.985;
            1e-7, -4.985;
            1e-6, -4.985;
            2e-6, -4.985;
            3e-6, -4.985;
            4e-6, -4.985;
            5e-6, -4.985;
            6e-6, -4.985;
            7e-6, -4.985;
            8e-6, -5.00;
            9e-6, -5.005;
        ];

        NFL2_FL = [
            -1e-10, -13.405;
            -1e-9, -11.885;
            -1e-8, -10.655;
            -1e-7, -9.42;
            -1e-6, -8.185;
            1e-10, -13.405;
            1e-9, -11.885;
            1e-8, -10.655;
            1e-7, -9.42;
            1e-6, -8.185;
            2e-6, -7.815;
            3e-6, -7.605;
            4e-6, -7.45;
            5e-6, -7.335;
            6e-6, -7.24;
            7e-6, -7.16;
            8e-6, -7.085;
            9e-6, -7.025;
        ];

        xl = [1e-10, 1];
        yl = [1e-16, 1e-2];

        Unscreened_NFL1(abs(Unscreened_NFL1(:,1))<xl(1),:) = [];
        NFL1_FL(abs(NFL1_FL(:,1))<xl(1),:) = [];
        NFL1_NFL2(abs(NFL1_NFL2(:,1))<xl(1),:) = [];

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

        ylim(yl);
        xlim([-1,1]);

        Xzerowidth = 0.08;

        NFL1_X = [Intsec(1,1), Intsec(1,2), xl(1), xl(1), -xl(1), -xl(1), Intsec(1,1)]; ...
        NFL1_Y = [Intsec(2,1), Intsec(2,2), power(10, polyval(A2,log(xl(1))/log(10)) ), yl(1), yl(1),...
                    power(10, polyval(A3,log(xl(1))/log(10)) ), Intsec(2,1)];
        FL1_X = [-xl(2), -xl(2), -power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), -xl(1), -xl(1), -xl(2)];
        FL1_Y = [yl(1), yl(2), yl(2), power(10, polyval(A3,log(xl(1))/log(10)) ), yl(1), yl(1)];
        FL2_X = [xl(2), xl(2), power(10, (log(yl(2))/log(10)-A2(2))/A2(1)), xl(1), xl(1), xl(2)];
        FL2_Y = [yl(1), yl(2), yl(2), power(10, polyval(A2,log(xl(1))/log(10)) ), yl(1), yl(1)];
        Unscreened_X = [Intsec(1,1), -power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), ...
                                power(10, (log(yl(2))/log(10)-A3(2))/A3(1)), Intsec(1,2), Intsec(1,1)];
        Unscreened_Y = [Intsec(2,1), yl(2), yl(2), Intsec(2,2), Intsec(2,1)];

        [NFL2_X, Idx] = sort(NFL1_NFL2(:,1).', 'ascend');
        NFL2_Y = power(10, NFL1_NFL2(:,2).');
        NFL2_Y = NFL2_Y(Idx);
        [tmp_X,Idx] = sort(NFL2_FL(:,1).', 'descend');
        tmp_X_pos = tmp_X(tmp_X>0);
        tmp_X_neg = tmp_X(tmp_X<0);
        
        tmp_Y = power(10, NFL2_FL(:,2).');
        tmp_Y = tmp_Y(Idx);

        NFL2_X = [NFL2_X, tmp_X, NFL2_X(1)];
        NFL2_Y = [NFL2_Y, tmp_Y, NFL2_Y(1)];

        [axes, lin2sym_X, sym2lin_X] = SLplot(X1, Y1, 'YScale', 'log', 'XexpLim', [log(xl(1))/log(10), log(xl(2))/log(10)], 'xzerowidth', Xzerowidth);


        NFL1_X = lin2sym_X(NFL1_X);
        FL1_X = lin2sym_X(FL1_X);
        FL2_X = lin2sym_X(FL2_X);
        Unscreened_X = lin2sym_X(Unscreened_X);
        NFL2_X = lin2sym_X(NFL2_X);
        %disp(NFL2_X);
        %disp(log(NFL2_Y)/log(10));

        FL1_X = [-1,-1,FL1_X(3:5),-1];
        FL2_X = [1,1,FL2_X(3:5),1];

        axes.XAxis.FontSize = 11;
        axes.YAxis.FontSize = 11;
        axes.XAxis.MinorTick = 'on';
        axes.Layer = 'top';
        title('$J_{0} = 0.2, K_{0} = 0.3$', 'interpreter', 'latex', 'fontsize', 21);
        xlabel(axes,'$I_{0}$', 'Interpreter', 'latex', 'fontsize', 16);
        ylabel(axes,'$\mathrm{T}$', 'Interpreter', 'latex', 'fontsize', 16);
            
        %NFL1_h = patch(NFL1_X, NFL1_Y,  [1,.718,.2], 'FaceAlpha', 0.7, 'linestyle', 'none');
        FL_h = patch(FL1_X, FL1_Y, [1,.753,.796], 'FaceAlpha', 0.7, 'linestyle', 'none');
        patch(FL2_X, FL2_Y, [1,.753,.796], 'FaceAlpha', 0.7, 'linestyle', 'none');
        Unscreened_h = patch(Unscreened_X, Unscreened_Y, [0,0.9,0.05], 'FaceAlpha', 0.7, 'linestyle', 'none');
        NFL2_h = patch(NFL2_X, NFL2_Y, [.576,.710,.949], 'FaceAlpha', 0.5, 'linestyle', 'none');


        patch([-Xzerowidth*[1,1], -0.1*Xzerowidth*[1,1], -Xzerowidth], ...
                [yl(1), yl(2)*[1,1], yl(1)*[1,1]], [.5,.5,.5], 'FaceAlpha', 0.5, 'linestyle', 'none');

        patch([0.1*Xzerowidth*[1,1], Xzerowidth*[1,1], 0.1*Xzerowidth], ...
                [yl(1), yl(2)*[1,1], yl(1)*[1,1]], [.5,.5,.5], 'FaceAlpha', 0.5, 'linestyle', 'none');

        %patch([-0.1*Xzerowidth*[1,1], 0.1*Xzerowidth*[1,1], -0.1*Xzerowidth], [yl(1), power(10, NFL1_NFL2(1,2))*[1,1], yl(1)*[1,1]], ...
        %        [.235 .000 .392], 'FaceAlpha', 0.5, 'linestyle', 'none');

        % [.988,.769,.8] : pink

        X1 = lin2sym_X(X1);
        X2 = lin2sym_X(X2);
        X3 = lin2sym_X(X3);
        Unscreened_NFL1_X = lin2sym_X(Unscreened_NFL1(:,1));
        NFL1_FL_X = lin2sym_X(NFL1_FL(:,1));

        Unscreened_NFL1_Y = power(10, Unscreened_NFL1(:,2));
        NFL1_FL_Y = power(10, NFL1_FL(:,2));

        NFL1_NFL2_X = lin2sym_X(NFL1_NFL2(:,1));
        NFL1_NFL2_Y = power(10, NFL1_NFL2(:,2));

        NFL2_FL_X = lin2sym_X(NFL2_FL(:,1));
        NFL2_FL_Y = power(10, NFL2_FL(:,2));

        %plot(X1,Y1,'--','color','black');
        %plot(X2,Y2,'--','color','black');
        %plot(X3,Y3,'--','color','black');
        plot(Unscreened_NFL1_X, Unscreened_NFL1_Y,'.','color',[.85, .325, .098],'MarkerSize',15);
        plot(NFL1_FL_X, NFL1_FL_Y,'.','color','blue','MarkerSize',15);
        plot(NFL1_NFL2_X, NFL1_NFL2_Y,'.','color',[.494 .184 .556],'MarkerSize',15);
        plot(NFL2_FL_X, NFL2_FL_Y,'.','color','black','MarkerSize',15);

        %handle = [NFL2_h, FL_h, Unscreened_h];
        handle = [NFL1_h, NFL2_h, FL_h, Unscreened_h];
        legends = {'$\mathrm{NFL}_{1}$', '$\mathrm{NFL}_{2}$', '$\mathrm{FL}$', '$\mathrm{Unscreened}$'};
        legend(handle, legends,'Interpreter','latex','Location','southeast','FontSize',15, 'AutoUpdate', 'off');

            
        hold off;

        saveas(gcf,fullfile('/home/hyunsung/MyWork/Figures','Red_Mag_Phase_Diagram_TI_J0=0.2_K0=0.3.png'),'png');
    end

end