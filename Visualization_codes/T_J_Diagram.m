function T_J_Diagram(I)

    if I == 0
        Blue_x = [0.2,0.1,0.09,0.08,0.07,0.0625,0.055,0.05,0.045,0.044,0.043,0.041,0.039,0.037,0.035,0.033];
        Blue_y = [-4.98, -9.98, -11.065, -12.405, -14.125,...
                    -15.77, -17.835, -19.56, -21.665, -22.165, -22.645, -23.715, -24.895, -26.2, -27.65, -29.275];
        x = Blue_x(2:end);
        y = Blue_y(2:end);
        f = polyfit(x, y, 1);

        yfit = polyval(f, x);          % Estimated  Regression Line
        SStot = sum((y-mean(y)).^2);                    % Total Sum-Of-Squares
        SSres = sum((y-yfit).^2);                       % Residual Sum-Of-Squares
        Rsq = 1- SSres/SStot;                            % R^2
        %disp(['R^2 : ',sprintf('%.15g',Rsq)]);

        Blue_y = power(10, Blue_y);

        Orange_x = [0.2,0.1,0.09,0.08,0.07,0.0625,0.055,0.05,0.045,0.044,0.043,0.041,0.04,0.039,0.037,0.035,0.033,0.01];
        Orange_y = [-3.065, -3.02, -3.02, -3.02, -3.015, -3.025, ...
                    -3.01, -3.01, -3.01, -3.015, -3.01, -3.01, -3.015, -3.01, -3.01, -3.01, -3.01, -3.02];
        Orange_y = power(10, Orange_y);

        figure;
        hold on;
        %plot(x, power(10,polyval(f,x)),'-','LineWidth',1,'Color',[.7, .7, .7]);
        %plot(Blue_x(2:end), power(9.03562*Blue_x(2:end), 23.9436),'-','LineWidth',1,'Color','Black');
        plot(Blue_x,Blue_y,'.','MarkerSize',15,'Color','blue');
        plot(Orange_x,Orange_y,'.','MarkerSize',15,'Color','#EDB120');
        set(gca,'XScale','linear','YScale','log','fontsize',20);
        xlabel('$J_{0}$','Interpreter','latex','FontSize',25);
        ylabel('$T$','Interpreter','latex','FontSize',25);
        title('$K_{0} = 0.3,\ I_{0} = 0$','Interpreter','latex','FontSize',25);
        %xlim([0, 0.22]);
        x = [0.412, 0.412];
        y = [0.9, 0.15];
        %annotation('arrow',x,y,'Linewidth',0.8);
        saveas(gcf,fullfile('/home/hyunsung/MyWork/Figures','T_I_J=0.0625.png'),'png');
        hold off;
    end
end