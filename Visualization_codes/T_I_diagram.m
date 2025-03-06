function T_I_diagram(J)
    
    if J == 0.0625
        Blue_x = [-10.^(-14:-11), 10.^(-14:-11)];
        Blue_y = -15.76*ones(1,8);
        Blue_y = power(10,Blue_y);
        
        Orange_x = [-10.^(-14:-3), 10.^(-14:-3)];
        Orange_y = -3.015*ones(1,24);
        Orange_y = power(10,Orange_y);
        
        Green_x = [-10.^(-6:-3), 10.^(-6:-3)];
        Green_y = [-15, -12.915, -10.68, -8.455];
        Green_y = [Green_y, Green_y];
        Green_y = power(10, Green_y);
        
        Red_x = [-10.^(-10:-3), 10.^(-10:-3)];
        Red_y = [-17.66, -15.74, -14.475, -13.165, -11.765, -10.27, -8.69, -6.975];
        Red_y = [Red_y, Red_y];
        Red_y = power(10,Red_y);
        
        Magenta_x = [-10.^(-10:-2), 10.^(-10:-2)];
        Magenta_y = [-16.575, -14.915, -13.725, -12.435, -11.045, -9.55, -7.975, -6.345, -4.655];
        Magenta_y = [Magenta_y, Magenta_y];
        Magenta_y = power(10, Magenta_y);
        
        X = {Blue_x, Orange_x, Green_x, Red_x, Magenta_x};
        Y = {Blue_y, Orange_y, Green_y, Red_y, Magenta_y};
        color = {'blue', '#EDB120', 'green', 'red', 'magenta'};
        figure;
        hold on;
        SymLogPlot(X, Y, color, 'xzerowidth', 0.07, 'YScale', 'log');
        xlabel('$I_{0}$','Interpreter','latex','FontSize',25);
        ylabel('$T$','Interpreter','latex','FontSize',25);
        title('$J_{0} = 0.0625,\ K_{0} = 0.3$','Interpreter','latex','FontSize',25);
        text(-0.05,1e-9,'NFL1','FontSize',15);
        text(-0.05,1e-17,'NFL2','FontSize',15);
        %text([-0.68, 0.58], [1e-13, 1e-13], 'NFL3?','FontSize',12);
        text([-0.85, 0.8], [1e-14, 1e-14], 'FL','FontSize',15);
        %x = [0.618, 0.618];
        %x = [0.788, 0.788];
        x = [0.764, 0.764];
        y = [0.9, 0.15];
        %annotation('arrow',x,y,'Linewidth',0.8);

        hold off;
    end


    if J == 0.1
        Blue_x = [-10.^(-14:-9), 10.^(-14:-9)];
        Blue_y = [-9.98*ones(1,4), -9.985, -9.99];
        Blue_y = [Blue_y, Blue_y];
        Blue_y = power(10,Blue_y);
        
        Orange_x = [-10.^(-14:-6), -10.^(-5:-3), 10.^(-14:-6), 10.^(-5:-3)];
        Orange_y = [-3.02*ones(1,18), -3.025*ones(1,6)];
        Orange_y = power(10,Orange_y);
        
        Green_x = [-10.^(-9:-3), 10.^(-9:-3)];
        Green_y = [-15.3, -14.12,-12.955, -12.075, -10.82, -9.295, -7.585];
        Green_y = [Green_y, Green_y];
        Green_y = power(10, Green_y);
        
        Red_x = [-10.^(-10:-3), 10.^(-10:-3)];
        Red_y = [-16.045, -14.335, -13.01, -11.665, -10.35, -9.085, -7.77, -6.27];
        Red_y = [Red_y, Red_y];
        Red_y = power(10,Red_y);
        
        Magenta_x = [-10.^(-8:-2), 10.^(-8:-2)];
        Magenta_y = [-11.8, -10.615, -9.485, -8.325, -7.07, -5.705, -4.245];
        Magenta_y = [Magenta_y, Magenta_y];
        Magenta_y = power(10, Magenta_y);

        Black_x = [-10.^(-10:-9), 10.^(-10:-9)];
        Black_y = [-14.695, -13.015];
        Black_y = [Black_y, Black_y];
        Black_y = power(10,Black_y);

        Gray_x = [-1e-2, 1e-2];
        Gray_y = -4.49*ones(1,2);
        Gray_y = power(10, Gray_y);


        X = {Blue_x, Orange_x, Green_x, Red_x, Magenta_x, Gray_x, Black_x};
        Y = {Blue_y, Orange_y, Green_y, Red_y, Magenta_y, Gray_y, Black_y};
        color = {'blue', '#EDB120', 'green', 'red', 'magenta', [.7, .7, .7], 'black'};
        figure;
        hold on;
        SymLogPlot(X, Y, color, 'xzerowidth', 0.07, 'YScale', 'log');
        xlabel('$I_{0}$','Interpreter','latex','FontSize',25);
        ylabel('$T$','Interpreter','latex','FontSize',25);
        title('$J_{0} = 0.1,\ K_{0} = 0.3$','Interpreter','latex','FontSize',25);
        text(-0.05,1e-6,'NFL1','FontSize',15);
        text(-0.05,1e-14,'NFL2','FontSize',15);
        %text([-0.68, 0.58], [1e-13, 1e-13], 'NFL3?','FontSize',12);
        text([-0.85, 0.8], [1e-14, 1e-14], 'FL','FontSize',15);
        %x = [0.618, 0.618];
        %x = [0.788, 0.788];
        %x = [0.764, 0.764];
        %y = [0.9, 0.15];
        %annotation('arrow',x,y,'Linewidth',0.8);

        hold off;
    end


    if J == 0.2
        Blue_x = 10.^(-14:-6);
        Blue_y = -4.98*ones(1,9);
        %Blue_y = [Blue_y, Blue_y];
        Blue_y = power(10,Blue_y);
        
        Orange_x = 10.^(-14:-4);
        Orange_y = -3.065*ones(1,11);
        Orange_y = power(10,Orange_y);
        
        Green_x = 10.^(-10:-4);
        Green_y = [-14.595, -13.13, -11.95, -10.76, -9.55, -8.34, -7.135];
        %Green_y = [Green_y, Green_y];
        Green_y = power(10, Green_y);
        
        Red_x = 10.^(-10:-4);
        Red_y = [-13.785, -12.295, -11.085, -9.87, -8.635, -7.395, -6.16];
        %Red_y = [Red_y, Red_y];
        Red_y = power(10,Red_y);
        
        Magenta_x = 10.^(-5:-2);
        Magenta_y = [-6.97, -5.775, -4.64, -3.49];
        %Magenta_y = [Magenta_y, Magenta_y];
        Magenta_y = power(10, Magenta_y);

        Black_x = 10.^(-10:-6);
        Black_y = [-13.405, -11.885, -10.655, -9.42, -8.185];
        Black_y = power(10,Black_y);

        Gray_x = 10.^(-3:-2);
        Gray_y = [-4.92, -3.065];
        Gray_y = power(10, Gray_y);


        X = {Blue_x, Orange_x, Green_x, Red_x, Magenta_x, Gray_x, Black_x};
        Y = {Blue_y, Orange_y, Green_y, Red_y, Magenta_y, Gray_y, Black_y};
        color = {'blue', '#EDB120', 'green', 'red', 'magenta', [.7, .7, .7], 'black'};
        figure;
        hold on;
        SymLogPlot(X, Y, color, 'xzerowidth', 0.07, 'YScale', 'log');
        xlabel('$I_{0}$','Interpreter','latex','FontSize',25);
        ylabel('$T$','Interpreter','latex','FontSize',25);
        title('$J_{0} = 0.2,\ K_{0} = 0.3$','Interpreter','latex','FontSize',25);
        text(-0.05,1e-4,'NFL1','FontSize',15);
        text(-0.05,1e-10,'NFL2','FontSize',15);
        %text([-0.68, 0.58], [1e-13, 1e-13], 'NFL3?','FontSize',12);
        text([0.75], [1e-13], 'FL','FontSize',15);
        %x = [0.618, 0.618];
        %x = [0.788, 0.788];
        %x = [0.764, 0.764];
        %y = [0.9, 0.15];
        %annotation('arrow',x,y,'Linewidth',0.8);

        hold off;
    end
end