clear;
J = [-0.01,-0.02,-0.05,-0.1,-0.15,-0.2,-0.25,-0.3];
K = -0.3*ones(1,8);
I = zeros(1,8);

J = [J,-0.05*ones(1,32)];
K = [K,-0.3*ones(1,32)];
I = [I,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7];
I = [I,-0.01,-0.03,-0.05,-0.07,-0.1,-0.3,-0.5,-0.7];
I = [I,10.^(-14:-7)];
I = [I,-10.^(-14:-7)];

J = [J,-0.1*ones(1,32)];
K = [K,-0.3*ones(1,32)];
I = [I,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7];
I = [I,-0.01,-0.03,-0.05,-0.07,-0.1,-0.3,-0.5,-0.7];
I = [I,10.^(-14:-7)];
I = [I,-10.^(-14:-7)];

J = [J,-0.2*ones(1,32)];
K = [K,-0.3*ones(1,32)];
I = [I,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7];
I = [I,-0.01,-0.03,-0.05,-0.07,-0.1,-0.3,-0.5,-0.7];
I = [I,10.^(-14:-7)];
I = [I,-10.^(-14:-7)];

figure;
hold on;
SymLogPlot({I,[0.7,-0.7]},{J,[-0.05,-0.05]},{[0,0,0],'red'},'style',{'.','*'},'XScale','symlog','YScale','linear','MarkerSize',10);
xlabel('$I_{0}$','Interpreter','latex','FontSize',25);
ylabel('$J_{0}$','Interpreter','latex','FontSize',25);
title('$K_{0} = -0.3$','Interpreter','latex','FontSize',25)
hold off;
