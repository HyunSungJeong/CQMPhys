clear;
order = 3;
J0 = 0.01;
K0 = 0.3;
I0 = 1e-1;
[J1,K1,I1,log_D] = RG_PoorMan(J0,K0,0,'order',order);
[J2,K2,I2,~] = RG_PoorMan(J0,K0,I0,'order',order);
figure;
hold on;
legend('AutoUpdate','on');
plot(log_D,J1,'Color','red','LineWidth',1,'LineStyle','-');
plot(log_D,K1,'Color','green','LineWidth',1,'LineStyle','-.');
plot(log_D,J2,'Color','red','LineWidth',1,'LineStyle','--');
plot(log_D,K2,'Color','green','LineWidth',1,'LineStyle','--');
plot(log_D,I2,'Color','blue','LineWidth',1,'LineStyle','--');
    
legends = {'$J: I_{0}=0$', '$K: I_{0}=0$', ['$J: I_{0}=',sprintf('%.15g',I0),'$'], ...
            ['$K: I_{0}=',sprintf('%.15g',I0),'$'], ['$I: I_{0}=',sprintf('%.15g',I0),'$']};
legend(legends,'Interpreter','latex','Location','northeast','FontSize',10);
set(gca,'XScale','linear','YScale','log','fontsize',20);
xlabel('$\log_{10} \frac{D}{D_{0}}$','Interpreter','latex','FontSize',25);
%ylabel('Parameter Values','Interpreter','latex','FontSize',25);
title(['RG flow by poor man''s scaling($J_{0}=',sprintf('%.15g',J0),',K_{0}=',sprintf('%.15g',K0),')$']...
            ,'Interpreter','latex','FontSize',15);