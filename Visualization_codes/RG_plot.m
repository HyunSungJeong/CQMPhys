clear;
order = 3;
J0 = 0.1;
K0 = 0.3;
I0 = 1e-6;
[J1,K1,I1,log_D] = RG_PoorMan(J0,K0,0,'order',order);
[J2,K2,I2,log_D] = RG_PoorMan(J0,K0,I0,'order',order);
figure;
hold on;
legend('AutoUpdate','on');
%{}
%plot(log_D/log(10),J1,'Color','red','LineWidth',1,'LineStyle','-');
%plot(log_D/log(10),K1,'Color','green','LineWidth',1,'LineStyle','-');
%plot(log_D/log(10),J2,'Color','red','LineWidth',1,'LineStyle','--');
%plot(log_D/log(10),K2,'Color','green','LineWidth',1,'LineStyle','--');
%plot(log_D/log(10),I2,'Color','blue','LineWidth',1,'LineStyle','--');
set(gca,'XScale','linear','YScale','log','fontsize',25);
%ylim([-0.5,5]);
%ylim([1e-2,1e+1]);
%}
%
X = {log_D/log(10), log_D/log(10), log_D/log(10), log_D/log(10), log_D/log(10)};
Y = {J1,K1,J2,K2,I2};
%Y = {nan,nan,J2,K2,I2};
color = {'red', 'green', 'red', 'green', 'blue'};
style = {'-','-.','-','-','-'};
%style = {'.','.','.','.','.'};
num = 5;
%SymLogPlot(X(3:num),Y(3:num),color(3:num),'XScale','linear','YScale','symlog','LineWidth',1,'style',style(3:end));
SymLogPlot(X(1:2),Y(1:2),color(1:2),'XScale','linear','YScale','symlog','LineWidth',1,'style',style(1:2));
%
legends = {'$J: I_{0}=0$', '$K: I_{0}=0$', ['$J: I_{0}=',sprintf('%.15g',I0),'$'], ...
            ['$K: I_{0}=',sprintf('%.15g',I0),'$'], ['$I: I_{0}=',sprintf('%.15g',I0),'$']};
%legend(legends(3:5),'Interpreter','latex','Location','best','FontSize',20);
legend(legends(1:2),'Interpreter','latex','Location','best','FontSize',20);

xlabel('$\log_{10} \frac{D}{D_{0}}$','Interpreter','latex','FontSize',25);
%ylabel('Parameter Values','Interpreter','latex','FontSize',25);
title(['RG flow by poor man''s scaling($J_{0}=',sprintf('%.15g',J0),',K_{0}=',sprintf('%.15g',K0),')$']...
            ,'Interpreter','latex','FontSize',15);