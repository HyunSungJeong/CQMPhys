clear;
figure;
hold on;
annotation('arrow',[0.17,0.88],[0.355,0.355]);
annotation('arrow',[0.517,0.517],[0.3,0.75]);
text(1.2,-0.03,'$\omega$','Interpreter','latex','FontSize',20);
text(0.02,0.48,'$\Delta(\omega)$','Interpreter','latex','FontSize',20);
plot([-1,-1],[0,0.3],'Color','blue','LineWidth',2);
plot([-1,1],[0.3,0.3],'Color','blue','LineWidth',2);
plot([1,1],[0.3,0],'Color','blue','LineWidth',2);
text(-1.03,-0.05,'-1','Interpreter','latex','FontSize',18);
text(0.99,-0.05,'1','Interpreter','latex','FontSize',18);
xlim([-1.3,1.3]);
ylim([-0.3,0.7]);

Lambda = 1.5;
N1 = 5;
N2 = 7;
sft_pos = 0.02;
sft_neg = 0.05;
for it = (1:N2)
    w_disc = power(Lambda,-it);
    plot([w_disc, w_disc], [0, 0.35],'red','LineWidth',1.5);
    plot([-w_disc, -w_disc], [0, 0.35], 'red','LineWidth',1.5);
    if it < N1
        text(w_disc-sft_pos,-0.05,['$\Lambda^{',sprintf('%d',-it),'}$'],'Interpreter','latex','FontSize',15);
        text(-w_disc-sft_neg,-0.05,['$-\Lambda^{',sprintf('%d',-it),'}$'],'Interpreter','latex','FontSize',15);
    end
end
text(-0.034,0.15,'...','Interpreter','latex','FontSize',20);
text(-0.025,-0.05,'...','Interpreter','latex','FontSize',14);