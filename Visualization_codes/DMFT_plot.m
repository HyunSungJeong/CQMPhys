% code to plot from calculated DMFT data
% type the following command to WSL terminal to syncronize data with SNU server
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/DMFT_QI_2CK_lat /mnt/c/Users/82104/Documents/Physics/Research/data

path = 'C:\Users\82104\Documents\Physics\Research\data\DMFT_QI_2CK_lat';

FileInfo = dir(path);

U = 2;
V = 0.5;
t_0 = 0.5;
phi_div_pi = 0.13;
T = 1e-10;
ndfix = 5;
Nfit = 27;
Nkeep = 3000;
path = [path,'\U=',sprintf('%.3g',U),'_V=',sprintf('%.3g',V),'_t_0=',sprintf('%.3g',t_0), ...
        '_phi_div_pi=',sprintf('%.3g',phi_div_pi),'_T=',sprintf('%.3g',T),'_ndfix=',sprintf('%.3g',ndfix), ...
            '_Nkeep=',sprintf('%d',Nkeep)];

tmp = load([path,'\ocont.mat']);
ocont = tmp.ocont;
tmp = load([path,'\RhoV2out.mat']);
RhoV2out = tmp.RhoV2out;
tmp = load([path,'\RhoV2in.mat']);
RhoV2in = tmp.RhoV2in;
tmp = load([path,'\SE_even.mat']);
SE_even = tmp.SE_save_e;
tmp = load([path,'\SE_odd.mat']);
SE_odd = tmp.SE_save_o;
tmp = load([path,'\Etots.mat']);
Etots = tmp.Etots;
tmp = load([path,'\Qtots.mat']);
Qtots = tmp.Qtots;
tmp = load([path,'\Qdiffs.mat']);
Qdiffs = tmp.Qdiffs;
tmp = load([path,'\ffs.mat']);
ffs = tmp.ffs;
tmp = load([path,'\ggs.mat']);
ggs = tmp.ggs;
tmp = load([path,'\dffs.mat']);
dffs = tmp.dffs;
tmp = load([path,'\dggs.mat']);
dggs = tmp.dggs;

Ndmft = 5;

plotE(Etots{Ndmft},Qtots{Ndmft},'title','$U(1)_{c} \times SU(2)_{sp}$','Emax',2,'legmax',13);

%% Self-energy: even sector
%{
figure;
[~,lin2sym_X,~,lin2sym_Y,~] = SLplot(ocont,-imag(SE_even(:,3,1,Ndmft)),'XScale','symlog','YScale','symlog','YexpLim',[-5,0]);
ocont_sym = lin2sym_X(ocont);
hold on;
legends = cell(1,9);
orbitals = {'1e', '2e', '3'};
linestyles = {'-','--','-.',':'};
for it1 = 1:3
    for it2 = 1:3
        plot(ocont_sym, -lin2sym_Y(imag(SE_even(:,it1,it2,Ndmft))),linestyles{rem(2*(it1-1)+it2,4)+1},'LineWidth',1);
        legends{3*(it1-1)+it2} = ['$\left( ',orbitals{it1},',',orbitals{it2},'\right)$'];
    end
end

set(gca,'fontsize',17);
xlabel('$\omega$','Interpreter','latex','FontSize',25);
ylabel('$\Sigma_{ij}^{\mathrm{even}} \left( \omega \right)$','Interpreter','latex','FontSize',25);
legend(legends,'Interpreter','latex','Location','eastoutside','FontSize',15);
title('$-\Im \left( \mathbf{\Sigma} \right) \mathrm{ : \ even \ sector}$','Interpreter','latex','FontSize',35);
hold off;
%}

%{
figure;
[~,lin2sym_X,~,lin2sym_Y,~] = SLplot(ocont,RhoV2out(:,1,1),'XScale','symlog','YScale','symlog');
ocont_sym = lin2sym_X(ocont);
hold on;
legends = cell(1,9);
orbitals = {'1e', '2e', '3'};
for it1 = 1:3
    for it2 = 1:3
        plot(ocont_sym, lin2sym_Y(real(SE_even(:,it1,it2,6))),'LineWidth',0.8);
        legends{3*(it1-1)+it2} = ['$\left( ',orbitals{it1},',',orbitals{it2},'\right)$'];
    end
end

set(gca,'fontsize',17);
xlabel('$\omega$','Interpreter','latex','FontSize',25);
ylabel('$\Sigma_{ij}^{\mathrm{even}} \left( \omega \right)$','Interpreter','latex','FontSize',25);
legend(legends,'Interpreter','latex','Location','eastoutside','FontSize',15);
title('$\Re \left( \mathbf{\Sigma} \right) \mathrm{ : \ even \ sector}$','Interpreter','latex','FontSize',35);
hold off;
%}

%{
figure;
hold on;
xlim([-15,15]);
legends = cell(1,9);
orbitals = {'1e', '2e', '3'};
for it1 = 1:3
    for it2 = 1:3
        plot(ocont, -imag(SE_even(:,it1,it2,Ndmft)),'LineWidth',0.8);
        legends{3*(it1-1)+it2} = ['$\left( ',orbitals{it1},',',orbitals{it2},'\right)$'];
    end
end

set(gca,'fontsize',17);
xlabel('$\omega$','Interpreter','latex','FontSize',25);
ylabel('$\Sigma_{ij}^{\mathrm{even}} \left( \omega \right)$','Interpreter','latex','FontSize',25);
legend(legends,'Interpreter','latex','Location','eastoutside','FontSize',15);
title('$-\Im \left( \mathbf{\Sigma} \right) \mathrm{ : \ even \ sector}$','Interpreter','latex','FontSize',35);
hold off;
%}

%{
figure;
hold on;
xlim([-15,15]);
legends = cell(1,9);
orbitals = {'1e', '2e', '3'};
for it1 = 1:3
    for it2 = 1:3
        plot(ocont, real(SE_even(:,it1,it2,Ndmft)),'LineWidth',0.8);
        legends{3*(it1-1)+it2} = ['$\left( ',orbitals{it1},',',orbitals{it2},'\right)$'];
    end
end

set(gca,'fontsize',17);
xlabel('$\omega$','Interpreter','latex','FontSize',25);
ylabel('$\Sigma_{ij}^{\mathrm{even}} \left( \omega \right)$','Interpreter','latex','FontSize',25);
legend(legends,'Interpreter','latex','Location','eastoutside','FontSize',15);
title('$\Re \left( \mathbf{\Sigma} \right) \mathrm{ : \ even \ sector}$','Interpreter','latex','FontSize',35);
hold off;
%}

%% Self-energy: odd sector

%{
figure;
[~,lin2sym_X,~,lin2sym_Y,~] = SLplot(ocont,-imag(SE_odd(:,2,1,Ndmft)),'XScale','symlog','YScale','symlog','YexpLim',[-5,0]);
ocont_sym = lin2sym_X(ocont);
hold on;
legends = cell(1,4);
orbitals = {'1o', '2o'};
linestyles = {'-','--','-.',':'};
for it1 = 1:2
    for it2 = 1:2
        plot(ocont_sym, -lin2sym_Y(imag(SE_odd(:,it1,it2,Ndmft))),linestyles{2*(it1-1)+it2},'LineWidth',1);
        legends{2*(it1-1)+it2} = ['$\left( ',orbitals{it1},',',orbitals{it2},'\right)$'];
    end
end

set(gca,'fontsize',17);
xlabel('$\omega$','Interpreter','latex','FontSize',25);
ylabel('$\Sigma_{ij}^{\mathrm{odd}} \left( \omega \right)$','Interpreter','latex','FontSize',25);
legend(legends,'Interpreter','latex','Location','eastoutside','FontSize',15);
title('$-\Im \left( \mathbf{\Sigma} \right) \mathrm{ : \ odd \ sector}$','Interpreter','latex','FontSize',35);
hold off;
%}

%{
figure;
[~,lin2sym_X,~,lin2sym_Y] = SLplot(ocont,RhoV2out(:,1,1),'XScale','symlog','YScale','symlog');
ocont_sym = lin2sym_X(ocont);
hold on;
legends = cell(1,4);
orbitals = {'1o', '2o'};
for it1 = 1:2
    for it2 = 1:2
        plot(ocont_sym, lin2sym_Y(real(SE_odd(:,it1,it2,6))),'LineWidth',0.8);
        legends{2*(it1-1)+it2} = ['$\left( ',orbitals{it1},',',orbitals{it2},'\right)$'];
    end
end

set(gca,'fontsize',17);
xlabel('$\omega$','Interpreter','latex','FontSize',25);
ylabel('$\Sigma_{ij}^{\mathrm{odd}} \left( \omega \right)$','Interpreter','latex','FontSize',25);
legend(legends,'Interpreter','latex','Location','eastoutside','FontSize',15);
title('$\Re \left( \mathbf{\Sigma} \right) \mathrm{ : \ odd \ sector}$','Interpreter','latex','FontSize',35);
hold off;
%}

%{
figure;
hold on;
xlim([-15,15]);
legends = cell(1,4);
orbitals = {'1o', '2o'};
for it1 = 1:2
    for it2 = 1:2
        plot(ocont, -imag(SE_odd(:,it1,it2,Ndmft)),'LineWidth',0.8);
        legends{2*(it1-1)+it2} = ['$\left( ',orbitals{it1},',',orbitals{it2},'\right)$'];
    end
end

set(gca,'fontsize',17);
xlabel('$\omega$','Interpreter','latex','FontSize',25);
ylabel('$\Sigma_{ij}^{\mathrm{odd}} \left( \omega \right)$','Interpreter','latex','FontSize',25);
legend(legends,'Interpreter','latex','Location','eastoutside','FontSize',15);
title('$-\Im \left( \mathbf{\Sigma} \right) \mathrm{ : \ odd \ sector}$','Interpreter','latex','FontSize',35);
hold off;
%}

%{
figure;
hold on;
xlim([-15,15]);
legends = cell(1,4);
orbitals = {'1o', '2o'};
for it1 = 1:2
    for it2 = 1:2
        plot(ocont, real(SE_odd(:,it1,it2,Ndmft)),'LineWidth',0.8);
        legends{2*(it1-1)+it2} = ['$\left( ',orbitals{it1},',',orbitals{it2},'\right)$'];
    end
end

set(gca,'fontsize',17);
xlabel('$\omega$','Interpreter','latex','FontSize',25);
ylabel('$\Sigma_{ij}^{\mathrm{odd}} \left( \omega \right)$','Interpreter','latex','FontSize',25);
legend(legends,'Interpreter','latex','Location','eastoutside','FontSize',15);
title('$\Re \left( \mathbf{\Sigma} \right) \mathrm{ : \ odd \ sector}$','Interpreter','latex','FontSize',35);
hold off;
%}

%%

%{
for itD = 2:4
    figure;
    [~,lin2sym_X,sym2lin_X] = SLplot(ocont,RhoV2out(:,1,1),'XScale','symlog','YScale','linear');
    ocont_sym = lin2sym_X(ocont);
    hold on;
    plot(ocont_sym, -imag(SE_even(:,1,1,itD)),'color','blue');
    plot(ocont_sym, -imag(SE_odd(:,1,1,itD)),'--','color','red');
    hold off;
end
%}

%{
figure;
hold on;
plot(ffs{4}{1}(:,1), 'linewidth', 1);
plot(ffs{4}{1}(:,2),'linewidth',1);
set(gca, 'XScale', 'linear', 'YScale', 'log', 'Fontsize', 20);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{f\,f}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{Intersite \ hoppings(ff)}$','Interpreter','latex','FontSize',35);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;
%}

%{
figure;
[~,lin2sym_X,sym2lin_X] = SLplot(ocont,RhoV2out(:,1,1),'XScale','symlog','YScale','log');
ocont_sym = lin2sym_X(ocont);
hold on;
set(gca, 'fontsize', 17);
plot(ocont_sym,RhoV2in(:,1,1),'-','color','blue','linewidth',1);
plot(ocont_sym,RhoV2in(:,2,1),'-','color','red','linewidth',1);
title('$\mathrm{Input \ hybridization \ function}$', 'Interpreter', 'latex', 'fontsize', 30);
xlabel('$\omega$', 'Interpreter', 'latex', 'fontsize', 25);
ylabel('$\Delta'''' \left( \omega \right)$', 'Interpreter', 'latex', 'fontsize', 25);
legend('$\mathrm{even \ sector}$', '$\mathrm{odd \ sector}$', 'Interpreter', 'latex', 'fontsize', 20, 'location', 'south');
hold off;
%}

%{
figure;
hold on;
plot(ffs{1}{1}(:,1), 'linewidth', 1);
plot(ffs{1}{1}(:,2),'linewidth',1);
set(gca, 'Fontsize', 20);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{f\,f}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{Intersite \ hoppings(ff)}$','Interpreter','latex','FontSize',35);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;


figure;
hold on;
plot(ggs{1}{1}(:,1), 'linewidth', 1);
plot(ggs{1}{1}(:,2),'linewidth',1);
set(gca, 'Fontsize', 20);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{g\,g}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{gg}$','Interpreter','latex','FontSize',35);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;


figure;
hold on;
plot(dffs{1}{1}(:,1), 'linewidth', 1);
plot(dffs{1}{1}(:,2),'linewidth',1);
set(gca, 'Fontsize', 20);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{d\,f\,f}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{dff}$','Interpreter','latex','FontSize',40);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;


figure;
hold on;
plot(dggs{1}{1}(:,1), 'linewidth', 1);
plot(dggs{1}{1}(:,2),'linewidth',1);
set(gca, 'Fontsize', 20);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{d\,g\,g}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{dgg}$','Interpreter','latex','FontSize',40);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;
%}


%{
figure;
hold on;
plot(ffs{1}{1}(:,1), 'linewidth', 1);
plot(ffs{1}{1}(:,2),'linewidth',1);
set(gca, 'XScale', 'linear', 'YScale', 'log', 'Fontsize', 20);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{f\,f}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{Intersite \ hoppings(ff)}$','Interpreter','latex','FontSize',35);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;


figure;
hold on;
[~, lin2sym_Y, sym2lin] = SLplot((1:numel(ggs{1}{1}(:,1))), ggs{1}{1}(:,1), 'XScale', 'linear', 'YScale', 'symlog');
gg_sym_e = lin2sym_Y(ggs{1}{1}(:,1));
gg_sym_o = lin2sym_Y(ggs{1}{1}(:,2));
set(gca, 'Fontsize', 15);
plot(gg_sym_e, 'linewidth', 1);
plot(gg_sym_o, 'linewidth', 1);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{g\,g}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{on-site \ energies(gg)}$','Interpreter','latex','FontSize',35);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;


figure;
hold on;
[~, lin2sym_Y, sym2lin] = SLplot((1:numel(dffs{1}{1}(:,1))), dffs{1}{1}(:,1), 'XScale', 'linear', 'YScale', 'symlog');
dff_sym_e = lin2sym_Y(dffs{1}{1}(:,1));
dff_sym_o = lin2sym_Y(dffs{1}{1}(:,2));
set(gca, 'Fontsize', 15);
plot(dff_sym_e, 'linewidth', 1);
plot(dff_sym_o, 'linewidth', 1);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{d\,f\,f}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{d\,f\,f}$','Interpreter','latex','FontSize',40);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;


figure;
hold on;
[~, lin2sym_Y, sym2lin] = SLplot((1:numel(dggs{1}{1}(:,1))), dggs{1}{1}(:,1), 'XScale', 'linear', 'YScale', 'symlog');
dgg_sym_e = lin2sym_Y(dggs{1}{1}(:,1));
dgg_sym_o = lin2sym_Y(dggs{1}{1}(:,2));
set(gca, 'Fontsize', 15);
plot(dgg_sym_e, 'linewidth', 1);
plot(dgg_sym_o, 'linewidth', 1);
xlabel('$N$', 'Interpreter', 'latex', 'Fontsize', 25);
ylabel('$\mathrm{d\,g\,g}$', 'Interpreter', 'latex', 'Fontsize', 25);
title('$\mathrm{d\,g\,g}$','Interpreter','latex','FontSize',40);
legend({'$\mathrm{even}$', '$\mathrm{odd}$'},'Interpreter','latex','Location','northeast','FontSize',25);
hold off;
%}


