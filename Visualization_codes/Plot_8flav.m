clear;
Hyb = 1;
U = 5;
J = 0;
N0 = 2;
T = 1e-8;
mu = 0;
Nkeep = 7000;
Lambda = 6;
dataFolder = ['Hyb=',sprintf('%.15g',Hyb),'_U=',sprintf('%.15g',U),'_J=',sprintf('%.15g',J),'_N0=',sprintf('%.15g',N0), ...
                '_T=',sprintf('%.15g',T),'_mu=',sprintf('%.15g',mu),'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];

STG_path = ['C:\Users\hsjun\OneDrive\Physics\Research\data\8flav\',dataFolder];

plotEflow = true;
plotSusc = true;
plotEnt = false;

if plotEflow
    Etot = load([STG_path,'\Etot.mat']);
    Qtot = load([STG_path,'\Qtot.mat']);
    nu = load([STG_path,'\nu.mat']);
    Etot = Etot.Etot;   Qtot = Qtot.Qtot;   nu = nu.nu;

    disp(nu);
    disp(sum(nu));
    plotE(Etot,Qtot);
end

if plotSusc
    ocont = load([STG_path,'\ocont.mat']);
    ocont = ocont.ocont;
    
    %{
    Charge_susc = cell(1,4);
    tmp = load([STG_path,'\NRG_Op=Charge1m.mat']);
    Charge_susc{1} = tmp.temp;
    tmp = load([STG_path,'\NRG_Op=Charge1p.mat']);
    Charge_susc{2} = tmp.temp;
    tmp = load([STG_path,'\NRG_Op=Charge2m.mat']);
    Charge_susc{3} = tmp.temp;
    tmp = load([STG_path,'\NRG_Op=Charge2p.mat']);
    Charge_susc{4} = tmp.temp;
    %}
    
    Sp_susc = cell(1,1);
    tmp = load([STG_path,'\NRG_Op=ImpSp1p.mat']);
    Sp_susc{1} = tmp.temp;
    %tmp = load([STG_path,'\NRG_Op=ImpSp1m.mat']);
    %Sp_susc{2} = tmp.temp;
    %tmp = load([STG_path,'\NRG_Op=ImpSp2p.mat']);
    %Sp_susc{3} = tmp.temp;
    %%tmp = load([STG_path,'\NRG_Op=ImpSp2m.mat']);
    %Sp_susc{4} = tmp.temp;
    

    Lambda_susc = cell(1,2);
    tmp = load([STG_path,'\NRG_Op=Lambda_p.mat']);
    Lambda_susc{1} = tmp.temp;
    tmp = load([STG_path,'\NRG_Op=Lambda_z.mat']);
    Lambda_susc{2} = tmp.temp;
    
    figure;
    hold on;
    set(gca,'XScale','log','YScale','log');
    
    legends = cell(1,3);
    flav = {'1+', '1-', '2+', '2-'};
    linestyle = {'-', '--', '-.', ':'};
 
    %{
    for itC = 1:4
        plot(ocont, Charge_susc{itC}, linestyle{itC}, 'LineWidth', 1.5);
        legends{itC} = ['$\chi_{\mathrm{c_{',flav{itC},'}}}^{\mathrm{imp}}$'];
    end
    %}
    
    for itS = 1:1
        plot(ocont(ocont>0), Sp_susc{itS}(ocont>0), linestyle{itS}, 'LineWidth', 1.5);
        legends{itS} = ['$\chi_{\mathrm{sp_{',flav{itS},'}}}^{\mathrm{imp}}$'];
    end

    for itS = 1:2
        plot(ocont(ocont>0), Lambda_susc{itS}(ocont>0), linestyle{itS}, 'LineWidth', 1.5);
    end

    legends{2} = '$\chi_{\Lambda +}$';
    legends{3} = '$\chi_{\Lambda -}$';
    
    legend(legends,'Interpreter','latex','FontSize',18, 'Location', 'southeast');
    xlim([1e-8, 1e2]);
    ylim([1e-15, 1]);
    hold off;
end


if plotEnt
    EntData = load([STG_path,'\EntData.mat']);
    EntData = EntData.EntData;
    
    beta = 1;
    idx = find(EntData.beta == beta);

    figure;
    hold on;
    set(gca,'XScale','log','YScale','linear');
    plot(EntData.Temps{idx}, EntData.S_imp{idx}, 'LineWidth', 2);
    hold off;
end
