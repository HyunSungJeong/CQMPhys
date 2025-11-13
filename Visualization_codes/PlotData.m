% code to plot from calculated data
% type the following command to WSL terminal to syncronize data with SNU server
% (To sync TsoK data)
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/TsoK /mnt/c/Users/hsjun/OneDrive/Physics/Research/data
% (To sync TsoK_Aniso data)
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/TsoK_Aniso /mnt/c/Users/hsjun/OneDrive/Physics/Research/data
% (To sync TCK_Aniso data)
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/TCK_Aniso /mnt/c/Users/hsjun/OneDrive/Physics/Research/data
% (To sync Kondo_Aniso data)
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/Kondo_Aniso /mnt/c/Users/hsjun/OneDrive/Physics/Research/data
% (To sync ThsoK data)
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/ThsoK /mnt/c/Users/hsjun/OneDrive/Physics/Research/data
% (To sync Quartic data)
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/Quartic /mnt/c/Users/hsjun/OneDrive/Physics/Research/data
% (To sync 8flav data)
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/8flav /mnt/c/Users/hsjun/OneDrive/Physics/Research/data
% (To sync LineSearch data)
% rsync -avz --delete --progress -e 'ssh -p 1018' hyunsung@147.46.44.239:/data/hyunsung/LineSearch /mnt/c/Users/hsjun/OneDrive/Physics/Research/data

%clear;
%% Choose Calculation Type
strtmp = cell(9,1);
strtmp{1} = '1: TsoK_NRG';
strtmp{2} = '2: TsoK_Aniso_NRG';
strtmp{3} = '3: TCK_Aniso_NRG';
strtmp{4} = '4: Kondo_Aniso_NRG';
strtmp{5} = '5: ThsoK_NRG';
strtmp{6} = '6: Quartic_NRG';
strtmp{7} = '7: Anderson_8flav_NRG';
strtmp{8} = '8: 2oAH_NRG';
strtmp{9} = '9: LineSearch';
dispbox('-width',75,strtmp{:});
fprintf('Choose index of calculation type to plot\n');
intmp = input('>>> ');
path = [];

while isempty(path)
    if intmp == 1
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK';
    elseif intmp == 2
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK_Aniso';
    elseif intmp == 3
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TCK_Aniso';
    elseif intmp == 4
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\Kondo_Aniso';
    elseif intmp == 5
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\ThsoK';
    elseif intmp == 6
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\Quartic';
    elseif intmp == 7
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\8flav';
    elseif intmp == 8
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\ToAH';    
    elseif intmp == 9
        path = 'C:\Users\hsjun\OneDrive\Physics\Research\data\LineSearch';
    else
        fprintf('WRN: Invalid input\n');
        intmp = input('>>> ');
    end
end

%%
if intmp == 1       % TsoK_NRG

    Lambda = 2.5;

    FileInfo = dir(path);
    strtmp = cell(0,0);
    cnt = 1;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1) == 'J' || DirName(1) == 'T'
            strtmp = cat(2,strtmp,[sprintf('%.15g',cnt),': ',DirName]);
            cnt = cnt + 1;
        end
    end
    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;

    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    if isequal(strtmp{idx}(numel(num2str(idx))+3 : numel(num2str(idx))+6), 'TEST')
        tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': TEST_Nk=%f_J0=%f_K0=%f_I0=%f_T=%f']);
        Nkeep = tmp(1);
        J0 = tmp(2);
        K0 = tmp(3);
        I0 = tmp(4);
        T = tmp(5);

        tmp = strtmp{idx};
        path = [path, filesep, tmp(numel(sprintf('%d',idx))+3:end)];
    else
        tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': J0=%f_K0=%f_I0=%f_T=%f']);
        if isempty(tmp)
        end
        J0 = tmp(1);
        K0 = tmp(2);
        I0 = tmp(3);
        T = tmp(4);
        
        tmp = strtmp{idx};
        path = [path, filesep, tmp(numel(sprintf('%d',idx))+3:end)];
    end

    FileInfo = dir(path);
    EraseIdx = [];                     % indices of unecessary file infos to be erased
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')         % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')    % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName(1:4),'TsoK')     % Erase log files
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];       % Erase selected unecessary file infos

    Eflow = cell(2,1);      % RG flow data. 1: Etot, 2: Qtot
    ImpDyn = cell(3,1);     % Impurity dynamical susceptibilities. 1: ImpSp, 2: ImpOrb, 3: ImpSpOrb
    BathDyn = cell(3,1);    % Bath dynamical susceptibilities. 1: BathSp, 2: BathOrb, 3: BathSpOrb
    avail = false(1,9);     
    % RGflow, ImpDyn, Second Derivatives of ImpDyn, BathDyn, Second Derivatives of BathDyn,
    % Entropy, spin-spin correlators(inner product), spin-spin correlators(square product)

    for it = (1:numel(FileInfo))
        tmp = load([FileInfo(it).folder,filesep,FileInfo(it).name]);
        field = fieldnames(tmp);
        
        switch FileInfo(it).name
            case 'Etot.mat'
                Eflow{1} = getfield(tmp,field{1});
                avail(1) = true;
            case 'Qtot.mat'
                Eflow{2} = getfield(tmp,field{1}); 
                avail(1) = true;
            case 'ocont.mat'
                ocont = getfield(tmp,field{1});
            case 'NRG_Op=ImpSp.mat'
                ImpDyn{1} = getfield(tmp,field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpOrb.mat'
                ImpDyn{2} = getfield(tmp,field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpSpOrb.mat'
                ImpDyn{3} = getfield(tmp,field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=BathSp.mat'
                BathDyn{1} = getfield(tmp,field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathOrb.mat'
                BathDyn{2} = getfield(tmp,field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathSpOrb.mat'
                BathDyn{3} = getfield(tmp,field{1});
                avail(4) = true;
                avail(5) = true;
            case 'Sent_imp.mat'
                Sent_imp = getfield(tmp,field{1});
                avail(6) = true;
            case 'Temps.mat'
                Temps = getfield(tmp,field{1});
            case 'Spin_correlators.mat'
                Sp_corr = getfield(tmp,field{1});
                avail(7) = true;
            case 'Orbital_correlators.mat'
                Orb_corr = getfield(tmp,field{1});
                avail(8) = true;
            case 'SpOrb_correlators.mat'
                SpOrb_corr = getfield(tmp,field{1});
                avail(9) = true;    
            otherwise
                fprintf('WRN: unknown data type');
        end

    end

    strtmp = cell(1,0);
    options = [];
    if avail(1)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': RG flow diagram']});
        options = [options,1];
    end
    if avail(2)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': Impurity dynamic susceptibilities']});
        options = [options,2];
    end
    if avail(3)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': Second derivatives of impurity dynamic susceptibilities']});
        options = [options,3];
    end
    if avail(4)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': Bath dynamic susceptibilities']});
        options = [options,4];
    end
    if avail(5)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': Second derivatives of bath dynamic susceptibilities']});
        options = [options,5];
    end
    if avail(6)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': Impurity contribution to entropy']});
        options = [options,6];
    end
    if avail(7)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': Spin correlators']});
        options = [options,7];
    end
    if avail(8)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': Orbital correlators']});
        options = [options,8];
    end
    if avail(9)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': SpOrb correlators']});
        options = [options,9];
    end
    strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': All of the above']});
    dispbox('-width',130,strtmp{:});

    fprintf('Which one do you want to plot?\n')
    intmp = input('>>> ');

    if intmp == numel(options)+1   % plot all
        chosen = avail;
    else
        chosen = false(1,9);
        chosen(options(intmp)) = true;
        while ~isempty(intmp) && ~isequal(chosen, avail)
            fprintf('Type the index of data you want to plot\n');
            fprintf('(Press enter to finish)\n');
            intmp = input('>>> ');
            if intmp == numel(options)+1
                chosen = avail;
            else
                chosen(options(intmp)) = true;
            end
        end
    end

    if avail(1) && chosen(1)
        plotE(Eflow{1},Eflow{2},'title',['J0=',sprintf('%.15g',J0),' K0=',sprintf('%.15g',K0),' I0=',sprintf('%.15g',I0),' T=',sprintf('%.15g',T)],'Emax',3,'legmax',13,'Qdiff',[0,0,0]);
    end
    
    if avail(2) && avail(3)
        
        names = {'sp', 'orb', 'sp-orb'};
        legends = cell(0,0);
        for it = (1:3)
            if isempty(ImpDyn{it})
                ImpDyn(it) = [];
            else
                legends = cat(2,legends,names{it});
            end
        end

        
        num_ImpDyn = numel(ImpDyn);         % Number of impurity dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % Log temperatures
        log_ImpDyn = cell(num_ImpDyn,1);            % Log impurity dynamic susceptibilities
        log_ImpDyn_1stDer = cell(num_ImpDyn,1);     % First derivatives
        log_T_1stDer = cell(num_ImpDyn,1);          % Log Temperatures for fisrt derivatives
        log_ImpDyn_2ndDer = cell(num_ImpDyn,1);     % Second derivatives
        log_T_2ndDer = cell(num_ImpDyn,1);          % Log Temperatures for second derivatives

        for it = (1:num_ImpDyn)
            tmp = ImpDyn{it};
            log_ImpDyn{it} = log(tmp(ocont>0))./log(10);

            log_ImpDyn_1stDer{it} = diff(log_ImpDyn{it},1)./diff(log_T,1);                  % First derivative
            tmp = movmean(log_T, [0,1]);
            log_T_1stDer{it} = tmp(1:end-1);                                                % Log Temperatures for fisrt derivatives

            log_ImpDyn_2ndDer{it} = diff(log_ImpDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % Second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                % Log Temperatures for second derivatives
        end

        if chosen(2)

            logT_grid = log10(ocont(ocont>0));      % log temperature grid
            log_Susc = cell(1,3);
            log_Susc{1} = log10(ImpDyn{1}(ocont>0));
            log_Susc{2} = log10(ImpDyn{2}(ocont>0));
            log_Susc{3} = log10(ImpDyn{2}(ocont>0));
    
            log_Sp_Susc_1stDer = diff(log_Susc{1},1)./diff(logT_grid,1);      % first derivative
            log_Orb_Susc_1stDer = diff(log_Susc{2},1)./diff(logT_grid,1);
            log_SpOrb_Susc_1stDer = diff(log_Susc{3},1)./diff(logT_grid,1);
            logT_grid_1stDer = movmean(logT_grid, [0,1]);           % log temperature grid for first derivative
            logT_grid_1stDer = logT_grid_1stDer(1:end-1);

            figure;
            hold on;
            ylim([-1.5,1]);
            plot(logT_grid_1stDer, log_Sp_Susc_1stDer,'LineWidth',2);
            plot(logT_grid_1stDer, log_Orb_Susc_1stDer,'LineWidth',2);
            plot(logT_grid_1stDer, log_SpOrb_Susc_1stDer, '--', 'LineWidth',2);
            hold off;
    
            figure;
            hold on;
            xlim([1e-12,1]);
            linestyle = {'-', '-', '-'};
            for it = 1:num_ImpDyn
                plot(ocont(ocont < 1),ImpDyn{it}(ocont < 1),'Linewidth',2,'LineStyle',linestyle{it});
            end
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Location','best','FontSize',25);

            legend('AutoUpdate','off');
            fit_range = [0, -3; -3.5, -5.5; -3.5, -5.5];

            if num_ImpDyn < 3
                log_ImpDyn = cat(1,log_ImpDyn,log_ImpDyn{1});
            end

            [a1,Rsq1,a2,Rsq2,a3,Rsq3] = Insert(log_T, log_ImpDyn,fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2 - 1;
            text_y = polyval(a1,text_x) - 2.2;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1) - 0.3;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            plot(x1,y1,'-','Color',[0,0,0],'LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',15);

            x2 = fit_range(2,:);
            text_x = (x2(1)+x2(2))/2 - 1;
            text_y = polyval(a2,text_x) - 1.1;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y2 = polyval(a2,x2) - 0.5;
            x2 = power(10,x2);
            y2 = power(10,y2);
            text2 = ['$\omega^{',sprintf('%.2f',a2(1)),'}$'];
            plot(x2,y2,'-','Color',[0,0,0],'LineWidth',1);
            text(text_x, text_y, text2,'Interpreter','latex','FontSize',15);
            
            x3 = fit_range(3,:);
            text_x = (x3(1)+x3(2))/2;
            text_y = polyval(a3,text_x) + 1;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y3 = polyval(a3,x3) + 0.5;
            x3 = power(10,x3);
            y3 = power(10,y3);
            text3 = ['$\omega^{',sprintf('%.2f',a3(1)),'}$'];
            plot(x3,y3,'-','Color',[0,0,0],'LineWidth',1);
            text(text_x, text_y, text3,'Interpreter','latex','FontSize',15);
            legend('AutoUpdate','on');

            log_ImpDyn(3) = [];

            set(gca,'XScale','log','YScale','log','fontsize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi^{\mathrm{imp}} (\omega)$','Interpreter','latex','FontSize',25);
            title(['Impurity Dynamic Susceptibilities (J0,K0,I0) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=',sprintf('%.15g',T)],'FontSize',20);
            hold off;
        end

        if chosen(3)
            figure;
            hold on;
            linestyle = {'-', '-.', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_ImpDyn)
                X = log_T_2ndDer{it};
                Y = log_ImpDyn_2ndDer{it};
                plot(X, Y,'LineWidth',2,'LineStyle',linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_ImpDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([log(T)/log(10)-1,3]);
            ylim([-3,4]);
            legend(legends,'Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['2nd Derivatives of Impurity Dynamic Susceptibilities (J0,K0,I0) = (',sprintf('%.15g',J0), ...
                        ', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=',sprintf('%.15g',T)],'FontSize',15);
            hold off;
        end

    end

    if avail(4) && avail(5)

        names = {'Bath Spin','Bath Orb','Bath SpOrb'};
        legends = cell(0,0);
        for it = (1:3)
            if isempty(BathDyn{it})
                BathDyn(it) = [];
            else
                legends = cat(2,legends,names{it});
            end
        end

        
        num_BathDyn = numel(BathDyn);         % Number of bath dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % Log temperatures
        log_BathDyn = cell(num_BathDyn,1);          % Log bath dynamic susceptibilities
        log_BathDyn_1stDer = cell(num_BathDyn,1);   % First derivatives
        log_T_1stDer = cell(num_BathDyn,1);         % Log Temperatures for fisrt derivatives
        log_BathDyn_2ndDer = cell(num_BathDyn,1);   % Second derivatives
        log_T_2ndDer = cell(num_BathDyn,1);         % Log Temperatures for second derivatives

        for it = (1:num_BathDyn)
            tmp = BathDyn{it};
            log_BathDyn{it} = log(tmp(ocont>0))./log(10);

            log_BathDyn_1stDer{it} = diff(log_BathDyn{it},1)./diff(log_T,1);                  % First derivative
            tmp = movmean(log_T, [0,1]);
            log_T_1stDer{it} = tmp(1:end-1);                                                % Log Temperatures for fisrt derivatives

            log_BathDyn_2ndDer{it} = diff(log_BathDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % Second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                % Log Temperatures for second derivatives
        end
        
        if chosen(4)

            %{}
            logT_grid = log10(ocont(ocont>0));      % log temperature grid
            log_Susc = cell(1,2);
            log_Susc{1} = log10(BathDyn{1}(ocont>0));
            log_Susc{2} = log10(BathDyn{2}(ocont>0));
    
            log_Sp_Susc_1stDer = diff(log_Susc{1},1)./diff(logT_grid,1);      % first derivative
            log_Orb_Susc_1stDer = diff(log_Susc{2},1)./diff(logT_grid,1);
            logT_grid_1stDer = movmean(logT_grid, [0,1]);           % log temperature grid for first derivative
            logT_grid_1stDer = logT_grid_1stDer(1:end-1);

            figure;
            hold on;
            ylim([-1.5,1]);
            plot(logT_grid_1stDer, log_Sp_Susc_1stDer,'LineWidth',2);
            plot(logT_grid_1stDer, log_Orb_Susc_1stDer,'LineWidth',2);
            plot([-25,0],[-1,-1],'--','Color','black','LineWidth',2);
            plot([-25,0],[-1.1,-1.1],'--','Color','black','LineWidth',2);
            plot([-25,0],[-1.2,-1.2],'--','Color','black','LineWidth',2);
            hold off;
            %}

            figure;
            hold on;
            for it = (1:num_BathDyn)
                plot(ocont(ocont < 1),BathDyn{it}(ocont < 1),'Linewidth',2);
            end
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Location','northeast','FontSize',25);
            set(gca,'XScale','log','YScale','log','fontsize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi'''' (\omega)$','Interpreter','latex','FontSize',25);
            title(['Bath Dynamic Susceptibilities (J0,K0,I0) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=',sprintf('%.15g',T)],'FontSize',20);

            
            hold off;
        end

        if chosen(5)
            figure;
            hold on;
            linestyle = {'-', '-', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_BathDyn)
                plot(log_T_2ndDer{it}, log_BathDyn_2ndDer{it}, 'LineWidth', 2, linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_BathDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([-23,3]);
            ylim([-3,4]);
            legend(legends,'Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['2nd Derivatives of Bath Dynamic Susceptibilities (J0,K0,I0) = (',sprintf('%.15g',J0), ...
                        ', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=',sprintf('%.15g',T)],'FontSize',15);
            hold off;
        end

    end

    if avail(6) && chosen(6)

        figure('Position', [100, 100, 800, 350]);
        hold on;

        set(gca, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 18);
        set(gca, 'XTick', 10.^(-15:5:0));
        set(gca, 'FontSize', 24);
        set(gca, 'LineWidth', 1);  % make axis lines (incl. ticks) bold
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
        set(gca, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

        plot(Temps,exp(Sent_imp),'Linewidth',2);
        xlim([1e-18,1]);
        ylim([1,4.1])
        plot([1e-22,1],[sqrt(8),sqrt(8)],'--','LineWidth',1,'color',[.7,.7,.7]);
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'tex';
        set(gca, 'TickLabelInterpreter', 'latex');
        ticks = {'1','2','$2 \times \sqrt{2}$','4'};
        yticklabels(ticks);
        yticks([1,2,sqrt(8),4]);

        ax = gca;
        xlabel('$T$','Interpreter','latex','FontSize',27);
        ylabel('$\mathrm{exp}(S_{\mathrm{imp}})$','Interpreter','latex','FontSize',27);
        %title(['Impurity contribution to entropy (J0,K0,I0)= (', ...
        %            sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=',sprintf('%.15g',T)],'FontSize',20);
        hold off;
    end
    
    if avail(7) && chosen(7)

        Sp_odd = Sp_corr(1:2:end);
        Sp_even = Sp_corr(2:2:end);
        X = (1:numel(Sp_corr));
        X_odd = (1:2:numel(Sp_corr));
        X_even = (2:2:numel(Sp_corr));
        X = Lambda.^(-X/2);
        X_odd = Lambda.^(-X_odd/2);
        X_even = Lambda.^(-X_even/2);
        
        Sp_avg = nan(1,numel(Sp_corr)-2);
        for it = 2:numel(Sp_corr)-1
            Sp_avg(it-1) = (2*Sp_corr(it) + Sp_corr(it+1) + Sp_corr(it-1)) / 4 ;
        end

        disp(sum(Sp_corr));

        figure;
        hold on;
        %{
        sgtitle(['$(J_{0},K_{0},I_{0}) = (',sprintf('%.15g',J0),', ',sprintf('%.15g',K0),...
                    ', ',SciNot(I0),')$'],'Interpreter','latex','FontSize',30);
        %}
        sgtitle(['$(J_{0},K_{0},I_{0}) = (',sprintf('%.15g',J0),', ',SciNot(K0),...
                    ', ',SciNot(I0),')$'],'Interpreter','latex','FontSize',30);

        %sgtitle(['$J_{0} = ',SciNot(J0),'$'],'Interpreter','latex','FontSize',30);
        % subplot 1
        ax1 = subplot(2,1,1);
        hold on;
        ylim([1.1*min(Sp_corr(2:end)), 1.1*max(Sp_corr(2:end))]);
        xlim([1e-20,X(1)]);
        set(gca,'XScale','log','YScale','linear');

        plot(X_odd,Sp_odd,'-','LineWidth',1.5);
        plot(X_even,Sp_even,'-.','LineWidth',1.5);
        plot(X(2:end-1),Sp_avg,'-.','LineWidth',1.5);
        plot(Lambda.^(1:numel(Sp_corr)), zeros(1,numel(Sp_corr)), '--', 'Color', [.5,.5,.5]);
        
        ax=gca;
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;   
        set(ax1, 'XMinorTick', 'off', 'YMinorTick', 'off');
        xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',20);
        ylabel('$\langle \vec{S} \cdot \vec{J}^{m}_{\mathrm{sp}} \rangle$', 'Interpreter','latex','FontSize',20);
        legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'}, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 20);

        % subplot 2
        ax2 = subplot(2,1,2);
        hold on;
        YexpLim = [-10,1];
        [~,lin2sym_Y,sym2lin_Y] = SLplot(X, Sp_corr, 'XScale', 'log', 'YScale', 'symlog', 'YexpLim', YexpLim, 'yzerowidth', 0.08);

        Sp_odd_SL = lin2sym_Y(Sp_odd);
        Sp_even_SL = lin2sym_Y(Sp_even);
        Sp_avg_SL = lin2sym_Y(Sp_avg);

        if abs(min(Sp_corr)) <= 10^YexpLim(1) || abs(max(Sp_corr)) <= 10^YexpLim(1)
            ylim_SL = lin2sym_Y([-10,10]);
        else
            ylim_SL = lin2sym_Y([min(Sp_corr)*5, max(Sp_corr)*5]);
        end

        Yexplim_SL = lin2sym_Y([-10^YexpLim(1), 10^YexpLim(1)]);

        xlim([1e-20,X(1)]);
        ylim(ylim_SL);

        plot(X_odd,Sp_odd_SL,'-','LineWidth',1.5);
        plot(X_even,Sp_even_SL,'-.','LineWidth',1.5);
        plot(X(2:end-1),Sp_avg_SL,'-.','LineWidth',1.5);
        plot(Lambda.^(1:numel(Sp_corr)), zeros(1,numel(Sp_corr)), '--', 'Color', [.5,.5,.5]);

        patch([1e-20, X(1), X(1), 1e-20], [repmat(Yexplim_SL(1), [1,2]), repmat(Yexplim_SL(2), [1,2])], [.7 .7 .7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        
        ax=gca;
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;     
        set(ax2, 'XMinorTick', 'off', 'YMinorTick', 'off');
        xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',20);
        ylabel('$\langle \vec{S} \cdot \vec{J}^{m}_{\mathrm{sp}} \rangle$', 'Interpreter','latex','FontSize',20);
        legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'}, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 20);
        hold off;
    end

    if avail(8) && chosen(8)

        Orb_odd = Orb_corr(1:2:end);
        Orb_even = Orb_corr(2:2:end);
        X = (1:numel(Orb_corr));
        X_odd = (1:2:numel(Orb_corr));
        X_even = (2:2:numel(Orb_corr));
        X = Lambda.^(-X/2);
        X_odd = Lambda.^(-X_odd/2);
        X_even = Lambda.^(-X_even/2);

        Orb_avg = nan(1,numel(Orb_corr)-2);
        for it = 2:numel(Orb_corr)-1
            Orb_avg(it-1) = (2*Orb_corr(it) + Orb_corr(it+1) + Orb_corr(it-1)) / 4 ;
        end

        disp(sum(Orb_corr));

        figure;
        hold on;
        %{
        sgtitle(['$(J_{0},K_{0},I_{0}) = (',sprintf('%.15g',J0),', ',sprintf('%.15g',K0),...
                    ', ',SciNot(I0),')$'],'Interpreter','latex','FontSize',30);
        %}
        sgtitle(['$(J_{0},K_{0},I_{0}) = (',sprintf('%.15g',J0),', ',SciNot(K0),...
                    ', ',SciNot(I0),')$'],'Interpreter','latex','FontSize',30);

        % subplot 1
        ax1 = subplot(2,1,1);
        hold on;
        ylim([-0.07, 0.01]);
        xlim([1e-20,X(1)]);
        set(gca,'XScale','log','YScale','linear');

        plot(X_odd,Orb_odd,'-','LineWidth',1.5);
        plot(X_even,Orb_even,'-.','LineWidth',1.5);
        plot(X(2:end-1),Orb_avg,'-.','LineWidth',1.5);
        Diff = interp1(X_odd,Orb_odd,X(2:end-1),'linear','extrap') - interp1(X_even,Orb_even,X(2:end-1),'linear','extrap');
        %plot(X(2:end-1),Diff,'-.','LineWidth',1.5);
        plot(Lambda.^(1:numel(Orb_corr)), zeros(1,numel(Orb_corr)), '--', 'Color', [.5,.5,.5]);
        
        ax=gca;
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;   
        set(ax1, 'XMinorTick', 'off', 'YMinorTick', 'off');
        xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',20);
        ylabel('$\langle \vec{T} \cdot \vec{J}^{m}_{\mathrm{orb}} \rangle$', 'Interpreter','latex','FontSize',20);
        legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'}, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 20);

        % subplot 2
        ax2 = subplot(2,1,2);
        hold on;
        YexpLim = [-16,1];
        [~,lin2sym_Y,sym2lin_Y] = SLplot(X, Orb_corr, 'XScale', 'log', 'YScale', 'symlog', 'YexpLim', YexpLim, 'yzerowidth', 0.08);

        Orb_odd_SL = lin2sym_Y(Orb_odd);
        Orb_even_SL = lin2sym_Y(Orb_even);
        Orb_avg_SL = lin2sym_Y(Orb_avg);
        Orb_diff_SL = lin2sym_Y(Diff);
        
        if abs(min(Orb_corr)) <= 10^YexpLim(1) || abs(max(Orb_corr)) <= 10^YexpLim(1)
            ylim_SL = lin2sym_Y([-10,10]);
        else
            ylim_SL = lin2sym_Y([min(Orb_corr)*5, max(Orb_corr)*5]);
        end

        Yexplim_SL = lin2sym_Y([-10^YexpLim(1), 10^YexpLim(1)]);

        xlim([1e-20,1]);
        ylim(ylim_SL);

        plot(X_odd,Orb_odd_SL,'-','LineWidth',2);
        plot(X_even,Orb_even_SL,'-.','LineWidth',2);
        plot(X(2:end-1),Orb_avg_SL,'-.','LineWidth',2);
        %plot(X(2:end-1),Orb_diff_SL,'-.','LineWidth',2);
        plot(Lambda.^(1:numel(Orb_corr)), zeros(1,numel(Orb_corr)), '--', 'Color', [.5,.5,.5]);

        patch([1e-20, X(1), X(1), 1e-20], [repmat(Yexplim_SL(1), [1,2]), repmat(Yexplim_SL(2), [1,2])], [.7 .7 .7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        
        ax=gca;
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;     
        set(ax2, 'XMinorTick', 'off', 'YMinorTick', 'off');
        xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',20);
        ylabel('$\langle \vec{T} \cdot \vec{J}^{m}_{\mathrm{orb}} \rangle$', 'Interpreter','latex','FontSize',20);
        legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'}, 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 20);
        hold off;
    end

    if avail(9) && chosen(9)
        
        SpOrb_odd = SpOrb_corr(1:2:end);
        SpOrb_even = SpOrb_corr(2:2:end);
        X = (1:numel(SpOrb_corr));
        X_odd = (1:2:numel(SpOrb_corr));
        X_even = (2:2:numel(SpOrb_corr));
        X = Lambda.^(-X/2);
        X_odd = Lambda.^(-X_odd/2);
        X_even = Lambda.^(-X_even/2);

        SpOrb_avg = nan(1,numel(SpOrb_corr)-2);
        for it = 2:numel(SpOrb_corr)-1
            SpOrb_avg(it-1) = (2*SpOrb_corr(it) + SpOrb_corr(it+1) + SpOrb_corr(it-1)) / 4 ;
        end

        disp(sum(SpOrb_corr));

        figure;
        hold on;
        %{
        sgtitle(['$(J_{0},K_{0},I_{0}) = (',sprintf('%.15g',J0),', ',sprintf('%.15g',K0),...
                    ', ',SciNot(I0),')$'],'Interpreter','latex','FontSize',30);
        %}
        sgtitle(['$(J_{0},K_{0},I_{0}) = (',sprintf('%.15g',J0),', ',SciNot(K0),...
                    ', ',SciNot(I0),')$'],'Interpreter','latex','FontSize',30);

        % subplot 1
        ax1 = subplot(2,1,1);
        hold on;
        ylim([-0.07, 0.01]);
        xlim([1e-20,X(1)]);
        set(gca,'XScale','log','YScale','linear');

        plot(X_odd,SpOrb_odd,'-','LineWidth',1.5);
        plot(X_even,SpOrb_even,'-.','LineWidth',1.5);
        plot(X(2:end-1),SpOrb_avg,'-.','LineWidth',1.5);
        plot(Lambda.^(1:numel(SpOrb_corr)), zeros(1,numel(SpOrb_corr)), '--', 'Color', [.5,.5,.5]);
        
        ax=gca;
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;   
        set(ax1, 'XMinorTick', 'off', 'YMinorTick', 'off');
        xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',20);
        ylabel('$\langle \vec{S} \cdot \vec{J}^{m}_{\mathrm{sp-orb}} \cdot \vec{T} \rangle$', 'Interpreter','latex','FontSize',20);
        legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'}, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 20);

        % subplot 2
        ax2 = subplot(2,1,2);
        hold on;
        YexpLim = [-7,1];
        [~,lin2sym_Y,sym2lin_Y] = SLplot(X, SpOrb_corr, 'XScale', 'log', 'YScale', 'symlog', 'YexpLim', YexpLim, 'yzerowidth', 0.08);

        SpOrb_odd_SL = lin2sym_Y(SpOrb_odd);
        SpOrb_even_SL = lin2sym_Y(SpOrb_even);
        SpOrb_avg_SL = lin2sym_Y(SpOrb_avg);

        if abs(min(SpOrb_corr)) <= 10^YexpLim(1) || abs(max(SpOrb_corr)) <= 10^YexpLim(1)
            ylim_SL = lin2sym_Y([-10,10]);
        else
            ylim_SL = lin2sym_Y([min(SpOrb_corr)*5, max(SpOrb_corr)*5]);
        end

        Yexplim_SL = lin2sym_Y([-10^YexpLim(1), 10^YexpLim(1)]);

        xlim([1e-20,X(1)]);
        ylim(ylim_SL);

        plot(X_odd,SpOrb_odd_SL,'-','LineWidth',1.5);
        plot(X_even,SpOrb_even_SL,'-.','LineWidth',1.5);
        plot(X(2:end-1),SpOrb_avg_SL,'-.','LineWidth',1.5);
        plot(Lambda.^(1:numel(SpOrb_corr)), zeros(1,numel(SpOrb_corr)), '--', 'Color', [.5,.5,.5]);

        patch([1e-20, X(1), X(1), 1e-20], [repmat(Yexplim_SL(1), [1,2]), repmat(Yexplim_SL(2), [1,2])], [.7 .7 .7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        
        ax=gca;
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;     
        set(ax2, 'XMinorTick', 'off', 'YMinorTick', 'off');
        xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',20);
        ylabel('$\langle \vec{S} \cdot \vec{J}^{m}_{\mathrm{sp-orb}} \cdot \vec{T} \rangle$', 'Interpreter','latex','FontSize',20);
        legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'}, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 20);
        hold off;
    end

elseif intmp == 2   % TsoK_Aniso_NRG
    
    FileInfo = dir(path);
    strtmp = cell(0,0);
    cnt = 1;
    count_iso = 0; count_aniso = 0;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1) == 'J' || DirName(1) == 'T'
            strtmp = cat(2, strtmp, [sprintf('%.15g',cnt),': ',DirName]);
            cnt = cnt + 1;

            tmp = sscanf(strtmp{end}, [sprintf('%d',numel(strtmp)),': J0=%f_K_perp=%f_K_z=%f_I0=%f_T=%f']);
            if numel(tmp) == 5
                if tmp(4) == 0
                    count_iso = count_iso + 1;
                else
                    count_aniso = count_aniso + 1;
                    disp(strtmp{end});
                end
            end

        end
    end
    disp(count_iso);
    disp(count_aniso);

    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;

    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    % extract parameters from folder name
    tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': J0=%f_K_perp=%f_K_z=%f_I0=%f_T=%f']);
    
    if numel(tmp) == 5
        J0 = tmp(1);
        K_perp = tmp(2);
        K_z = tmp(3);
        I0 = tmp(4);
        T = tmp(5);
        I_aniso = false;
    else
        tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': J0=%f_K_perp=%f_K_z=%f_I_perp=%f_I_z=%f_T=%f']);
        J0 = tmp(1);
        K_perp = tmp(2);
        K_z = tmp(3);
        I_perp = tmp(4);
        I_z = tmp(5);
        T = tmp(6);
        I_aniso = true;
    end

    tmp = strtmp{idx};
    path = [path, filesep, tmp(numel(num2str(idx))+3:end)];

    % make a list of unnecessary file infos
    FileInfo = dir(path);
    EraseIdx = [];
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')                 % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')            % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'TsoK')          % Erase log files
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];    % Erase selected unnecessary file infos

    Eflow = cell(2,1);      % RG flow data. 1: Etot, 2: Qtot
    ImpDyn = cell(3,1);     % Impurity dynamical susceptibilities. 1: ImpSp, 2: ImpOrb_plus, 3: ImpOrb_z
    BathDyn = cell(3,1);    % Bath dynamical susceptibilities. 1: BathSp, 2: BathOrb_plus, 3: BathOrb_z
    avail = false*ones(1,10);
    % RGflow, ImpDyn, Second derivatives of ImpDyn, BathDyn, Second derivatives of BathDyn,
    % Spin-Spin correlators, Orbital-Orbital correlators, 
    % Cumulated spin correlators, Cumulated orbital correlators

    EntData_exist = false;
    % Check availabel and unavailable data
    for it = (1:numel(FileInfo))
        tmp = load([FileInfo(it).folder, filesep, FileInfo(it).name]);
        field = fieldnames(tmp);

        switch FileInfo(it).name
            case 'Etot.mat'
                Eflow{1} = getfield(tmp, field{1});
                avail(1) = true;
            case 'Qtot.mat'
                Eflow{2} = getfield(tmp, field{1});
                avail(1) = true;
            case 'ocont.mat'
                ocont = getfield(tmp, field{1});
            case 'NRG_Op=ImpSp.mat'
                ImpDyn{1} = getfield(tmp, field{1});
                avail(2) = true;
                avial(3) = true;
            case 'NRG_Op=ImpOrb_plus.mat'
                ImpDyn{2} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpOrb_z.mat'
                ImpDyn{3} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=BathSp.mat'
                BathDyn{1} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathOrb_plus.mat'
                BathDyn{2} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathOrb_z.mat'
                BathDyn{3} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'Sent_imp.mat'
                Sent_imp = getfield(tmp, field{1});
                avail(6) = true;
            case 'Sent_imp_beta=2.mat'
                Sent_imp = getfield(tmp, field{1});
                avail(6) = true;
            case 'Temps.mat'
                Temps = getfield(tmp, field{1});
            case 'EntData.mat'
                EntData_exist = true;
                EntData = getfield(tmp, field{1});
                avail(6) = true;
            case 'spin_spin_correlators.mat'
                Sp_corr = getfield(tmp,field{1});
                avail(7) = true;
                avail(9) = true;
            case 'orbital_orbital_correlators.mat'
                Orb_corr = getfield(tmp,field{1});
                avail(8) = true;
                avail(10) = true;
            otherwise
                fprintf('WRN: unknown data type');
        end
    end

    beta = 1.2;
    if EntData_exist
        idx = find(abs(EntData.beta - beta) < 1e-3);
        Temps = EntData.Temps{idx};
        Sent_imp = EntData.S_imp{idx};
    end

    strtmp = cell(1,0);
    options = [];
    if avail(1)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': RG flow diagram']});
        options = [options,1];
    end
    if avail(2)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity dynamic susceptibilities']});
        options = [options,2];
    end
    if avail(3)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of impurity dynamic susceptibilities']});
        options = [options,3];
    end
    if avail(4)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Bath dynamic susceptibilities']});
        options = [options,4];
    end
    if avail(5)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of bath dynamic susceptibilities']});
        options = [options,5];
    end
    if avail(6)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity contribution to entropy']});
        options = [options,6];
    end
    if avail(7)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': sp-sp correlators']});
        options = [options,7];
    end
    if avail(8)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': orbz-orbz correlators']});
        options = [options,8];
    end
    if avail(9)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': cumulated spin-spin correlators']});
        options = [options,9];
    end
    if avail(10)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': cumulated orbz-orbz correlators']});
        options = [options,10];
    end
    strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': All of the above']});
    dispbox('-width',130,strtmp{:});

    fprintf('Which one do you want to plot?\n');
    intmp = input('>>> ');

    if intmp == numel(options)+1    % plot all
        chosen = avail;
    else
        chosen = false*ones(1,10);
        chosen(options(intmp)) = true;
        while ~isempty(intmp) && ~isequal(chosen, avail)
            fprintf('Type the index of data you want to plot\n');
            fprintf('(Press enter to finish)\n');
            intmp = input('>>> ');
            if intmp == numel(options) + 1
                chosen = avail;
            else
                chosen(options(intmp)) = true;
            end
        end
    end

    if avail(1) && chosen(1)
        %{
        plotE(Eflow{1}, Eflow{2}, 'title', '$ \mathrm{ RGflow: \ U(1)_{c1} \times U(1)_{c2} \times SU(2)_{sp}} $', 'FontSize', 15, ...
                                                'Emax',1.3,'legmax',8,'Qdiff',[0,0,0]);
        %}

        %{}
        if I_aniso
            plotE(Eflow{1}, Eflow{2}, 'title', ['$J_{0}=',sprintf('%.15g',J0),'\, K_{\perp}=',sprintf('%.15g',K_perp),'\, K_{z}=',sprintf('%.15g',K_z),'\, I_{\perp}=',sprintf('%.15g',I_perp), ...
                                                    '\, I_{z}=',sprintf('%.15g',I_z), '$'], 'Emax',2,'legmax',15,'Qdiff',[0,0,0]);
        else
            plotE(Eflow{1}, Eflow{2}, 'title', ['$J_{0}=',sprintf('%.15g',J0),'\, K_{\perp}=',sprintf('%.15g',K_perp),'\, K_{z}=',sprintf('%.15g',K_z),'\, I_{0}=',sprintf('%.15g',I0),'$'], ...
                                                    'Emax',2,'legmax',15,'Qdiff',[0,0,0]);
        end
        %}
    end

    if avail(2) && avail(3)

        names = {'$\vec{S}$', '$T^{+}$', '$T^{z}$'};
        legends = cell(0,0);
        for it = (1:3)
            if isempty(ImpDyn{it})
                ImpDyn{it} = [];
            else
                legends = cat(2,legends,names{it});
            end
        end
        
        num_ImpDyn = numel(ImpDyn);         % number of impurity dynamic susceptibilities

        log_T = log10(ocont(ocont>0));              % log temperatures
        log_ImpDyn = cell(num_ImpDyn,1);            % log impurity dynamic sysceptibilities
        log_ImpDyn_1stDer = cell(num_ImpDyn,1);     % first derivatives
        log_T_1stDer = cell(num_ImpDyn,1);          % log temperatures for first derivatives
        log_ImpDyn_2ndDer = cell(num_ImpDyn,1);     % second derivatives
        log_T_2ndDer = cell(num_ImpDyn,1);          % log temperatures for second derivatives

        for it = (1:num_ImpDyn)
            tmp = ImpDyn{it};
            log_ImpDyn{it} = log(tmp(ocont>0))./log(10);

            log_ImpDyn_1stDer{it} = diff(log_ImpDyn{it},1)./diff(log_T,1);      % first derivative
            tmp = movmean(log_T, [0,1]);                                        
            log_T_1stDer{it} = tmp(1:end-1);                                    % log temperatures for first derivatives

            log_ImpDyn_2ndDer{it} = diff(log_ImpDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                    % log temperatures for second derivatives
        end

        if chosen(2)

            %%%%%
            %{}
            figure;
            hold on;
            %}
            %%%%%
            legend('AutoUpdate','off');
            linestyle = {'-', '--', ':'};
            color = {[0, .4470, .7410], [.8500 .3250 .0980], [.9290 .6940 .1250]};
            for it = (1:num_ImpDyn)
                plot(ocont(ocont<1), ImpDyn{it}(ocont<1), 'LineWidth',2,'LineStyle',linestyle{it},'color',color{it});
            end
            xlim([1e-16, 1]);
            ylim([1e-6, 1e12]);
            
            ax = gca;
            legend(legends,'Interpreter','latex','Location','best','FontSize',25);

            set(gca,'XScale','log','YScale','log','FontSize',18);
            xlabel('$\omega$','Interpreter','latex','FontSize',24);
            ylabel('$\chi_{\mathrm{imp}} (\omega)$','Interpreter','latex','FontSize',24);

            if I_aniso
                %title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{\perp}, I_{z}) = (', ...
                %            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I_perp),', ',sprintf('%.15g',I_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
                title(['$(J_{0}, K_{\perp}, K_z, I_{\perp}, I_{z}) = (', sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp), ...
                            ', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I_perp),', ',sprintf('%.15g',I_z),')$'],'Interpreter','latex','FontSize',20);
            else
                title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            end

            %{}
            legend('AutoUpdate','off');
            fit_range = {[-9, -13]};
            [Coeff,Rsq] = LinearFit(log_T, log_ImpDyn(1),fit_range);
            colors = orderedcolors("gem");

            textShifts = {[-0.7, 2.5]};
            lineShifts = [0.7];

            for it = 1:numel(Coeff)
                text_x = mean(fit_range{it}) + textShifts{it}(1);               
                text_y = polyval(Coeff{it},text_x) + textShifts{it}(2);     

                text_x = power(10, text_x);
                text_y = power(10, text_y);
                
                y = polyval(Coeff{it},fit_range{it}) + lineShifts(it);
                x = power(10,fit_range{it});
                y = power(10,y);

                power_txt = ['$\omega^{',sprintf('%.1f',Coeff{it}(1)),'}$'];
                plot(x,y,'-','Color','black','LineWidth',1);
                text(text_x, text_y, power_txt,'Interpreter','latex','FontSize',20);
                
            end

            legend('AutoUpdate','on');
            hold off;
        end

        if chosen(3)
            figure;
            hold on;
            linestyle = {'-', '-.', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_ImpDyn)
                X = log_T_2ndDer{it};
                Y = log_ImpDyn_2ndDer{it};
                plot(X, Y,'LineWidth',2,'LineStyle',linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_ImpDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([log(T)/log(10)-1,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            if I_aniso
                title(['$\mathrm{2nd \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{\perp}, I_{z}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I_perp),', ',sprintf('%.15g',I_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            else
                title(['$\mathrm{2nd \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            end
            hold off;
        end

    end

    if avail(4) && avail(5)

        names = {'$\vec{J}_{\mathrm{sp}}$','$J_{\mathrm{orb}}^{+}$','$J_{\mathrm{orb}}^{z}$'};
        legends = cell(0,0);
        for it = (1:3)
            if isempty(BathDyn{it})
                BathDyn(it) = [];
            else
                legends = cat(2,legends,names{it});
            end
        end

        
        num_BathDyn = numel(BathDyn);         % Number of bath dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % Log temperatures
        log_BathDyn = cell(num_BathDyn,1);          % Log bath dynamic susceptibilities
        log_BathDyn_1stDer = cell(num_BathDyn,1);   % First derivatives
        log_T_1stDer = cell(num_BathDyn,1);         % Log Temperatures for fisrt derivatives
        log_BathDyn_2ndDer = cell(num_BathDyn,1);   % Second derivatives
        log_T_2ndDer = cell(num_BathDyn,1);         % Log Temperatures for second derivatives

        for it = (1:num_BathDyn)
            tmp = BathDyn{it};
            log_BathDyn{it} = log(tmp(ocont>0))./log(10);

            log_BathDyn_1stDer{it} = diff(log_BathDyn{it},1)./diff(log_T,1);                  % First derivative
            tmp = movmean(log_T, [0,1]);
            log_T_1stDer{it} = tmp(1:end-1);                                                % Log Temperatures for fisrt derivatives

            log_BathDyn_2ndDer{it} = diff(log_BathDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % Second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                % Log Temperatures for second derivatives
        end
        
        if chosen(4)

            figure;
            hold on;
            for it = (1:num_BathDyn)
                plot(ocont(ocont < 1),BathDyn{it}(ocont < 1),'Linewidth',2);
            end
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','log','YScale','log','fontsize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi'''' (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{Bath \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(5)
            figure;
            hold on;
            linestyle = {'-', '-', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_BathDyn)
                plot(log_T_2ndDer{it}, log_BathDyn_2ndDer{it}, 'LineWidth', 2, 'LineStyle', linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_BathDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([-23,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Bath \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(6) && chosen(6)

        figure('Position', [100, 100, 800, 350]);
        hold on;

        set(gca, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 18);
        set(gca, 'XTick', 10.^(-15:5:0));
        set(gca, 'FontSize', 24);
        set(gca, 'LineWidth', 2);  % make axis lines (incl. ticks) bold
        set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
        set(gca, 'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.3);

        plot(Temps,exp(Sent_imp),'Linewidth',2);
        xlim([1e-16,1]);
        ylim([0.97,4.1])
        plot([1e-22,1],[sqrt(2),sqrt(2)],'--','LineWidth',1,'color',[.7,.7,.7]);
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'tex';
        set(gca, 'TickLabelInterpreter', 'latex');
        ticks = {'1','$\sqrt{2}$','2','$2 \times \sqrt{2}$','4'};
        yticklabels(ticks);
        yticks([1,sqrt(2),2,sqrt(8),4]);

        xlabel('T','Interpreter','latex','FontSize',27);
        ylabel('$\mathrm{exp}(S_{\mathrm{imp}})$','Interpreter','latex','FontSize',27);
        %title('$\mathrm{Impurity \ contribution \ to \ entropy}$','Interpreter','latex','FontSize',30);
        %{}
        if I_aniso
            title(['$ J_{0}=', sprintf('%.15g',J0),',K_{\perp}=',sprintf('%.15g',K_perp),',K_{z}=',sprintf('%.15g',K_z),'I_{\perp}=',sprintf('%.15g',I_perp),'I_{z}=',sprintf('%.15g',I_z),'$'],'Interpreter','latex','FontSize',20);
        else
            title(['$ J_{0}=', sprintf('%.15g',J0),',K_{\perp}=',sprintf('%.15g',K_perp),',K_{z}=',sprintf('%.15g',K_z),'I_{0}=',sprintf('%.15g',I0),'$'],'Interpreter','latex','FontSize',20);
        end
        %}
        hold off;
    end

    if avail(7)

        Sp_corr = cell2mat(Sp_corr);
        Sp_odd = Sp_corr(1:2:end);
        Sp_even = Sp_corr(2:2:end);
        Sp_avg = [];
        Sp_cum = [];
        Sp_cum_avg = [];

        for it = (2:numel(Sp_corr)-1)
            Sp_avg = [Sp_avg, (2*Sp_corr(it) + Sp_corr(it+1) + Sp_corr(it-1))/4];
        end

        for it = (1:numel(Sp_corr))
            Sp_cum(it) = sum(Sp_corr(1:it));
        end

        for it = (2:numel(Sp_corr)-1)
            Sp_cum_avg = [Sp_cum_avg, (2*Sp_cum(it) + Sp_cum(it+1) + Sp_cum(it-1))/4];
        end

        X = (1:numel(Sp_corr));
        X_odd = (1:2:numel(Sp_corr));
        X_even = (2:2:numel(Sp_corr));
        X = (2.5).^(X/2);
        X_odd = (2.5).^(X_odd/2);
        X_even = (2.5).^(X_even/2);

        if chosen(7)

            figure;
            hold on;
            
            plot(X_odd,abs(Sp_odd),'--','LineWidth',1.5);
            plot(X_even,abs(Sp_even),'--','LineWidth',1.5);
            plot(X(2:end-1),abs(Sp_avg),'-','LineWidth',1.5);
    
            xlim([X(1),X(end)]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 15;
            ax.YAxis.FontSize = 15;
            set(gca,'XScale','log','YScale','log','Xdir','reverse');
            xlabel('$\Lambda^{m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle$','Interpreter','latex','FontSize',25);

            legend({'$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle \mathrm{odd}$', '$\rangle \vec{S}_{imp} \cdot \vec{S_{m}} \rangle \mathrm{even}$', ...
                            '$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle \mathrm{average}$'},'Location','southeast','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Spin-Spin \ Correlators} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(9)

            figure;
            hold on;
            
            plot(X_odd,Sp_cum(1:2:end),'--','LineWidth',1.5);
            plot(X_even,Sp_cum(2:2:end),'--','LineWidth',1.5);
            plot(X(2:end-1),Sp_cum_avg,'-','LineWidth',1.5);
    
            xlim([X(1),X(end)]);
            plot([X(end),X(1)],[-0.75,-0.75],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle \vec{S}_{\mathrm{imp}} \cdot \left( \sum_{i=1}^{m} \vec{S_{i}} \right) \rangle$','Interpreter','latex','FontSize',25);

            set(gca, 'TickLabelInterpreter', 'latex');
            ticks = {'-0.75','-0.6','-0.4','-0.2','0'};
            yticklabels(ticks);
            yticks([-0.75,-0.6,-0.4,-0.2,0]);

            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Cumulated \ Spin-Spin \ Correlators} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

    end

    if avail(8)

        Orb_corr = cell2mat(Orb_corr);
        Orb_odd = Orb_corr(1:2:end);
        Orb_even = Orb_corr(2:2:end);
        Orb_avg = [];
        Orb_cum = [];
        Orb_cum_avg = [];

        for it = (2:numel(Orb_corr)-1)
            Orb_avg = [Orb_avg, (2*Orb_corr(it) + Orb_corr(it+1) + Orb_corr(it-1))/4];
            %Orb_cum(it-1) = sum(Orb_avg(1:it-1));
        end

        for it = (1:numel(Orb_corr))
            Orb_cum(it) = sum(Orb_corr(1:it));
        end

        for it = (2:numel(Orb_corr)-1)
            Orb_cum_avg = [Orb_cum_avg, (2*Orb_cum(it) + Orb_cum(it+1) + Orb_cum(it-1))/4];
        end

        X = (1:numel(Orb_corr));
        X_odd = (1:2:numel(Orb_corr));
        X_even = (2:2:numel(Orb_corr));
        X = 4.^(X/2);
        X_odd = 4.^(X_odd/2);
        X_even = 4.^(X_even/2);

        if chosen(8)

            figure;
            hold on;
    
            plot(X_odd,abs(Orb_odd),'--','LineWidth',1.5);
            plot(X_even,abs(Orb_even),'--','LineWidth',1.5);
            plot(X(2:end-1),abs(Orb_avg),'-','LineWidth',1.5);
    
            xlim([X(1),X(end)]);
            %xlim([1e-16,X(1)]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 15;
            ax.YAxis.FontSize = 15;
            set(gca,'XScale','log','YScale','log','Xdir','reverse');
            xlabel('$\Lambda^{m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$< T_{\mathrm{imp}}^{z} T_{\mathrm{imp}}^{z} >$','Interpreter','latex','FontSize',25);
            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Orbital-Orbital \ Correlators} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(10)

            figure;
            hold on;
            
            plot(X_odd,Orb_cum(1:2:end),'--','LineWidth',1.5);
            plot(X_even,Orb_cum(2:2:end),'--','LineWidth',1.5);
            plot(X(2:end-1),Orb_cum_avg,'-','LineWidth',1.5);
    
            xlim([X(1),X(end)]);
            %xlim([1e-10,X(1)]);
            %ylim([-0.26,0.01]);
            ylim([-0.25,0.01]);
            %plot([X(end),X(1)],[-0.25,-0.25],'--','LineWidth',1,'color',[.7,.7,.7]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 25;
            ax.YAxis.FontSize = 25;
            %set(gca,'XScale','linear','YScale','linear');
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle T_{\mathrm{imp}} \left( \sum_{i=1}^{m} T_{i} \right) \rangle$','Interpreter','latex','FontSize',25);

            set(gca, 'TickLabelInterpreter', 'latex');
            %{}
            ticks = {'-0.25','-0.2','-0.1','0'};
            yticklabels(ticks);
            yticks([-0.25,-0.2,-0.1,0]);
            %}

            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Cumulated \ Orbz-Orbz \ Correlators} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end
    end

elseif intmp == 3           %% TCK_Aniso

    FileInfo = dir(path);
    strtmp = cell(0,0);
    cnt = 1;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1) == 'J' || DirName(1) == 'T'
            strtmp = cat(2, strtmp, [sprintf('%.15g',cnt),': ',DirName]);
            cnt = cnt + 1;
        end
    end
    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;

    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    % extract parameters from folder name
    tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': J_perp=%f_J_z=%f_T=%f']);
    J_perp = tmp(1);
    J_z = tmp(2);
    T = tmp(3);

    tmp = strtmp{idx};
    path = [path, filesep, tmp(numel(num2str(idx))+3:end)];

    % make a list of unnecessary file infos
    FileInfo = dir(path);
    EraseIdx = [];
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')                 % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')            % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'TCK')          % Erase log files
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];    % Erase selected unnecessary file infos

    Eflow = cell(2,1);      % RG flow data. 1: Etot, 2: Qtot
    ImpDyn = cell(2,1);     % Impurity dynamical susceptibilities. 1: ImpSp_plus, 2: ImpSp_z
    BathDyn = cell(2,1);    % Bath dynamical susceptibilities. 1: BathSp_plus, 2: BathSp_z
    avail = [false, false, false, false, false, false];
    % RGflow, ImpDyn, Second derivatives of ImpDyn, BathDyn, Second derivatives of BathDyn

    % Check availabel and unavailable data
    for it = (1:numel(FileInfo))
        tmp = load([FileInfo(it).folder, filesep, FileInfo(it).name]);
        field = fieldnames(tmp);

        switch FileInfo(it).name
            case 'Etot.mat'
                Eflow{1} = getfield(tmp, field{1});
                avail(1) = true;
            case 'Qtot.mat'
                Eflow{2} = getfield(tmp, field{1});
                avail(1) = true;
            case 'ocont.mat'
                ocont = getfield(tmp, field{1});
            case 'NRG_Op=ImpSp_plus.mat'
                ImpDyn{1} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpSp_z.mat'
                ImpDyn{2} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=BathSp_plus.mat'
                BathDyn{1} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathSp_z.mat'
                BathDyn{2} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'Sent_imp.mat'
                Sent_imp = getfield(tmp, field{1});
                avail(6) = true;
            case 'Temps.mat'
                Temps = getfield(tmp, field{1});
            otherwise
                fprintf('WRN: unknown data type');
        end
    end

    strtmp = cell(1,0);
    options = [];
    if avail(1)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': RG flow diagram']});
        options = [options,1];
    end
    if avail(2)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity dynamic susceptibilities']});
        options = [options,2];
    end
    if avail(3)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of impurity dynamic susceptibilities']});
        options = [options,3];
    end
    if avail(4)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Bath dynamic susceptibilities']});
        options = [options,4];
    end
    if avail(5)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of dynamic susceptibilities']});
        options = [options,5];
    end
    if avail(6)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity contribution to entropy']});
        options = [options,6];
    end
    strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': All of the above']});
    dispbox('-width',130,strtmp{:});

    fprintf('Which one do you want to plot?\n');
    intmp = input('>>> ');

    if intmp == numel(options)+1    % plot all
        chosen = avail;
    else
        chosen = [false, false, false, false, false, false];
        chosen(options(intmp)) = true;
        while ~isempty(intmp) && ~isequal(chosen, avail)
            fprintf('Type the index of data you want to plot\n');
            fprintf('(Press enter to finish)\n');
            intmp = input('>>> ');
            if intmp == numel(options) + 1
                chosen = avail;
            else
                chosen(options(intmp)) = true;
            end
        end
    end

    if avail(1) && chosen(1)
        plotE(Eflow{1}, Eflow{2}, 'title', ['$J_{\perp}=',sprintf('%.15g',J_perp),'\, J_{z}=',sprintf('%.15g',J_z),'$'], ...
                                                'Emax',3,'legmax',13,'Qdiff',[0,0,0]);
    end

    if avail(2) && avail(3)

        names = {'$S^{+}$', '$S^{z}$'};
        legends = cell(0,0);
        for it = (1:2)
            if isempty(ImpDyn{it})
                ImpDyn{it} = [];
            else
                legends = cat(2,legends,names{it});
            end
        end
        
        num_ImpDyn = numel(ImpDyn);         % number of impurity dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % log temperatures
        log_ImpDyn = cell(num_ImpDyn,1);            % log impurity dynamic sysceptibilities
        log_ImpDyn_1stDer = cell(num_ImpDyn,1);     % first derivatives
        log_T_1stDer = cell(num_ImpDyn,1);          % log temperatures for first derivatives
        log_ImpDyn_2ndDer = cell(num_ImpDyn,1);     % second derivatives
        log_T_2ndDer = cell(num_ImpDyn,1);          % log temperatures for second derivatives

        for it = (1:num_ImpDyn)
            tmp = ImpDyn{it};
            log_ImpDyn{it} = log(tmp(ocont>0))./log(10);

            log_ImpDyn_1stDer{it} = diff(log_ImpDyn{it},1)./diff(log_T,1);      % first derivative
            tmp = movmean(log_T, [0,1]);                                        
            log_T_1stDer{it} = tmp(1:end-1);                                    % log temperatures for first derivatives

            log_ImpDyn_2ndDer{it} = diff(log_ImpDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                    % log temperatures for second derivatives
        end

        if chosen(2)

            figure;
            
            hold on;
            linestyle = {'-', '-', '-'};
            for it = (1:num_ImpDyn)
                plot(ocont(ocont<1), ImpDyn{it}(ocont<1), 'LineWidth',2,'LineStyle',linestyle{it});
            end
            
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','best','FontSize',25);

            set(gca,'XScale','log','YScale','log','FontSize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi_{\mathrm{imp}} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities} \ (J_{\perp}, J_z) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);

            %{
            legend('AutoUpdate','off');
            fit_range = [-7, -12; -13, -18; -13, -18];
            [a1,Rsq1,a2,Rsq2,a3,Rsq3] = Insert(log_T, log_ImpDyn,fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2;
            text_y = polyval(a1,text_x) + 1;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1)+0.5;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            plot(x1,y1,'-','Color',[0,0,0],'LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',15);

            x2 = fit_range(2,:);
            text_x = (x2(1)+x2(2))/2;
            text_y = polyval(a2,text_x) + 1;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y2 = polyval(a2,x2)+0.5;
            x2 = power(10,x2);
            y2 = power(10,y2);
            text2 = ['$\omega^{',sprintf('%.2f',a2(1)),'}$'];
            plot(x2,y2,'-','Color',[0,0,0],'LineWidth',1);
            text(text_x, text_y, text2,'Interpreter','latex','FontSize',15);

            x3 = fit_range(3,:);
            text_x = (x3(1)+x3(2))/2;
            text_y = polyval(a3,text_x) - 2.3;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y3 = polyval(a3,x3) - 0.5;
            x3 = power(10,x3);
            y3 = power(10,y3);
            text3 = ['$\omega^{',sprintf('%.2f',a3(1)),'}$'];
            plot(x3,y3,'-','Color',[0,0,0],'LineWidth',1);
            text(text_x, text_y, text3,'Interpreter','latex','FontSize',15);
            legend('AutoUpdate','on');
            %}

            hold off;
        end

        if chosen(3)
            figure;
            hold on;
            linestyle = {'-', '-.', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_ImpDyn)
                X = log_T_2ndDer{it};
                Y = log_ImpDyn_2ndDer{it};
                plot(X, Y,'LineWidth',2,'LineStyle',linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_ImpDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([log(T)/log(10)-1,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities} \ (J_{\perp}, J_{z}) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(4) && avail(5)

        names = {'$J_{\mathrm{sp}}^{+}$','$J_{\mathrm{sp}}^{z}$'};
        legends = cell(0,0);
        for it = (1:2)
            if isempty(BathDyn{it})
                BathDyn(it) = [];
            else
                legends = cat(2,legends,names{it});
            end
        end

        
        num_BathDyn = numel(BathDyn);         % Number of bath dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % Log temperatures
        log_BathDyn = cell(num_BathDyn,1);          % Log bath dynamic susceptibilities
        log_BathDyn_1stDer = cell(num_BathDyn,1);   % First derivatives
        log_T_1stDer = cell(num_BathDyn,1);         % Log Temperatures for fisrt derivatives
        log_BathDyn_2ndDer = cell(num_BathDyn,1);   % Second derivatives
        log_T_2ndDer = cell(num_BathDyn,1);         % Log Temperatures for second derivatives

        for it = (1:num_BathDyn)
            tmp = BathDyn{it};
            log_BathDyn{it} = log(tmp(ocont>0))./log(10);

            log_BathDyn_1stDer{it} = diff(log_BathDyn{it},1)./diff(log_T,1);                  % First derivative
            tmp = movmean(log_T, [0,1]);
            log_T_1stDer{it} = tmp(1:end-1);                                                % Log Temperatures for fisrt derivatives

            log_BathDyn_2ndDer{it} = diff(log_BathDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % Second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                % Log Temperatures for second derivatives
        end
        
        if chosen(4)

            figure;
            hold on;
            for it = (1:num_BathDyn)
                plot(ocont(ocont < 1),BathDyn{it}(ocont < 1),'Linewidth',2);
            end
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','log','YScale','log','fontsize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi'''' (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{Bath \ Dynamic \ Susceptibilities} \ (J_{\perp}, J_{z}) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(5)
            figure;
            hold on;
            linestyle = {'-', '-', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_BathDyn)
                plot(log_T_2ndDer{it}, log_BathDyn_2ndDer{it}, 'LineWidth', 2, 'LineStyle', linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_BathDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([-23,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Bath \ Dynamic \ Susceptibilities} \ (J_{\perp}, J_{z}) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(6) && chosen(6)

        figure;
        hold on;
        plot(Temps,exp(Sent_imp),'Linewidth',1.5);
        xlim([1e-22,1]);
        ylim([1,4.1])
        plot([1e-22,1],[sqrt(8),sqrt(8)],'--','LineWidth',1,'color',[.7,.7,.7]);
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'tex';
        set(gca, 'TickLabelInterpreter', 'latex');
        ticks = {'1','2','$2 \times \sqrt{2}$','4'};
        yticklabels(ticks);
        yticks([1,2,sqrt(8),4]);

        ax = gca;
        ax.XAxis.FontSize = 5;
        ax.YAxis.FontSize = 5;
        set(gca,'XScale','log','YScale','linear','fontsize',20);
        xlabel('T','Interpreter','latex','FontSize',25);
        ylabel('$\mathrm{exp}(S_{\mathrm{imp}})$','Interpreter','latex','FontSize',25);
        title(['$\mathrm{Impurity \ contribution \ to \ entropy} \ (J_{0}, K_{\perp}, K_{z}, I_{0}) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
        hold off;
    end


elseif intmp == 4       %% Kondo_Aniso_NRG


    FileInfo = dir(path);
    strtmp = cell(0,0);
    cnt = 1;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1) == 'J' || DirName(1) == 'T'
            strtmp = cat(2, strtmp, [sprintf('%.15g',cnt),': ',DirName]);
            cnt = cnt + 1;
        end
    end
    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;

    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    % extract parameters from folder name
    tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': J_perp=%f_J_z=%f_T=%f']);
    J_perp = tmp(1);
    J_z = tmp(2);
    T = tmp(3);

    tmp = strtmp{idx};
    path = [path, filesep, tmp(numel(num2str(idx))+3:end)];

    % make a list of unnecessary file infos
    FileInfo = dir(path);
    EraseIdx = [];
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')                 % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')            % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName(1:5),'Kondo')          % Erase log files
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];    % Erase selected unnecessary file infos

    Eflow = cell(2,1);      % RG flow data. 1: Etot, 2: Qtot
    ImpDyn = cell(2,1);     % Impurity dynamical susceptibilities. 1: ImpSp_plus, 2: ImpSp_z
    BathDyn = cell(2,1);    % Bath dynamical susceptibilities. 1: BathSp_plus, 2: BathSp_z
    avail = [false, false, false, false, false, false];
    % RGflow, ImpDyn, Second derivatives of ImpDyn, BathDyn, Second derivatives of BathDyn

    % Check availabel and unavailable data
    for it = (1:numel(FileInfo))
        tmp = load([FileInfo(it).folder, filesep, FileInfo(it).name]);
        field = fieldnames(tmp);

        switch FileInfo(it).name
            case 'Etot.mat'
                Eflow{1} = getfield(tmp, field{1});
                avail(1) = true;
            case 'Qtot.mat'
                Eflow{2} = getfield(tmp, field{1});
                avail(1) = true;
            case 'ocont.mat'
                ocont = getfield(tmp, field{1});
            case 'NRG_Op=ImpSp_plus.mat'
                ImpDyn{1} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpSp_z.mat'
                ImpDyn{2} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=BathSp_plus.mat'
                BathDyn{1} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathSp_z.mat'
                BathDyn{2} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'Sent_imp.mat'
                Sent_imp = getfield(tmp, field{1});
                avail(6) = true;
            case 'Temps.mat'
                Temps = getfield(tmp, field{1});
            otherwise
                fprintf('WRN: unknown data type');
        end
    end

    strtmp = cell(1,0);
    options = [];
    if avail(1)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': RG flow diagram']});
        options = [options,1];
    end
    if avail(2)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity dynamic susceptibilities']});
        options = [options,2];
    end
    if avail(3)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of impurity dynamic susceptibilities']});
        options = [options,3];
    end
    if avail(4)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Bath dynamic susceptibilities']});
        options = [options,4];
    end
    if avail(5)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of dynamic susceptibilities']});
        options = [options,5];
    end
    if avail(6)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity contribution to entropy']});
        options = [options,6];
    end
    strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': All of the above']});
    dispbox('-width',130,strtmp{:});

    fprintf('Which one do you want to plot?\n');
    intmp = input('>>> ');

    if intmp == numel(options)+1    % plot all
        chosen = avail;
    else
        chosen = [false, false, false, false, false, false];
        chosen(options(intmp)) = true;
        while ~isempty(intmp) && ~isequal(chosen, avail)
            fprintf('Type the index of data you want to plot\n');
            fprintf('(Press enter to finish)\n');
            intmp = input('>>> ');
            if intmp == numel(options) + 1
                chosen = avail;
            else
                chosen(options(intmp)) = true;
            end
        end
    end

    if avail(1) && chosen(1)
        plotE(Eflow{1}, Eflow{2}, 'title', ['$J_{\perp}=',sprintf('%.15g',J_perp),'\, J_{z}=',sprintf('%.15g',J_z),'$'], ...
                                                'Emax',3,'legmax',13,'Qdiff',[0,0]);
    end

    if avail(2) && avail(3)

        names = {'$S^{+}$', '$S^{z}$'};
        legends = cell(0,0);
        for it = (1:2)
            if isempty(ImpDyn{it})
                ImpDyn{it} = [];
            else
                legends = cat(2,legends,names{it});
            end
        end
        
        num_ImpDyn = numel(ImpDyn);         % number of impurity dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % log temperatures
        log_ImpDyn = cell(num_ImpDyn,1);            % log impurity dynamic sysceptibilities
        log_ImpDyn_1stDer = cell(num_ImpDyn,1);     % first derivatives
        log_T_1stDer = cell(num_ImpDyn,1);          % log temperatures for first derivatives
        log_ImpDyn_2ndDer = cell(num_ImpDyn,1);     % second derivatives
        log_T_2ndDer = cell(num_ImpDyn,1);          % log temperatures for second derivatives

        for it = (1:num_ImpDyn)
            tmp = ImpDyn{it};
            log_ImpDyn{it} = log(tmp(ocont>0))./log(10);

            log_ImpDyn_1stDer{it} = diff(log_ImpDyn{it},1)./diff(log_T,1);      % first derivative
            tmp = movmean(log_T, [0,1]);                                        
            log_T_1stDer{it} = tmp(1:end-1);                                    % log temperatures for first derivatives

            log_ImpDyn_2ndDer{it} = diff(log_ImpDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                    % log temperatures for second derivatives
        end

        if chosen(2)

            figure;
            
            hold on;
            linestyle = {'-', '-', '-'};
            for it = (1:num_ImpDyn)
                plot(ocont(ocont<1), ImpDyn{it}(ocont<1), 'LineWidth',2,'LineStyle',linestyle{it});
            end
            
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','best','FontSize',25);

            set(gca,'XScale','log','YScale','log','FontSize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi^{\mathrm{imp}} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities} \ (J_{\perp}, J_z) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);

            %{}
            legend('AutoUpdate','off');
            log_ImpDyn = cat(1,log_ImpDyn,log_ImpDyn(1));
            fit_range = [-4, -7; -4, -7; 0, 0];
            %fit_range = [-7, -12; -7, -12; 0, 0];
            [a1,Rsq1,a2,Rsq2,a3,Rsq3] = Insert(log_T, log_ImpDyn,fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2;
            text_y = polyval(a1,text_x) + 1.3;
            %text_x = (x1(1)+x1(2))/2 - 0.3;
            %text_y = polyval(a1,text_x) + 1.2;

            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1) + 0.5;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            plot(x1,y1,'-','Color',[0 0.4470 0.7410],'LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',15);

            x2 = fit_range(2,:);
            text_x = (x2(1)+x2(2))/2;
            text_y = polyval(a2,text_x) - 1;
            %text_x = (x2(1)+x2(2))/2;
            %text_y = polyval(a2,text_x) - 0.8;

            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y2 = polyval(a2,x2) - 0.5;
            x2 = power(10,x2);
            y2 = power(10,y2);
            text2 = ['$\omega^{',sprintf('%.2f',a2(1)),'}$'];
            plot(x2,y2,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1);
            text(text_x, text_y, text2,'Interpreter','latex','FontSize',15);

            log_ImpDyn(3) = [];
            %}

            hold off;
        end

        if chosen(3)
            figure;
            hold on;
            linestyle = {'-', '-.', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_ImpDyn)
                X = log_T_2ndDer{it};
                Y = log_ImpDyn_2ndDer{it};
                plot(X, Y,'LineWidth',2,'LineStyle',linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_ImpDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([log(T)/log(10)-1,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','best','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities} \ (J_{\perp}, J_{z}) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(4) && avail(5)

        names = {'$J_{\mathrm{sp}}^{+}$','$J_{\mathrm{sp}}^{z}$'};
        legends = cell(0,0);
        for it = (1:2)
            if isempty(BathDyn{it})
                BathDyn(it) = [];
            else
                legends = cat(2,legends,names{it});
            end
        end

        
        num_BathDyn = numel(BathDyn);         % Number of bath dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % Log temperatures
        log_BathDyn = cell(num_BathDyn,1);          % Log bath dynamic susceptibilities
        log_BathDyn_1stDer = cell(num_BathDyn,1);   % First derivatives
        log_T_1stDer = cell(num_BathDyn,1);         % Log Temperatures for fisrt derivatives
        log_BathDyn_2ndDer = cell(num_BathDyn,1);   % Second derivatives
        log_T_2ndDer = cell(num_BathDyn,1);         % Log Temperatures for second derivatives

        for it = (1:num_BathDyn)
            tmp = BathDyn{it};
            log_BathDyn{it} = log(tmp(ocont>0))./log(10);

            log_BathDyn_1stDer{it} = diff(log_BathDyn{it},1)./diff(log_T,1);                  % First derivative
            tmp = movmean(log_T, [0,1]);
            log_T_1stDer{it} = tmp(1:end-1);                                                % Log Temperatures for fisrt derivatives

            log_BathDyn_2ndDer{it} = diff(log_BathDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % Second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                % Log Temperatures for second derivatives
        end
        
        if chosen(4)

            figure;
            hold on;
            for it = (1:num_BathDyn)
                plot(ocont(ocont < 1),BathDyn{it}(ocont < 1),'Linewidth',2);
            end
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','best','FontSize',25);
            set(gca,'XScale','log','YScale','log','fontsize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi^{\mathrm{bath}} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{Bath \ Dynamic \ Susceptibilities} \ (J_{\perp}, J_{z}) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);

            legend('AutoUpdate','off');
            log_BathDyn = cat(1,log_BathDyn,log_BathDyn(1));
            fit_range = [-4, -7; -4, -7; 0, 0];
            %fit_range = [-7, -12; -7, -12; 0, 0];
            [a1,Rsq1,a2,Rsq2,a3,Rsq3] = Insert(log_T, log_BathDyn,fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2;
            text_y = polyval(a1,text_x) + 1.3;
            %text_x = (x1(1)+x1(2))/2 - 0.3;
            %text_y = polyval(a1,text_x) + 1.2;

            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1) + 0.5;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            plot(x1,y1,'-','Color',[0 0.4470 0.7410],'LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',15);

            x2 = fit_range(2,:);
            text_x = (x2(1)+x2(2))/2;
            text_y = polyval(a2,text_x) - 1;
            %text_x = (x2(1)+x2(2))/2;
            %text_y = polyval(a2,text_x) - 0.8;

            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y2 = polyval(a2,x2) - 0.5;
            x2 = power(10,x2);
            y2 = power(10,y2);
            text2 = ['$\omega^{',sprintf('%.2f',a2(1)),'}$'];
            plot(x2,y2,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1);
            text(text_x, text_y, text2,'Interpreter','latex','FontSize',15);

            log_BathDyn(3) = [];

            hold off;
        end

        if chosen(5)
            figure;
            hold on;
            linestyle = {'-', '-', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_BathDyn)
                plot(log_T_2ndDer{it}, log_BathDyn_2ndDer{it}, 'LineWidth', 2, 'LineStyle', linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_BathDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([-23,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Bath \ Dynamic \ Susceptibilities} \ (J_{\perp}, J_{z}) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(6) && chosen(6)

        figure;
        hold on;
        plot(Temps,exp(Sent_imp),'Linewidth',1.5);
        xlim([1e-16,1]);
        ylim([1,2.5]);
        plot([1e-16,1],[2,2],'--','LineWidth',1,'color',[.7,.7,.7]);

        ax = gca;
        ax.XAxis.FontSize = 5;
        ax.YAxis.FontSize = 5;
        set(gca,'XScale','log','YScale','linear','fontsize',20);
        xlabel('T','Interpreter','latex','FontSize',25);
        ylabel('$\mathrm{exp}(S_{\mathrm{imp}})$','Interpreter','latex','FontSize',25);
        title(['$\mathrm{Impurity \ contribution \ to \ entropy} \ (J_{0}, K_{\perp}, K_{z}, I_{0}) = (', ...
                        sprintf('%.15g',J_perp),', ',sprintf('%.15g',J_z),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
        hold off;
    end

elseif intmp == 5   % ThsoK_NRG
    
    FileInfo = dir(path);
    strtmp = cell(0,0);
    cnt = 1;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1) == 'J' || DirName(1) == 'T'
            strtmp = cat(2, strtmp, [sprintf('%.15g',cnt),': ',DirName]);
            cnt = cnt + 1;
        end
    end
    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;

    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    % extract parameters from folder name
    tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': J0=%f_K0=%f_I0=%f_T=%f']);
    J0 = tmp(1);
    K0 = tmp(2);
    I0 = tmp(3);
    T = tmp(4);

    tmp = strtmp{idx};
    path = [path, filesep, tmp(numel(num2str(idx))+3:end)];

    % make a list of unnecessary file infos
    FileInfo = dir(path);
    EraseIdx = [];
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')                 % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')            % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'TsoK')          % Erase log files
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];    % Erase selected unnecessary file infos

    Eflow = cell(2,1);      % RG flow data. 1: Etot, 2: Qtot
    ImpDyn = cell(3,1);     % Impurity dynamical susceptibilities. 1: ImpSp, 2: ImpOrb_plus, 3: ImpOrb_z
    BathDyn = cell(3,1);    % Bath dynamical susceptibilities. 1: BathSp, 2: BathOrb_plus, 3: BathOrb_z
    avail = false*ones(1,10);
    % RGflow, ImpDyn, Second derivatives of ImpDyn, BathDyn, Second derivatives of BathDyn,
    % Spin-Spin correlators, Orbital-Orbital correlators, 
    % Cumulated spin correlators, Cumulated orbital correlators

    % Check availabel and unavailable data
    for it = (1:numel(FileInfo))
        tmp = load([FileInfo(it).folder, filesep, FileInfo(it).name]);
        field = fieldnames(tmp);

        switch FileInfo(it).name
            case 'Etot.mat'
                Eflow{1} = getfield(tmp, field{1});
                avail(1) = true;
            case 'Qtot.mat'
                Eflow{2} = getfield(tmp, field{1});
                avail(1) = true;
            case 'ocont.mat'
                ocont = getfield(tmp, field{1});
            case 'NRG_Op=ImpSp.mat'
                ImpDyn{1} = getfield(tmp, field{1});
                avail(2) = true;
                avial(3) = true;
            case 'NRG_Op=ImpOrb_plus.mat'
                ImpDyn{2} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpOrb_z.mat'
                ImpDyn{3} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=BathSp.mat'
                BathDyn{1} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathOrb_plus.mat'
                BathDyn{2} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathOrb_z.mat'
                BathDyn{3} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'Sent_imp.mat'
                Sent_imp = getfield(tmp, field{1});
                avail(6) = true;
            case 'Temps.mat'
                Temps = getfield(tmp, field{1});
            case 'spin_spin_correlators.mat'
                Sp_corr = getfield(tmp,field{1});
                avail(7) = true;
                avail(9) = true;
            case 'orbital_orbital_correlators.mat'
                Orb_corr = getfield(tmp,field{1});
                avail(8) = true;
                avail(10) = true;
            otherwise
                fprintf('WRN: unknown data type');
        end
    end

    strtmp = cell(1,0);
    options = [];
    if avail(1)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': RG flow diagram']});
        options = [options,1];
    end
    if avail(2)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity dynamic susceptibilities']});
        options = [options,2];
    end
    if avail(3)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of impurity dynamic susceptibilities']});
        options = [options,3];
    end
    if avail(4)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Bath dynamic susceptibilities']});
        options = [options,4];
    end
    if avail(5)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of bath dynamic susceptibilities']});
        options = [options,5];
    end
    if avail(6)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity contribution to entropy']});
        options = [options,6];
    end
    if avail(7)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': sp-sp correlators']});
        options = [options,7];
    end
    if avail(8)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': orbz-orbz correlators']});
        options = [options,8];
    end
    if avail(9)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': cumulated spin-spin correlators']});
        options = [options,9];
    end
    if avail(10)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': cumulated orbz-orbz correlators']});
        options = [options,10];
    end
    strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': All of the above']});
    dispbox('-width',130,strtmp{:});

    fprintf('Which one do you want to plot?\n');
    intmp = input('>>> ');

    if intmp == numel(options)+1    % plot all
        chosen = avail;
    else
        chosen = false*ones(1,10);
        chosen(options(intmp)) = true;
        while ~isempty(intmp) && ~isequal(chosen, avail)
            fprintf('Type the index of data you want to plot\n');
            fprintf('(Press enter to finish)\n');
            intmp = input('>>> ');
            if intmp == numel(options) + 1
                chosen = avail;
            else
                chosen(options(intmp)) = true;
            end
        end
    end

    if avail(1) && chosen(1)
        plotE(Eflow{1}, Eflow{2}, 'title', ['$J_{0}=',sprintf('%.15g',J0),'\, K_{0}=',sprintf('%.15g',K0),'\, I_{0}=',sprintf('%.15g',I0),'$'], ...
                                                'Emax',3,'legmax',13,'Qdiff',[0,0,0,0]);
    end

    if avail(2) && avail(3)

        names = {'$\vec{S}$', '$T^{+}$', '$T^{z}$'};
        legends = cell(0,0);
        for it = (1:3)
            if isempty(ImpDyn{it})
                ImpDyn{it} = [];
            else
                legends = cat(2,legends,names{it});
            end
        end
        
        num_ImpDyn = numel(ImpDyn);         % number of impurity dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % log temperatures
        log_ImpDyn = cell(num_ImpDyn,1);            % log impurity dynamic sysceptibilities
        log_ImpDyn_1stDer = cell(num_ImpDyn,1);     % first derivatives
        log_T_1stDer = cell(num_ImpDyn,1);          % log temperatures for first derivatives
        log_ImpDyn_2ndDer = cell(num_ImpDyn,1);     % second derivatives
        log_T_2ndDer = cell(num_ImpDyn,1);          % log temperatures for second derivatives

        for it = (1:num_ImpDyn)
            tmp = ImpDyn{it};
            log_ImpDyn{it} = log(tmp(ocont>0))./log(10);

            log_ImpDyn_1stDer{it} = diff(log_ImpDyn{it},1)./diff(log_T,1);      % first derivative
            tmp = movmean(log_T, [0,1]);                                        
            log_T_1stDer{it} = tmp(1:end-1);                                    % log temperatures for first derivatives

            log_ImpDyn_2ndDer{it} = diff(log_ImpDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                    % log temperatures for second derivatives
        end

        if chosen(2)

            figure;
            
            hold on;
            linestyle = {'-', '-', '-'};
            for it = (1:num_ImpDyn)
                plot(ocont(ocont<1), ImpDyn{it}(ocont<1), 'LineWidth',2,'LineStyle',linestyle{it});
            end
            
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','best','FontSize',25);

            set(gca,'XScale','log','YScale','log','FontSize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi_{\mathrm{imp}} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities} \ (J_{0}, K_{0}, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);

            %{}
            legend('AutoUpdate','off');
            fit_range = [-3, -7; -9, -14; -9, -14];
            [a1,Rsq1,a2,Rsq2,a3,Rsq3] = Insert(log_T, log_ImpDyn,fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2;
            text_y = polyval(a1,text_x) + 1;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1)+0.5;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            plot(x1,y1,'-','Color',[0, 0.447,0.741],'LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',15);

            x2 = fit_range(2,:);
            text_x = (x2(1)+x2(2))/2;
            text_y = polyval(a2,text_x) + 1;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y2 = polyval(a2,x2)+0.5;
            x2 = power(10,x2);
            y2 = power(10,y2);
            text2 = ['$\omega^{',sprintf('%.2f',a2(1)),'}$'];
            plot(x2,y2,'-','Color',[0.85,0.325,0.098],'LineWidth',1);
            text(text_x, text_y, text2,'Interpreter','latex','FontSize',15);

            x3 = fit_range(3,:);
            text_x = (x3(1)+x3(2))/2;
            text_y = polyval(a3,text_x) - 2.3;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y3 = polyval(a3,x3) - 0.5;
            x3 = power(10,x3);
            y3 = power(10,y3);
            text3 = ['$\omega^{',sprintf('%.2f',a3(1)),'}$'];
            plot(x3,y3,'-','Color',[0.929,0.694,0.125],'LineWidth',1);
            text(text_x, text_y, text3,'Interpreter','latex','FontSize',15);
            legend('AutoUpdate','on');
            %}

            hold off;
        end

        if chosen(3)
            figure;
            hold on;
            linestyle = {'-', '-.', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_ImpDyn)
                X = log_T_2ndDer{it};
                Y = log_ImpDyn_2ndDer{it};
                plot(X, Y,'LineWidth',2,'LineStyle',linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_ImpDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([log(T)/log(10)-1,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities} \ (J_{0}, K_{0}, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(4) && avail(5)

        names = {'$\vec{J}_{\mathrm{sp}}$','$J_{\mathrm{orb}}^{+}$','$J_{\mathrm{orb}}^{z}$'};
        legends = cell(0,0);
        for it = (1:3)
            if isempty(BathDyn{it})
                BathDyn(it) = [];
            else
                legends = cat(2,legends,names{it});
            end
        end

        
        num_BathDyn = numel(BathDyn);         % Number of bath dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % Log temperatures
        log_BathDyn = cell(num_BathDyn,1);          % Log bath dynamic susceptibilities
        log_BathDyn_1stDer = cell(num_BathDyn,1);   % First derivatives
        log_T_1stDer = cell(num_BathDyn,1);         % Log Temperatures for fisrt derivatives
        log_BathDyn_2ndDer = cell(num_BathDyn,1);   % Second derivatives
        log_T_2ndDer = cell(num_BathDyn,1);         % Log Temperatures for second derivatives

        for it = (1:num_BathDyn)
            tmp = BathDyn{it};
            log_BathDyn{it} = log(tmp(ocont>0))./log(10);

            log_BathDyn_1stDer{it} = diff(log_BathDyn{it},1)./diff(log_T,1);                  % First derivative
            tmp = movmean(log_T, [0,1]);
            log_T_1stDer{it} = tmp(1:end-1);                                                % Log Temperatures for fisrt derivatives

            log_BathDyn_2ndDer{it} = diff(log_BathDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % Second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                % Log Temperatures for second derivatives
        end
        
        if chosen(4)

            figure;
            hold on;
            for it = (1:num_BathDyn)
                plot(ocont(ocont < 1),BathDyn{it}(ocont < 1),'Linewidth',2);
            end
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','log','YScale','log','fontsize',20);
            xlabel('$\omega$','Interpreter','latex','FontSize',25);
            ylabel('$\chi'''' (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{Bath \ Dynamic \ Susceptibilities} \ (J_{0}, K_{0}, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(5)
            figure;
            hold on;
            linestyle = {'-', '-', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_BathDyn)
                plot(log_T_2ndDer{it}, log_BathDyn_2ndDer{it}, 'LineWidth', 2, 'LineStyle', linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_BathDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([-23,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Bath \ Dynamic \ Susceptibilities} \ (J_{0}, K_{0}, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(6) && chosen(6)

        figure;
        hold on;
        plot(Temps,exp(Sent_imp),'Linewidth',1.5);
        xlim([1e-22,1]);
        ylim([1,4.1])
        plot([1e-22,1],[sqrt(8),sqrt(8)],'--','LineWidth',1,'color',[.7,.7,.7]);
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'tex';
        set(gca, 'TickLabelInterpreter', 'latex');
        ticks = {'1','2','$2 \times \sqrt{2}$','4'};
        yticklabels(ticks);
        yticks([1,2,sqrt(8),4]);

        ax = gca;
        ax.XAxis.FontSize = 5;
        ax.YAxis.FontSize = 5;
        set(gca,'XScale','log','YScale','linear','fontsize',20);
        xlabel('T','Interpreter','latex','FontSize',25);
        ylabel('$\mathrm{exp}(S_{\mathrm{imp}})$','Interpreter','latex','FontSize',25);
        title(['$\mathrm{Impurity \ contribution \ to \ entropy} \ (J_{0}, K_{0}, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
        hold off;
    end

    if avail(7)

        Sp_corr = cell2mat(Sp_corr);
        Sp_odd = Sp_corr(1:2:end);
        Sp_even = Sp_corr(2:2:end);
        Sp_avg = [];
        Sp_cum = [];
        Sp_cum_avg = [];

        for it = (2:numel(Sp_corr)-1)
            Sp_avg = [Sp_avg, (2*Sp_corr(it) + Sp_corr(it+1) + Sp_corr(it-1))/4];
        end

        for it = (1:numel(Sp_corr))
            Sp_cum(it) = sum(Sp_corr(1:it));
        end

        for it = (2:numel(Sp_corr)-1)
            Sp_cum_avg = [Sp_cum_avg, (2*Sp_cum(it) + Sp_cum(it+1) + Sp_cum(it-1))/4];
        end

        X = (1:numel(Sp_corr));
        X_odd = (1:2:numel(Sp_corr));
        X_even = (2:2:numel(Sp_corr));
        X = 4.^(-X/2);
        X_odd = 4.^(-X_odd/2);
        X_even = 4.^(-X_even/2);

        if chosen(7)

            figure;
            hold on;
            
            plot(X_odd,Sp_odd,'--','LineWidth',1.5);
            plot(X_even,Sp_even,'--','LineWidth',1.5);
            plot(X(2:end-1),Sp_avg,'-','LineWidth',1.5);
    
            xlim([X(end),X(1)]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 15;
            ax.YAxis.FontSize = 15;
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle$','Interpreter','latex','FontSize',25);

            legend({'$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle \mathrm{odd}$', '$\rangle \vec{S}_{imp} \cdot \vec{S_{m}} \rangle \mathrm{even}$', ...
                            '$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle \mathrm{average}$'},'Location','southeast','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Spin-Spin \ Correlators} \ (J_{0}, K_{0}, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(9)

            figure;
            hold on;
            
            plot(X_odd,Sp_cum(1:2:end),'--','LineWidth',1.5);
            plot(X_even,Sp_cum(2:2:end),'--','LineWidth',1.5);
            plot(X(2:end-1),Sp_cum_avg,'-','LineWidth',1.5);
    
            xlim([X(end),X(1)]);
            plot([X(end),X(1)],[-1,-1],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle \vec{S}_{\mathrm{imp}} \cdot \left( \sum_{i=1}^{m} \vec{S_{i}} \right) \rangle$','Interpreter','latex','FontSize',25);

            %{
            set(gca, 'TickLabelInterpreter', 'latex');
            ticks = {'-1','-0.8','-0.6','-0.4','-0.2','0'};
            yticklabels(ticks);
            yticks([-1,-0.8,-0.6,-0.4,-0.2,0]);
            %}

            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Cumulated \ Spin-Spin \ Correlators} \ (J_{0}, K_{0}, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

    end

    X = false;
    if avail(8)

        Orb_corr = cell2mat(Orb_corr);
        Orb_odd = Orb_corr(1:2:end);
        Orb_even = Orb_corr(2:2:end);
        Orb_avg = [];
        Orb_cum = [];
        Orb_cum_avg = [];

        for it = (2:numel(Orb_corr)-1)
            Orb_avg = [Orb_avg, (2*Orb_corr(it) + Orb_corr(it+1) + Orb_corr(it-1))/4];
            %Orb_cum(it-1) = sum(Orb_avg(1:it-1));
        end

        for it = (1:numel(Orb_corr))
            Orb_cum(it) = sum(Orb_corr(1:it));
        end

        for it = (2:numel(Orb_corr)-1)
            Orb_cum_avg = [Orb_cum_avg, (2*Orb_cum(it) + Orb_cum(it+1) + Orb_cum(it-1))/4];
        end

        X = (1:numel(Orb_corr));
        X_odd = (1:2:numel(Orb_corr));
        X_even = (2:2:numel(Orb_corr));
        X = 4.^(-X/2);
        X_odd = 4.^(-X_odd/2);
        X_even = 4.^(-X_even/2);

        if chosen(8)

            figure;
            hold on;
    
            plot(X_odd,Orb_odd,'--','LineWidth',1.5);
            plot(X_even,Orb_even,'--','LineWidth',1.5);
            plot(X(2:end-1),Orb_avg,'-','LineWidth',1.5);
    
            xlim([X(end),X(1)]);
            %xlim([1e-16,X(1)]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 15;
            ax.YAxis.FontSize = 15;
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$< T_{\mathrm{imp}}^{z} T_{\mathrm{imp}}^{z} >$','Interpreter','latex','FontSize',25);
            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Orbital-Orbital \ Correlators} \ (J_{0}, K_{0}, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(10)

            figure;
            hold on;
            
            plot(X_odd,Orb_cum(1:2:end),'--','LineWidth',1.5);
            plot(X_even,Orb_cum(2:2:end),'--','LineWidth',1.5);
            plot(X(2:end-1),Orb_cum_avg,'-','LineWidth',1.5);
    
            xlim([X(end),X(1)]);
            %xlim([1e-10,X(1)]);
            %ylim([-0.26,0.01]);
            ylim([-0.25,0.01]);
            %plot([X(end),X(1)],[-0.25,-0.25],'--','LineWidth',1,'color',[.7,.7,.7]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 25;
            ax.YAxis.FontSize = 25;
            %set(gca,'XScale','linear','YScale','linear');
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle T_{\mathrm{imp}} \left( \sum_{i=1}^{m} T_{i} \right) \rangle$','Interpreter','latex','FontSize',25);

            set(gca, 'TickLabelInterpreter', 'latex');
            %{}
            ticks = {'-0.25','-0.2','-0.1','0'};
            yticklabels(ticks);
            yticks([-0.25,-0.2,-0.1,0]);
            %}

            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Cumulated \ Orbz-Orbz \ Correlators} \ (J_{0}, K_{0}, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K0),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end
    end

elseif intmp == 6   % Quartic_NRG
    
    FileInfo = dir(path);
    strtmp = cell(0,0);
    cnt = 1;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1) == 'K' || DirName(1) == 'T'
            strtmp = cat(2, strtmp, [sprintf('%.15g',cnt),': ',DirName]);
            cnt = cnt + 1;
        end
    end
    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;

    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    % extract parameters from folder name
    tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': K_z=%f_Q=%f_T=%f']);
    K_z = tmp(1);
    Q = tmp(2);
    T = tmp(3);

    tmp = strtmp{idx};
    path = [path, filesep, tmp(numel(num2str(idx))+3:end)];

    % make a list of unnecessary file infos
    FileInfo = dir(path);
    EraseIdx = [];
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')                 % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')            % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'Quartic')          % Erase log files
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];    % Erase selected unnecessary file infos

    Eflow = cell(2,1);      % RG flow data. 1: Etot, 2: Qtot
    ImpDyn = cell(1,1);     % Impurity dynamical susceptibilities. 1: ImpOrb_z
    BathDyn = cell(4,1);    % Bath dynamical susceptibilities. 1: BathSp1, 2: BathSp2, 3: BathSptot, 4: BathOrb_z
    avail = false*ones(1,10);
    % RGflow, ImpDyn, Second derivatives of ImpDyn, BathDyn, Second derivatives of BathDyn,
    % Spin-Spin correlators, Orbital-Orbital correlators, 
    % Cumulated spin correlators, Cumulated orbital correlators

    % Check availabel and unavailable data
    for it = (1:numel(FileInfo))
        tmp = load([FileInfo(it).folder, filesep, FileInfo(it).name]);
        field = fieldnames(tmp);

        switch FileInfo(it).name
            case 'Etot.mat'
                Eflow{1} = getfield(tmp, field{1});
                avail(1) = true;
            case 'Qtot.mat'
                Eflow{2} = getfield(tmp, field{1});
                avail(1) = true;
            case 'ocont.mat'
                ocont = getfield(tmp, field{1});
            case 'NRG_Op=ImpOrb_z.mat'
                ImpDyn{1} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=BathSp1.mat'
                BathDyn{1} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathSp2.mat'
                BathDyn{2} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathSptot.mat'
                BathDyn{3} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'NRG_Op=BathOrb_z.mat'
                BathDyn{4} = getfield(tmp, field{1});
                avail(4) = true;
                avail(5) = true;
            case 'Sent_imp.mat'
                Sent_imp = getfield(tmp, field{1});
                avail(6) = true;
            case 'Temps.mat'
                Temps = getfield(tmp, field{1});
            case 'spin_spin_correlators.mat'
                Sp_corr = getfield(tmp,field{1});
                avail(7) = true;
                avail(9) = true;
            case 'orbital_orbital_correlators.mat'
                Orb_corr = getfield(tmp,field{1});
                avail(8) = true;
                avail(10) = true;
            otherwise
                fprintf('WRN: unknown data type');
        end
    end

    strtmp = cell(1,0);
    options = [];
    if avail(1)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': RG flow diagram']});
        options = [options,1];
    end
    if avail(2)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity dynamic susceptibilities']});
        options = [options,2];
    end
    if avail(3)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of impurity dynamic susceptibilities']});
        options = [options,3];
    end
    if avail(4)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Bath dynamic susceptibilities']});
        options = [options,4];
    end
    if avail(5)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of bath dynamic susceptibilities']});
        options = [options,5];
    end
    if avail(6)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity contribution to entropy']});
        options = [options,6];
    end
    if avail(7)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': sp-sp correlators']});
        options = [options,7];
    end
    if avail(8)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': orbz-orbz correlators']});
        options = [options,8];
    end
    if avail(9)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': cumulated spin-spin correlators']});
        options = [options,9];
    end
    if avail(10)
        strtmp = cat(2,strtmp,{[sprintf('%d',numel(strtmp)+1),': cumulated orbz-orbz correlators']});
        options = [options,10];
    end
    strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': All of the above']});
    dispbox('-width',130,strtmp{:});

    fprintf('Which one do you want to plot?\n');
    intmp = input('>>> ');

    if intmp == numel(options)+1    % plot all
        chosen = avail;
    else
        chosen = false*ones(1,10);
        chosen(options(intmp)) = true;
        while ~isempty(intmp) && ~isequal(chosen, avail)
            fprintf('Type the index of data you want to plot\n');
            fprintf('(Press enter to finish)\n');
            intmp = input('>>> ');
            if intmp == numel(options) + 1
                chosen = avail;
            else
                chosen(options(intmp)) = true;
            end
        end
    end

    if avail(1) && chosen(1)
        %{
        plotE(Eflow{1}, Eflow{2}, 'title', '$ \mathrm{ RGflow: \ U(1)_{c1} \times U(1)_{c2} \times SU(2)_{sp}} $', 'FontSize', 15, ...
                                                'Emax',1.3,'legmax',8,'Qdiff',[0,0,0]);
        %}
        %{}
        plotE(Eflow{1}, Eflow{2}, 'title', '$ \mathrm{ RGflow: \ U(1)_{c1} \times U(1)_{c2} \times SU(2)_{sp1} \times SU(2)_{sp2}} $', 'FontSize', 13, ...
                                                'Emax',3,'legmax',13,'Qdiff',[0,0,0,0]);
        %}
    end

    if avail(2) && avail(3)

        names = {'$T^{z}$'};
        legends = cell(0,0);
        for it = (1:1)
            if isempty(ImpDyn{it})
                ImpDyn{it} = [];
            else
                legends = cat(2,legends,names{it});
            end
        end
        
        num_ImpDyn = numel(ImpDyn);         % number of impurity dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % log temperatures
        log_ImpDyn = cell(num_ImpDyn,1);            % log impurity dynamic sysceptibilities
        log_ImpDyn_1stDer = cell(num_ImpDyn,1);     % first derivatives
        log_T_1stDer = cell(num_ImpDyn,1);          % log temperatures for first derivatives
        log_ImpDyn_2ndDer = cell(num_ImpDyn,1);     % second derivatives
        log_T_2ndDer = cell(num_ImpDyn,1);          % log temperatures for second derivatives

        for it = (1:num_ImpDyn)
            tmp = ImpDyn{it};
            log_ImpDyn{it} = log(tmp(ocont>0))./log(10);

            log_ImpDyn_1stDer{it} = diff(log_ImpDyn{it},1)./diff(log_T,1);      % first derivative
            tmp = movmean(log_T, [0,1]);                                        
            log_T_1stDer{it} = tmp(1:end-1);                                    % log temperatures for first derivatives

            log_ImpDyn_2ndDer{it} = diff(log_ImpDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                    % log temperatures for second derivatives
        end

        if chosen(2)

            figure;
            
            hold on;
            linestyle = {'-', '-', '-'};
            for it = (1:num_ImpDyn)
                plot(ocont(ocont<1), ImpDyn{it}(ocont<1), 'LineWidth',2,'LineStyle',linestyle{it});
            end
            
            %xlim([1e-25,1]);
            %ylim([1e-7,1]);
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','best','FontSize',30);

            set(gca,'XScale','log','YScale','log','FontSize',25);
            xlabel('$\omega$','Interpreter','latex','FontSize',30);
            ylabel('$\chi_{\mathrm{imp}} (\omega)$','Interpreter','latex','FontSize',30);
            %title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities}$'],'Interpreter','latex','FontSize',30);
            %{}
            title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities} \ (\lambda_{z}, \lambda_{x}) = (', ...
                        sprintf('%.15g',K_z/2),', ',sprintf('%.15g',Q),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',25);
            %}

            %{}
            legend('AutoUpdate','off');
            log_ImpDyn = cat(1,log_ImpDyn,{ zeros(1,numel(log_T)) ; zeros(1,numel(log_T)) } );
            %log_ImpDyn = cat(1,log_ImpDyn,{log_ImpDyn{1}; zeros(1,numel(log_T))} );
            fit_range = [-6, -11; 0, 0; 0, 0];
            [a1,Rsq1,a2,Rsq2,~,~] = Insert(log_T, log_ImpDyn,fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2 - 0.5;
            text_y = polyval(a1,text_x) + 0.5;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1) + 0.15;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            %plot(x1,y1,'-','Color',[0, 0.447,0.741],'LineWidth',1);
            plot(x1,y1,'-','Color','black','LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',20);

            %{
            x2 = fit_range(2,:);
            text_x = (x2(1)+x2(2))/2 - 0.7;
            text_y = polyval(a2,text_x) + 1.8;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y2 = polyval(a2,x2)+0.5;
            x2 = power(10,x2);
            y2 = power(10,y2);
            text2 = ['$\omega^{',sprintf('%.2f',a2(1)),'}$'];
            plot(x2,y2,'-','Color','black','LineWidth',1);
            text(text_x, text_y, text2,'Interpreter','latex','FontSize',20);
            %}
            %}

            hold off;
        end

        if chosen(3)
            figure;
            hold on;
            linestyle = {'-', '-.', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_ImpDyn)
                X = log_T_2ndDer{it};
                Y = log_ImpDyn_2ndDer{it};
                plot(X, Y,'LineWidth',2,'LineStyle',linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_ImpDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([log(T)/log(10)-1,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(4) && avail(5)

        names = {'$\vec{J}_{\mathrm{sp1}}$','$\vec{J}_{\mathrm{sp2}}$','$\vec{J}_{\mathrm{Sptot}}$','$J_{\mathrm{orb}}^{z}$'};
        legends = cell(0,0);
        for it = (1:4)
            if isempty(BathDyn{it})
                BathDyn(it) = [];
            else
                legends = cat(2,legends,names{it});
            end
        end

        
        num_BathDyn = numel(BathDyn);         % Number of bath dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % Log temperatures
        log_BathDyn = cell(num_BathDyn,1);          % Log bath dynamic susceptibilities
        log_BathDyn_1stDer = cell(num_BathDyn,1);   % First derivatives
        log_T_1stDer = cell(num_BathDyn,1);         % Log Temperatures for fisrt derivatives
        log_BathDyn_2ndDer = cell(num_BathDyn,1);   % Second derivatives
        log_T_2ndDer = cell(num_BathDyn,1);         % Log Temperatures for second derivatives

        for it = (1:num_BathDyn)
            tmp = BathDyn{it};
            log_BathDyn{it} = log(tmp(ocont>0))./log(10);

            log_BathDyn_1stDer{it} = diff(log_BathDyn{it},1)./diff(log_T,1);                  % First derivative
            tmp = movmean(log_T, [0,1]);
            log_T_1stDer{it} = tmp(1:end-1);                                                % Log Temperatures for fisrt derivatives

            log_BathDyn_2ndDer{it} = diff(log_BathDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % Second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                % Log Temperatures for second derivatives
        end
        
        if chosen(4)

            figure;
            hold on;
            linestyle = {'-','-.','--',':'};
            for it = (1:num_BathDyn)
                plot(ocont(ocont < 1),BathDyn{it}(ocont < 1),linestyle{it},'Linewidth',2);
            end
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',30);
            set(gca,'XScale','log','YScale','log','fontsize',25);
            xlabel('$\omega$','Interpreter','latex','FontSize',30);
            ylabel('$\chi'''' (\omega)$','Interpreter','latex','FontSize',30);
            title(['$\mathrm{Bath \ Dynamic \ Susceptibilities} \ (\lambda_{z}, \lambda_{x}) = (', ...
                        sprintf('%.15g',K_z/2),', ',sprintf('%.15g',Q),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',25);

            %xlim([1e-25,1]);
            %ylim([1e-20,1e5]);
            %{}
            legend('AutoUpdate','off');
            fit_range = [-5, -10; -2.5, -4.5; -8, -13];
            [a1,Rsq1,a2,Rsq2,a3,Rsq3] = Insert(log_T, log_BathDyn(2:4),fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2 - 0.4;
            text_y = polyval(a1,text_x) - 1.5;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1) - 0.5;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            %plot(x1,y1,'-','Color',[0, 0.447,0.741],'LineWidth',1);
            plot(x1,y1,'-','Color','black','LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',20);

            x2 = fit_range(2,:);
            text_x = (x2(1)+x2(2))/2 - 0.7;
            text_y = polyval(a2,text_x) + 2.3;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y2 = polyval(a2,x2)+0.5;
            x2 = power(10,x2);
            y2 = power(10,y2);
            text2 = ['$\omega^{',sprintf('%.2f',a2(1)),'}$'];
            %plot(x2,y2,'-','Color','black','LineWidth',1);
            %plot(x2,y2,'-','Color',[0.85,0.325,0.098],'LineWidth',1);
            %text(text_x, text_y, text2,'Interpreter','latex','FontSize',20);

            x3 = fit_range(3,:);
            text_x = (x3(1)+x3(2))/2 - 0.3;
            text_y = polyval(a3,text_x) + 1.3;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y3 = polyval(a3,x3) + 0.5;
            x3 = power(10,x3);
            y3 = power(10,y3);
            text3 = ['$\omega^{',sprintf('%.2f',a3(1)),'}$'];
            %plot(x3,y3,'-','Color',[0.929,0.694,0.125],'LineWidth',1);
            plot(x3,y3,'-','Color','black','LineWidth',1);
            text(text_x, text_y, text3,'Interpreter','latex','FontSize',20);
            legend('AutoUpdate','on','location','southwest');
            %}

            hold off;
        end

        if chosen(5)
            figure;
            hold on;
            linestyle = {'-', '-', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_BathDyn)
                plot(log_T_2ndDer{it}, log_BathDyn_2ndDer{it}, 'LineWidth', 2, 'LineStyle', linestyle{it});
            end
            legend('AutoUpdate','off');
            for it = (1:num_BathDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_BathDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([-23,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Bath \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end

    if avail(6) && chosen(6)

        figure;
        hold on;
        Sent_imp = exp(Sent_imp);
        plot(Temps,Sent_imp,'Linewidth',2);
        xlim([1e-24,1]);
        ylim([1,2.1])
        %plot([1e-22,1],[sqrt(8),sqrt(8)],'--','LineWidth',1.5,'color',[.7,.7,.7]);
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'tex';
        set(gca, 'TickLabelInterpreter', 'latex');
        ticks = {'1','2'};
        yticklabels(ticks);
        yticks([1,2]);

        ax = gca;
        ax.XAxis.FontSize = 5;
        ax.YAxis.FontSize = 5;
        set(gca,'XScale','log','YScale','linear','fontsize',30);
        xlabel('T','Interpreter','latex','FontSize',30);
        ylabel('$\mathrm{exp}(S_{\mathrm{imp}})$','Interpreter','latex','FontSize',30);
        title('$\mathrm{Impurity \ contribution \ to \ entropy}$','Interpreter','latex','FontSize',35);
        %{
        title(['$\mathrm{Impurity \ contribution \ to \ entropy} \ (K_{z}, Q) = (', ...
                        sprintf('%.15g',K_z),', ',sprintf('%.15g',Q),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
        %}
        hold off;
    end

    if avail(7)

        Sp_corr = cell2mat(Sp_corr);
        Sp_odd = Sp_corr(1:2:end);
        Sp_even = Sp_corr(2:2:end);
        Sp_avg = [];
        Sp_cum = [];
        Sp_cum_avg = [];

        for it = (2:numel(Sp_corr)-1)
            Sp_avg = [Sp_avg, (2*Sp_corr(it) + Sp_corr(it+1) + Sp_corr(it-1))/4];
        end

        for it = (1:numel(Sp_corr))
            Sp_cum(it) = sum(Sp_corr(1:it));
        end

        for it = (2:numel(Sp_corr)-1)
            Sp_cum_avg = [Sp_cum_avg, (2*Sp_cum(it) + Sp_cum(it+1) + Sp_cum(it-1))/4];
        end

        X = (1:numel(Sp_corr));
        X_odd = (1:2:numel(Sp_corr));
        X_even = (2:2:numel(Sp_corr));
        X = 4.^(-X/2);
        X_odd = 4.^(-X_odd/2);
        X_even = 4.^(-X_even/2);

        if chosen(7)

            figure;
            hold on;
            
            plot(X_odd,Sp_odd,'--','LineWidth',1.5);
            plot(X_even,Sp_even,'--','LineWidth',1.5);
            plot(X(2:end-1),Sp_avg,'-','LineWidth',1.5);
    
            xlim([X(end),X(1)]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 15;
            ax.YAxis.FontSize = 15;
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle$','Interpreter','latex','FontSize',25);

            legend({'$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle \mathrm{odd}$', '$\rangle \vec{S}_{imp} \cdot \vec{S_{m}} \rangle \mathrm{even}$', ...
                            '$\langle \vec{S}_{\mathrm{imp}} \cdot \vec{S_{\mathrm{imp}}} \rangle \mathrm{average}$'},'Location','southeast','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Spin-Spin \ Correlators} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(9)

            figure;
            hold on;
            
            plot(X_odd,Sp_cum(1:2:end),'--','LineWidth',1.5);
            plot(X_even,Sp_cum(2:2:end),'--','LineWidth',1.5);
            plot(X(2:end-1),Sp_cum_avg,'-','LineWidth',1.5);
    
            xlim([X(end),X(1)]);
            plot([X(end),X(1)],[-0.75,-0.75],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle \vec{S}_{\mathrm{imp}} \cdot \left( \sum_{i=1}^{m} \vec{S_{i}} \right) \rangle$','Interpreter','latex','FontSize',25);

            set(gca, 'TickLabelInterpreter', 'latex');
            ticks = {'-0.75','-0.6','-0.4','-0.2','0'};
            yticklabels(ticks);
            yticks([-0.75,-0.6,-0.4,-0.2,0]);

            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Cumulated \ Spin-Spin \ Correlators} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

    end

    if avail(8)

        Orb_corr = cell2mat(Orb_corr);
        Orb_odd = Orb_corr(1:2:end);
        Orb_even = Orb_corr(2:2:end);
        Orb_avg = [];
        Orb_cum = [];
        Orb_cum_avg = [];

        for it = (2:numel(Orb_corr)-1)
            Orb_avg = [Orb_avg, (2*Orb_corr(it) + Orb_corr(it+1) + Orb_corr(it-1))/4];
            %Orb_cum(it-1) = sum(Orb_avg(1:it-1));
        end

        for it = (1:numel(Orb_corr))
            Orb_cum(it) = sum(Orb_corr(1:it));
        end

        for it = (2:numel(Orb_corr)-1)
            Orb_cum_avg = [Orb_cum_avg, (2*Orb_cum(it) + Orb_cum(it+1) + Orb_cum(it-1))/4];
        end

        X = (1:numel(Orb_corr));
        X_odd = (1:2:numel(Orb_corr));
        X_even = (2:2:numel(Orb_corr));
        X = 4.^(-X/2);
        X_odd = 4.^(-X_odd/2);
        X_even = 4.^(-X_even/2);

        if chosen(8)

            figure;
            hold on;
    
            plot(X_odd,Orb_odd,'--','LineWidth',1.5);
            plot(X_even,Orb_even,'--','LineWidth',1.5);
            plot(X(2:end-1),Orb_avg,'-','LineWidth',1.5);
    
            xlim([X(end),X(1)]);
            %xlim([1e-16,X(1)]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 15;
            ax.YAxis.FontSize = 15;
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$< T_{\mathrm{imp}}^{z} T_{\mathrm{imp}}^{z} >$','Interpreter','latex','FontSize',25);
            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Orbital-Orbital \ Correlators} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end

        if chosen(10)

            figure;
            hold on;
            
            plot(X_odd,Orb_cum(1:2:end),'--','LineWidth',1.5);
            plot(X_even,Orb_cum(2:2:end),'--','LineWidth',1.5);
            plot(X(2:end-1),Orb_cum_avg,'-','LineWidth',1.5);
    
            xlim([X(end),X(1)]);
            %xlim([1e-10,X(1)]);
            %ylim([-0.26,0.01]);
            ylim([-0.25,0.01]);
            %plot([X(end),X(1)],[-0.25,-0.25],'--','LineWidth',1,'color',[.7,.7,.7]);
            plot([X(end),X(1)],[0,0],'--','LineWidth',1,'color',[.7,.7,.7]);
    
            ax=gca;
            ax.XAxis.FontSize = 25;
            ax.YAxis.FontSize = 25;
            %set(gca,'XScale','linear','YScale','linear');
            set(gca,'XScale','log','YScale','linear','Xdir','reverse');
            xlabel('$\Lambda^{-m/2} D$','Interpreter','latex','FontSize',25);
            ylabel('$\langle T_{\mathrm{imp}} \left( \sum_{i=1}^{m} T_{i} \right) \rangle$','Interpreter','latex','FontSize',25);

            set(gca, 'TickLabelInterpreter', 'latex');
            %{}
            ticks = {'-0.25','-0.2','-0.1','0'};
            yticklabels(ticks);
            yticks([-0.25,-0.2,-0.1,0]);
            %}

            legend({'$\mathrm{odd}$', '$\mathrm{even}$', '$\mathrm{average}$'},'Location','best','Interpreter','latex','FontSize',25);
            title(['$\mathrm{ Cumulated \ Orbz-Orbz \ Correlators} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                            sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
            hold off;
        end
    end


elseif intmp == 7   % Anderson_8flav_NRG
    
    FileInfo = dir(path);
    strtmp = cell(0,0);
    cnt = 1;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1) == 'H' || DirName(1) == 'T'
            strtmp = cat(2, strtmp, [sprintf('%.15g',cnt),': ',DirName]);
            cnt = cnt + 1;
        end
    end
    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;

    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    % extract parameters from folder name
    tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': Hyb=%f_U=%f_J=%f_N0=%f_T=%f_Nkeep=%f']);
    Hyb = tmp(1);
    U = tmp(2);
    J = tmp(3);
    N0 = tmp(4);
    T = tmp(5);

    tmp = strtmp{idx};
    path = [path, filesep, tmp(numel(num2str(idx))+3:end)];

    % make a list of unnecessary file infos
    FileInfo = dir(path);
    EraseIdx = [];
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')                 % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')            % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'8flav')         % Erase log files
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];    % Erase selected unnecessary file infos

    Eflow = cell(2,1);      % RG flow data. 1: Etot, 2: Qtot
    ImpDyn = cell(8,1);     % Impurity dynamical susceptibilities. 1: ImpSp1m, 2: ImpSp1p, 3: ImpSp2m, 4: ImpSp2p, 5: Charge1m, 6: Charge1p, 7: Charge2m, 8: Charge2p
    BathDyn = cell(0,1);    % Bath dynamic susceptibilities. 
    avail = false*ones(1,4);
    % RGflow, ImpDyn, Second derivatives of ImpDyn, Impurity contribution to entropy

    % Check availabel and unavailable data
    for it = (1:numel(FileInfo))
        tmp = load([FileInfo(it).folder, filesep, FileInfo(it).name]);
        field = fieldnames(tmp);

        switch FileInfo(it).name
            case 'Etot.mat'
                Eflow{1} = getfield(tmp, field{1});
                avail(1) = true;
            case 'Qtot.mat'
                Eflow{2} = getfield(tmp, field{1});
                avail(1) = true;
            case 'ocont.mat'
                ocont = getfield(tmp, field{1});
            case 'NRG_Op=ImpSp1m.mat'
                ImpDyn{1} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpSp1p.mat'
                ImpDyn{2} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpSp2m.mat'
                ImpDyn{3} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpSp2p.mat'
                ImpDyn{4} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=Charge1m.mat'
                ImpDyn{5} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=Charge1p.mat'
                ImpDyn{6} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=Charge2m.mat'
                ImpDyn{7} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=Charge2p.mat'
                ImpDyn{8} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'Sent_imp.mat'
                Sent_imp = getfield(tmp, field{1});
                avail(4) = true;
            case 'Temps.mat'
                Temps = getfield(tmp, field{1});
            otherwise
                fprintf('WRN: unknown data type');
        end
    end

    strtmp = cell(1,0);
    options = [];
    if avail(1)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': RG flow diagram']});
        options = [options,1];
    end
    if avail(2)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity dynamic susceptibilities']});
        options = [options,2];
    end
    if avail(3)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of impurity dynamic susceptibilities']});
        options = [options,3];
    end
    if avail(4)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity contribution to entropy']});
        options = [options,4];
    end
    strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': All of the above']});
    dispbox('-width',130,strtmp{:});

    fprintf('Which one do you want to plot?\n');
    intmp = input('>>> ');

    if intmp == numel(options)+1    % plot all
        chosen = avail;
    else
        chosen = false*ones(1,4);
        chosen(options(intmp)) = true;
        while ~isempty(intmp) && ~isequal(chosen, avail)
            fprintf('Type the index of data you want to plot\n');
            fprintf('(Press enter to finish)\n');
            intmp = input('>>> ');
            if intmp == numel(options) + 1
                chosen = avail;
            else
                chosen(options(intmp)) = true;
            end
        end
    end

    if avail(1) && chosen(1)
        %{
        plotE(Eflow{1}, Eflow{2}, 'title', '$ \mathrm{ RGflow: \ U(1)_{c1} \times U(1)_{c2} \times SU(2)_{sp}} $', 'FontSize', 15, ...
                                                'Emax',1.3,'legmax',8,'Qdiff',[0,0,0]);
        %}
        %{}
        plotE(Eflow{1}, Eflow{2}, 'title', '$ \mathrm{ RGflow: \ U(1)_{c1+} \times U(1)_{c1-} \times U(1)_{c2+} \times U(1)_{c2-} \times SU(2)_{sp1} \times SU(2)_{sp2}} $', 'FontSize', 13, ...
                                                'Emax',3,'legmax',13,'Qdiff',[0,0,0,0,0,0]);
        %}
    end

    if avail(2) && avail(3)

        names = {'$\chi_{\mathrm{sp1+}}$', '$\chi_{\mathrm{sp1-}}$', '$\chi_{\mathrm{sp2+}}$', '$\chi_{\mathrm{sp2-}}$', ...
                    '$\chi_{\mathrm{c1+}}$', '$\chi_{\mathrm{c1-}}$', '$\chi_{\mathrm{c2+}}$', '$\chi_{\mathrm{c2-}}$'};
        legends = cell(0,0);
        for it = (1:numel(names))
            if isempty(ImpDyn{it})
                ImpDyn{it} = [];
            else
                legends = cat(2,legends,names{it});
            end
        end
        
        num_ImpDyn = numel(ImpDyn);         % number of impurity dynamic susceptibilities

        log_T = log(ocont(ocont>0))./log(10);       % log temperatures
        log_ImpDyn = cell(num_ImpDyn,1);            % log impurity dynamic sysceptibilities
        log_ImpDyn_1stDer = cell(num_ImpDyn,1);     % first derivatives
        log_T_1stDer = cell(num_ImpDyn,1);          % log temperatures for first derivatives
        log_ImpDyn_2ndDer = cell(num_ImpDyn,1);     % second derivatives
        log_T_2ndDer = cell(num_ImpDyn,1);          % log temperatures for second derivatives

        disp(ImpDyn);
        disp(size(ocont));
        for it = (1:num_ImpDyn)
            tmp = ImpDyn{it};
            log_ImpDyn{it} = log(tmp(ocont>0))./log(10);

            log_ImpDyn_1stDer{it} = diff(log_ImpDyn{it},1)./diff(log_T,1);      % first derivative
            tmp = movmean(log_T, [0,1]);                                        
            log_T_1stDer{it} = tmp(1:end-1);                                    % log temperatures for first derivatives

            log_ImpDyn_2ndDer{it} = diff(log_ImpDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                    % log temperatures for second derivatives
        end

        if chosen(2)

            figure;
            
            hold on;
            linestyle = {'-', '-.', '--'};
            for it = (1:num_ImpDyn)
                plot(ocont(ocont<1), ImpDyn{it}(ocont<1), 'LineWidth',2,'LineStyle',linestyle{rem(it,3)+1});
            end
            
            %xlim([1e-25,1]);
            %ylim([1e-7,1]);
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','best','FontSize',30);

            set(gca,'XScale','log','YScale','log','FontSize',25);
            xlabel('$\omega$','Interpreter','latex','FontSize',30);
            ylabel('$\chi_{\mathrm{imp}} (\omega)$','Interpreter','latex','FontSize',30);
            %title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities}$'],'Interpreter','latex','FontSize',30);
            %{}
            title(['$ \mathrm{Impurity \ Dynamic \ Susceptibilities} \ (\lambda_{z}, \lambda_{x}) = (', ...
                        sprintf('%.15g',K_z/2),', ',sprintf('%.15g',Q),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',25);
            %}

            %{}
            legend('AutoUpdate','off');
            log_ImpDyn = cat(1,log_ImpDyn,{ zeros(1,numel(log_T)) ; zeros(1,numel(log_T)) } );
            %log_ImpDyn = cat(1,log_ImpDyn,{log_ImpDyn{1}; zeros(1,numel(log_T))} );
            fit_range = [-6, -11; 0, 0; 0, 0];
            [a1,Rsq1,a2,Rsq2,~,~] = Insert(log_T, log_ImpDyn,fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2 - 0.5;
            text_y = polyval(a1,text_x) + 0.5;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1) + 0.15;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            %plot(x1,y1,'-','Color',[0, 0.447,0.741],'LineWidth',1);
            plot(x1,y1,'-','Color','black','LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',20);

            %{
            x2 = fit_range(2,:);
            text_x = (x2(1)+x2(2))/2 - 0.7;
            text_y = polyval(a2,text_x) + 1.8;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y2 = polyval(a2,x2)+0.5;
            x2 = power(10,x2);
            y2 = power(10,y2);
            text2 = ['$\omega^{',sprintf('%.2f',a2(1)),'}$'];
            plot(x2,y2,'-','Color','black','LineWidth',1);
            text(text_x, text_y, text2,'Interpreter','latex','FontSize',20);
            %}
            %}

            hold off;
        end

        if chosen(3)
            figure;
            hold on;
            linestyle = {'-', '-.', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_ImpDyn)
                X = log_T_2ndDer{it};
                Y = log_ImpDyn_2ndDer{it};
                plot(X, Y,'LineWidth',2,'LineStyle',linestyle{rem(it,3)+1});
            end
            legend('AutoUpdate','off');
            for it = (1:num_ImpDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([log(T)/log(10)-1,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',25);
            set(gca,'XScale','linear','YScale','linear','fontsize',20);
            xlabel('$\log T$','Interpreter','latex','FontSize',25);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',25);
            title(['$\mathrm{2nd \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities} \ (J_{0}, K_{\perp}, K_z, I_{0}) = (', ...
                        sprintf('%.15g',J0),', ',sprintf('%.15g',K_perp),', ',sprintf('%.15g',K_z),', ',sprintf('%.15g',I0),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',15);
            hold off;
        end

    end
    
    if avail(6) && chosen(6)

        figure;
        hold on;
        Sent_imp = exp(Sent_imp);
        plot(Temps,Sent_imp,'Linewidth',2);
        xlim([1e-24,1]);
        ylim([1,2.1])
        %plot([1e-22,1],[sqrt(8),sqrt(8)],'--','LineWidth',1.5,'color',[.7,.7,.7]);
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'tex';
        set(gca, 'TickLabelInterpreter', 'latex');
        ticks = {'1','2'};
        yticklabels(ticks);
        yticks([1,2]);

        ax = gca;
        ax.XAxis.FontSize = 5;
        ax.YAxis.FontSize = 5;
        set(gca,'XScale','log','YScale','linear','fontsize',30);
        xlabel('T','Interpreter','latex','FontSize',30);
        ylabel('$\mathrm{exp}(S_{\mathrm{imp}})$','Interpreter','latex','FontSize',30);
        title('$\mathrm{Impurity \ contribution \ to \ entropy}$','Interpreter','latex','FontSize',35);
        %{
        title(['$\mathrm{Impurity \ contribution \ to \ entropy} \ (K_{z}, Q) = (', ...
                        sprintf('%.15g',K_z),', ',sprintf('%.15g',Q),'), T=10^{',sprintf('%d',round(log(T)/log(10))),'}$'],'Interpreter','latex','FontSize',20);
        %}
        hold off;
    end


elseif intmp == 8   % ToAH_NRG
    
    FileInfo = dir(path);
    strtmp = cell(0,0);
    cnt = 1;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if numel(DirName) >=3
            if isequal(DirName(1:3), 'Hyb')
                strtmp = cat(2, strtmp, [sprintf('%.15g',cnt),': ',DirName]);
                cnt = cnt + 1;
            end
        end
    end
    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;

    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    % extract parameters from folder name
    tmp = sscanf(strtmp{idx}, [sprintf('%d',idx),': Hyb=%f_U=%f_J=%f_J_L=%f_T=%f_mu=%f_Nkeep=%f_Lambda=%f']);
    Hyb = tmp(1);
    U = tmp(2);
    J = tmp(3);
    J_L = tmp(4);
    T = tmp(5);
    mu = tmp(6);
    Nkeep = tmp(7);
    Lambda = tmp(8);

    tmp = strtmp{idx};
    path = [path, filesep, tmp(numel(num2str(idx))+3:end)];

    % make a list of unnecessary file infos
    FileInfo = dir(path);
    EraseIdx = [];
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')                 % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')            % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif numel(DirName) >= 7              % Erase log files
            if isequal(DirName(1:4), 'Susc') || isequal(DirName(1:7), 'MuSweep') || isequal(DirName(1:7), 'Entropy')
                EraseIdx = cat(2,EraseIdx,it);
            end
        end
    end
    FileInfo(EraseIdx) = [];    % Erase selected unnecessary file infos

    Eflow = cell(2,1);      % RG flow data. 1: Etot, 2: Qtot
    ImpDyn = cell(6,1);     % Impurity dynamical susceptibilities. 
                            % 1: ImpCharge_minus, 2: ImpCharge_plus, 3: ImpCharge_tot, 4: ImpSp, 5: ImpOrb_plus, 6: ImpOrb_plus
    BathDyn = cell(0,1);    % Bath dynamic susceptibilities. 
    avail = false(1,4);     % RGflow, ImpDyn, Second derivatives of ImpDyn, Impurity contribution to entropy

    % Check availabel and unavailable data
    for it = (1:numel(FileInfo))
        tmp = load([FileInfo(it).folder, filesep, FileInfo(it).name]);
        field = fieldnames(tmp);

        switch FileInfo(it).name
            case 'Etot.mat'
                Eflow{1} = getfield(tmp, field{1});
                avail(1) = true;
            case 'Qtot.mat'
                Eflow{2} = getfield(tmp, field{1});
                avail(1) = true;
            case 'ocont.mat'
                ocont = getfield(tmp, field{1});
            case 'NRG_Op=ImpCharge_minus.mat'
                ImpDyn{1} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpCharge_plus.mat'
                ImpDyn{2} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpCharge_tot.mat'
                ImpDyn{3} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpSp.mat'
                ImpDyn{4} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpOrb_plus.mat'
                ImpDyn{5} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'NRG_Op=ImpOrb_z.mat'
                ImpDyn{6} = getfield(tmp, field{1});
                avail(2) = true;
                avail(3) = true;
            case 'EntData.mat'
                EntData = getfield(tmp, field{1});
                avail(4) = true;
            case 'nu.mat'
                nu = getfield(tmp, field{1});
                fprintf('Filling = %.3f', sum(nu));
            otherwise
                fprintf('WRN: unknown data type');
        end
    end

    strtmp = cell(1,0);
    options = [];
    if avail(1)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': RG flow diagram']});
        options = [options,1];
    end
    if avail(2)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity dynamic susceptibilities']});
        options = [options,2];
    end
    if avail(3)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Second derivatives of impurity dynamic susceptibilities']});
        options = [options,3];
    end
    if avail(4)
        strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': Impurity contribution to entropy']});
        options = [options,4];
    end
    strtmp = cat(2, strtmp, {[sprintf('%d',numel(strtmp)+1),': All of the above']});
    dispbox('-width',130,strtmp{:});

    fprintf('Which one do you want to plot?\n');
    intmp = input('>>> ');

    if intmp == numel(options)+1    % plot all
        chosen = avail;
    else
        chosen = false*ones(1,4);
        chosen(options(intmp)) = true;
        while ~isempty(intmp) && ~isequal(chosen, avail)
            fprintf('Type the index of data you want to plot\n');
            fprintf('(Press enter to finish)\n');
            intmp = input('>>> ');
            if intmp == numel(options) + 1
                chosen = avail;
            else
                chosen(options(intmp)) = true;
            end
        end
    end

    if avail(1) && chosen(1)
        plotE(Eflow{1}, Eflow{2}, 'title', '$ \mathrm{ RGflow: \ U(1)_{c+} \times U(1)_{c-} \times SU(2)_{sp}} $', ...
                                                'FontSize', 13, 'Emax',3,'legmax',13,'Qdiff',[0,0,0]);
    end

    if avail(2) && avail(3)

        %names = {'$\chi^{-}_{\mathrm{c}}$', '$\chi^{+}_{\mathrm{c}}$', '$\chi^{\mathrm{tot}}_{\mathrm{c}}$', ...
        %            '$\chi_{\mathrm{sp}}$', '$\chi^{+}_{\mathrm{orb}}$', '$\chi^{\mathrm{z}}_{\mathrm{orb}}$'};
        names = {'$\chi_{\mathrm{c}}$', '$\chi_{\mathrm{sp}}$', '$\chi_{\mathrm{orb}}$'};

        legends = cell(0,0);
        for it = (1:numel(names))
            if isempty(ImpDyn{it})
                ImpDyn{it} = [];
            else
                legends = cat(2,legends,names{it});
            end
        end
        
        num_ImpDyn = numel(ImpDyn);         % number of impurity dynamic susceptibilities

        log_T = log10(ocont(ocont>0));       % log temperatures
        log_ImpDyn = cell(num_ImpDyn,1);            % log impurity dynamic sysceptibilities
        log_ImpDyn_1stDer = cell(num_ImpDyn,1);     % first derivatives
        log_T_1stDer = cell(num_ImpDyn,1);          % log temperatures for first derivatives
        log_ImpDyn_2ndDer = cell(num_ImpDyn,1);     % second derivatives
        log_T_2ndDer = cell(num_ImpDyn,1);          % log temperatures for second derivatives

        for it = 1:num_ImpDyn
            tmp = ImpDyn{it};
            log_ImpDyn{it} = log10(tmp(ocont>0));

            log_ImpDyn_1stDer{it} = diff(log_ImpDyn{it},1)./diff(log_T,1);      % first derivative
            tmp = movmean(log_T, [0,1]);                                        
            log_T_1stDer{it} = tmp(1:end-1);                                    % log temperatures for first derivatives

            log_ImpDyn_2ndDer{it} = diff(log_ImpDyn_1stDer{it},1)./diff(log_T_1stDer{it},1);    % second derivative
            tmp = movmean(log_T_1stDer{it}, [0,1]);
            log_T_2ndDer{it} = tmp(1:end-1);                                                    % log temperatures for second derivatives
        end

        if chosen(2)

            figure;
            hold on;
            xlim([1e-4, 1e2]);

            linestyle = {'--', '-', ':'};
            for it = 3:5 %(1:num_ImpDyn)
                plot(ocont, ImpDyn{it}, 'LineWidth',2,'LineStyle',linestyle{rem(it,3)+1});
            end
            
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            legend(legends,'Interpreter','latex','Location','northwest','FontSize',24);

            set(gca,'XScale','log','YScale','log','FontSize',18);
            xlabel('$\omega$','Interpreter','latex','FontSize',24);
            ylabel('$\chi_{\mathrm{imp}} (\omega)$','Interpreter','latex','FontSize',24);
            title(['$(U, J, J_{L}, \Delta'''') = (', ...
                        sprintf('%.15g',U),', ',sprintf('%.15g',J),', ',sprintf('%.15g',J_L),', ',sprintf('%.15g',Hyb), ...
                            '), T=10^{',sprintf('%d',round(log10(T))),'}$'],'Interpreter','latex','FontSize',21);
         
            legend('AutoUpdate','off');

            log_ImpDyn = cat(1,log_ImpDyn,{ zeros(1,numel(log_T)) ; zeros(1,numel(log_T)) } );
            fit_range = [-6, -11; -6, -11; -6, -11];
            [a1,Rsq1,a2,Rsq2,~,~] = Insert(log_T, log_ImpDyn, fit_range);

            x1 = fit_range(1,:);
            text_x = (x1(1)+x1(2))/2 - 0.5;
            text_y = polyval(a1,text_x) + 1.2;
            text_x = power(10, text_x);
            text_y = power(10, text_y);
            y1 = polyval(a1,x1) + 0.15;
            x1 = power(10,x1);
            y1 = power(10,y1);
            text1 = ['$\omega^{',sprintf('%.2f',a1(1)),'}$'];
            plot(x1,y1,'-','Color','black','LineWidth',1);
            text(text_x, text_y, text1,'Interpreter','latex','FontSize',18);

            hold off;
        end

        if chosen(3)

            figure;
            hold on;
            linestyle = {'-', '-.', '--'};
            legend('AutoUpdate','on');
            for it = (1:num_ImpDyn)
                X = log_T_2ndDer{it};
                Y = log_ImpDyn_2ndDer{it};
                plot(X, Y,'LineWidth',2,'LineStyle',linestyle{rem(it,3)+1});
            end
            legend('AutoUpdate','off');
            for it = (1:num_ImpDyn)
                [Minima, ~, MinPos] = LocMin(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                [Maxima, ~, MaxPos] = LocMax(log_T_2ndDer{it}, log_ImpDyn_2ndDer{it});
                plot(MinPos, Minima, 'o', 'Color', 'blue', 'LineWidth', 2);
                plot(MaxPos, Maxima, 'o', 'Color', 'red', 'LineWidth', 2);
            end
            legend('AutoUpdate','on');
            ax = gca;
            ax.XAxis.FontSize = 5;
            ax.YAxis.FontSize = 5;
            xlim([log(T)/log(10)-1,3]);
            ylim([-3,4]);
            legend(legends,'Interpreter','latex','Location','northeast','FontSize',24);
            set(gca,'XScale','linear','YScale','linear','fontsize',18);
            xlabel('$\log T$','Interpreter','latex','FontSize',24);
            ylabel('$\frac{d^{2} \log \chi''''}{d^{2} \log T} (\omega)$','Interpreter','latex','FontSize',24);
            title(['$\mathrm{2nd \ Derivatives \ of \ Impurity \ Dynamic \ Susceptibilities} \ (U, J, J_{L}, \Delta'''') = (', ...
                        sprintf('%.15g',U),', ',sprintf('%.15g',J),', ',sprintf('%.15g',J_L),', ',sprintf('%.15g',Hyb), ...
                            '), T=10^{',sprintf('%d',round(log10(T))),'}$'],'Interpreter','latex','FontSize',24);
            hold off;
        end

    end
    
    if avail(4) && chosen(4)

        fprintf('Choose beta for plotting impurity contribution to entropy');
        strtmp = arrayfun(@(x,i) sprintf('%d: %g', i, x), EntData.beta, 1:numel(EntData.beta), 'UniformOutput', false);
        dispbox('-width',130,strtmp{:});
        idx = input('>>> ');

        Temps = EntData.Temps{idx};
        Sent_imp = EntData.S_imp{idx};

        figure;
        hold on;
        Sent_imp = exp(Sent_imp);
        plot(Temps,Sent_imp,'Linewidth',2);
        xlim([1e-24,1]);
        ylim([1,4.1])
        %plot([1e-22,1],[sqrt(8),sqrt(8)],'--','LineWidth',1.5,'color',[.7,.7,.7]);
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'tex';
        set(gca, 'TickLabelInterpreter', 'latex');
        ticks = {'1','2','3','4'};
        yticklabels(ticks);
        yticks([1,2,3,4]);

        ax = gca;
        ax.XAxis.FontSize = 5;
        ax.YAxis.FontSize = 5;
        set(gca,'XScale','log','YScale','linear','fontsize',30);
        xlabel('T','Interpreter','latex','FontSize',30);
        ylabel('$\mathrm{exp}(S_{\mathrm{imp}})$','Interpreter','latex','FontSize',30);
        title(['$S_{\mathrm{imp}} , \ (U, J, J_{L}, \Delta'''') = (', ...
                        sprintf('%.15g',U),', ',sprintf('%.15g',J),', ',sprintf('%.15g',J_L),', ',sprintf('%.15g',Hyb), ...
                            '), T=10^{',sprintf('%d',round(log10(T))),'}$'],'Interpreter','latex','FontSize',21);
        hold off;
    end


else    % LineSearch

    %% Choose Folder
    fprintf('Choose the index of a folder\n')
    strtmp = {'1: J0_search', '2: K0_search', '3: I0_search'}; 
    dispbox('-width',130,strtmp{:});
    
    ValidIdx = false;
    while ~ValidIdx
        idx = input('>>> \n');
        if ismember(idx,[1,2,3])
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end

    switch idx
        case 1
            path = [path,'\J0_search'];
            Bi_Par = 'J0';                  % Binary search parameter
        case 2
            path = [path,'\K0_search'];
            Bi_Par = 'K0';
        case 3
            path = [path,'\I0_search'];
            Bi_Par = 'I0';
    end

    %% Choose Binary Search Dataset
    FileInfo = dir(path);
    EraseIdx = [];                     % indices of unecessary file infos to be erased
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')         % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')    % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];       % Erase selected unecessary file info

    strtmp = cell(0,0);
    cnt = 1;
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1:2) == Bi_Par
            strtmp = cat(2,strtmp,[sprintf('%.15g',cnt),': ',DirName]);
            cnt = cnt + 1;
        end
    end

    dispbox('-width',130,strtmp{:});
    fprintf('Choose the index of the data set to plot\n');

    ValidIdx = false;
    while ~ValidIdx
        idx = input('>>> ');
        if ismember(idx,(1:numel(strtmp)))
            ValidIdx = true;
        else
            fprintf('WRN: Invalid input\n');
        end
    end
    chosenfd = erase(strtmp{idx},[sprintf('%d',idx),': ']);
    path = [path,'\',chosenfd];

    switch Bi_Par       % get values of fixed parameters
        case 'J0'
            tmp = sscanf(chosenfd, 'J0_search_%f_to_%f_K0=%f_I0=%f');
            K0 = tmp(1);
            I0 = tmp(2);
        case 'K0'
            tmp = sscanf(chosenfd, 'K0_search_%f_to_%f_J0=%f_I0=%f');
            J0 = tmp(1);
            I0 = tmp(2);
        case 'I0'
            tmp = sscanf(chosenfd, 'I0_search_%f_to_%f_J0=%f_K0=%f');
            J0 = tmp(1);
            K0 = tmp(2);
    end

    %%  Choose Final RG Flow Data to Plot
    FileInfo = dir(path);
    EraseIdx = [];                     % indices of unecessary file infos to be erased
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if isequal(DirName,'.')         % Erase '.'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName,'..')    % Erase '..'
            EraseIdx = cat(2,EraseIdx,it);
        elseif isequal(DirName(1:2),Bi_Par)     % Erase log files
            EraseIdx = cat(2,EraseIdx,it);
        end
    end
    FileInfo(EraseIdx) = [];       % Erase selected unecessary file infos

    strtmp = cell(0,0);
    Etot = cell(0,0);
    Qtot = cell(0,0);

    switch Bi_Par       % Choose RG flow data to plot. Treat Bi_Par = 'J0', 'K0', 'I0' cases separately
        case 'J0'       % Bi_Par = 'J0'

            J0 = [];
            for it = (1:numel(FileInfo))        % Read all available RG flow data
                FileName = FileInfo(it).name;
                if isequal(FileName(1:4), 'Etot')
                    tmp = sscanf(FileName, 'Etot_J0=%f.mat');
                    J0 = [J0,tmp];
                    strtmp = cat(2, strtmp, [sprintf('%d',it), ': J0 = ', sprintf('%.15g',tmp)]);
                end
            end
            strtmp = cat(2, strtmp, [sprintf('%d', numel(J0)+1), ': All of above']);
            dispbox('-width',130,strtmp{:});

            StopAsk = false;
            cnt = 0;
            while ~StopAsk && cnt < numel(J0)       % Choose RG flow data to plot
                fprintf(['Type the index of the RG flow data you want to plot. (K0, I0) = (', ...
                        sprintf('%.15g',K0), ', ', sprintf('%.15g',I0),')\n']);
                fprintf('(Press enter to finish)\n');
                idx = input('>>> ');

                if idx == numel(J0)+1     % Plotting all data
                    StopAsk = true;
                    cnt = numel(J0);

                    for it = (1:cnt)        % Loading all the data in folder

                        tmp = load([path, '\Etot_J0=', sprintf('%.15g',J0(it)), '.mat']);
                        field = fieldnames(tmp);
                        Etot = cat(2, Etot, getfield(tmp, field{1}));
                        tmp = load([path, '\Qtot_J0=', sprintf('%.15g',J0(it)), '.mat']);
                        field = fieldnames(tmp);
                        Qtot = cat(2, Qtot, getfield(tmp, field{1}));
                    end

                elseif ismember(idx, (1:numel(strtmp)-1))   % Certain data chosen
                    cnt = cnt + 1;
                    tmp = load([path, '\Etot_J0=', sprintf('%.15g',J0(idx)), '.mat']);
                    field = fieldnames(tmp);
                    Etot = cat(2, Etot, getfield(tmp, field{1}));
                    tmp = load([path, '\Qtot_J0=', sprintf('%.15g',J0(idx)), '.mat']);
                    field = fieldnames(tmp);
                    Qtot = cat(2, Qtot, getfield(tmp, field{1}));

                elseif isempty(idx)                         % Pressed enter to finish
                    StopAsk = true;
                else
                    fprintf('WRN: Invalid Input\n');
                end
            end

        case 'K0'       % Bi_Par = 'K0'

            K0 = [];
            for it = (1:numel(FileInfo))        % Read all available RG flow data
                FileName = FileInfo(it).name;
                if isequal(FileName(1:4), 'Etot')
                    tmp = sscanf(FileName, 'Etot_K0=%f.mat');
                    K0 = [K0,tmp];
                    strtmp = cat(2, strtmp, [sprintf('%d',it), ': K0 = ', sprintf('%.15g',tmp)]);
                end
            end
            strtmp = cat(2, strtmp, [sprintf('%d', numel(K0)+1), ': All of above']);
            dispbox('-width',130,strtmp{:});

            StopAsk = false;
            cnt = 0;
            while ~StopAsk && cnt < numel(K0)       % Choose RG flow data to plot
                fprintf(['Type the index of the RG flow data you want to plot. (J0, I0) = (', ...
                        sprintf('%.15g',J0), ', ', sprintf('%.15g',I0),')\n']);
                fprintf('(Press enter to finish)\n');
                idx = input('>>> ');

                if idx == numel(K0) + 1     % Plotting all data
                    StopAsk = true;
                    cnt = numel(K0);

                    for it = (1:cnt)        % Loading all the data in folder

                        tmp = load([path, '\Etot_K0=', sprintf('%.15g',K0(it)), '.mat']);
                        field = fieldnames(tmp);
                        Etot = cat(2, Etot, getfield(tmp, field{1}));
                        tmp = load([path, '\Qtot_K0=', sprintf('%.15g',K0(it)), '.mat']);
                        field = fieldnames(tmp);
                        Qtot = cat(2, Qtot, getfield(tmp, field{1}));
                    end

                elseif ismember(idx, (1:numel(strtmp)-1))   % Certain data chosen
                    cnt = cnt + 1;
                    tmp = load([path, '\Etot_K0=', sprintf('%.15g',K0(idx)), '.mat']);
                    field = fieldnames(tmp);
                    Etot = cat(2, Etot, getfield(tmp, field{1}));
                    tmp = load([path, '\Qtot_K0=', sprintf('%.15g',K0(idx)), '.mat']);
                    field = fieldnames(tmp);
                    Qtot = cat(2, Qtot, getfield(tmp, field{1}));

                elseif isempty(idx)                         % Pressed enter to finish
                    StopAsk = true;
                else
                    fprintf('WRN: Invalid Input\n');
                end
            end

        case 'I0'       % Bi_Par = 'I0'

            I0 = [];
            for it = (1:numel(FileInfo))        % Read all available RG flow data
                FileName = FileInfo(it).name;
                if isequal(FileName(1:4), 'Etot')
                    tmp = sscanf(FileName, 'Etot_I0=%f.mat');
                    I0 = [I0,tmp];
                    strtmp = cat(2, strtmp, [sprintf('%d',it), ': I0 = ', sprintf('%.15g',tmp)]);
                end
            end
            strtmp = cat(2, strtmp, [sprintf('%d', numel(I0)+1), ': All of above']);
            dispbox('-width',130,strtmp{:});

            StopAsk = false;
            cnt = 0;
            while ~StopAsk && cnt < numel(I0)       % Choose RG flow data to plot
                fprintf(['Type the index of the RG flow data you want to plot. (J0, K0) = (', ...
                        sprintf('%.15g',J0), ', ', sprintf('%.15g',K0),')\n']);
                fprintf('(Press enter to finish)\n');
                idx = input('>>> ');

                if idx == numel(I0) + 1     % Plotting all data
                    StopAsk = true;
                    cnt = numel(I0);

                    for it = (1:cnt)        % Loading all the data in folder

                        tmp = load([path, '\Etot_I0=', sprintf('%.15g',I0(it)), '.mat']);
                        field = fieldnames(tmp);
                        Etot = cat(2, Etot, getfield(tmp, field{1}));
                        tmp = load([path, '\Qtot_I0=', sprintf('%.15g',I0(it)), '.mat']);
                        field = fieldnames(tmp);
                        Qtot = cat(2, Qtot, getfield(tmp, field{1}));
                    end

                elseif ismember(idx, (1:numel(strtmp)-1))   % Certain data chosen
                    cnt = cnt + 1;
                    tmp = load([path, '\Etot_I0=', sprintf('%.15g',I0(idx)), '.mat']);
                    field = fieldnames(tmp);
                    Etot = cat(2, Etot, getfield(tmp, field{1}));
                    tmp = load([path, '\Qtot_I0=', sprintf('%.15g',I0(idx)), '.mat']);
                    field = fieldnames(tmp);
                    Qtot = cat(2, Qtot, getfield(tmp, field{1}));

                elseif isempty(idx)                         % Pressed enter to finish
                    StopAsk = true;
                else
                    fprintf('WRN: Invalid Input\n');
                end
            end

    end

    for it = (1:cnt)        % Plot all chosen RG flow data
        switch Bi_Par
            case 'J0'
                plotE(Etot(:,it),Qtot(:,it),'title',['J0=',sprintf('%.15g',J0(it)),' K0=',...
                      sprintf('%.15g',K0),' I0=',sprintf('%.15g',I0)],'Emax',3,'legmax',10,'Qdiff',[0,0,0]);
            case 'K0'
                plotE(Etot(:,it),Qtot(:,it),'title',['J0=',sprintf('%.15g',J0),' K0=',...
                    sprintf('%.15g',K0(it)),' I0=',sprintf('%.15g',I0)],'Emax',3,'legmax',10,'Qdiff',[0,0,0]);
            case 'I0'
                plotE(Etot(:,it),Qtot(:,it),'title',['J0=',sprintf('%.15g',J0),' K0=',...
                    sprintf('%.15g',K0),' I0=',sprintf('%.15g',I0(it))],'Emax',3,'legmax',10,'Qdiff',[0,0,0]);
        end
    end


end