function AllParams = TsoK_selected_data()

    path = 'C:\Users\hsjun\OneDrive\Physics\Research\data';
    dataTypes = {'Etot.mat', 'Qtot.mat', 'Sent_imp.mat', 'Temps.mat', 'NRG_Op=ImpOrb.mat', 'NRG_Op=ImpSp.mat', 'ocont.mat'};
    AllParams = zeros(0,4);
    
    %% TI diagram, (J0, K0) = (0.1, 0.3)
    J0 = 0.1; K0 = 0.3;
    Params = ...        % (I0, Nkeep)
        [1e-14, 5000;
         1e-13, 5000;
         1e-12, 5000;
         1e-11, 5000;
         1e-10, 5000;
         1e-9, 5000;
         2e-9, 5000;
         3e-9, 5000;
         4e-9, 5000;
         5e-9, 5000;
         6e-9, 5000;
         7e-9, 5000;
         8e-9, 5000;
         9e-9, 5000;
         1e-8, 5000;
         1e-7, 5000;
         1e-6, 5000;
         1e-5, 5000;
         1e-4, 5000;
         1e-3, 5000;
         1e-2, 5000;
         1e-1, 5000;
         -1e-14, 5000;
         -1e-13, 5000;
         -1e-12, 5000;
         -1e-11, 5000;
         -1e-10, 5000;
         -1e-9, 5000;
         -2e-9, 5000;
         -3e-9, 5000;
         -4e-9, 5000;
         -5e-9, 5000;
         -6e-9, 5000;
         -7e-9, 5000;
         -8e-9, 5000;
         -9e-9, 5000;
         -1e-8, 5000;
         -1e-7, 5000;
         -1e-6, 5000;
         -1e-5, 5000;
         -1e-4, 5000;
         -1e-3, 5000;
         -1e-2, 5000;
         -1e-1, 5000;
         ];
    
    for it = 1:size(Params,1)
    
        pathSelected = [path,filesep,'TsoK_selected',filesep,'J0=',sprintf('%.15g',J0),'_K0=',sprintf('%.15g',K0), ...
                                '_I0=',sprintf('%.15g',Params(it,1)),'_T=1e-24_Nkeep=',sprintf('%.15g',Params(it,2))];
  
        pathTsoK = [path,filesep,'TsoK',filesep,'J0=',sprintf('%.15g',J0),'_K0=',sprintf('%.15g',K0), ...
                        '_I0=',sprintf('%.15g',Params(it,1)),'_T=1e-24_Nkeep=',sprintf('%.15g',Params(it,2))];
        if ~exist(pathSelected,'dir')
            mkdir(pathSelected);
        end
    
        TsoKInfo = dir(pathTsoK);
        SelectedInfo = dir(pathSelected);
        DirSelected = {SelectedInfo.name};
        for itx = 1:numel(TsoKInfo)
            
            DirTsoK = TsoKInfo(itx).name;
            if ismember(DirTsoK, dataTypes) && ~ismember(DirTsoK, DirSelected)
                copyfile([pathTsoK,filesep,DirTsoK], pathSelected);
                disp(['Copied ','J0=',sprintf('%.15g',J0),'_K0=',sprintf('%.15g',K0),'_I0=',sprintf('%.15g',Params(it,1)), ...
                        '_T=1e-24_Nkeep=',sprintf('%.15g',Params(it,2)),filesep,DirTsoK,' to ''TsoK_selected'' folder']);
            end
        end % itx
        

        AllParams = [AllParams; [J0, K0, Params(it,:)]];
    end


    %% TI diagram, (J0, K0) = (0.2, 0.3)
    J0 = 0.2; K0 = 0.3;
    Params = ...        % (I0, Nkeep)
        [1e-14, 5000;
         1e-13, 5000;
         1e-12, 5000;
         1e-11, 5000;
         1e-10, 5000;
         1e-9, 5000;
         1e-8, 5000;
         1e-7, 5000;
         1e-6, 5000;
         2e-6, 5000;
         3e-6, 5000;
         4e-6, 5000;
         5e-6, 5000;
         6e-6, 5000;
         7e-6, 5000;
         8e-6, 5000;
         9e-6, 5000;
         1e-5, 5000;
         1e-4, 5000;
         1e-3, 5000;
         1e-2, 5000;
         1e-1, 5000;
         -1e-14, 5000;
         -1e-13, 5000;
         -1e-12, 5000;
         -1e-11, 5000;
         -1e-10, 5000;
         -1e-9, 5000;
         -1e-8, 5000;
         -1e-7, 5000;
         -1e-6, 5000;
         -2e-6, 5000;
         -3e-6, 5000;
         -4e-6, 5000;
         -5e-6, 5000;
         -6e-6, 5000;
         -7e-6, 5000;
         -8e-6, 5000;
         -9e-6, 5000;
         -1e-5, 5000;
         -1e-4, 5000;
         -1e-3, 5000;
         -1e-2, 5000;
         -1e-1, 5000;
         ];
    
    for it = 1:size(Params,1)
    
        pathSelected = [path,filesep,'TsoK_selected',filesep,'J0=',sprintf('%.15g',J0),'_K0=',sprintf('%.15g',K0), ...
                                '_I0=',sprintf('%.15g',Params(it,1)),'_T=1e-24_Nkeep=',sprintf('%.15g',Params(it,2))];
  
        pathTsoK = [path,filesep,'TsoK',filesep,'J0=',sprintf('%.15g',J0),'_K0=',sprintf('%.15g',K0), ...
                        '_I0=',sprintf('%.15g',Params(it,1)),'_T=1e-24_Nkeep=',sprintf('%.15g',Params(it,2))];
        if ~exist(pathSelected,'dir')
            mkdir(pathSelected);
        end
    
        TsoKInfo = dir(pathTsoK);
        SelectedInfo = dir(pathSelected);
        DirSelected = {SelectedInfo.name};
        for itx = 1:numel(TsoKInfo)

            DirTsoK = TsoKInfo(itx).name;
            if ismember(DirTsoK, dataTypes) && ~ismember(DirTsoK, DirSelected)
                copyfile([pathTsoK,filesep,DirTsoK], pathSelected);
                disp(['Copied ','J0=',sprintf('%.15g',J0),'_K0=',sprintf('%.15g',K0),'_I0=',sprintf('%.15g',Params(it,1)), ...
                        '_T=1e-24_Nkeep=',sprintf('%.15g',Params(it,2)),filesep,DirTsoK,' to ''TsoK_selected'' folder']);
            end
            
        end % itx

        AllParams = [AllParams; [J0, K0, Params(it,:)]];
    end




    AllParams = table(AllParams(:,1), AllParams(:,2), AllParams(:,3), AllParams(:,4));
    AllParams.Properties.VariableNames = {'J0', 'K0', 'I0', 'Nkeep'};
end

