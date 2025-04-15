clear;
Datapath = 'C:\Users\hsjun\OneDrive\Physics\Research\data\Quartic';
NewPath = 'C:\Users\hsjun\OneDrive\Physics\Research\CQMPhys\Visualization_codes\PKData_Espectrum';
FileInfo = dir(Datapath);

cnt = 1;
lambda_x = [];
lambda_z = [];
for it = (1:numel(FileInfo))
    DirName = FileInfo(it).name;
    if DirName(1) == 'K'
        tmp = sscanf(DirName, 'K_z=%f_Q=%f_T=%f_Nkeep=%f');
        K_z = tmp(1);
        Q = tmp(2);
        T = tmp(3);
        Nkeep = tmp(4);
        if Nkeep == 3000
            LoadPath = [Datapath,filesep,DirName];
            NewFolder = [NewPath,filesep,'lambda_z=',sprintf('%.15g',K_z/2),'_lambda_x=',sprintf('%.15g',Q),'_T=',sprintf('%.15g',T),'_Nkeep=3000'];

            if ~exist(NewFolder)
                Etot = load([LoadPath,filesep,'Etot.mat']);
                Etot = Etot.Etot;
                Qtot = load([LoadPath,filesep,'Qtot.mat']);
                Qtot = Qtot.Qtot;
    
                mkdir(NewFolder);
                save([NewFolder,filesep,'Etot.mat'], 'Etot');
                save([NewFolder,filesep,'Qtot.mat'], 'Qtot');
    
                disp(['Moved ',sprintf('%d',cnt),' espectrum data']);
                cnt = cnt + 1;
            else
                lambda_z = [lambda_z; K_z/2];
                lambda_x = [lambda_x; Q];
            end
            
        end
    end
end % it

if cnt == 1
    disp('All data are up to date');
    Avail = table(lambda_z, lambda_x);
    disp(Avail);
end

