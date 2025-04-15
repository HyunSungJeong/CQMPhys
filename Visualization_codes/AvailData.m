function Avail = AvailData()
    % <Description>
    %
    % <Output>
    % Table showing the parameters for which energy spectrum data are available

    path = [fileparts(mfilename('fullpath')),filesep,'PKData_Espectrum'];
    FileInfo = dir(path);

    lambda_x = [];
    lambda_z = [];
    for it = (1:numel(FileInfo))
        DirName = FileInfo(it).name;
        if DirName(1) == 'l'
            tmp = sscanf(DirName, 'lambda_z=%f_lambda_x=%f_T=%f_Nkeep=%f');
            lambda_z = [lambda_z; tmp(1)];
            lambda_x = [lambda_x; tmp(2)];
        end
    end

    Avail = table(lambda_z, lambda_x);
end