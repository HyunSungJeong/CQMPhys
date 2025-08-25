function TsoK_Aniso_selected_data()
    
    % Parameters
    baseOriginal = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK_Aniso';
    baseTarget = 'C:\Users\hsjun\OneDrive\Physics\Research\data\TsoK_Aniso_selected';
    suffix = '_Nkeep=6000_Lambda=2.5';
    
    % Get all subfolders in baseOriginal
    allFolders = dir(baseOriginal);
    allFolders = allFolders([allFolders.isdir]);  % keep only directories
    allFolders = allFolders(~ismember({allFolders.name}, {'.', '..'}));  % exclude . and ..
    
    % Filter folders whose names end with the specified suffix
    isMatch = cellfun(@(name) endsWith(name, suffix), {allFolders.name});
    matchedFolders = allFolders(isMatch);
    
    for k = 1:length(matchedFolders)
        origFolderName = matchedFolders(k).name;
        origFolderPath = fullfile(baseOriginal, origFolderName);
        
        % Define corresponding target folder
        targetFolderPath = fullfile(baseTarget, origFolderName);
        
        % Create target folder if it doesn't exist
        if ~exist(targetFolderPath, 'dir')
            mkdir(targetFolderPath);
        end
    
        % List files in original folder (non-recursive)
        origFiles = dir(fullfile(origFolderPath, '*'));
        origFiles = origFiles(~[origFiles.isdir]);  % exclude subfolders
    
        if ~ismember('Temps.mat', {origFiles.name})
            disp(origFolderPath);
        end

        if isempty({origFiles.name})
            %disp(['empty : ', origFolderPath]);
        end

        % Get existing files in target folder
        if exist(targetFolderPath, 'dir')
            targetFiles = dir(fullfile(targetFolderPath, '*'));
            targetFileNames = {targetFiles(~[targetFiles.isdir]).name};
        else
            targetFileNames = {};
        end
    
        % Copy only missing files
        for i = 1:length(origFiles)
            if ~ismember(origFiles(i).name, targetFileNames)
                srcFile = fullfile(origFolderPath, origFiles(i).name);
                destFile = fullfile(targetFolderPath, origFiles(i).name);
                copyfile(srcFile, destFile);
                fprintf('Copied %s to %s\n', srcFile, destFile);
            end
        end
    end
    
    disp('Synchronization complete.');


end

