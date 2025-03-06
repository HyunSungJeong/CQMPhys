function disp2 (varargin)
% < Description >
%
% disp2 (str1, [str2, ...])
%
% Display the input char array(s) on the screen if the environment variable
% $LOG_FILE is not set. Otherwise, print on the file whose path is
% specified by $LOG_FILE, without showing on the screen. This function uses
% fprintf for both ways.
%
% < Input >
% str1, str2, ... : [char array] Each line of message to be displayed or
%       printed to a file. After each input, a line change '\n' is
%       inserted.
%
% Written by S.Lee (Feb.06,2022)
% Updated by S.Lee (Feb.21,2023): to escape special characters '%' and '\'.

for it1 = (1:numel(varargin))
    if isempty(varargin{it1})
        varargin{it1} = '';
    end
end

if ~isempty(getenv('LOG_FILE')) % print to a file
    fileID = fopen(getenv('LOG_FILE'),'A');
    for it1 = (1:numel(varargin))
        if ischar(varargin{it1})
            fprintf(fileID,[varargin{it1},'\n']);
        else
            fprintf(fileID,[evalc('disp(varargin{it1})'),'\n']);
        end
    end
    fclose(fileID);
    
else % display on a screen
    for it1 = (1:numel(varargin))
        if ischar(varargin{it1})
            % "escape" special characters '%' and '\'
            varargin{it1} = strrep(varargin{it1},'%','%%');
            varargin{it1} = strrep(varargin{it1},'\','\\');
            fprintf([varargin{it1},'\n']);
        else
            disp(varargin{it1});
        end
    end
end

end