function varargout = dispbox (varargin)
% < Description >
%
% [res =] dispbox (['-width', ..,] input1 [, input2, input3, ...]  )
%
% Convert input strings to boxed character array. Line breaks between
% inputs.
%
% < Input >
% input* : [cell array] Each input* corresponds to a line. That is, there
%       is a line change between different input*.
%       If input* is
%       - char array: show as it is.
%       - numeric: converted to string with sprintf/fprintf format
%           '%.4g, '. The last remaining ', ' is truncated.
%
% < Option >
% '-width', .. : (integer) Width of box. If > 0, the width is fixed by the
%       value. If = 0, the width is determined by the maximum length of
%       input lines. If < 0, the width is lower bounded by the absolute
%       value and can be larger.
%       NOTE : This option needs to be the first input.
%       (Default: 88)
%
% < Output >
% res : [char array] Char array that describes the content of input and box
%       boundary. If the output is suppressed (i.e. not explicitly 
%       requested), the char array is displayed on screen. If the ouput is
%       requested, it is not displayed on screen.
%
% Written by S.Lee (2015)
% Updated by S.Lee (Jan.15,2019): Changed Unicode characters for the box
%       boundary.
% Updated by S.Lee (Jan.16,2019): Rewrote for clarity and performance.
% Updated by S.Lee (Mar.18,2019): Now it checks locale. If the locale is
%       UTF-8, then it uses Unicode box-drawing characters to draw the box
%       boundary. Otherwise, it uses '+', '-', and '|' to draw the box
%       boundary.
% Updated by S.Lee (Feb.05,2022): Use 'disp2' that can print either on
%       screen or on a log file.
% Updated by S.Lee (Feb.21,2023): escaping special characters '%' and '\'
%       is now done by 'disp2', not within this routine.

% try

% check locale
isUTF = false; % default
loc_info = feature('locale');
if contains(loc_info.encoding,'UTF') || contains(loc_info.encoding,'utf')
    isUTF = true;
end

maxlen = 88; % default

if ischar(varargin{1}) && strcmp(varargin{1},'-width')
    maxlen = varargin{2};
    varargin(1:2) = [];
end

res = cell(numel(varargin),1);

wmax = 0; % temporary variable that measures the actual width

for it1 = (1:numel(varargin))
    if ischar(varargin{it1})
        restmp = varargin{it1};
    elseif isnumeric(varargin{it1})
        restmp = sprintf('%.4g, ',varargin{it1});
        restmp(end-1:end) = [];
    else
        error(['ERR: Data type (',class(varargin{it1}),') of input ', ...
            sprintf('%i',it1),' is not supported.']);
    end
    
    res2 = cell(size(restmp,1),1);
    
    for it2 = (1:size(restmp,1))
        res_part = restmp(it2,:); % row vector
        
        if maxlen > 0
            res_part_tmp = cell(max(ceil(numel(res_part)/maxlen),1),1);
            for it3 = (1:(numel(res_part_tmp)-1))
                res_part_tmp{it3} = res_part(((it3-1)*maxlen+1):(it3*maxlen));
            end
            res_part_tmp{end} = res_part(((numel(res_part_tmp)-1)*maxlen+1):end);
            if numel(res_part_tmp{end}) < maxlen
                res_part_tmp{end} = [res_part_tmp{end},repmat(' ',[1 maxlen-numel(res_part_tmp{end})])];
            end
            res2{it2} = cell2mat(res_part_tmp);
        elseif (maxlen < 0) && (numel(res_part) < abs(maxlen))
            res2{it2} = [res_part,repmat(' ',[1 abs(maxlen)-numel(res_part)])];            
        else
            res2{it2} = res_part;
        end
        
        if wmax < size(res2{it2},2)
            wmax = size(res2{it2},2);
        end
    end
    
    res{it1} = res2;
end

for it1 = (1:numel(res))
    for it2 = (1:numel(res{it1}))
        if size(res{it1}{it2},2) < wmax
            res{it1}{it2} = [res{it1}{it2}, repmat(' ',[size(res{it1}{it2},1) wmax-size(res{it1}{it2},2)])];
        end
    end
    res{it1} = cell2mat(res{it1});
end

res = cell2mat(res);

if isUTF % UTF : use box-drawing unicode characters to describe the box boundary
    hltmp = repmat(char(9472),[1 size(res,2)]); % horizontal line
    res = [repmat([char(9474),' '],[size(res,1) 1]),res,repmat([' ',char(9474)],[size(res,1) 1])]; % left and right borders
    res = [[char([9484 9472]),hltmp,char([9472 9488])];res;[char([9492 9472]),hltmp,char([9472 9496])]]; % top and bottom borders
else % use '+', '-', and '|' to describe the box boundary
    hltmp = repmat('-',[1 size(res,2)]); % horizontal line
    res = [repmat('| ',[size(res,1) 1]),res,repmat(' |',[size(res,1) 1])]; % left and right borders
    res = [['+-',hltmp,'-+'];res;['+-',hltmp,'-+']]; % top and bottom borders
end

if nargout > 0
    varargout{1} = res;
else
    res2 = num2cell(res,2); 
    disp2(res2{:});
end

% catch e
%     disp2(getReport(e));
%     keyboard;
%     rethrow(e);
% end

end