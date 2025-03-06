function [Boundaries, Num_DegMuls] = N_Ary_Search(Bi_Par,T,J0,K0,I0,lo,hi,MaxSep,Thres,varargin)
    
    % <Description>
    %
    % N-ary search with for phase of 2soK model at zero temperature with two of (J0,K0,I0) fixed
    %
    % <Input>
    % Bi_Par : [char] Variable to be searched. One of J0, K0, I0.
    % T : [numeric] Temperature.
    % J0 : [numeric] 2soK model parameter. if Bi_Par == 'J0', the value of this input is not used (just leave it empty: J0 = []).
    % K0 : [numeric] 2soK model parameter. if Bi_Par == 'K0', the value of this input is not used (just leave it empty: K0 = []).
    % I0 : [numeric] 2soK model parameter. if Bi_Par == 'I0', the value of this input is not used (just leave it empty: I0 = []).
    % hi : [numeric] Upper bound of Bi_Par.
    % lo : [numeric] Lower bound of Bi_Par.
    % MaxSep : [numeric] Maximum separation of energy levels that are regarded as degenerate.
    % Thres : [numeric] Minimum difference of hi and lo in binary search.
    %
    % <Option>
    % 'NumComp',.. : [numeric] Number of lowest degenerate multiplets to be compared. (Default: 20).
    % 'NRGdata',... : [char] Path to temporarily store NRG data.
    %
    % <Output>
    % Boundaries : [numeric vector] Phase boundaries of Bi_Par
    % Num_DegMuls : [cell array of numeric vectors] Multiplet degeneracies for each phase divided by Boundaries.
    %           Each cell elements are numeric vectors consisting of multiplet degeneracies for each phase divied by Boundaries.
    %           In order of increasing value of Bi_Par.
    
        NumComp = 20;   % Deflault value of NumComp
    while ~isempty(varargin)
        switch varargin{1}
            case 'NumComp'
                NumComp = varargin{2};
                varargin(1:2) = [];
            case 'NRGdata'
                NRGdata = varargin{2};
                varargin(1:2) = [];
                
            otherwise
                if ischar(varargin{10})
                    error(['ERR: Unknown input ',varargin{10}]);
                else
                    disp2(varargin{10});
                    error('ERR: Unknown input');
                end
        end
    end

    output = BiSearch(Bi_Par,T,J0,K0,I0,hi,lo,MaxSep,Thres,'NumComp',NumComp);
    Binary = output{1};

    if ~Binary
        lo_new = output{2};
        hi_new = output{3};
        mid = output{4};
        [BND_lower, NDM_lower] = N_Ary_Search(Bi_Par,T,J0,K0,I0,lo_new,mid,MaxSep,Thres,'NumComp',NumComp);

        disp2('');
        disptime('');
        disp2('################################################################');
        strtmp = cell(3,1);
        strtmp{1} = ['Binary Search for interval ',Bi_Par,' = [',sprintf('%.15g',lo_new),', ',sprintf('%.15g',mid),']',' is complete!'];
        strtmp{2} = 'Boundaries are at :';
        tmp = cellfun(@(x) [sprintf('%.15g',BND_lower(x)),' '],num2cell(1:numel(BND_lower)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{3} = [Bi_Par,' = [ ',tmp,']'];
        dispbox('-width',130,strtmp{:});
        disp2('################################################################');
        disp2('');

        [BND_upper, NDM_upper] = N_Ary_Search(Bi_Par,T,J0,K0,I0,mid,hi_new,MaxSep,Thres,'NumComp',NumComp);

        disp2('');
        disptime('');
        disp2('################################################################');
        strtmp = cell(3,1);
        strtmp{1} = ['Binary Search for interval ',Bi_Par,' = [',sprintf('%.15g',mid),', ',sprintf('%.15g',hi_new),']',' is complete!'];
        strtmp{2} = 'Boundaries are at :';
        tmp = cellfun(@(x) [sprintf('%.15g',BND_upper(x)),' '],num2cell(1:numel(BND_upper)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{3} = [Bi_Par,' = [ ',tmp,']'];
        dispbox('-width',130,strtmp{:});
        disp2('################################################################');
        disp2('');

        Boundaries = cat(2,BND_lower,BND_upper);
        Num_DegMuls = cat(2,NDM_lower(1:end),NDM_upper(2:end));
    else
        Boundary = output{2};
        Num_DegMuls = cell(1,2);

        Boundaries = Boundary;
        Num_DegMuls{1} = output{3};
        Num_DegMuls{2} = output{4};
    end

end