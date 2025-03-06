function BiSearch_result = BiSearch(Bi_Par,T,J0,K0,I0,lo,hi,DegMul_lo,DegMul_hi,MaxSep,Thres,varargin)

    % <Description>
    %
    % Binary search with for phase of 2soK model at zero temperature with two of (J0,K0,I0) fixed
    %
    % <Input>
    % Bi_Par : [char] Variable to be binary-searched. One of J0, K0, I0.
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
    % 'NumComp',.. : [numeric] Number of lowest degenerate multiplets to be compared. (Default: 5).
    %
    % <Output>
    % Binary : [Logical] true if the phase is binary in interval [hi,lo].
    %
    % CASE I : Binary == true
    %
    %   Boundary : [numeric] The boundary value of Bi_Par.
    % 
    % CASE II : Binary == false
    % 
    %   lo : [numeric] A value of Bi_Par. The system is in same phase if Bi_Par is between this value and the input hi.
    %   hi : [numeric] A value of Bi_Par. The system is in same phase if Bi_Par is between the input lo and this value.
    %   mid : [numeric] Intermediate value of Bi_Par, which differs in phase from Bi_Par = hi or Bi_Par = lo.
    %   DegMul_mid : [numeric vector] Multiplet degeneracies at Bi_Par = mid.
    %   Etot : Energy flow data for Bi_Par = mid.
    %   Qtot : Energy flow quantum number data for Bi_Par = mid.
    %   lo < mid < hi is always true.

    hi_orig = hi;
    lo_orig = lo;

    NumComp = 5;   % Deflault value of NumComp
    if ~isempty(varargin)
        if varargin{1} == 'NumComp'
            NumComp = varargin{2};
        else
            disp2('ERR: invalid option');
            error('ERR: invalid option');
        end
    end

    if numel(DegMul_lo) ~= NumComp
        disp2(['Length of DegMul_lo (',sprintf('%.15g',numel(DegMul_lo)),') does not match NumComp (',sprintf('%.15g',NumComp),')']);
        error(['Length of DegMul_lo (',sprintf('%.15g',numel(DegMul_lo)),') does not match NumComp (',sprintf('%.15g',NumComp),')']);
    end

    if numel(DegMul_hi) ~= NumComp
        disp2(['Length of DegMul_hi (',sprintf('%.15g',numel(DegMul_hi)),') does not match NumComp (',sprintf('%.15g',NumComp),')']);
        error(['Length of DegMul_hi (',sprintf('%.15g',numel(DegMul_hi)),') does not match NumComp (',sprintf('%.15g',NumComp),')']);
    end

    % NRG parameters
    % ***************************************************
    %% Lambda and Nkeep must be syncronized with LineSearch
    % ***************************************************
    Lambda = 4;
    N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
    Nkeep = 5000;

    [ff, gg] = doZLD([-1;1],[1;1],Lambda,N,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));
    
    BiSearch_result{1} = true;

    if Bi_Par == 'J0'     

        disp2('');
        disptime('Binary Search Start');
        
        strtmp = cell(9,1);
        strtmp{1} = ['Searching Interval ',Bi_Par,' = [',sprintf('%.15g',lo),',',sprintf('%.15g',hi),']'];
        strtmp{2} = ['K0 = ',sprintf('%.15g',K0),', I0 = ',sprintf('%.15g',I0)];
        strtmp{3} = ['MaxSep = ',sprintf('%.15g',MaxSep),', Thres = ',sprintf('%.15g',Thres)];

        tmp = cellfun(@(x) [sprintf('%.15g',DegMul_lo(x)),' '],num2cell(1:numel(DegMul_lo)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{4} = ' ';
        strtmp{5} = ['Numbers of degenerate multiplets at ',Bi_Par,' = ',sprintf('%.15g',lo),' :'];
        strtmp{6} = ['[ ',tmp,']'];
        strtmp{7} = ' ';
        tmp = cellfun(@(x) [sprintf('%.15g',DegMul_hi(x)),' '],num2cell(1:numel(DegMul_hi)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{8} = ['Numbers of degenerate multiplets at ',Bi_Par,' = ',sprintf('%.15g',hi),' :'];
        strtmp{9} = ['[ ',tmp,']'];
        dispbox('-width',130,strtmp{:});    
        disp2('');    

        while abs(hi - lo) > Thres
            
            mid = (hi+lo)/2;
            [DegMul_mid, Etot, Qtot] = DegMul(mid,K0,I0,ff,Nkeep,Lambda,MaxSep);
            DegMul_mid = DegMul_mid(1:NumComp);
            
            if isequal(DegMul_hi,DegMul_mid)
                
                hi = mid;
                disp2('');
                disptime('hi = mid');
            elseif isequal(DegMul_lo,DegMul_mid)

                lo = mid;
                disp2('');
                disptime('lo = mid');
            else        % Not binary

                BiSearch_result{1} = false;
                BiSearch_result{2} = lo;
                BiSearch_result{3} = hi;
                BiSearch_result{4} = mid;
                BiSearch_result{5} = DegMul_mid;
                hi = lo;
            end

            if BiSearch_result{1}
                disp2(['Boundary is somewhere inside interval ',Bi_Par,' = [',sprintf('%.15g',lo),', ',sprintf('%.15g',hi),']']);
                disp2('');
            end

        end

    elseif Bi_Par == 'K0'

        disp2('');
        disptime('Binary Search Start');
        
        strtmp = cell(9,1);
        strtmp{1} = ['Searching Interval ',Bi_Par,' = [',sprintf('%.15g',lo),',',sprintf('%.15g',hi),']'];
        strtmp{2} = ['J0 = ',sprintf('%.15g',J0),', I0 = ',sprintf('%.15g',I0)];
        strtmp{3} = ['MaxSep = ',sprintf('%.15g',MaxSep),', Thres = ',sprintf('%.15g',Thres)];

        tmp = cellfun(@(x) [sprintf('%.15g',DegMul_lo(x)),' '],num2cell(1:numel(DegMul_lo)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{4} = ' ';
        strtmp{5} = ['Numbers of degenerate multiplets at ',Bi_Par,' = ',sprintf('%.15g',lo),' :'];
        strtmp{6} = ['[ ',tmp,']'];
        strtmp{7} = ' ';
        tmp = cellfun(@(x) [sprintf('%.15g',DegMul_hi(x)),' '],num2cell(1:numel(DegMul_hi)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{8} = ['Numbers of degenerate multiplets at ',Bi_Par,' = ',sprintf('%.15g',hi),' :'];
        strtmp{9} = ['[ ',tmp,']'];
        dispbox('-width',130,strtmp{:});    
        disp2('');

        while abs(hi - lo) > Thres  

            mid = (hi+lo)/2;
            [DegMul_mid, Etot, Qtot] = DegMul(J0,mid,I0,ff,Nkeep,Lambda,MaxSep);
            DegMul_mid = DegMul_mid(1:NumComp);
            
            if isequal(DegMul_hi,DegMul_mid)

                disp2('');
                disptime('hi = mid');
                hi = mid;
            elseif isequal(DegMul_lo,DegMul_mid)

                disp2('');
                disptime('lo = mid');
                lo = mid;
            else        % Not Binary

                BiSearch_result{1} = false;
                BiSearch_result{2} = lo;
                BiSearch_result{3} = hi;
                BiSearch_result{4} = mid;
                BiSearch_result{5} = DegMul_mid;
                hi = lo;
            end

            if BiSearch_result{1}
                disp2(['Boundary is somewhere inside interval ',Bi_Par,' = [',sprintf('%.15g',lo),', ',sprintf('%.15g',hi),']']);
                disp2('');
            end

        end

    elseif Bi_Par == 'I0'

        disp2('');
        disptime('Binary Search Start');
        
        strtmp = cell(9,1);
        strtmp{1} = ['Searching Interval ',Bi_Par,' = [',sprintf('%.15g',lo),',',sprintf('%.15g',hi),']'];
        strtmp{2} = ['J0 = ',sprintf('%.15g',J0),', K0 = ',sprintf('%.15g',K0)];
        strtmp{3} = ['MaxSep = ',sprintf('%.15g',MaxSep),', Thres = ',sprintf('%.15g',Thres)];

        tmp = cellfun(@(x) [sprintf('%.15g',DegMul_lo(x)),' '],num2cell(1:numel(DegMul_lo)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{4} = ' ';
        strtmp{5} = ['Numbers of degenerate multiplets at ',Bi_Par,' = ',sprintf('%.15g',lo),' :'];
        strtmp{6} = ['[ ',tmp,']'];
        strtmp{7} = ' ';
        tmp = cellfun(@(x) [sprintf('%.15g',DegMul_hi(x)),' '],num2cell(1:numel(DegMul_hi)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{8} = ['Numbers of degenerate multiplets at ',Bi_Par,' = ',sprintf('%.15g',hi),' :'];
        strtmp{9} = ['[ ',tmp,']'];

        disp2('');
        dispbox('-width',130,strtmp{:});    
        disp2('');

        while abs(hi - lo) > Thres 

            mid = (hi+lo)/2;
            [DegMul_mid, Etot, Qtot] = DegMul(J0,K0,mid,ff,Nkeep,Lambda,MaxSep);
            DegMul_mid = DegMul_mid(1:NumComp);

            if isequal(DegMul_hi,DegMul_mid)

                disp2('');
                disptime('hi = mid');
                hi = mid;
            elseif isequal(DegMul_lo,DegMul_mid)

                disp2('');
                disptime('lo = mid');
                lo = mid;
            else        % Not Binary

                BiSearch_result{1} = false;
                BiSearch_result{2} = lo;
                BiSearch_result{3} = hi;
                BiSearch_result{4} = mid;
                BiSearch_result{5} = DegMul_mid;
                hi = lo;
            end

            if BiSearch_result{1}
                disp2(['Boundary is somewhere inside interval ',Bi_Par,' = [',sprintf('%.15g',lo),', ',sprintf('%.15g',hi),']']);
                disp2('');
            end

        end

    else 
        disp2('ERR: Bi_Par must be one of J0, K0, I0');
        error('ERR: Bi_Par must be one of J0, K0, I0');
    end

    if BiSearch_result{1}
        BiSearch_result{2} = (hi+lo)/2;   % Boundary
    else
        BiSearch_result{6} = Etot;
        BiSearch_result{7} = Qtot;

        disp2('');
        disptime(['Search Fork at ',Bi_Par,' = ',sprintf('%.15g',mid)]);
        strtmp = cell(7,1);
        strtmp{1} = ['J0 = [',sprintf('%.15g',lo_orig),',',sprintf('%.15g',hi_orig),'] ---> [', ...
                        sprintf('%.15g',BiSearch_result{2}),',',sprintf('%.15g',BiSearch_result{3}),'] ===> [', ...
                            sprintf('%.15g',BiSearch_result{2}),',',sprintf('%.15g',mid),'] / [', ...
                                sprintf('%.15g',mid),',',sprintf('%.15g',BiSearch_result{3}),']'];

        tmp = cellfun(@(x) [sprintf('%.15g',DegMul_mid(x)),' '],num2cell(1:numel(DegMul_mid)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{2} = ' ';
        strtmp{3} = ['Numbers of degenerate multiplets at ',Bi_Par,' = ',sprintf('%.15g',mid),' :'];
        strtmp{4} = ['[ ',tmp,']'];
        strtmp{5} = ' ';

        if Bi_Par == 'J0'
            strtmp{6} = ['K0 = ',sprintf('%.15g',K0),', I0 = ',sprintf('%.15g',I0)];
        elseif Bi_Par == 'K0'
            strtmp{6} = ['J0 = ',sprintf('%.15g',J0),', I0 = ',sprintf('%.15g',I0)];
        elseif Bi_Par == 'I0'
            strtmp{6} = ['J0 = ',sprintf('%.15g',J0),', K0 = ',sprintf('%.15g',K0)];
        end

        strtmp{7} = ['MaxSep = ',sprintf('%.15g',MaxSep),', Thres = ',sprintf('%.15g',Thres)];
        dispbox('-width',130,strtmp{:});
        disp2('');
    end
end