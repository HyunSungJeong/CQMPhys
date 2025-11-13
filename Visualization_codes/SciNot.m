function sci = SciNot(Num,varargin)
    % <Description>
    % Generates a string representing the input number in scientific notation
    %
    % <Input>
    % Num : [numeric] input number
    %
    % <Option>
    % 'Signif', ... : [numeric] number of significant figures to be represented in the output
    %                       (Default: 1)
    %
    % <Output>
    % sci : [char] input represented in scientific notation

    %% Parse input
    if ~isnumeric(Num)
        error('ERR: ''Num'' must be a real number');
    end

    %% Parse options

    NumSig = 1;

    while ~isempty(varargin)
        switch varargin{1}
            case 'Signif'
                if ~isnumeric(varargin{2})
                    error('ERR: ''Signif'' must be a natural number');
                else
                    NumSig = varargin{2};
                    varargin(1:2) = [];
                end

            otherwise
                if ischar(varargin{1})
                    error(['ERR: unknown input ',varargin{1}]);
                else
                    error('ERR: unknown input');
                end
        end % switch - case
    end

    if Num == 0
        sci = '0';

    else
        Expon = floor(log10(abs(Num)));
        Coeff = round(Num/power(10,Expon), NumSig, 'significant');
    
        if Coeff == 1
            sci = ['10^{',sprintf('%.15g',Expon),'}'];
        elseif Coeff == -1
            sci = ['-10^{',sprintf('%.15g',Expon),'}'];
        else
            sci = [sprintf('%.15g',Coeff),' \times 10^{',sprintf('%.15g',Expon),'}'];
        end
        
    end
end