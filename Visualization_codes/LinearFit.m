function [Coeff, Rsq] = LinearFit(log_T,log_ImpDyn,fit_range)

    Coeff = cell(1,numel(log_ImpDyn));
    Rsq = zeros(1,numel(log_ImpDyn));

    for it = 1:numel(log_ImpDyn)

        x = log_T(log_T < fit_range{it}(1));
        y = log_ImpDyn{it}(log_T < fit_range{it}(1));
        y = y(x > fit_range{it}(2));
        x = x(x > fit_range{it}(2));
        Coeff{it} = polyfit(x,y,1);
        yfit = polyval(Coeff{it}, x);       % Estimated  Regression Line
        SStot = sum((y-mean(y)).^2);    % Total Sum-Of-Squares
        SSres = sum((y-yfit).^2);      % Residual Sum-Of-Squares
        Rsq(it) = 1-SSres/SStot;           % R^2
    end

end