function [a1,Rsq1,a2,Rsq2,a3,Rsq3] = Insert(log_T,log_ImpDyn,fit_range)

    x = log_T;
    y = log_ImpDyn{1};
    x1 = x(x < fit_range(1,1));
    y1 = y(x < fit_range(1,1));
    y1 = y1(x1 > fit_range(1,2));
    x1 = x1(x1 > fit_range(1,2));
    a1 = polyfit(x1,y1,1);
    yfit1 = polyval(a1, x1);          % Estimated  Regression Line
    SStot = sum((y1-mean(y1)).^2);                    % Total Sum-Of-Squares
    SSres = sum((y1-yfit1).^2);                       % Residual Sum-Of-Squares
    Rsq1 = 1-SSres/SStot;                            % R^2


    x = log_T;
    y = log_ImpDyn{2};
    x2 = x(x < fit_range(2,1));
    y2 = y(x < fit_range(2,1));
    y2 = y2(x2 > fit_range(2,2));
    x2 = x2(x2 > fit_range(2,2));
    a2 = polyfit(x2,y2,1);
    yfit2 = polyval(a2, x2);          % Estimated  Regression Line
    SStot = sum((y2-mean(y2)).^2);                    % Total Sum-Of-Squares
    SSres = sum((y2-yfit2).^2);                       % Residual Sum-Of-Squares
    Rsq2 = 1- SSres/SStot;                            % R^2

    
    x = log_T;
    y = log_ImpDyn{3};
    x3 = x(x < fit_range(3,1));
    y3 = y(x < fit_range(3,1));
    y3 = y3(x3 > fit_range(3,2));
    x3 = x3(x3 > fit_range(3,2));
    a3 = polyfit(x3,y3,1);
    yfit3 = polyval(a3, x3);          % Estimated  Regression Line
    SStot = sum((y3-mean(y3)).^2);                    % Total Sum-Of-Squares
    SSres = sum((y3-yfit3).^2);                       % Residual Sum-Of-Squares
    Rsq3 = 1-SSres/SStot;                            % R^2
end