function [slope, intercept, correlation] = LinReg(X,Y)

    % check that X and Y are the same length
    if length(X) ~= length(Y)
        error('LinReg: Vectors are not same length.');
    end
    
    Xbar = mean(X);
    Ybar = mean(Y);
    
    % calculate the covariance and variance of the two variables
    CovArray = (X - Xbar).*(Y - Ybar);
    Cov = (1/(length(X) - 1))*sum(CovArray);
    
    XVarArray = (Xbar - X).^2;
    XVar = (1/(length(X) - 1))*sum(XVarArray);
    YVarArray = (Ybar - Y).^2;
    YVar = (1/(length(Y) - 1))*sum(YVarArray);
    
    % compute the slope and intercept of the best-fit line
    slope = Cov/XVar;
    intercept = Ybar - slope*Xbar;
    
    % compute the correlation coefficient r
    correlation = Cov/((XVar*YVar)^(1/2));

end
