function [val] = SSD(a,k,x,y)

% a function to return the sum of squared deviations

    val = sum((a*exp(k*x)-y).^2);

end