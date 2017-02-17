function [G,H] = GradHess(a,k,x,y)

% a function to return the gradient and Hessian for nonlinear least
% squares optimization
% inputs:
    % x - x values
    % y - y values
    % a and k - parameters of exponential function
% outputs:
    % G - gradient vector with respect to a and k
    % H - Hessian matrix with respect to a and k

    % difference between actual and expected
    f = @(x) a*exp(k*x);
    diff = arrayfun(f,x) - y;
    
    % functions for the partial derivatives
    apartial = @(x) exp(k*x);
    kpartial = @(x) x*a*exp(k*x);
    da = arrayfun(apartial,x);
    dk = arrayfun(kpartial,x);
    
    % compute the output gradient G and Hessian H
    G = [sum(diff.*da); sum(diff.*dk)];
    H = [sum(da.^2) sum(da.*dk); sum(da.*dk) sum(dk.^2)];

        
end