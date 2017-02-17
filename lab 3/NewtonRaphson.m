function [x,niter] = NewtonRaphson(x0,b,A,maxiter,tol)

% a function for the Newton-Raphson method
% inputs:
    % x0 - initial value
    % b - a vector
    % A - a matrix
    % maxiter - max iterations
    % tol - a tolerance level
% outputs:
    % x - best solution vector
    % niter - the number of iterations
    
    niter = 0;

    % move in the first step
    x_old = x0;
    x = x_old - A\Gradient(x_old,b,A);
    
    while and(niter < maxiter, norm(x - x_old) > tol)
        x_new = x - A\Gradient(x,b,A);
        x_old = x;
        x = x_new;
        niter = niter + 1;
    end
    
end