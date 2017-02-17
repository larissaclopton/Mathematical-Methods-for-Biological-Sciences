function [x,niter] = GradDescent(x0,b,A,maxiter,tol)

% a function for the gradient descent method
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
    r = 1; % affects size of step
    
    % move in the first step
    x_old = x0;
    x = x_old - r*Gradient(x_old,b,A);
    
    while and(niter < maxiter, norm(x - x_old) > tol)
        x_new = x - r*Gradient(x,b,A);
        % if updated solution better, accept
        if QuadFuncVal(x_new,b,A) < QuadFuncVal(x,b,A)
            x_old = x;
            x = x_new;
        % else, halve r and try again
        else
            r = r/2;
        end
        niter = niter+1;
    end
    
end
