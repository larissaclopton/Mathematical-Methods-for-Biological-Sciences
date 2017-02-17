function [x,niter] = LevenbergMarquardt(x0,b,A,maxiter,tol,S)

% a function for the Levenberg-Marquardt method
% inputs:
    % x0 - initial value
    % b - a vector
    % A - a matrix
    % maxiter - max iterations
    % tol - a tolerance level
    % S - a scaling factor of lambda (between 2-10)
% outputs:
    % x - best solution vector
    % niter - the number of iterations

    lambda = 1000;
    niter = 0;

    % move in the first step
    x_old = x0;
    x = x_old - (A + lambda*diag(A))\Gradient(x_old,b,A);
    
    while and(niter < maxiter, norm(x - x_old) > tol)
        x_new = x - (A + lambda*diag(A))\Gradient(x,b,A);
        % if updated solution better, accept
        if QuadFuncVal(x_new,b,A) < QuadFuncVal(x,b,A)
            x_old = x;
            x = x_new;
            lambda = lambda/S;
        else
            lambda = lambda*S;
        end
        niter = niter + 1;
    end
end