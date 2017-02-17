function [a,k,niter] = NonlinearLSO(a0,k0,x,y,maxiter,tol,S)

% a function to return the best-fit parameters a and k
% for the exponential function f(x) = a*exp(k*x)
% using the Levenberg-Marquardt method

% inputs:
    % a0 - initial guess for a
    % k0 - initial guess for k
    % x - vector for estimating f(x) = a*exp(k*x)
    % y - vector of actual values for the exponential function 
    % maxiter - max iterations
    % tol - a tolerance level
    % S - a scaling factor of lambda (between 2-10)
% outputs:
    % a - best-fit parameter a
    % k - best-fit parameter k
    % niter - the number of iterations

    lambda = 1000;
    niter = 0;

    % bundle the parameters into a column vector
    p_old = [a0; k0];
    
    % move in the first step
    [G,H] = GradHess(p_old(1),p_old(2),x,y);
    p = p_old - (H + lambda*diag(H))\G;
    
     while and(niter < maxiter, norm(p - p_old) > tol)
         [G,H] = GradHess(p(1),p(2),x,y);
         p_new = p - (H + lambda*diag(H))\G;
         % if updated solution better, accept
         if SSD(p_new(1),p_new(2),x,y) < SSD(p(1),p(2),x,y)
            p_old = p;
            p = p_new;
            lambda = lambda/S;
        else
            lambda = lambda*S;
        end
        niter = niter + 1;
     end
     
     % unbundle the parameters
     a = p(1);
     k = p(2);
end