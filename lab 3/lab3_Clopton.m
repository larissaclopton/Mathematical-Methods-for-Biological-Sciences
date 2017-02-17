%% Lab 3: Nonlinear optimization and data fitting
% BIOS 26211 winter 2017
% Larissa Clopton

%% Part 1: Linear regression with log transform
% You will find a data file called KaiC_data.mat which contains the data
% of concentration of a circadian protein KaiC (which was transformed to decay 
% to 0 instead of approaching a positive asymptote). Read in the file and
% you will see two variables: Time (in hours) and KaiC (in arbitrary
% units). Plot the time series and you will see that it looks like
% exponential decay. However, the log-transform of the dependent variable,
% has a linear relationship with time:
%
%    KaiC=a*exp(k*Time)
%
%    log(KaiC) = log(a)+k*Time
% 
% 1.1 Use the linear regression function from lab 1 and the log-transformed
% data for KaiC to find the best-fit parameters A and k. Plot the data and
% overlay the fit using the exponential function with the parameters you
% found. Report the sum of squared errors of the fit.

clear all; close all;
load('KaiC_data.mat');

% plot the time series
figure; hold on;
plot(Time,KaiC,'k');
xlabel('Time (hours)'); ylabel('KaiC concentration');
title('Time vs KaiC concentration: Linear fit');

% perform linear regression on the log transform to obtain 
% best-fit parameters a and k of the exponential function
log_KaiC = log(KaiC);
[k,log_a,corr] = LinReg(Time,log_KaiC);
a = exp(log_a);

% overlay the model exponential function
model_function = a*exp(k*Time);
plot(Time,model_function,'b');
legend('data','linear regression');

% calculate sum of squared errors
SSE = sum((model_function-KaiC).^2);
display(SSE);

% The sum of squared errors is 782.6108. As shown by the time series,
% linear fitting techniques do not do a great job of representing the data
% as nonlinear factors are at work.

%% Part 2: Multidimensional optimization 
% The following are exercises in implementing multidimensional optimization
% algorithms in the case of a quadratic function with known linear and
% quardatic terms, as follows:
%
%   y = - b*x + 1/2*x'*A*x
%
% (where x and b are column vectors and A is a matrix)
%
% Use the following values for testing your code:
b = [5 -10]; % made into a row vector
A1 = [2 1; 1 2];
A2 = [1 -2; -2 20];

% 2.1 Write a function to return the value of the quadratic function y for
% a given vector x, vector b, and matrix A. The output should be a scalar
% number (should be one line). Then write another function to return the
% gradient of the quadratic function y for a given x, b, and A (should
% be one line). It should return a vector of the two partial derivatives of
% y w.r.t. the two components x1 and x2. 

disp('testing QuadFuncVal and Gradient');

x = [2; 4]; % Note: this is arbitrary

% return the value of the quadratic function for given x, b, and A
val1 = QuadFuncVal(x,b,A1);
val2 = QuadFuncVal(x,b,A2);
display(val1); display(val2);

% return the gradient of the quadratic function for given x, b, and A
% taking partial derivatives of x(1) and x(2)
g1 = Gradient(x,b,A1);
g2 = Gradient(x,b,A2);
display(g1); display(g2);

% 2.2 Implement the *gradient descent* method in a function. The tricky
% step is determining how far to move down the gradient - it can be done by
% implementing a golden section search, or any other one-dimensional
% optimization algorithm. I recommend an easy way out - multiply the
% gradient by some arbitrary number r (e.g. 1); if the updated solution
% is better, accept it, and if not, that means you overshot, so divide r by
% some factor (e.g. 2) and try again.
%
% * Input: initial value x0, vector b, matrix A, max number of steps,
% tolerance.
% * DO: Implement the gradient descent method; call the gradient function
% from 2.1 to obtain the gradient at each value and the quadratic function
% for the value of y.
% * Output: best solution vector x, number of steps 
%
% Test the performance of your function using both Hessian
% matrices A1 and A2. Try different starting values, use the tolerance of
% 1e-5, and report the number of steps it took to converge.

disp('Gradient descent method');

tol = 1*10^(-5); % tolerance
maxiter = 500; % max number of iterations

% test performance on A1 and A2 for multiple starting values
x0 = [-10 10 50; 20 -15 10];
for i = 1:3
    [x1,niter_A1] = GradDescent(x0(:,i),b,A1,maxiter,tol);
    [x2,niter_A2] = GradDescent(x0(:,i),b,A2,maxiter,tol);
    display(x1); display(niter_A1);
    display(x2); display(niter_A2);
end

% For a number of starting values, the gradient descent method converges
% quicker for the A1 Hessian matrix than the A2 Hessian matrix. In
% addition, the gradient descent with the A1 Hessian matrix converges to 
% the exact same value while the gradient descent with the A2 Hessian matrix 
% converges to approximately the same value.

% 2.3 Implement the *Newton-Raphson* method 
%
% * Input: initial value x0, vector b, matrix A, max number of steps,
% tolerance.
% * DO: Implement the Newton-Raphson method; call the gradient function
% from 2.1 to obtain the gradient at each value and the quadratic function
% for the value of y.
% * Output: best solution vector x, number of steps 
%
% Test the function using both Hessian matrices A1 and A2. Try different
% starting values and use a tolerance of 1e-5, and report the number of
% steps it took to converge. 

disp('Newton-Raphson method');

% test performance on A1 and A2 for multiple starting values
for i = 1:3
    [x1,niter_A1] = NewtonRaphson(x0(:,i),b,A1,maxiter,tol);
    [x2,niter_A2] = NewtonRaphson(x0(:,i),b,A2,maxiter,tol);
    display(x1); display(niter_A1);
    display(x2); display(niter_A2);
end
    
% For different starting values and both Hessian matrixes A1 and A2, the
% Newton-Raphson method converges in just one iteration to the exact same
% value, performing quite efficiently.

% 2.4 Implement the *Levenberg-Marquardt* method 
%
% * Input: initial value x0, vector b, matrix A, max number of steps,
% tolerance, scaling factor x.
% * DO: Implement the Levenberg-Marquardt method; call the gradient function
% from 2.1 to obtain the gradient at each value and the quadratic function
% for the value of y.
% * Output: best solution vector x, number of steps 
%
% Test the function using both Hessian matrices A1 and A2. Try different
% starting values and use a tolerance of 1e-5, and report the number of
% steps it took to converge. 

disp('Levenberg-Marquardt method');

% test performance on A1 and A2 for multiple starting values
S = 5; % a scaling factor
for i = 1:3
    [x1,niter_A1] = LevenbergMarquardt(x0(:,i),b,A1,maxiter,tol,S);
    [x2,niter_A2] = LevenbergMarquardt(x0(:,i),b,A2,maxiter,tol,S);
    display(x1); display(niter_A1);
    display(x2); display(niter_A2);
end

% For different starting values and both Hessian matrixes A1 and A2, the 
% Levenberg-Marquardt method appears take around the same number of
% iterations for both Hessians, with the number of iterations for A2 just
% slightly higher than that of A1. Across all starting values for a given
% Hessian, the Levenberg-Marquardt method converges reliably to the same
% minimum value.

%% Part 3: Nonlinear least squares optimization
%
% In this part you will use the Levenberg-Marquardt method to perform a
% nonlinear least squares fit, using the exponential decay function:
% 
%   f(x,a,k) = a*exp(k*x)
%
% 3.1 Write a function to compute the sum of square deviations (the objective
% function).
%
% * Input: initial values a0 and k0, vector x, vector y.
% * DO: add up all the square differences between predicted and actual
% values of y (should be one line).
% * Output: the sum of square differences - a single scalar. 

% Write a function to read in a data set of x and y values and
% calculate the gradient (G) and Hessian (H) for the nonlinear least squares
% optimization. 
% * Input: values a and k, vector x, vector y.
% * DO: The gradient  component in the direction of parameter a is the sum 
% over all data points x(i), y(i) of the following term:
%   G(1) = sum[(f(x(i)-y(i))*d(f(x(i))/da]
% The gradient component in the direction of parameter k is the sum 
% over all data points x(i), y(i) of the following term:
%   G(2) = sum[(f(x(i)-y(i))*d(f(x(i))/dk]
%
% The Hessian matrix is then computed as follows:
%   H(1,1) = sum[(d(f(x(i))/da)^2]
%   H(2,2) = sum[(d(f(x(i))/dk)^2]
%   H(1,2) = H(2,1)= sum[(d(f(x(i))/da)*(d(f(x(i))/dk)]
%
% You'll have to do some work on paper to calculate the partial
% derivatives of the function with respect to a and k.
%
% * Output: vector G (1 by 2) and matrix H (2 by 2).

% 3.2  Modify the function from 2.4 to implement the Levenberg-Marquardt
% method for finding best-fit parameters a and k for the data set KaiC.
%
% * Input: initial values a0, k0, vector x, vector y, max number of steps,
% tolerance, scaling factor x.
% * DO: Perform L-M optimization, use the functions you wrote in 3.1 to 
% * Output: best fit parameters a and k, and the number of steps it took to
% converge.
%
% Try different starting values of a and k, and report how
% quickly the optimization converges. Compare the fit obtained using
% log-transform and linear regression in Part 1 and the fit from the L-M
% optimization by inputing both sets of values a nd k into the sum of
% squared errors function. Plot the data (as circles) and overlay plots 
% of two exponential functions: one with parameters from 1.1, and the other
% from Levenberg-Marquardt optimization. Report the sum of squared
% differences for the nonlinear fit and compare it with the log-transform
% fit.

clear all; close all;
load('KaiC_data.mat');

% obtain the linear regression parameters
log_KaiC = log(KaiC);
[linreg_k,log_a,corr] = LinReg(Time,log_KaiC);
linreg_a = exp(log_a);

lin_function = linreg_a*exp(linreg_k*Time);
SSE_lin = sum((lin_function-KaiC).^2);
display(SSE_lin);

% define matrix of starting parameter values
p0 = [25 20 30; -0.6 -0.7 -0.5]; % first row a0, second row k0
tol = 1*10^(-5); % tolerance
S = 5; % a scaling factor

for i = 1:3
    [a,k,niter] = NonlinearLSO(p0(1,i),p0(2,i),Time,KaiC,500,tol,S);
    display(a); display(k); display(niter);
    
    % compare linear and nonlinear fits by their SSE
    nonlin_function = a*exp(k*Time);
    SSE_nonlin = sum((nonlin_function-KaiC).^2); 
    display(SSE_nonlin);
 
    % plot data and overlay the two fitted functions
    figure; hold on;
    scatter(Time,KaiC,'k','filled')
    plot(Time,lin_function,'b');
    plot(Time,nonlin_function,'g');
    legend('Data','Linear fit', 'Nonlinear fit');
    xlabel('Time (hours)'); ylabel('KaiC concentration');
    title('Time versus KaiC concentration: Linear and Nonlinear fits');
    
end

% We see that the output parameters from the nonlinear least squares
% optimization fits the data better (green line) than the output parameters
% from the linear regression (blue line). This is reflected in a
% much lower sum of squared errors for the nonlinear versus linear
% fit.

% Total points: 30
% Point breakdown:
% 1.1 4 pts (2 pts for code, 1 for plot, 1 for answer)
% 2.1 3 pts (1 pt for objective function, 2 for gradient function)
% 2.2 4 pts (2 for function, 2 for minimization results of A1 and A2)
% 2.3 4 pts (2 for function, 2 for minimization results of A1 and A2)
% 2.4 4 pts (2 for function, 2 for minimization results of A1 and A2)
% 3.1 5 pts (2 for objective function, 3 for gradient and Hessian)
% 3.2 6 pts (3 pts for function, 2 for results and experimentation, 1 for plot)
