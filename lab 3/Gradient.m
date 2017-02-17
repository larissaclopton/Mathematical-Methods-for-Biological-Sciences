function [gradient] = Gradient(x,b,A) 

% a function returning the gradient of a quadratic form

    gradient = [(-b(1) + A(1,1)*x(1) + A(1,2)*x(2)); ...
                (-b(2) + A(2,1)*x(1) + A(2,2)*x(2))];

end
