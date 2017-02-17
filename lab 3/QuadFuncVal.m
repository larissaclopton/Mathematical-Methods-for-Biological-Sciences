function [val] = QuadFuncVal(x,b,A)

% a function returning the value of a quadratic form

    val = -b*x + (1/2)*x'*A*x;
end