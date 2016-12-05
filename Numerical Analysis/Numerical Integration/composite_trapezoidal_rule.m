%% Composite Trapezoidal Rule
% Numerical Integration Approximation

% Input:       f - Function of interest defined over interval [a,b]
%              a - Leftpoint of interval of interest
%              b - Endpoint of interval of interest
%              n - Number of points used to determine mesh size over [a,b]
% Output:      value - Approximate value of the integral of f over [a,b]

function [ value ] = composite_trapezoidal_rule( f, a, b, n )

% Mesh size
h = ( b - a ) / n;

% Get row-vector of x-values between [a,b]
% Distances between each pair of points given by mesh size h
x = a : h : b;

% Initialize approximate solution
value = 0;

% Composite trapezoidal rule
value = value + h * 0.5 * feval(f,a);
for i = 2:n
    value = value + h * feval(f,x(i));
end
value = value + h * 0.5 * feval(f,b);

end