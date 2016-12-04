%% Simpson's Rule
% Numerical Integration Approximation

% Input:       f - Function of interest defined over interval [a,b]
%              a - Leftpoint of interval of interest
%              b - Endpoint of interval of interest
% Output:      value - Approximate value of the integral of f over [a,b]

function [ value ] = simpsons_rule( f, a, b )

% Compute the midpoint between points a & b
c = ( a + b ) / 2;

% Simpson's rule on [a,b]
h = ( b - a ) / 6;
value = h * ( feval(f,a) + 4 * feval(f,c) + feval(f,b) );

end