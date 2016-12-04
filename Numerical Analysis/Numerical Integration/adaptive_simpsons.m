%% Adaptive Simpson's approximation of integral of f
% Numerical Integration Approximation

% Input:       f - Function of interest defined over interval [a,b]
%              a - Leftpoint of interval of interest
%              b - Endpoint of interval of interest
%              epsilon - Tolerance level upon which to determine method
%              convergence
% Output:      value - Approximate value of the integral of f over [a,b]

% Example: 
% Given the tolerance epsilon = 10^(-8), report the approximate value of
% the integral over [1,pi] of function f(x) = x^2 * sin(x)

function [ value ] = adaptive_simpsons( f, a, b, epsilon )

% Apply Simpon's rule on the whole interval [a,b]
whole = simpsons_rule( f, a, b );

% Apply Adaptive Simpson's Rule on [a,b]
value = recursive_simpsons( f, a, b, epsilon, whole );

end