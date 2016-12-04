%% Recursive Simpson's Rule
% Numerical Integration Approximation

% Input:       f - Function of interest defined over interval [a,b]
%              a - Leftpoint of interval of interest
%              b - Endpoint of interval of interest
%              epsilon - Tolerance level upon which to determine method
%              convergence
%              whole = simpsons_rule( f, a, b )
% Output:      value - Approximate value of the integral of f over [a,b]

function [ value ] = recursive_simpsons( f, a, b, epsilon, whole )

% Compute the midpoint between points a & b
c = ( a + b ) / 2;

% Compute for left and right intervals
left = simpsons_rule( f, a, c );
right = simpsons_rule( f, c, b );

% Estimate/check error tolerance
% Constant 15 computed by considering error on each interval
if abs(left + right - whole) <= ( 15 * (b - a) * epsilon ) 
    value = left + right;
else
value = recursive_simpsons( f, a, c, epsilon, left ) + ...
    recursive_simpsons( f, c, b, epsilon, right );
end

end