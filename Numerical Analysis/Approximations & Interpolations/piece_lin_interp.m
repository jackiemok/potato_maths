%% Piecewise Linear Polynomial Interpolation Example

% Consider the function f(x) = 1 / ( 1 + x^2 ) on [ -5, 5 ].
% Approximate f using piecewise linear polynomial interpolation 
% with equally distributed points.

% Input:     interpX - (1 x 1000) vector containing equidistant points
%            within the interval [ -5, 5 ]
% Output:    L - (1 x 1000) linear interpolator for our function f(x) with
%            approximate values to f at the corresponding points
%            error - Normed error between the values given by the 
%            linear interpolator & the true values given by function f(x)

% To run, let:
% interpX = linspace( -5, 5, 1000 );
% More generally, let interpX = linspace( -5, 5, [Any integer] ).

function [ L, error ] = piece_lin_interp( interpX )

% Using 11 equidistant points over [ -5, 5], get
% row-vector of x-values and row-vector of f(x)-values.
x = -5 : 5;                     
y = 1 ./ ( 1 + x.^2 );

% Initialize linear interpolator L
L = zeros( size( interpX ) ); 
 
% For each entry j in interpX, 
% check which interval I_i = [ x_{i-1}, x_i ] contains current entry.
% If entry is in I_i, compute the linear interpolation at given value.
for j = 1:length(interpX);
    for i = 2:length(x)

        % Check if current entry is in current interval:
        if interpX(j) >= x(i-1) && interpX(j) <= x(i)
            % Get approximate value to f at given point
            L1 = ( ( x(i) - interpX(j) ) / ( x(i) - x(i-1) ) ) * y(i-1);
            L2 = ( ( interpX(j) - x(i-1) ) / ( x(i) - x(i-1) ) ) * y(i);
            L(j) = L1 + L2; 
        end
    end
end

% Plot true function values along with approximation values.
hold on
plot( x, y, '*', 'Color', 'red' );                 % True values
plot( interpX, L, 'Color', 'magenta' );            % Approximate values

% Plot true values with 1000 equidistant points in [ -5, 5 ].
true_x = linspace( -5, 5, length(interpX) );
true_y = 1 ./ ( 1 + true_x.^2 );
plot( true_x, true_y, 'Color', 'blue' );

% Compute error: [ value( linear interpolator ) - value( y = f(x) ) ]
error = norm( L - true_y );

end
 