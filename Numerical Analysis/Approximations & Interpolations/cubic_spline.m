%% Clamped Cubic Spline Approximation Example

% Consider the function f(x) = 1 / ( 1 + x^2 ) on [ -5, 5 ].
% Approximate f using the clamped cubic spline method with
% equally distributed points of mesh size h.

% Input:     interpX - (1 x 1000) vector containing equidistant points
%            within the interval [ -5, 5 ]
% Output:    p - (1 x 1000) linear interpolator for our function f(x) with
%            approximate values to f at the corresponding points
%            error - Normed error between the values given by the 
%            linear interpolator & the true values given by function f(x)

% To run, let:
% interpX = linspace( -5, 5, 1000 );
% More generally, let interpX = linspace( -5, 5, [Any integer] ).

function [ p, error ] = cubic_spline( interpX )

% Use 11 equidistant points within the interval of interest.
startpt = -5; endpt = 5;  
x = linspace( startpt, endpt, 11 );
num_pts = length(x);

% Function f(x) = 1 / ( 1 + x^2 )
% Get row-vector of funtion values for desired x-values.
y = 1 ./ ( 1 + x.^2 );

% Initialize and define mesh sizes h:
% As h decreases in size, we expect a better approximation to f.
h = zeros( size(x) );
for i = 2:num_pts
 h(i-1) = abs( x(i) - x(i-1) );
end
h(num_pts) = abs( x(num_pts) - x(num_pts-1) );      % Handle last entry

% Set up linear system:
% Initialize and define (n x n) matrix A
A = zeros( num_pts, num_pts );
% Rows 2 to (n-1) of A
for i = 2:(num_pts - 1)
 A(i,i-1) = h(i-1);                                 % Subdiagonal terms
 A(i,i) = 2 * ( h(i-1) + h(i) );                    % Diagonal terms
 A(i,i+1) = h(i);                                   % Upper-diagonal terms
end
% Rows 1 and n of A
A(1,1) = 2 * h(1);
A(1,2) = h(1);
A(num_pts,num_pts) = 2 * h(num_pts);
A(num_pts,num_pts-1) = h(num_pts);

% Initialize and define (n x 1) vector b
b = zeros( num_pts, 1 );
% Rows 2 to (n-1) of b
for i = 2:(num_pts - 1)
b(i) = ( 1/h(i) * y(i-1) ) - ( 1/h(i) + 1/h(i+1) ) * y(i) + ( 1/h(i+1) * y(i+1) );
end
% Rows 1 and n of b:
% f'(x) = ( - 2x ) / ( x^2 + 1 )^2
b(1) = (2 * startpt) / (1 + startpt^2)^2 - ( 1/h(1) * y(1)) + ( 1/h(1) * y(2));
b(num_pts) = - (2 * endpt) / (1 + endpt^2)^2 + ( 1/h(num_pts) * y(num_pts-1) ) - ( 1/h(num_pts) * y(num_pts) );
b = 6 * b; 

% Solve the linear system for vector alpha
alpha = A \ b;

% Initialize and compute p
p = zeros( size(interpX) ); 
constA = zeros(num_pts,1);
constB = zeros(num_pts,1);

% For each entry j in interpX, 
% check which interval I_i = [ x_{i-1}, x_i ] contains current entry.
% If entry is in I_i, compute the approximation to f at given value.
for j = 1:length(interpX);
 for i = 2:length(x)
     
     % Check if current entry is in current interval:
     if interpX(j) >= x(i-1) && interpX(j) <= x(i)
         constA(i) = y(i-1) - alpha(i-1) * (h(i)^2) / 6;
         constB(i) = y(i) - alpha(i) * (h(i)^2) / 6;
         term1 = alpha(i-1) * (x(i) - interpX(j))^3 / (6 * h(i));
         term2 = alpha(i) * (interpX(j) - x(i-1))^3 / (6 * h(i));
         term3 = constA(i) * (x(i) - interpX(j)) / h(i);
         term4 = constB(i) * (interpX(j) - x(i-1)) / h(i);
         
         % Get approximate value to f at given point
         p(j) = term1 + term2 + term3 + term4; 
     end
 end
end

% Plot true function values along with approximation values.
hold on
plot( x, y, '*', 'Color', 'red' );                 % True values
plot( interpX, p, 'Color', 'magenta' );            % Approximate values

% Plot true values with 1000 equidistant points in [ -5, 5 ].
true_x = linspace( startpt, endpt, length(interpX) );
true_y = 1 ./ ( 1 + true_x.^2 );
plot( true_x, true_y, 'Color', 'blue' );

% Compute error: [ value( linear interpolator ) - value( y = f(x) ) ]
error = norm( p - true_y );

end
