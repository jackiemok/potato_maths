%% Broyden I Method Example
% Solve F(x*) = 0 satisfying G(x*) = x*

% Example: Consider the following system.
% x1^2 + x2^2 = 1
% (x1 - 1)^2 + x2^2 = 1

% Input:  x_0 - (n x k) initial guess for the solution approximating x
%         num_iterations - Number of allowed iterations
% Output: sol - (n x k) matrix/vector approximation to true x solving Ax = b
%         num_iterations - Number of iterations required for convergence

function[ sol , num_iterations ] = broyden1( x_0, num_iterations )

% F(x) = [ f1(x1,x2) ; f2(x1,x2) ]
F = @(x1,x2) [ (x1^2 + x2^2 - 1); ((x1 - 1)^2 + x2^2 - 1) ];

sol = x_0;                   % Initialize approximate solution
B_k = eye(2);                % Initialize B0 ~= F'(x) (Jacobian)
 
% Run until convergence reached
for k = 1:num_iterations

    % Solve linear system
    s = B_k \ ( -1 * F(sol(1),sol(2)) ); 
    
    % Updating steps
    x_k = sol + s;
    s = x_k - sol;
    v = F(x_k(1),x_k(2)) - F(sol(1),sol(2));
    B_k = B_k + 1 / (transpose(s)*s) * ((v - B_k*s) * transpose(s));
  
    % Check convergence => Perform next iteration?
    if ( norm( x_k - sol )/norm(x_k ) ) < 10^(-6)
        sol = x_k;
        num_iterations = k;
        break;
    else 
        sol = x_k;
    end
end