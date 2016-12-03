%% Newton's Method Example

% Example: Consider the following system.
% x1^2 + x2^2 = 1
% (x1 - 1)^2 + x2^2 = 1

% Newton's Method: 
% Solve F(x*) = 0 satisfying G(x*) = x*

% Iterative step: 
% x_{k+1} = x_k - [ J(x_k)_inv ] * F(x_k)

% Input:  x_0 - (n x k) initial guess for the solution approximating x
%         num_iterations - Number of allowed iterations
% Output: sol - (n x k) matrix/vector approximation to true x solving Ax = b
%         num_iterations - Number of iterations required for convergence

function[ sol, num_iterations ] = newton( x_0, num_iterations )

% F(x) = [ f1(x1,x2) ; f2(x1,x2) ]
F = @(x1,x2) [ (x1^2 + x2^2 - 1); ((x1 - 1)^2 + x2^2 - 1) ];

% Jacobian J = F'(x)
J = @(x1,x2) [ 2*x1 2*x2; (2*x1 - 2)  2*x2 ];

% Initialize approximate solution
sol = x_0; 

% Run until convergence reached
for k = 1:num_iterations

    % Solve linear system
    delta_xk = J(sol(1),sol(2)) \ ( -1 * F(sol(1),sol(2))  );          
    
    % Updating step
    x_k = sol + delta_xk;          
  
    % Check convergence => Perform next iteration?
    if ( norm(x_k - sol) / norm(x_k) ) < 10^(-6)
        sol = x_k;
        num_iterations = k;
        break;
    else 
        sol = x_k;
    end
    
end