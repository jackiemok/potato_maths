%% Jacobi Linear Iterative Method
% Approximately solve Ax = b, where A = [ D + ( L + U ) ] =: [ D + R ]

% Iterative step: 
% The k-th iteration is the k-th approximation of x solving Ax = b
% x_k = D_inv( b - R * x_{k-1} ) = [ (D_inv * b) - (D_inv * R) * x_{k-1} ]
% =: [ C - T * x_{k-1} ]

% Input:  A - (n x n) (strictly or irreducibly) diagonally dominant matrix
%         b - (n x k) matrix/vector equivalent to A*x
%         x_0 - (n x k) initial guess for the solution approximating x
%         tol - Tolerance level upon which to determine method convergence
%         iteration_lim - Upper bound for number of iterations in method
% Output: sol - (n x k) matrix/vector approximation to true x solving Ax = b
%         num_iterations - Number of iterations required for convergence

function[ sol, num_iterations ] = jacobi( A, b, x_0, tol, iteration_lim )                         

% Initialize matrices in regular splitting of A = D + R
D = eye(n);                            
R = A;                                 

% Create diagonal matrix D and matrix R = L + U
for i = 1:size(A,1)
    D(i,i) = A(i,i);                   % Insert diagonal values of A
    R(i,i) = 0;                        % Insert zeros along diagonal of R
end

sol = x_0;                             % Initialize approximate solution

% Compute C = D_inv*b and T = D_inv*R for iterative step
C = D \ b;                           
T = D \ R;  

% Main loop: Run until residual < tolerance level
% Residual: r = b - A * x_(k-1)
for k = 1:iteration_lim
    
    % Obtain better and better approximate solutions to Ax = b
    % x_k = C - T * x_{k-1} 
    sol = ( C - T * sol );
    
    % Check residual => Perform next iteration?
    if ( norm( b - A * sol ) / norm(b) ) < tol
        num_iterations = k;
        break;   
    end  
end

end