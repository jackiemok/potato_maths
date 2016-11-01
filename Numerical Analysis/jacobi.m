%% Jacobi Linear Iterative Method
% Approximately solve Ax = b, where A = [ D + ( L + U ) ] =: [ D + R ]
% Iterative step: 
% The k-th iteration is the k-th approximation of x
% x_k = D_inv( b - R*x_{k-1} ) = [ D_inv*b - D_inv*R*x_{k-1} ]
% =: [ C - T*x_{k-1} ]

% Input:  A - (n x n) (strictly or irreducibly) diagonally dominant matrix
%         b - (n x k) matrix/vector equivalent to A*x
%         x_0 - (n x k) initial guess for the solution approximating x
%         tol - Tolerance level upon which to determine method convergence
%         iteration_lim - Upper bound for number of iterations in method
% Output: sol - Approximate solution to true x solving Ax = b
%         num_iterations - Number of iterations required for convergence

function[ sol, num_iterations ] = jacobi( A, b, x_0, tol, iteration_lim )                         

% Find regular splitting of A = D + R
n = size(A,1);
D = eye(n);                            % Initialize D = I
R = A;                                 % Initialize R = A
sol = x_0;                             % Initialize approximate solution

% Create diagonal matrix D
for i = 1:n
    D(i,i) = A(i,i);                   % Insert diagonal values of A
end

% Create matrix R = L + U
for i_1 = 1:n 
    R(i_1,i_1) = 0;                    % Insert zeros along diagonal of R
end

% Compute C = D_inv*b and T = D_inv*R for iterative step
C = D \ b;                           
T = D \ R;  

% Main loop: Run until residual < tolerance level
% Residual: r = [ b - A*x_(k-1) ]
for k = 1:iteration_lim
    
    % Obtain better and better approximate solutions to Ax = b
    % x_k = [ C - T*x_{k-1} ]
    x_k = ( C - T*sol );
    sol = x_k;
    
    % Check residual => Perform next iteration?
    if ( norm( b - A * sol ) / norm(b) ) < tol
        num_iterations = k;
        break;   
    end  
end

end