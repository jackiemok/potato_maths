%% Gauss-Seidel Linear Iterative Method
% Approximately solve Ax = b, where A = [ ( D + L ) + U ] =: [ R + U ].

% Iterative step: 
% The k-th iteration is the k-th approximation of x solving Ax = b
% x_k = [ R_inv( b - U*x_{k-1} ) ]
% = [ R_inv*b - R_inv*U*x_{k-1} ] 
% =: [ C - T * x_{k-1} ]

% Input:  A - (n x n) (strictly or irreducibly) diagonally dominant matrix
%         b - (n x k) matrix/vector equivalent to A*x
%         x_0 - (n x k) initial guess for the solution approximating x
%         tol - Tolerance level upon which to determine method convergence
%         iteration_lim - Upper bound for number of iterations in method
% Output: sol - (n x k) matrix/vector approximation to true x solving Ax = b
%         num_iterations - Number of iterations required for convergence

function[ sol, num_iterations ] = gauss_seidel( A, b, x_0, tol, iteration_lim )

% Initialize matrices in regular splitting of A = D + R + U
U = zeros( size(A,1) );                          
R = A;                                 

% Create lower-triangular matrix R & upper-triangular matrix U
for j = 2:size(A,1)                    % Loop over columns
    for i = 1:j-1                      % Loop over rows
    R(i,j) = 0;                        % Insert zeros in upper-right matrix
    U(i,j) = A(i,j);                   % Insert values of A in upper-right matrix
    end
end

sol = x_0;                             % Initialize approximate solution

% Compute C = R_inv*b & T = R_inv*U for iterative step
C = R \ b;                           
T = R \ U;  

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