%% Richardson Linear Iterative Method
% Approximately solve Ax = b, where A = [ D + ( L + U ) ] =: [ D + R ]

% Iterative step: 
% The k-th iteration is the k-th approximation of x solving Ax = b
% x_k = x_{k-1} + [ wI( b - A*x_{k-1} ) ]
% =: x_{k-1} + W( b - A*x_{k-1} )
% =: x_{k-1} + Wb - WA*x_{k-1}

% Input:  A - (n x n) (strictly or irreducibly) diagonally dominant matrix
%         b - (n x k) matrix/vector equivalent to A*x
%         x_0 - (n x k) initial guess for the solution approximating x
%         tol - Tolerance level upon which to determine method convergence
%         iteration_lim - Upper bound for number of iterations in method
% Output: sol - (n x k) matrix/vector approximation to true x solving Ax = b
%         num_iterations - Number of iterations required for convergence

function[ sol, num_iterations ] = richardson( A, b, x_0, tol, iteration_lim )

% Initialize approximate solution
sol = x_0;                          

% Compute optimal value of omega: [ 2 / ( lambda_min(A) + lambda_max(A) ) ]
w = 2 / ( min(eig(A)) + max(eig(A)) );

% Test performance of less optimal value of omega
% w = 2 / max(eig(A));

% Initialize W, Wb, WA
W = w * eye( size(A,1) );
Wb = W * b;
WA = W * A;

% Main loop: Run until residual < tolerance level
% Residual: r = b - A * x_(k-1)
for k = 1:iteration_lim
    
    % Obtain better and better approximate solutions to Ax = b
    % x_k = x_{k-1} + Wb - WA*x_{k-1}
    sol = ( sol + Wb - WA * sol );
    
    % Check residual => Perform next iteration?
    if ( norm( b - A * sol ) / norm(b) ) < tol
        num_iterations = k;
        break;   
    end  
end

end