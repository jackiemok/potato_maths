%% Successive Over-Relaxation Method (SOR)
% Approximately solve Ax = b, where A = D + L + U  

% Iterative step: 
% The k-th iteration is the k-th approximation of x solving Ax = b
% x_k = x_{k-1} + inv( w^(-1)*D + L ) * [ b - A*x_{k-1} ]
% =: x_{k-1} + B( b - A*x_{k-1} ) 
% =: x_{k-1} + Bb - BA*x_{k-1} 

% Input:  A - (n x n) (strictly or irreducibly) diagonally dominant matrix
%         b - (n x k) matrix/vector equivalent to A*x
%         x_0 - (n x k) initial guess for the solution approximating x
%         tol - Tolerance level upon which to determine method convergence
%         iteration_lim - Upper bound for number of iterations in method
% Output: sol - (n x k) matrix/vector approximation to true x solving Ax = b
%         num_iterations - Number of iterations required for convergence

function[ sol, num_iterations ] = sor( A, b, x_0, tol, iteration_lim )

% Initialize approximate solution
sol = x_0;  

% Initialize matrices in regular splitting of A = D + L + U
n = size(A,1);
D = eye(n);                            
L = A;                                 

% Create diagonal matrix D and lower-triangular matrix L
for i = 1:n;
    D(i,i) = A(i,i);                   % Insert diagonal values of A
    L(i,i) = 0;                        % Insert zeros along diagonal of L
end

% Clean up upper-triangular portion of L
for j = 2:n                            % Loop over columns
    for i = 1:j-1                      % Loop over rows
    L(i,j) = 0;                        % Insert zeros in upper-right matrix
    end
end

% Compute optimal omega: w = 2 / [ 1 + sqrt( 1 - rho( I - inv(D)*A )^2 ) ]
rho = max( eig( eye(n) - D \ A ) );
w = 2 / ( 1 + sqrt( 1 - rho^2 ) );

% Test performance of less optimal value of omega
% w = 1.5;

% Initialize B, Bb, BA
B_inv = ( w^(-1) * D + L );
Bb = B_inv \ b;
BA = B_inv \ A;

% Main loop: Run until residual < tolerance level
% Residual: r = b - A * x_(k-1)
for k = 1:iteration_lim
    
    % Obtain better and better approximate solutions to Ax = b
    % x_k = x_{k-1} + B*b - B*A*x_{k-1}
    sol = ( sol + Bb - BA * sol );
    
    % Check residual => Perform next iteration?
    if ( norm( b - A * sol ) / norm(b) ) < tol
        num_iterations = k;
        break;   
    end  
end

end