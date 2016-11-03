%% Test Linear Iterative Methods:
% Approximately solve Ax = b, where A = D + L + U, using:
% (1) Richardson
% (2) Jacobi
% (3) Gauss-Seidel
% (4) SOR

% Let n in { 2^k | k = 4,5,6,7,8 } be the sizes of our linear systems.
% Set the tolerance level: || b - A*x^k ||_2 / ||b||_2 < 10^(-6)

for k = 4:8
    
   % Initialize matrix A and vector b
   n = 2^k;                         
   A = zeros(n);                     
   b = zeros(n,1);
   
   % Create matrix A:
   % A is a tridiagonal matrix with 2's along the diagonal & 
   % (-1)'s on the upper and lower bands
   for i = 1:n-1
       A(i+1,i) = -1;               % Insert (-1)'s in the upper band
       A(i,i+1) = -1;               % Insert (-1)'s in the lower band
       A(i,i) = 2;                  % Insert 2's along the diagonal
   end
   A(n,n) = 2;

   % Create vector b = (1,0,...,0,-1)
   b(1) = 1;
   b(n) = -1;
    
   % Find a good initial guess x_0
   % x = A \ b;   
   % y = zeros(n,1);
   
   % for j = 1:n
       % y(j) = 0.00001;
       % y(j) = 0.000001;
   % end
   
   % x_0 = x + y;
   
   % Initialize other input values
   x_0 = zeros(n,1);
   tol = 10^(-6);
   iteration_lim = 100000;
   
   %[ sol, num_iterations ] = richardson( A, b, x_0, tol, iteration_lim );
   %[ sol, num_iterations ] = jacobi( A, b, x_0, tol, iteration_lim );
   %[ sol, num_iterations ] = gauss_seidel( A, b, x_0, tol, iteration_lim );
   [ sol, num_iterations ] = sor( A, b, x_0, tol, iteration_lim );
   
   disp(num_iterations);
     
end