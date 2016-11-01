%% LU Decomposition without Pivoting
% Compute the L & U factors of a random matrix A of size (n x n)
% using Gaussian Elimination without pivoting.

% Input:   n - The size of the linear system ( Random A in R^{n x n} )
% Output:  A - (n x n) random matrix equivalent to the product of L*U
%          L - (n x n) lower factor of A
%          U - (n x n) upper factor of A

function [ A, L, U ] = LU_decomp_nopivot( n )

% Initialize matrices
A = randn(n);  
L = eye(n);
U = A;

% To demonstrate the weakness of Gaussian Elimination without pivoting,
% guarantee at least one of the pivots is small.
% A(1,1) = 50 * eps * A(1,1);

% Gaussian Elimination without Pivoting
for k = 1:n-1
    for j = k+1:n
        
        % Check for zeros along the diagonal
        if U(k,k) == 0
            disp( 'A has no unpivoted LU factorization.' )
            return
        end
        
        % Compute multiplier of row j
        L(j,k) = U(j,k) / U(k,k);
        
        % Row j <-- Row j - multiplier * (pivot row)
        U(j,k:n) = U(j,k:n) - L(j,k) * U(k,k:n);  
    end 
end 

%% Error Analysis
% Pivoting decreases the backward error and norms of L and U,
% but the performance of Gaussian Elimination without pivoting is
% usually not conspicuously bad and is more or less stable.

% Print infinity norms of L, U, and the backward error E = ( LU - A )
fprintf( 'The infinity norm of L is %e.\n', norm( L, inf ) );
fprintf( 'The infinity norm of U is %e.\n', norm( U, inf ) );
fprintf( 'The infinity norm of backward error ||E|| = ||LU - A|| = %e.\n', ...
    norm( ( L*U - A ), inf ) );

end