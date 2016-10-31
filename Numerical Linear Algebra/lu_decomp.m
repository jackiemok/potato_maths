%% LU Decomposition
% Gaussian Elimination without Pivoting

% Input:   n - The size of the linear system ( Random A in R^{nxn} )
% Output:  The random matrix A and its factors L and U

function [ A, L, U ] = lu_decomp( n )

% Initialize matrices
A = randn(n);  
L = eye(n);
U = A;

% To demonstrate the weakness of Gaussian Elimination without pivoting,
% guarantee at least one of the pivots is small.
A(1,1) = 50*eps*A(1,1);

% Compute LU Factorization
for k = 1:n-1
    for j = k+1:n
        
        % Check for zeros along the diagonal
        if U(k,k) == 0
            disp( 'A has no unpivoted LU factorization.' )
            return
        end
        
        L(j,k) = U(j,k) / U(k,k);
        
        % Row j <-- Row j - multiplier * (pivot row)
        U(j,k:n) = U(j,k:n) - L(j,k) * U(k,k:n);  
    end 
end 

%% Error Analysis
% Pivoting decreases the backward error and norms of L and U,
% but the performance of Gaussian Elimination without pivoting is
% usually not conspicuously bad and is more or less stable.
fprintf( 'The infinity norm of L is %e.\n', norm( L, inf ) );
fprintf( 'The infinity norm of U is %e.\n', norm( U, inf ) );
fprintf( 'The infinity norm of backward error ||E|| = ||LU - A|| = %e.\n', ...
    norm( ( L*U - A ), inf ) );

end