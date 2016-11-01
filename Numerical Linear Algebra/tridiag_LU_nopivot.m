%% Tridiagonal Matrix LU Decomposition without Pivoting
% Compute the LU factorization of a tridiagonal matrix A via
% Gaussian Elimination without pivoting.

% Input:    A - (n x n) tridiagonal matrix   
% Output:   L - (n x n) lower factor of A
%           U - (n x n) upper factor of A

function [ L, U ] = tridiag_LU_nopivot( A )

[m,n] = size(A);
if m ~= n 
    disp( 'Error: A is not square.' )
    return
end

% Initialize L & U factors of A
L = eye(n);
U = A; 

% Gaussian Elimination without Pivoting
for k = 1:n-1
    
    % Check for zeros along the diagonal
    if U(k,k) == 0
        disp( 'A has no unpivoted LU factorization' )
        return
    end
    
    % Compute the multiplier of row (k+1) in L
    L(k+1,k) = U(k+1,k) / U(k,k);                        
    
    % Update row (k+1) in U
    U(k+1,k:min(k+1,n)) = U(k+1,k:min(k+1,n)) - L(k+1,k) * U(k,min(k+1,n));
    % Clean up the lower band in U
    U(k+1,k) = 0;  
end

end