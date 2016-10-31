%% Tridiagonal Matrix LU Decomposition without Pivoting

% Input:    A - (n x n) tridiagonal matrix   
% Output:   L - Lower factor of A
%           U - Upper factor of A

function [ L, U ] = tridiag_LU_nopivot( A )

[m,n] = size(A);
if m ~= n 
    disp( 'Error: A is not square.' )
    return
end

% Initialize LU factors of A
L = eye(n);
U = A; 

for k = 2:n
    
    % Check for zeros along the diagonal
    if U(k,k) == 0
        disp( 'A has no unpivoted LU factorization' )
        return
    end
    
    % Multiplier in L
    L(k,k-1) = U(k,k-1) / U(k-1,k-1);   
    
    % Clean up lower band
    U(k,k-1) = 0;                          
    
    U(k,k) = U(k,k) - L(k,k-1) * U(k-1,min(k+1,n));
end

end