%% LU Solve without Pivoting
% Perform A = LU decomposition and solve Ax = (LU)x = b
% ( A = LU ) => ( Ax = b <=> (LU)x = b )

% Input:    A - (n x n) matrix
%           b - n-length vector equivalent to A*xsol
%           xsol - [Optional input] The exact solution used to obtain b
% Output:   L - (n x n) lower factor of A
%           U - (n x n) upper factor of A
%           x - The computed solution via forward & backward substitutions

function [ L, U, x ] = LU_solve_nopivot( A, b, xsol )

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
    for j = k+1:n
        
        % Check for zeros along the diagonal
        if U(k,k) == 0
            disp( 'Error: A has no unpivoted LU factorization.' )
            return
        end
        
        % Compute multiplier for row j
        L(j,k) = U(j,k) / U(k,k);
        
        % Row j <-- Row j - multiplier * (pivot row)
        U(j,k:n) = U(j,k:n) - L(j,k) * U(k,k:n);  
    end 
end 

%% Solve the system

% Forward substitution for intermediate solution: Ly = b
% ( A = LU ) => ( Ax = b <=> (LU)x = b )
% <=> L(Ux) = b, where y = (Ux)
y = L \ b;  
% Backward substitution for final solution: Ux = y
x = U \ y; 

%% Error Analysis

% Compute the residual
r = b - A*x;

% Print relative residual norm
rep = sprintf( 'The relative residual norm is: %0.5g', ...
    norm(r) / norm(b) );
disp(rep)

% Print relative accuracy of solution
if nargin == 3
    rep = sprintf( 'The relative accuracy is: %0.5g', ...
        norm(xsol - x) / norm(xsol) );
    disp(rep)
end

end