%% PA = LU Decomposition: Gaussian Elimination with Partial Pivoting
% Given matrix A, find permutation matrix P, lower-diagonal matrix L, and
% upper-diagonal matrix U, such that PA = LU.

% Input:      A - (n x n) Matrix to be decomposed into L & U factors
% Output:     P - (n x n) Permutation matrix
%             L - (n x n) Lower factor of A
%             U - (n x n) Upper factor of A

function[ P, L, U ] = PALU_decomp( A )

% Initialization of matrices
n = size(A,1);                      % Size of matrix A
P = eye(n);                         % Permutation matrix
L = eye(n);                         % Lower factor of A
U = A;                              % Upper factor of A

%% WORK ON PERMUTATION P   
for k = 1:(n - 1)                   % Loop over columns
% Choose pivot p such that | u_pk | = max| u_ik |
p = k;
M = U(p,k);                         % Maximal entry in column k
    
    % Check entries down column k for maximal entry
    for m = (k + 1):n
        if abs( U(m,k) ) > abs(M)
            M = U(m,k);             % Update maximal entry M
            p = m;                  % Update p
        else
            p = k;                  % If u_kk is already maximal entry
        end   
    end

% Permute p-th and k-th rows of U and P
temp = U(k,:);
U(k,:) = U(p,:);
U(p,:) = temp(1,:);

temp = P(k,:);
P(k,:) = P(p,:);
P(p,:) = temp(1,:);

% Permute p-th and k-th rows' entries of 1st (k - 1) columns of L
temp = L(k,:);
for q = 1:(k-1)
    L(k,q) = L(p,q);
    L(p,q) = temp(1,q);
end

%% WORK ON MATRICES L, U
    for i = (k + 1):n                            % Loop over rows
        % Compute l_ji
        L(i,k) = U(i,k) / U(k,k);
        U(i,k) = 0;                              % Clean up lower part of U
        
        for j= k+1:n                             % Loop over columns
        U(i,j) = U(i,j) - L(i,k) * U(k,j);       % Update U
        end
    end    
end

%% Error Analysis
fprintf( 'The infinity-norm of L is %e.\n', norm( L, inf ) );
fprintf( 'The infinity-norm of U is %e.\n', norm( U, inf ) );
fprintf( 'The infinity-norm of backward error ||E|| = ||LU - A|| = %e.\n', ...
    norm( ( L * U - A ), inf ) );
end