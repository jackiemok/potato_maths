%% QR-Factorization of Random (m x n) Matrix
% Create a random (m x n) matrix and compute its QR-factorization.

% A method to compute the QR-factorization of a given matrix A is to
% orthogonalize the column vectors of A (Gram-Schmidt orthogonalization).

% Input:                  m - The number of rows in our random matrix A
%                         n - The number of columns in our random matrix A
% Output:                 A - (m x n) Random matrix
%                         Q_c - Obtained Q using CGS, such that A = QR
%                         R_c - Obtained R using CGS, such that A = QR
%                         Q_m - Obtained Q using MGS, such that A = QR
%                         R_m - Obtained R using MGS, such that A = QR

function [ A, Q_c, R_c , Q_m, R_m, error_c, error_m ] = qrhat( m, n )

format long e;

% Create random (m x n) matrix A
A = randn(m,n);

% Initialize Q
Q_c = zeros(size(A));
R_c = zeros(size(A));

%% Classical Graham Schmidt (CGS) Method
% CGS is column-oriented and unstable, so solutions obtained by CGS are 
% sensitive to perturbation.

for j = 1:n
    
   vj = A(:,j);
   
   for i = 1:j-1
     R_c(i,j) = Q_c(:,i)' * A(:,j);
     vj = vj - R_c(i,j) * Q_c(:,i);
   end
   
   R_c(j,j) = norm(vj);
   Q_c(:,j) = vj / R_c(j,j);
end

%% Modified Graham Schmidt (MGS) Method
% MGS is row-oriented and numerically more stable than CGS, so it is less
% sensitive to roundoff errors.

V = A;
Q_m = zeros(size(A));
R_m = zeros(size(A));

for i = 1:n
   
   R_m(i,i) = norm( V(:,i) );  
   Q_m(:,i) = V(:,i) / R_m(i,i); 
   
   for j = i+1:n
    R_m(i,j) = Q_m(:,i)' * V(:,j);
    V(:,j)= V(:,j) - R_m(i,j) * Q_m(:,i);
   end 
end

%% Error Analysis
error_c = norm(A - Q_c * R_c);
error_m = norm(A - Q_m * R_m);

end