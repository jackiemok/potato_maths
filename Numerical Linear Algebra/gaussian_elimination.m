%% Naive Gaussian Elimination
% Solve the linear system Ax = b.

% Input         A - (n x n) matrix
%               b - (n x 1) vector equivalent to A*x
% Output        x - (n x 1) vector such that A*x = b

function [ x ] = gaussian_elimination( A, b )

% Initialization of solution x
n = size(A,1);   
x = zeros(n,1);

%% Gaussian Elimination without Pivoting

for k = 1:(n - 1)                 % Loop over columns
    for i = (k + 1):n             % Loop over rows
        multiplier = A(i,k) / A(k,k);                        
        A(i,k+1:n) = A(i,k+1:n) - multiplier * A(k,k+1:n);       
        b(i) = b(i) - multiplier * b(k);                      
    end 
end

%% Back Substitution

for i = n:-1:1    
    b(i) = b(i) - A(i,i+1:n) * x(i+1:n);                    
    x(i) = b(i) / A(i,i);                       
end