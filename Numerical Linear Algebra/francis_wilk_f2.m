%% Francis Iteration on (n x n) Symmetric Tridiagonal Matrix of Degree One
% Exercise #5.6.3(f.2) from:
% 'Fundamentals of Matrix Computations, 3rd Ed.' by David S. Watkins

% Generate matrix t:
n = 6;
row = [ 2 1 zeros(1,n-2) ];
t = toeplitz(row);
a = t;

evals = eig(t);             % Eigenvalues of matrix t

% Iterate until | a(n,n-1) | < 10^(-16)
step = 1;
for n = 6:-1:2
while abs( a(n,n-1) ) >= 10^(-16)
    fprintf('Running Francis iteration with Wilkinson shift: Trial %d', ...
        step)
    francis_wilk;
    step = step + 1;
end
end

step = step - 1;