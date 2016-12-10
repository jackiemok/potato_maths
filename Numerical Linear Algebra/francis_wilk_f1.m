%% Francis Iteration on (n x n) Symmetric Tridiagonal Matrix of Degree One
% Exercise #5.6.3(f.1) from:
% 'Fundamentals of Matrix Computations, 3rd Ed.' by David S. Watkins

% Generate matrix t:
% n = 6;
% row = [ 2 1 zeros(1,n-2) ];
% t = toeplitz(row);
% a = t;

%% (i) Compute the Rayleigh quotient shift

shift = 0.753020396282533;

%% (ii) Once the shift has been chosen, we need a rotator transformation
% Create bulge

cs = a(1,1) - shift;
sn = a(2,1);
r = norm( [cs sn] );
cs = cs / r;
sn = sn / r;

% Normalizing to make cs^2 + sn^2 = 1
q0 = [ cs -sn; sn cs ];                    % Givens rotator
a(1:2,:) = q0' * a(1:2,:);                 % Left multiplication
a(:,1:2) = a(:,1:2) * q0;                  % Right multiplication

%% (iii) Bulge chasing
for k = 1:(n - 2)
    
    % Chase the bulge from position (k + 2, k)
    cs = a(k+1,k);
    sn = a(k+2,k);
    r = norm( [cs sn] );
    cs = cs / r;
    sn = sn / r;
    
    a(k+1,k) = r;
    a(k+2,k) = 0;
    qi = [ cs -sn; sn cs ];
    
    % Givens rotator to chase the bulge
    a(k+1:k+2,k+1:n) = qi' * a(k+1:k+2,k+1:n);
    a(:,k+1:k+2) = a(:,k+1:k+2) * qi;
end

%% (d) Print some matrix entries
% Our objective is to make the subdiagonal entries go to zero and the
% diagonal entries converge to the eigenvalues.

format short e;
subdiag = diag(a,-1);        % Subdiagonal entries of A exponential format
format long;
bottom_entry = a(n,n);       % Diagonal entries of A in long format
disp(shift);
disp(a);