%% Francis Iteration on (n x n) Symmetric Tridiagonal Matrix of Degree One
% Exercise #5.6.3(a)-(d) from:
% 'Fundamentals of Matrix Computations, 3rd Ed.' by David S. Watkins

% Generate matrix t:
% n = 6;
% row = [ 2 1 zeros(1,n-2) ];
% t = toeplitz(row);
% a = t;

%% (a) Compute the Wilkinson shift 

htr = .5 * ( a(n-1,n-1) + a(n,n) );                              % Half a (2 x 2) trace
dscr = sqrt( (.5 * (a(n-1,n-1) - a(n,n)))^2 + a(n,n-1)^2 );      % Discriminant

% Avoid cancellation
if htr < 0
    dscr = -dscr;
end

root1 = htr + dscr;                                              % Quadratic formula
if root1 == 0                                                    % Case almost never happens
    root2 = 0;
else                                                             % Case almost always happens
    % (2x2) determinant = product of roots
    det = a(n-1,n-1) * a(n,n) - a(n,n-1)^2;
    root2 = det / root1;
end

if abs( a(n,n) - root1 ) < abs( a(n,n) - root2 )
    shift = root1;
else
    shift = root2;
end

%% (b) Create bulge
% Once the shift has been chosen, we need a rotator transformation

cs = a(1,1) - shift;
sn = a(2,1);
r = norm( [cs sn] );
cs = cs / r;
sn = sn / r;

% Normalizing to make cs^2 + sn^2 = 1
q0 = [ cs -sn; sn cs ];                                 % Givens rotator
a(1:2,:) = q0' * a(1:2,:);                              % Left multiplication
a(:,1:2) = a(:,1:2) * q0;                               % Right multiplication

%% (c) Bulge chasing
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
