%% Demo: Condition Number & Residual Norm with Hilbert (n x n) Matrix 
% Exercise #2.6.6 from:
% 'Fundamentals of Matrix Computations, 3rd Ed.' by David S. Watkins

%% Consider the cases when n = 4,8,12,16
for i = 1:4
    
n = 4 * i;
fprintf( 'Considering the n = %d case: \n', n );

z = ones(n,1);               % z = 1 in R^n
H = hilb(n);                 % H = Hilbert(n x n) matrix

% Solve the system: H*x = b
b = H * z;
xhat = H \ b;
disp(xhat);

% If we solve the system H*x = b, in theory we should obtain x = z = ones(n,1)
% But what do we get in practice?
% As the size of our linear system n increases, notice how the
% condition number and error changes.
fprintf( 'The 2-norm condition number of H is %e.\n', cond(H) );
fprintf( 'The 2-norm difference ||xhat - z|| = %e.\n', norm(xhat - z) );
fprintf( 'The 2-norm residual ||b - H*xhat|| = %e.\n', norm(b - H*xhat) );

end

%% Exercise #2.6.7 from:
% 'Fundamentals of Matrix Computations, 3rd Ed.' by David S. Watkins

%% PART A
format long e
n = 12;
z = ones(n,1);               % z = 1 in R^n
H = hilb(n);                 % H = Hilbert(n x n) matrix

% Solve the system: H*y = b
y = H \ z;
disp(y);
fprintf( 'The 2-norm of y is %e.\n', norm(y) );

%% PART B
b = H * y;
disp(b);
fprintf( 'The 2-norm difference ||b - z|| = %e.\n', norm(b - z) );

%% PART C
x = H \ b;
disp(x);
fprintf( 'The 2-norm relative error ||x-y||/||y|| = %e.\n', norm(x - y) / norm(y) );

%% PART D
fprintf( 'The 2-norm residual ||b - H*x|| = %e.\n', norm(b - H * x) );
fprintf( 'The 2-norm condition number of H is %e.\n', cond(H) );
fprintf( 'The upper-bound for relative error is %e.', ( (cond(H) * norm(b - H*x)) / norm(b) ) );
