%% QR-Factorization & Solve Demo
% Exercise #3.3.10 from:
% 'Fundamentals of Matrix Computations, 3rd Ed.' by David S. Watkins

%% PART (a)

% Create a random (n x m) matrix A
n = 6; m = 3;
A = randn(n,m);
fprintf( 'Our random matrix A is \n' );
disp(A);

% Compute QR-factorization of A
[Q,R] = qr(A);
fprintf( 'We get A = QR, where Q = \n' )
disp(Q);
fprintf( 'and R = \n' );
disp(R);

fprintf( 'Transpose(Q)*Q = \n' );
disp(Q' * Q);

fprintf( 'Norm of (I - Transpose(Q)*Q) = \n');
disp( norm(eye(n) - Q'*Q) );
fprintf( 'Norm of (A - Q*R) = \n' );
disp( norm(A - Q*R) );

%% PART (b)
% A = QR => b = Ax = (QR)x = Q(Rx)
% [1] Solve Qc = b, where Rx = c, by multiplying: c = (Q^T)b
% [2] Solve Rx = c (for x) by computing: x = R\c

%% (i) Least Squares Problem for Exercise 3.3.9(a)
% Quadratic polynomial: f(t) = 1 + t + t^2

% Time A = [ t_1 ... t_7 ]
A = [ -1; -0.75; -0.5; 0; 0.25; 0.5; 0.75 ];

% Data Y = [ y_1 ... y_7 ]
Y = [ 1; 0.8125; 0.75; 1.00; 1.3125; 1.75; 2.3125 ];

B = zeros(7,3);
% Build B corresponding to basis vectors
for i = 1: 7
   B(i,1) = 1;
   B(i,2) = A(i,1);
   B(i,3) = A(i,1)^2;
end

% BX = Y
% QRX = Y
% Qc = Y
% c = Q'Y
% X = R\(Q'Y)
[Q,R] = qr(B);
X = R \ (Q.' * Y);

%% (ii) Least Squares Problem for Exercise 3.3.9(b)
% Linear polynomial with the basis vectors: 
% f1(t) = 50; f2(t) = t - 1065

% Time A = [ t_1 ... t_6 ]
A = [ 1000; 1050; 1060; 1080; 1110; 1130 ];

% Data Yi = [ y_1 ... y_6 ] for i = 1,2
Y1 = [ 6010; 6153; 6421; 6399; 6726; 6701 ];
Y2 = [ 9422; 9300; 9220; 9150; 9042; 8800 ];

B = zeros(6,2);
% Build B corresponding to basis vectors
for i = 1:6
   B(i,1) = 50;
   B(i,2) = A(i,1) - 1065;
end

% BX = Y
% QRX = Y
% Qc = Y
% c = Q'Y
% X = R\(Q'Y)
[Q,R] = qr(B);
X1 = R \ (Q.' * Y1);
X2 = R \ (Q.' * Y2);
