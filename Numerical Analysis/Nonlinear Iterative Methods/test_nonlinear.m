%% Test Nonlinear Iterative Methods: Newton, Broyden I, & Broyden II
% Solve F(x*) = 0 satisfying G(x*) = x*

% Example: Consider the following system.
% x1^2 + x2^2 = 1
% (x1 - 1)^2 + x2^2 = 1

% The system has two solutions: 
% ( 1/2, +sqrt(3)/2 ) & ( 1/2, -sqrt(3)/2 )

%% Newton's Method

num_iterations = 10; 

% Initial guess:
% Choose close to solution ( 1/2, +sqrt(3)/2 )
x_0 = [ 1 ; 1 ];

[ sol, num_iterations ] = newton( x_0, num_iterations );
disp(sol); disp(num_iterations);

% Initial guess:
% Choose close to solution ( 1/2, -sqrt(3)/2 )
x_0 = [ 1 ; -1 ];

[ sol, num_iterations ] = newton( x_0, num_iterations );
disp(sol); disp(num_iterations);

%% Broyden I Method

num_iterations = 10; 

% Initial guess:
% Choose close to solution ( 1/2, +sqrt(3)/2 )
x_0 = [ 1 ; 1 ];

[ sol , num_iterations ] = broyden1( x_0, num_iterations );
disp(sol); disp(num_iterations);

% Initial guess:
% Choose close to solution ( 1/2, -sqrt(3)/2 )
x_0 = [ 1 ; -1 ];

[ sol , num_iterations ] = broyden1( x_0, num_iterations );
disp(sol); disp(num_iterations);

%% Broyden II Method

num_iterations = 10; 

% Initial guess:
% Choose close to solution ( 1/2, +sqrt(3)/2 )
x_0 = [ 1 ; 1 ];

[ sol, num_iterations ] = broyden2( x_0, num_iterations );
disp(sol); disp(num_iterations);

% Initial guess:
% Choose close to solution ( 1/2, -sqrt(3)/2 )
x_0 = [ 1 ; -1 ];

[ sol, num_iterations ] = broyden2( x_0, num_iterations );
disp(sol); disp(num_iterations);
