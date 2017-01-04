%% Test Numerical ODE Solvers:
% Forward Euler, Backward Euler, Crank-Nicolson, & RK4 Methods

% Implement the above methods to approximately solve the given IVP problem:
% y'(t) = 6y(t) - 6y(t)^2,   for t in (0,20]
% y(0) = 0.5

%% TEST FORWARD EULER METHOD

% Test for step size h = 0.25:
fprintf( 'Testing Forward Euler Method for step size h = 0.25\n' );
h = 0.25;
[ y ] = euler_fwd( h );
disp(y);

% Test for step size h = 0.5:
fprintf( 'Testing Forward Euler Method for step size h = 0.5\n' );
h = 0.5;
[ y ] = euler_fwd( h );
disp(y);

% Test for step size h = 0.1:
fprintf( 'Testing Forward Euler Method for step size h = 0.1\n' );
h = 0.1;
[ y ] = euler_fwd( h );
disp(y);

%% TEST BACKWARD EULER METHOD

% Test for step size h = 0.25:
fprintf( 'Testing Backard Euler Method for step size h = 0.25\n' );
h = 0.25;
[ y ] = euler_bkwd( h );
disp(y);

% Test for step size h = 0.5:
fprintf( 'Testing Backward Euler Method for step size h = 0.5\n' );
h = 0.5;
[ y ] = euler_bkwd( h );
disp(y);

% Test for step size h = 0.1:
fprintf( 'Testing Backward Euler Method for step size h = 0.1\n' );
h = 0.1;
[ y ] = euler_bkwd( h );
disp(y);

%% TEST CRANK-NICOLSON METHOD

% Test for step size h = 0.25:
fprintf( 'Testing Crank-Nicolson Method for step size h = 0.25\n' );
h = 0.25;
[ y ] = crank_nicolson( h );
disp(y);

% Test for step size h = 0.5:
fprintf( 'Testing Crank-Nicolson Method for step size h = 0.5\n' );
h = 0.5;
[ y ] = crank_nicolson( h );
disp(y);

% Test for step size h = 1:
fprintf( 'Testing Crank-Nicolson Method for step size h = 0.1\n' );
h = 0.1;
[ y ] = crank_nicolson( h );
disp(y);

%% TEST RUNGE-KUTTA 4-STAGE (RK4) METHOD

% Test for step size h = 0.25:
fprintf( 'Testing RK4 Method for step size h = 0.25\n' );
h = 0.25;
[ y ] = runge_kutta( h );
disp(y);

% Test for step size h = 0.5:
fprintf( 'Testing RK4 Method for step size h = 0.5\n' );
h = 0.5;
[ y ] = runge_kutta( h );
disp(y);

% Test for step size h = 0.1:
fprintf( 'Testing RK4 Method for step size h = 0.1\n' );
h = 0.1;
[ y ] = runge_kutta( h );
disp(y);
