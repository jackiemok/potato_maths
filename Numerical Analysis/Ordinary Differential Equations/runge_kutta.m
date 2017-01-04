%% Runge-Kutta 4-Stage (RK4) Method Example

% Using step size h, approximate the solution to the following ODE:
% y'(t) = 6y(t) - 6y(t)^2, for t in (0,20]
% y(0) = 0.5

% Try using step sizes h = 0.1, 0.25, 0.5

function [ y ] = runge_kutta( h )

% Define function f = y'(t)
% Print initial conditions
f = @(t,y) 6 * (y - y^2);
fprintf( 'Step 0: t = 0, y(0) = 0.5\n' );

steps = 20 / h;
t = zeros(1,(steps + 1));
y = zeros(1,(steps + 1));
y(1) = 0.5;

 for i = 1:steps
     
    % Increment time step
    t(i+1) = t(i) + h;
     
    % Compute increments and slopes
    k1 = f( t(i), y(i) );
    k2 = f( t(i) + h/2, y(i) + h/2*k1 );
    k3 = f( t(i) + h/2, y(i) + h/2*k2 );
    k4 = f( t(i) + h, y(i) + h*k3 );
    
    % Update approximate solution and print current result
    y(i+1) = y(i) + h/6 * ( k1 + 2*k2 + 2*k3 + k4 );
    fprintf( 'Step %d: t = %15.15f, y%d = %15.15f\n', i, t(i+1), i, y(i+1) );
        
 end

% Plot results
hold on;
title('RK4 Method')
plot( t, y, 'Color', 'Blue' );
xlabel('t');
ylabel('y(t)');
 
end