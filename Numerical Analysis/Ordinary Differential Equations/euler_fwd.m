%% Forward Euler Method Example

% Using step size h, approximate the solution to the following ODE:
% y'(t) = 6y(t) - 6y(t)^2, for t in (0,20]
% y(0) = 0.5

% Try using step sizes h = 0.1, 0.25, 0.5

function euler_fwd( h )

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
    
    % Update approximate solution and print current result
    y(i+1) = y(i) + h * f( t(i), y(i) );
    fprintf( 'Step %d: t = %15.15f, y%d = %15.15f\n', i, t(i+1), i, y(i+1) );
    
end
 
% Plot results
hold on;
title('Forward Euler Method')
plot( t, y, 'Color', 'Blue' );
xlabel('t');
ylabel('y(t)'); 
    
end