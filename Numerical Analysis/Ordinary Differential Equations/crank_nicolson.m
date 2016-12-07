%% Crank-Nicolson Method Example

% Using step size h, approximate the solution to the following ODE:
% y'(t) = 6y(t) - 6y(t)^2, for t in (0,20]
% y(0) = 0.5

% Try using step sizes h = 0.1, 0.25, 0.5

function crank_nicolson( h )

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
    
    for j = 1:steps
        % Update approximate solution and print result
        y(j+1) = y(j) + h/2 * f( t(j), y(j) ) + f( t(j+1), y(j+1) );
    end
    
    % Print result for current step
    fprintf( 'Step %d: t = %15.15f, y%d = %15.15f\n', i, t(i+1), i, y(i+1) );
    
end
 
% Plot results
hold on;
title('Crank Nicolson Method')
plot( t, y, 'Color', 'Blue' );
xlabel('t');
ylabel('y(t)');

end