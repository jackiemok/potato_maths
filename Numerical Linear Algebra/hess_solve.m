%% Solve Upper-Hessenberg Linear System
% Solve the linear system Hx = b, where H is upper-Hessenberg.

% Input:              H - (n x n) Upper-Hessenberg matrix
%                     b - (n x 1) RHS vector equivalent to H*x
% Output:             b - (n x 1) Vector x stored as b, where H*x = b

function [ b ] = hess_solve( H, b )

n = length(b);

% Check dimensions of system
[ m1, m2 ] = size(H); 
if ( m1 ~= n || m2 ~= n )
    disp('Error: Dimensions are not correct!')
    return
end

for k = 1:(n - 1)
    x = H(k:k+1,k);
    
    v = x; 
    v(1) = v(1) + norm(x) * sign(x(1));
    v = v / norm(v);
    
    H(k:k+1,k:n) = H(k:k+1,k:n) - (2 * v) * (v' * H(k:k+1,k:n));
    b(k:k+1) = b(k:k+1) - (2 * v) * (v' * b(k:k+1));
end

% Back Substitution
for j = n:-1:1
    b(j) = b(j) - H(j,j+1:n) * b(j+1:n);
  
    % Is H nearly singular?
    if abs( H(j,j) ) < ( 10 * eps )
        disp( 'Warning: Matrix H may be singular.' )
    end
    
    b(j) = b(j) / H(j,j);
end
end