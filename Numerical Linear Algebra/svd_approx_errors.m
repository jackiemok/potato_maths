%% Jacqueline Mok

% Quick Demo: 
% Store and plot Frobenius norm errors corresponding to
% k-valued target rank matrix approximations to a given
% (random) matrix A.


% Select matrix dimensions
m = 500;
n = 1000;

% Create random matrix
A = randn(m, n);
save('500x1000_matrix.txt', 'A', '-ASCII');
r = rank(A);

% Compute svd of matrix A
[U, S, V] = svd(A);
save('500x1000_U.txt', 'U', '-ASCII');
save('500x1000_S.txt', 'S', '-ASCII');
save('500x1000_V.txt', 'V', '-ASCII');

% Store all singular values
sigmas = svds(A, r);

% Initialize vector storing approximation errors
rank_errors = zeros(r, 1);

% Find optimal k-rank-value for approximating A
for k = 1 : r
    
    % Calculate Frobenius norm error: ||A - A_k||_F
    temp = sigmas(k+1:r);
    temp = sum(temp .* temp);
    temp = sqrt(temp);
    
    rank_errors(k) = temp;
end

save('500x1000_rank_errors.txt', 'rank_errors', '-ASCII');
plot(rank_errors)