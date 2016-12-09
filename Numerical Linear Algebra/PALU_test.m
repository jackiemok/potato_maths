%% Partial Pivoting Decomposition: Test Cases
% Test cases/examples for PALU_decomp.m

% Matrices to test factorization method
A_1 = [2 1 2; 1 2 3; 4 1 2];
A_2 = [10 1 1; 1 10 1; 1 1 20];
A_3 = hilb(5);
A_4 = hilb(10);
A_5 = hilb(20);

%% PA = LU Decomposition for each A_i
% Check Frobenius norms

% A_1
[P,L,U] = PALU_decomp(A_1);
error_1 = norm(P*A_1 - L*U,'fro') / (norm(L, 'fro')*norm(U,'fro'));
error_2 = norm(P*A_1 - L*U,'fro') / norm(A_1,'fro');
disp(error_1); disp(error_2);

% A_2
[P,L,U] = PALU_decomp(A_2);
error_3 = norm(P*A_2 - L*U,'fro') / (norm(L, 'fro')*norm(U,'fro'));
error_4 = norm(P*A_2 - L*U,'fro') / norm(A_2,'fro');
disp(error_3); disp(error_4);

% A_3
[P,L,U] = PALU_decomp(A_3);
error_5 = norm(P*A_3 - L*U,'fro') / (norm(L, 'fro')*norm(U,'fro'));
error_6 = norm(P*A_3 - L*U,'fro') / norm(A_3,'fro');
disp(error_5); disp(error_6);

% A_4
[P,L,U] = PALU_decomp(A_4);
error_7 = norm(P*A_4 - L*U,'fro') / (norm(L, 'fro')*norm(U,'fro'));
error_8 = norm(P*A_4 - L*U,'fro') / norm(A_4,'fro');
disp(error_7); disp(error_8);

% A_5
[P,L,U] = PALU_decomp(A_5);
error_9 = norm(P*A_5 - L*U,'fro') / (norm(L, 'fro')*norm(U,'fro'));
error_10 = norm(P*A_5 - L*U,'fro') / norm(A_5,'fro');
disp(error_9); disp(error_10);