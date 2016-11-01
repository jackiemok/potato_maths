%% Test Cases: Banded Matrix LU Decomposition without Pivoting
% A few simple cases for testing the function 'banded_LU_solve_nopivot'

%% s_L = 1 = s_U:

A_1 = [1 2 0 0; 3 4 2 0; 0 5 7 1; 0 0 4 3];
A_2 = [7 3 0 0; 1 8 4 0; 0 9 2 5; 0 0 7 6];

b_1 = [1; 4; 3; 6];
b_2 = [2; 6; 9; 1];

%% s_L = 1; s_U = 2:

A_3 = [1 2 3 0 0 0; 5 8 9 2 0 0; 0 5 4 3 1 0; 0 0 1 3 5 9; 0 0 0 8 6 4; 0 0 0 0 2 6];
A_4 = [9 2 7 0 0 0; 1 5 3 4 0 0; 0 10 9 6 1 0; 0 0 12 3 9 19; 0 0 0 8 16 4; 0 0 0 0 2 16];

b_3 = [1; 3; 8; 5; 4; 6];
b_4 = [2; 5; 10; 4; 1; 3];

%% Test Function 'banded_LU_solve_nopivot'

[ L1, U1 ] = banded_LU_solve_nopivot( A_1, b_1, 1, 1 )
[ L2, U2 ] = banded_LU_solve_nopivot( A_2, b_2, 1, 1 )
[ L3, U3 ] = banded_LU_solve_nopivot( A_3, b_3, 2, 1 )
[ L4, U4 ] = banded_LU_solve_nopivot( A_4, b_4, 2, 1 )
