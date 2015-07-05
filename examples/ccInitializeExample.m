path('../bin/',path);

A = [0,1;0,0];
B = [1,0];
C = [0,1;1,1];
D = [1,1];

[fwd, bwd] = ccInitialize(A,B,C,D);