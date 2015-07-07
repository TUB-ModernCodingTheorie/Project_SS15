path('../bin/',path);

A = [0,1;0,0];
B = [1,0];
C = [0,1;1,1];
D = [1,1];

[fwd, bwd] = ccInitialize(A,B,C,D);

s0 = 1;

seq = [1,1];

[c, sN] = ccEncode(fwd,s0,seq);

Bmetric = zeros(4,2);
Bmetric(1,1) = 0;
Bmetric(1,2) = 1;
Bmetric(2,1) = 1;
Bmetric(2,2) = 0;
Bmetric(3,1) = 1;
Bmetric(3,2) = 0;
Bmetric(4,1) = 0;
Bmetric(4,2) = 1;

[m , ch] = ccDecode(bwd,c,Bmetric,s0,sN);