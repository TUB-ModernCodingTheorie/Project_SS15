path('../bin/',path);

forwardTrellis = uint64(zeros(2,2,2));

forwardTrellis(1,1,1) = 0; %s=0, u=0
forwardTrellis(1,1,2) = 1;

forwardTrellis(1,2,1) = 1; %s=0, u=1
forwardTrellis(1,2,2) = 1;

forwardTrellis(2,1,1) = 1; %s=1, u=0
forwardTrellis(2,1,2) = 0;

forwardTrellis(2,2,1) = 0; %s=1, u=1
forwardTrellis(2,2,2) = 0;

s0 = [1];

seq = uint64([0,1,0,1]);

[c, sN] = ccEncode(forwardTrellis,s0,seq);