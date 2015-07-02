path('../bin/',path);

forwardTrellis = uint64(zeros(2,2,2));

forwardTrellis(1,1,1) = 01; %s=0, u=0
forwardTrellis(1,1,2) = 1;

forwardTrellis(1,2,1) = 11; %s=0, u=1
forwardTrellis(1,2,2) = 1;

forwardTrellis(2,1,1) = 10; %s=1, u=0
forwardTrellis(2,1,2) = 0;

forwardTrellis(2,2,1) = 00; %s=1, u=1
forwardTrellis(2,2,2) = 0;

s0 = [1];

seq = [1,1,1,0];

[c, sN] = ccEncode(forwardTrellis,s0,seq,1,2);