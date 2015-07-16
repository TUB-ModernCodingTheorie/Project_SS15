clear all
path('../bin/',path);

p = 0.3

A = [0,1;0,0];
B = [1,0];
C = [0,1;1,1];
D = [1,1];

[fwd, bwd] = ccInitialize(A,B,C,D);

s0 = 0;

seq = round(rand(1,10000));

[c, sN] = ccEncode(fwd,seq,s0);

c_b = reshape(de2bi(c)',1,[]);
%% Channel
Bmetric = zeros(length(c)*2,2);
for i = 0:length(c)-1
    switch c(i+1)
        case 0 %00
            Bmetric(i*2+1,1) = log(1-p);
            Bmetric(i*2+1,2) = log(p);
            Bmetric(i*2+2,1) = log(1-p);
            Bmetric(i*2+2,2) = log(p);
        case 1 %10
            Bmetric(i*2+1,1) = log(p);
            Bmetric(i*2+1,2) = log(1-p);
            Bmetric(i*2+2,1) = log(1-p);
            Bmetric(i*2+2,2) = log(p);
        case 2 %01
            Bmetric(i*2+1,1) = log(1-p);
            Bmetric(i*2+1,2) = log(p);
            Bmetric(i*2+2,1) = log(p);
            Bmetric(i*2+2,2) = log(1-p);
        case 3 %11
            Bmetric(i*2+1,1) = log(p);
            Bmetric(i*2+1,2) = log(1-p);
            Bmetric(i*2+2,1) = log(p);
            Bmetric(i*2+2,2) = log(1-p);
    end
end
%% Decode
[m , ch] = ccDecode(bwd,length(c),Bmetric,s0,sN);

%% BER
display(['Result (message): ',num2str(length(find(m ~= seq))/length(m)) '% of error']);