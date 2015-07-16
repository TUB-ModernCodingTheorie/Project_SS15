clear all
path('../bin/',path);

p = 0:0.01:0.5;
N = 10^5;

A = [0,1;0,0];
B = [1,0];
C = [0,1;1,1];
D = [1,1];

[fwd, bwd] = ccInitialize(A,B,C,D);

s0 = 0;
seq = round(rand(1,N));

[c, sN] = ccEncode(fwd,seq,s0);

%% BSC
randomNoise = rand(1,N*2);
for idx = 1:numel(p)
    
    channel = zeros(1,N*2);
    channel(randomNoise < p(idx)) = 1;

    c_b = reshape(de2bi(c)',1,[]);
    tmp = c_b;

    c_b(randomNoise <= p(idx) & tmp == 1) = 0;
    c_b(randomNoise <= p(idx) & tmp == 0) = 1;
    %% Channel
    Bmetric = zeros(length(c_b),2);

    for i = 1:length(c_b)
        if c_b(i) == 1
            Bmetric(i,1) = log(p(idx));
            Bmetric(i,2) = log(1-p(idx));
        else
            Bmetric(i,1) = log(1-p(idx));
            Bmetric(i,2) = log(p(idx));
        end
    end
    %% Decode
    [m , ch] = ccDecode(bwd,length(c),Bmetric,s0,sN);

    %% BER
    BER(idx) = length(find(m ~= seq))/length(m);
end

semilogy(p,BER);
xlabel('Error probability');
ylabel('BER')
title('BSC Channel')