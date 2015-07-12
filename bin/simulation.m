function simulation()
A = [0 1;0 0];
B = [1 0];
C = [0 1;1 1];
D = [1 1];
m = 2;
n = 2;
k = 1;
SNR = 3; % 3dB
[fwd bwd] = ccInitialize(A,B,C,D);
% generating random bits
input = randi([0 1],1,10); 

disp('Input bits:')
disp(input)
% initial states
s0 = 0;
[c,sN] = ccEncode(fwd,s0,input);
disp('output code:')
disp(c)
disp('last state:')
disp(sN)

channelCode = zeros(1,length(c)*n);
for i=1:length(c)
    tmp = de2bi(c(i),n);
    channelCode(2*i-1) = tmp(1);
    channelCode(2*i) = tmp(2);
end
% output code after channel coding mapping into 2BSK
channelCode(channelCode==0)=-1;
disp('code after channel coding(2ASK):')
disp(channelCode)

% add white Gaussian noise
disp('code after AWGN channel:')
channelCodeNoise = awgn(channelCode,SNR);
disp(channelCodeNoise)










end

function    b=de2bi(d,K)

two=repmat(2,size(d));
for k=1:K
    b(:,k)=mod(d,two);
    d=floor(d./two);
end
end
