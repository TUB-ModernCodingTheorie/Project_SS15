A = [0 1;0 0];
B = [1 0];
C = [0 1;1 1];
D = [1 1];

display = 0;
SNR = 3; % in dB
Eb = 1; % Bit energy

[fwd bwd] = ccInitialize(A,B,C,D);

n = fwd.ldOutputs;
k = fwd.ldInputs;

%% generating random bits
input = randi([0 1],1,10); 

if display
    disp('Input bits:')
    disp(input)
end

%% Encoding
s0 = 0;
[c,sN] = ccEncode(fwd,input,s0);

if display
    disp('output code:')
    disp(c)
    disp('last state:')
    disp(sN)
end

%% 2BSK
channelCode_str = reshape(dec2bin(c,n),1,[]);
for i = 1:numel(channelCode_str)
    channelCode(i) = str2num(channelCode_str(i));
end

channelCode_2BSK = Eb*(2*channelCode(:)'-1);

if display
    disp('code after channel coding(2BSK):');
    disp(channelCode_2BSK);
end

%% AWG Channel
Y = channelCode_2BSK + randn(1,numel(channelCode));

if display
    disp('code after AWGN channel:');
    disp(Y);
end

%% define metric matrix
