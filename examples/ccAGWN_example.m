clear all
addpath('../bin');
N = 5;

A = [0 1;0 0];
B = [1 0];
C = [0 1;1 1];
D = [1 1];

display = 0;
SNR = 3; % in dB
Eb = 1; % Bit energy


%   X     X_coded    X_coded_BPSK     Y       X_hat
%  -->[enc]----->[BPSK]---->[Channel]--->[dec]-->

%% Initialization
[fwd, bwd] = ccInitialize(A,B,C,D);

n = fwd.ldOutputs;
k = fwd.ldInputs;

%% generating random bits
X = randi([0 1],1,N); 

if display
    disp('Input bits:')
    disp(X)
end

%% Encoding
s0 = 0;
[c,sN] = ccEncode(fwd,X,s0);

c_str = reshape(dec2bin(c,n),1,[]);
X_coded = zeros(1,numel(c_str));
for i = 1:numel(c_str)
    X_coded(i) = str2double(c_str(i));
end

if display
    disp('output code:')
    disp(X_coded)
    disp('last state:')
    disp(sN)
end

%% 2BSK
X_coded_BPSK = Eb*(2*X_coded(:)'-1);

if display
    disp('code after channel coding(BPSK):');
    disp(X_coded_BPSK);
end

%% AWG Channel
Y = X_coded_BPSK + randn(1,numel(X_coded_BPSK));

if display
    disp('code after AWGN channel:');
    disp(Y);
end

%% define metric matrix
metric(:,1) = normpdf(Y+Eb);    % prob that we the Y the encoded bit was 0
metric(:,2) = normpdf(Y-Eb);    % prob that we the Y the encoded bit was 1

%% Decode
[X_hat, Y_hat] = ccDecode(bwd, Y, metric, s0, sN);

disp(['BER = ' num2str(sum(Y_hat~=X_coded)/length(Y_hat))]);
disp(['BER = ' num2str(X_hat~=X)]);
