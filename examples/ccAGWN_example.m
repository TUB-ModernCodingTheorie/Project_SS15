clear all
addpath('../bin');

N = 1000;

A = [0 1;0 0];
B = [1 0];
C = [0 1;1 1];
D = [1 1];

display = 0;
SNR = 10; % in dB
Eb = 10^(SNR/10)/2; % Bit energy

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

%{
X_coded = zeros(1,length(c)*n);
for i = 1:length(c)
    switch c(i)
        case 0
            X_coded((i-1)*n+1) = 0;
            X_coded((i-1)*n+2) = 0;
        case 1
            X_coded((i-1)*n+1) = 1;
            X_coded((i-1)*n+2) = 0;
        case 2
            X_coded((i-1)*n+1) = 0;
            X_coded((i-1)*n+2) = 1;
        case 3
            X_coded((i-1)*n+1) = 1;
            X_coded((i-1)*n+2) = 1;
    end    
end
%}

c_str = reshape(fliplr(dec2bin(c)),1,[]);
X_coded = zeros(1,length(c_str)*n);
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
X_BPSK = sqrt(Eb)*(1 - 2*X_coded); % -sqrt(Eb) <=> 1 ; +sqrt(Eb) <=> 0

if display
    disp('code after channel coding(BPSK):');
    disp(X_BPSK);
end

%% AWG Channel
Y = X_BPSK + randn(1,numel(X_BPSK))/sqrt(2);

if display
    disp('code after AWGN channel:');
    disp(Y);
end

%% define metric matrix
metric = zeros(length(Y),2);
metric(:,1) = abs(Y(:) - sqrt(Eb)*(1 - 2*1)).^2;
metric(:,2) = abs(Y(:) - sqrt(Eb)*(1 - 2*0)).^2;

%% Decode
[X_hat, Y_hat] = ccDecode(bwd, length(c), metric, s0, sN);

disp(['BER = ' num2str(sum(Y_hat~=c)/length(Y_hat))]);
disp(['BER = ' num2str(sum(X_hat~=X)/length(X_hat))]);
