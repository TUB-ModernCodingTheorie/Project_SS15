addpath('../bin');

N = 10^6;

A = [0 1;0 0];
B = [1 0];
C = [0 1;1 1];
D = [1 1];

display = 0;

SNR_tab = -30:10; % in dB

for idx = 1:numel(SNR_tab)
    
    SNR = SNR_tab(idx);    
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

    X_coded = reshape(de2bi(c)',1,[]);

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
    metric(:,2) = -abs(Y - sqrt(Eb)*(1 - 2*1)).^2;
    metric(:,1) = -abs(Y - sqrt(Eb)*(1 - 2*0)).^2;

    %% Decode
    tic
    [X_hat, Y_hat] = ccDecode(bwd, length(c), metric, s0, sN);
    toc
    BER_tab(idx) = sum(X_hat ~= X)/length(X_hat);
end

semilogy(SNR_tab,BER_tab,'LineWidth',2);
grid on
xlabel('SNR (dB)');
ylabel('BER');
ylim([1e-5 1e0]);
title('Transmission over AWGN Channel')