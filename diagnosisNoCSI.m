%% Diagnosis without CSI
clear all
close all

rng(1) % random seed

%% Parameter Setting
H = 16;         % num of reflecting elements along the vertical direction
W = 16;         % num of reflecting elements along the horizontal direction
NantRX = 4;     % number of antennas at the RX
NantTX = 1;     % number of antennas at the TX
NrayRX = 4;     % num of sub-paths in IRS-RX channel
NrayTX = 4;     % num of sub-paths in TX-IRS channel
nFault = 3;     % num of faulty reflecting elements
Kfactor = 0.8;  % compression ratio
SNR = 20;       % dB
maxIter = 1000; % max num of iterations of ADMM algorithm
psSet = [1+1j,1-1j,-1+1j,-1-1j]/sqrt(2); % 2-bit phase shifting

%% Channel
[H_RX, H_TX] = channelGen(H, W, NantRX, NantTX, NrayRX, NrayTX);

%% Failure Masks
fMask = ones(H*W,1);
fIndex = randperm(H*W,nFault);
fMask(fIndex) = rand(nFault,1).*exp(1j*2*pi*rand(nFault,1));
fMask = diag(fMask);

%% Measurements
K = ceil(Kfactor*H*W);
F = psSet(randi(4,K,H*W));
y = zeros(K,NantRX);
for k=1:K
    
    Theta = diag(F(k,:));
    y(k,:) = H_RX*Theta*fMask*H_TX*ones(NantTX,1) + (10^(-SNR/20))*(randn(NantRX,1) + 1j*randn(NantRX,1))/sqrt(2);
    
end

%% Recovery
rho = 1; % penalty parameter
lambda = K*0.006/sqrt(NantRX); % regularization parameter of failure mask
tau = K*0.004*sqrt(NantRX); % regularization parameter of channel

tic
if(NantRX == 1)
    [h_rec, m_rec] = ANM_smv(H, W, rho, lambda, tau, maxIter, F, y);
else
    [Hchan_rec, m_rec] = ANM_mmv(H, W, NantRX, rho, lambda, tau, maxIter, F, y);
end
runningTime = toc;

%% Plot results
close all
fMaskMat = reshape(diag(fMask),H,W);

figure(2)
mesh(abs(fMaskMat - 1))
title('Failure masks')

figure(3)
mesh(reshape(abs(m_rec - 1),H,W))
title('Recovered masks')

NMSE = 10*log10(norm(m_rec - fMaskMat(:))^2/norm(fMaskMat(:))^2) % dB
runningTime