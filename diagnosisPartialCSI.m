%% Diagnosis with partial CSI
clear all
close all

rng(1) % random seed

%% Parameter Setting
H = 16;         % num of reflecting elements along the vertical direction
W = 16;         % num of reflecting elements along the horizontal direction
NrayTX = 4;     % num of sub-paths in TX-IRS channel
nFault = 3;     % num of faulty reflecting elements
Kfactor = 0.8;  % compression ratio
SNR = 20;       % dB
psSet = [1+1j,1-1j,-1+1j,-1-1j]/sqrt(2); % 2-bit phase shifting

%% Channel
% IRS-RX channel
thetaIRSout = rand*pi - pi/2;
phiIRSout = rand*pi - pi/2;

h_RX = kron(exp(1j*2*pi*0.5*(0:W-1)'*sin(thetaIRSout)*cos(phiIRSout)),...
    exp(1j*2*pi*0.5*(0:H-1)'*sin(thetaIRSout)*sin(phiIRSout)));

% TX-IRS channel
alphaTX = (randn(NrayTX,1) + 1j*randn(NrayTX,1))/sqrt(2);
thetaTXout = rand(NrayTX,1)*pi - pi/2;
thetaIRSin = rand(NrayTX,1)*pi - pi/2;
phiIRSin = rand(NrayTX,1)*pi - pi/2;
h_TX = 0;

for l=1:NrayTX
    
    h_TX = h_TX + alphaTX(l)*kron(exp(1j*2*pi*0.5*(0:W-1)'*sin(thetaIRSin(l))*cos(phiIRSin(l))),...
        exp(1j*2*pi*0.5*(0:H-1)'*sin(thetaIRSin(l))*sin(phiIRSin(l)))) ./ sqrt(NrayTX);
    
end

%% Failure Masks
fMask = ones(H*W,1);
fIndex = randperm(H*W,nFault);
fMask(fIndex) = rand(nFault,1).*exp(1j*2*pi*rand(nFault,1));
fMask = diag(fMask);

%% Measurements
K = ceil(Kfactor*H*W);
F = psSet(randi(4,K,H*W));
y = zeros(K,1);
for k=1:K
    
    Theta = diag(F(k,:));
    y(k,:) = h_RX.'*Theta*fMask*h_TX + (10^(-SNR/20))*(randn + 1j*randn)/sqrt(2);
    
end

%% Recovery
tic
[h_rec, m_rec] = cSLRMD(H, W, F*diag(h_RX), y, 0.35, sqrt(K)*10^(-SNR/20));
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