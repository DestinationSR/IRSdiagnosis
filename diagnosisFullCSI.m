%% Diagnosis with full CSI
clear all
close all

rng(1) % random seed

%% Parameter Setting
H = 16;         % num of reflecting elements along the vertical direction
W = 16;         % num of reflecting elements along the horizontal direction
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
thetaIRSin = rand*pi - pi/2;
phiIRSin = rand*pi - pi/2;

h_TX = kron(exp(1j*2*pi*0.5*(0:W-1)'*sin(thetaIRSin)*cos(phiIRSin)),...
    exp(1j*2*pi*0.5*(0:H-1)'*sin(thetaIRSin)*sin(phiIRSin)));


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
ideal_y = F*(h_TX.*h_RX);
tic
m_rec = lasso_ADMM(F*diag(h_TX.*h_RX), y - ideal_y, 0.65*K*(10^(-SNR/20))^2) + 1;
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