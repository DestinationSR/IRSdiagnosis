function [H_RX, H_TX] = channelGen(H, W, NantRX, NantTX, NrayRX, NrayTX)
%
% MmWave narrowband clustered channel generator
%
% INPUTS:
%   H       : number of reflecting elements along the vertical direction
%   W       : number of reflecting elements along the horizontal direction
%   NantRX  : number of antennas at the RX
%   NantTX  : number of antennas at the TX
%   NrayRX  : number of sub-paths in the IRS-RX channel
%   NrayTX  : number of sub-paths in the TX-IRS channel
%
% OUTPUTS:
%   H_RX    : IRS-RX channel matrix
%   H_TX    : TX-IRS channel matrix

%% IRS-RX channel
alphaRX = (randn(NrayRX,1) + 1j*randn(NrayRX,1))/sqrt(2); % sub-path gain
thetaRXin = rand(NrayRX,1)*pi - pi/2; % AoA
thetaIRSout = rand(NrayRX,1)*pi - pi/2; % elevation AoD
phiIRSout = rand(NrayRX,1)*pi - pi/2; % azimuth AoD

H_RX = 0;
for l=1:NrayRX
    
    H_RX = H_RX + alphaRX(l)*exp(1j*2*pi*0.5*(0:NantRX-1)'*sin(thetaRXin(l)))*...
        kron(exp(1j*2*pi*0.5*(0:W-1)'*sin(thetaIRSout(l))*cos(phiIRSout(l))),...
        exp(1j*2*pi*0.5*(0:H-1)'*sin(thetaIRSout(l))*sin(phiIRSout(l)))).' ./ sqrt(NrayRX);
    
end

%% TX-IRS channel
alphaTX = (randn(NrayTX,1) + 1j*randn(NrayTX,1))/sqrt(2); % sub-path gain
thetaTXout = rand(NrayTX,1)*pi - pi/2; % AoD
thetaIRSin = rand(NrayTX,1)*pi - pi/2; % elevation AoA
phiIRSin = rand(NrayTX,1)*pi - pi/2; % azimuth AoA

H_TX = 0;
for l=1:NrayTX
    
    H_TX = H_TX + alphaTX(l)*kron(exp(1j*2*pi*0.5*(0:W-1)'*sin(thetaIRSin(l))*cos(phiIRSin(l))),...
        exp(1j*2*pi*0.5*(0:H-1)'*sin(thetaIRSin(l))*sin(phiIRSin(l))))*...
        exp(1j*2*pi*0.5*(0:NantTX-1)'*sin(thetaTXout(l))).' ./ sqrt(NrayTX);
    
end

end