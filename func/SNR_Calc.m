function snr = SNR_Calc(I, In)
% SNR_Calc - Calculate the Signal-to-Noise Ratio (SNR) of a noisy speech signal.
%
% This function calculates the SNR for a given clean and noisy speech signal.
%
% Syntax: snr = SNR_Calc(I, In)
%
% Inputs:
%   I  - Clean speech signal
%   In - Noisy speech signal
%
% Outputs:
%   snr - Calculated SNR in decibels (dB)

% Ensure both signals are row vectors
I = I(:)';    
In = In(:)';

% Calculate the energy of the clean signal
Ps = sum((I - mean(I)).^2);          

% Calculate the energy of the noise
Pn = sum((I - In).^2);               

% Calculate the SNR (in dB)
snr = 10 * log10(Ps / Pn);           

end
