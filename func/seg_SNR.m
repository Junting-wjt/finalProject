function segSNR = seg_SNR(I, I_hat, N)
% seg_SNR - Calculate segmental SNR (Signal-to-Noise Ratio).
%
% This function calculates the segmental SNR for a given clean and estimated (or noisy) speech signal.
%
% Syntax: segSNR = seg_SNR(I, I_hat, N)
%
% Inputs:
%   I     - Clean speech signal
%   I_hat - Estimated or noisy speech signal
%   N     - Length of each segment
%
% Outputs:
%   segSNR - Calculated segmental SNR

% Ensure both signals are column vectors
I = I(:);    
I_hat = I_hat(:);

% Initialize the total segmental SNR sum
totalSegSNR = 0;
numSegments = floor(length(I) / N); % Calculate the number of complete segments

% Loop through each segment
for m = 0:numSegments-1
    startIdx = m * N + 1;
    endIdx = startIdx + N - 1;
    
    % Extract the current segment of the signals
    segment_I = I(startIdx:endIdx);
    segment_I_hat = I_hat(startIdx:endIdx);
    
    % Calculate the energy of each segment
    Ps = sum((segment_I - mean(segment_I)).^2);  % Energy of the clean signal
    Pn = sum((segment_I - segment_I_hat).^2);    % Energy of the noise
    
    % Avoid division by zero
    if Pn == 0
        currentSegSNR = Inf;  % If noise energy is zero, SNR is infinity
    else
        currentSegSNR = 10 * log10(Ps / Pn);  % Calculate SNR for the current segment
    end
    
    % Accumulate the segmental SNR
    totalSegSNR = totalSegSNR + currentSegSNR;
end

% Calculate the average segmental SNR
segSNR = totalSegSNR / numSegments;

end
