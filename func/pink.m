function pink_noise = pink(N)
% pink - Generates pink noise of length N.
%
% Syntax: pink_noise = pink(N)
%
% Inputs:
%   N - The length of the pink noise vector
%
% Outputs:
%   pink_noise - Generated pink noise

    % Generate white noise
    white_noise = randn(N, 1);

    % Use an FIR filter to change the spectrum, making it 1/f
    b = fir1(10, 0.1, 'high');  % High-pass filter design
    pink_noise = filter(b, 1, white_noise);

    % Normalize the noise
    pink_noise = pink_noise / std(pink_noise);
end
