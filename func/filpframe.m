function frameout = filpframe(x, win, inc)
% filpframe - Reconstructs a signal from framed data.
%
% Syntax: frameout = filpframe(x, win, inc)
%
% Inputs:
%   x   - Framed signal matrix, where each row represents a frame
%   win - Window function or a scalar (if no windowing was applied)
%   inc - Frame increment (hop size)
%
% Outputs:
%   frameout - Reconstructed signal from the framed data

% Get the number of frames (nf) and frame length (len)
[nf, len] = size(x);

% Calculate the length of the original signal
nx = (nf - 1) * inc + len;                 

% Initialize output signal
frameout = zeros(nx, 1);

% Get the length of the window function
nwin = length(win);                   

% Check if a window function was applied
if (nwin ~= 1)                        
    % Remove the effect of the window function
    winx = repmat(win', nf, 1);       % Replicate window vector to match frame size
    x = x ./ winx;                    % Divide frames by the window
    x(isinf(x)) = 0;                  % Remove infinities resulted from division by zero
end

% Initialize the signal for reconstruction
sig = zeros((nf - 1) * inc + len, 1);

% Reconstruct the signal from frames
for i = 1:nf
    start = (i - 1) * inc + 1;    
    xn = x(i, :)';
    sig(start:start + len - 1) = sig(start:start + len - 1) + xn;
end

% Output the reconstructed signal
frameout = sig;

end
