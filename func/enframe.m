function frameout = enframe(x, win, inc)

% Length of input signal
nx = length(x(:));

% Length of window function
nwin = length(win);

if (nwin == 1)  % If the window length is 1, it indicates no window function is set (win is just a number)
   len = win;   % Set frame length to win
else
   len = nwin;  % Otherwise, set frame length to window length
end

if (nargin < 3)  % If only two parameters are provided, set frame increment equal to frame length
   inc = len;
end

% Subtract one frame to avoid filling the last frame incompletely
nf = fix((nx - len + inc) / inc);  % Calculate the number of frames
frameout = zeros(nf, len);         % Initialize the output matrix

% Calculate the starting index of each frame
indf = inc * (0:(nf-1)).';  % Set the displacement positions of each frame in x

% Create an index array with a length equal to the frame length
inds = (1:len);  % Indices corresponding to each frame

% Extract each frame's data using indexing and store it in the output matrix frameout
frameout(:) = x(indf(:, ones(1, len)) + inds(ones(nf, 1), :));  % Frame data

% Apply window function (if available)
if (nwin > 1)  % If a window function is included in the parameters, multiply each frame by the window function
    w = win(:)';  % Convert win to a row vector
    frameout = frameout .* w(ones(nf, 1), :);  % Apply window function
end

end
