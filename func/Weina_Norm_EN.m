% --------------------------------------------------------------------------------------------------------------
% Weina_Norm_EN - Wiener filtering with improved noise estimation and adjustable spectral subtraction parameters
%
% This function performs Wiener filtering with improved noise estimation and spectral subtraction
% to enhance speech signals. 
%
% Syntax: enhanced = Weina_Norm_EN(x, wind, inc, NIS, alpha, beta)
%
% Inputs:
%   x        - Input speech signal
%   wind     - Window function or frame size
%   inc      - Frame overlap length
%   NIS      - Number of silent frames
%   alpha    - Suppression parameter for spectral subtraction
%   beta     - Suppression parameter for spectral subtraction
%
% Outputs:
%   enhanced - Enhanced speech signal
% ---------------------------------------------------------------------------------------------------------------

function enhanced = Weina_Norm_EN(x, wind, inc, NIS, alpha, beta)
    % Determine window length and frame size
    nwin = length(wind);
    if (nwin == 1)  % Check if window length is 1, indicating no window function is set
        framesize = wind;  % If so, frame size = wind
        wnd = hamming(framesize);  % Set default window function to Hamming
    else
        framesize = nwin;  % Otherwise, frame size = window length
        wnd = wind;
    end

    % Frame the input signal
    y = enframe(x, wnd, inc)'; % Framing
    framenum = size(y, 2); % Number of frames
    y_fft = fft(y); % FFT
    y_a = abs(y_fft); % Magnitude
    y_phase = angle(y_fft); % Phase angle
    y_a2 = y_a.^2; % Power spectrum

    % Initial noise estimation
    noise = mean(y_a2(:, 1:NIS), 2);

    % Decision-directed approach
    lambda_d = 0.99; % Smoothing factor, typically between 0.9 and 0.999

    % Process each frame
    for i = 1:framenum
        frame = y(:, i);
        y_fft = fft(frame);
        y_fft2 = abs(y_fft).^2;
        Mag_y = y_a(:, i);

        % Update noise power estimation using decision-directed approach
        if i > 1
            noise = lambda_d * noise + (1 - lambda_d) * min(y_fft2, noise);
        end

        signal = zeros(framesize, 1);
        % Spectral subtraction
        for k = 1:framesize
            if abs(y_fft2(k)) >= alpha * noise(k)
                signal(k) = y_fft2(k) - alpha * noise(k);
                if signal(k) < 0
                    signal(k) = 0; % Avoid negative power
                end
            else
                signal(k) = beta * noise(k);
            end
        end

        % Smoothing process - Simple moving average
        smoothed_signal = zeros(1, framesize);
        window_length = 5; % Smoothing window length
        for k = window_length + 1:framesize - window_length
            smoothed_signal(k) = mean(signal(k - window_length:k + window_length));
        end
        smoothed_signal(1:window_length) = mean(signal(1:window_length + 1));
        smoothed_signal(framesize - window_length + 1:end) = mean(signal(framesize - window_length:end));
        smoothed_signal = smoothed_signal';

        % Compute Wiener filter gain
        Hw = (smoothed_signal ./ (smoothed_signal + 1 * noise)).^1;

        % Apply Wiener filter to magnitude
        M = Mag_y .* Hw; % Magnitude after Wiener filtering

        % Reconstruct signal with original phase
        Mn = M .* exp(1i .* y_phase(:, i));
        yt(:, i) = real(ifft(Mn)); % Inverse FFT to time domain
    end

    % Reconstruct enhanced signal
    enhanced = filpframe(yt', wnd, inc);
end
