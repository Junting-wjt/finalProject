% -------------------------------------------------------------------- 
% Find the optimal parameters. (It takes some time to run)
% --------------------------------------------------------------------  
close all;
clc;
warning off;
addpath(genpath(pwd));

% Folder paths
cleanPath = 'audio\clean\clean';
noisyPath = 'audio\airport_0dB\0dB';

% File listings
cleanFiles = dir(fullfile(cleanPath, '*.wav'));
noisyFiles = dir(fullfile(noisyPath, '*.wav'));

% Initialization
IS = 0.25;              % Initial silence duration
wlen = 200;             % Frame length
inc = 80;               % Frame step
%SNR = 5;               % Set signal-to-noise ratio (SNR)

% Define ranges for alpha and beta
alpha_range = 0:0.5:3;
beta_range = 0:0.1:2;

% Matrix to store STOI values for each alpha and beta combination
stoi_grid = zeros(length(alpha_range), length(beta_range));

% Main loop over files
for fileIdx = 1:length(cleanFiles)
    cleanFile = fullfile(cleanFiles(fileIdx).folder, cleanFiles(fileIdx).name);
    noisyFile = fullfile(noisyFiles(fileIdx).folder, noisyFiles(fileIdx).name);

    % Read audio files
    [x, fs] = audioread(cleanFile);
    %[signal, ~] = awgn(x, SNR, 'measured', 'db');
    [signal, ~] = audioread(noisyFile);

    %%
    % pink_noise = pink(length(x));  % Generate pink noise
    % 
    % % Calculate the power of the original signal to set an appropriate noise level
    % signalPower = sum(x.^2) / length(x);
    % %noisePower = signalPower / 10^(5/10);  % Assume desired SNR is 5dB
    % noisePower = signalPower / 10^(0/10);  % Assume desired SNR is 0dB
    % 
    % % Adjust noise level
    % adjustedPinkNoise = pink_noise * sqrt(noisePower / (sum(pink_noise.^2) / length(pink_noise)));
    % 
    % % Add noise to signal
    % signal = x + adjustedPinkNoise;

    %%
    % Preprocess signal
    minLength = min(length(x), length(signal));
    x = x(1:minLength);
    signal = signal(1:minLength);
    x = x - mean(x);
    x = x / max(abs(x));
    signal = signal / max(abs(signal));

    NIS = fix((IS * fs - wlen) / inc + 1);  % Number of initial silence frames

    % Grid search over alpha and beta
    for i = 1:length(alpha_range)
        for j = 1:length(beta_range)
            alpha = alpha_range(i);
            beta = beta_range(j);
            
            % Wiener filtering
            %output = Weina_Norm_OPP(signal, wlen, inc, NIS, alpha, beta);
            output = Weina_Norm_ENfindhuan(signal, wlen, inc, NIS, alpha, beta);
            output = real(output / max(abs(output)));
            if length(output) < minLength
                output = [output; zeros(minLength - length(output), 1)];
            else
                output = output(1:minLength);
            end
            
            % Calculate STOI
            current_stoi = stoi(output, x, fs);
            stoi_grid(i, j) = stoi_grid(i, j) + current_stoi;
        end
    end
end

% Average STOI values over all files
stoi_grid = stoi_grid / length(cleanFiles);

% Find the maximum STOI value and corresponding alpha and beta
[max_stoi, idx] = max(stoi_grid(:));
[max_i, max_j] = ind2sub(size(stoi_grid), idx);
best_alpha = alpha_range(max_i);
best_beta = beta_range(max_j);

% Display results
disp(['Best Alpha: ', num2str(best_alpha)]);
disp(['Best Beta: ', num2str(best_beta)]);
disp(['Maximum STOI: ', num2str(max_stoi)]);

% Plot the results
figure;
surf(alpha_range, beta_range, stoi_grid');
title('STOI Index across Alpha and Beta Values');
xlabel('Alpha');
ylabel('Beta');
zlabel('Average STOI');
colorbar;
