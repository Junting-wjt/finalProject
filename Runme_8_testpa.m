% -------------------------------------------------------------------- 
% Evaluate different Wiener filtering algorithms with given optimal parameters
% --------------------------------------------------------------------  
clear;
close all;
clc;
warning off;
addpath(genpath(pwd));

% Assume there are two folders containing clean and noisy speech files
cleanPath = 'audio\clean\clean';
%noisyPath = 'audio\airport_5dB\5dB';
%noisyPath = 'audio\car_5dB\5dB';

% Get file list
cleanFiles = dir(fullfile(cleanPath, '*.wav'));
%noisyFiles = dir(fullfile(noisyPath, '*.wav'));

%% Initialization
stoiValues = zeros(1, length(cleanFiles));
enhancedstoiValues = zeros(1, length(cleanFiles));
enhancedstoiValues2 = zeros(1, length(cleanFiles));
enhancedstoiValues3 = zeros(1, length(cleanFiles));

noisySpeechpesqValues = zeros(1, length(cleanFiles));
enhanceSpeechpesqValues = zeros(1, length(cleanFiles));
enhanceSpeechpesqValues2 = zeros(1, length(cleanFiles));
enhanceSpeechpesqValues3 = zeros(1, length(cleanFiles));

snr1Values = zeros(1, length(cleanFiles));           
snr2Values = zeros(1, length(cleanFiles));   
snr3Values = zeros(1, length(cleanFiles));     
snr4Values = zeros(1, length(cleanFiles)); 

%% Set parameters
IS = 0.25;                             % Set leading silence length
wlen = 200;                            % Set frame length to 25ms
inc = 80;                              % Set frame shift to 10ms
%SNR = 5;                               % Set signal-to-noise ratio (SNR)

alpha = 0.5;
beta = 0.2;

%%
% Process each file in the loop
for i = 1:length(cleanFiles)
    cleanFile = fullfile(cleanFiles(i).folder, cleanFiles(i).name);
    %noisyFile = fullfile(noisyFiles(i).folder, noisyFiles(i).name);
    
    % Read audio files
    [x, fs] = audioread(cleanFile);
    %[signal, ~] = awgn(x, SNR, 'measured', 'db');  % Add noise
    %[signal, ~] = audioread(noisyFile);
    audiowrite("clean_pesq.wav", x, fs);

    %%
    pink_noise = pink(length(x));  % Generate pink noise

    % Calculate the power of the original signal to set an appropriate noise level
    signalPower = sum(x.^2) / length(x);
    noisePower = signalPower / 10^(5/10);  % Assume desired SNR is 5dB
    %noisePower = signalPower / 10^(0/10);  % Assume desired SNR is 0dB

    % Adjust noise level
    adjustedPinkNoise = pink_noise * sqrt(noisePower / (sum(pink_noise.^2) / length(pink_noise)));

    % Add noise to signal
    signal = x + adjustedPinkNoise;

    %%
    audiowrite("noise_pesq.wav", signal, fs)

    %%
    % Ensure signal length consistency
    minLength = min([length(x), length(signal)]);
    x = x(1:minLength);
    signal = signal(1:minLength);

    NIS = fix((IS * fs - wlen) / inc + 1);  % Calculate the number of silent frames
    % Speech preprocessing
    x = x - mean(x);
    x = x / max(abs(x));
    signal = signal / max(abs(signal));

    %%
    % Apply Wiener filtering
    output = Weina_Norm_ENhuan(signal, wlen, inc, NIS);  
    output2 = Weina_Norm_OPP(signal, wlen, inc, NIS, alpha, beta); 
    output3 = Weina_Norm_OPPES(signal, wlen, inc, NIS, alpha, beta); 

    %%
    output = real(output / max(abs(output)));
    output2 = real(output2 / max(abs(output2)));
    output3 = real(output3 / max(abs(output3)));

    % Ensure signal length consistency for output
    if length(output) < minLength
        warning('Output length is shorter than expected. Check Weina_Norm function.');
        output = [output; zeros(minLength - length(output), 1)];  % Pad with zeros to minLength
    else
        output = output(1:minLength);
    end

    % Ensure signal length consistency for output2
    if length(output2) < minLength
        warning('Output length is shorter than expected. Check Weina_Norm function.');
        output2 = [output2; zeros(minLength - length(output2), 1)];  % Pad with zeros to minLength
    else
        output2 = output2(1:minLength);
    end

    % Ensure signal length consistency for output3
    if length(output3) < minLength
        warning('Output length is shorter than expected. Check Weina_Norm function.');
        output3 = [output3; zeros(minLength - length(output3), 1)];  % Pad with zeros to minLength
    else
        output3 = output3(1:minLength);
    end

    %%
    audiowrite("pesq_output.wav", output, fs)
    audiowrite("pesq_output2.wav", output2, fs);
    audiowrite("pesq_output3.wav", output3, fs);

    %%
    % Signal-to-noise ratio
    snr1Values(i) = SNR_Calc(x, signal);  % Calculate initial SNR
    snr2Values(i) = SNR_Calc(x, output);  % Calculate SNR after noise reduction
    snr3Values(i) = SNR_Calc(x, output2); % Calculate SNR after noise reduction
    snr4Values(i) = SNR_Calc(x, output3); % Calculate SNR after noise reduction

    %%
    % Calculate STOI
    stoiValues(i) = stoi(signal, x, fs);
    enhancedstoiValues(i) = stoi(output, x, fs);
    enhancedstoiValues2(i) = stoi(output2, x, fs);
    enhancedstoiValues3(i) = stoi(output3, x, fs);

    %%
    pesqResult = pesq('clean_pesq.wav', 'noise_pesq.wav');
    noisySpeechpesqValues(i) = pesqResult(1);

    pesqResult = pesq('clean_pesq.wav', 'pesq_output.wav');
    enhanceSpeechpesqValues(i) = pesqResult(1);

    pesqResult = pesq('clean_pesq.wav', 'pesq_output2.wav');
    enhanceSpeechpesqValues2(i) = pesqResult(1);

    pesqResult = pesq('clean_pesq.wav', 'pesq_output3.wav');
    enhanceSpeechpesqValues3(i) = pesqResult(1);
end

% Calculate average STOI values
meanSTOI = mean(stoiValues);
disp(['Average STOI index: ', num2str(meanSTOI)]);
meanSTOI = mean(enhancedstoiValues);
disp(['Average STOI index 2: ', num2str(meanSTOI)]);
meanSTOI = mean(enhancedstoiValues2);
disp(['Average STOI index 3: ', num2str(meanSTOI)]);
meanSTOI = mean(enhancedstoiValues3);
disp(['Average STOI index 4: ', num2str(meanSTOI)]);

%%
noisySpeechpesq = mean(noisySpeechpesqValues);
enhanceSpeechpesq = mean(enhanceSpeechpesqValues);
enhanceSpeechpesq2 = mean(enhanceSpeechpesqValues2);
enhanceSpeechpesq3 = mean(enhanceSpeechpesqValues3);
disp(['Noise PESQ: ', num2str(noisySpeechpesq)]);
disp(['Enhanced PESQ: ', num2str(enhanceSpeechpesq)]);
disp(['Enhanced PESQ 2: ', num2str(enhanceSpeechpesq2)]);
disp(['Enhanced PESQ 3: ', num2str(enhanceSpeechpesq3)]);

%%
% Signal-to-noise ratio
snr1 = mean(snr1Values);  % Calculate initial SNR
snr2 = mean(snr2Values);  % Calculate SNR after noise reduction
snr3 = mean(snr3Values);  % Calculate SNR after noise reduction
snr4 = mean(snr4Values);  % Calculate SNR after noise reduction
fprintf('snr1=%5.4f   snr2=%5.4f   snr3=%5.4f   snr4=%5.4f\n', snr1, snr2, snr3, snr4);
