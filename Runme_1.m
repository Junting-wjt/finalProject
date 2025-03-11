% -------------------------------------------------------------------- 
% Basic Test, this code is modified with reference to the content of 
% the following website:
% https://blog.csdn.net/qq_42233059/article/details/
% 126502456?ops_request_misc=%257B%2522request%255Fid%
% 2522%253A%2522172445642616800172515220%2522%252C%2522scm%
% 2522%253A%252220140713.130102334.pc%255Fblog.%2522%257D
% &request_id=172445642616800172515220&biz_id=0&utm_medium
% =distribute.pc_search_result.none-task-blog-2~blog~
% first_rank_ecpm_v1~rank_v31_ecpm-1-126502456-null-null.
% nonecase&utm_term=%E7%BB%B4%E7%BA%B3%E6%BB%A4%E6%B3%A2&spm=
% 1018.2226.3001.4450
% -------------------------------------------------------------------- 
clc;
clear;
close all;
warning off;
addpath(genpath(pwd));

%[xx, fs] = audioread('dat.wav');           % Read data file
%[xx, fs] = audioread('lathe.wav');         % Read data file
[xx, fs] = audioread("audio\clean\clean\sp03.wav");

%[xx, fs] = audioread('handel_echo.wav'); 
%[xx, fs] = audioread('output_LMS.wav');   % Read data file
%xx = xx;
xx = xx - mean(xx);                         % Remove DC component
x = xx / max(abs(xx));                      % Amplitude normalization

audiowrite("clean_pesq.wav", x, fs);
IS = 0.25;                                  % Set leading silence length
wlen = 200;                                 % Set frame length to 25ms
inc = 80;                                   % Set frame shift to 10ms
SNR = 0;                                    % Set signal-to-noise ratio (SNR)
NIS = fix((IS * fs - wlen) / inc + 1);      % Calculate number of silent frames
alpha = 3;
beta = 0.01;

%%
signal = awgn(x, SNR, 'measured', 'db');    % Add noise
%audiowrite("noise_pesq.wav", signal, fs)

%%
%[signal, fs] = audioread("audio\car_0dB\0dB\sp03_car_sn5.wav");

%%
% pink_noise = pink(length(x));  % Generate pink noise
% 
% % Calculate the power of the original signal to set an appropriate noise level
% signalPower = sum(x.^2) / length(x);
% %noisePower = signalPower / 10^(5/10);  % Assume desired SNR is 5dB
% noisePower = signalPower / 10^(0/10);  % Assume desired SNR is 0dB
% 
% % Adjust the noise level
% adjustedPinkNoise = pink_noise * sqrt(noisePower / (sum(pink_noise.^2) / length(pink_noise)));
% 
% % Add noise to signal
% signal = x + adjustedPinkNoise;

%%
audiowrite("noise_pesq.wav", signal, fs)

%%
% Add pink noise to the signal
%signal = x + pink_noise;

%======== Limiter test
%noise = randn(size(x));
%amplifiedNoise = noise;
%signal = x + amplifiedNoise;
%-------------------------------------

%%
output = Weina_Norm(signal, wlen, inc, NIS, alpha, beta);

%-------------------- Improvement
output2 = Weina_Norm(signal, wlen, inc, NIS, alpha, beta);
output3 = Weina_Norm(signal, wlen, inc, NIS, alpha, beta);
output4 = Weina_Norm(signal, wlen, inc, NIS, alpha, beta);
output5 = Weina_Norm(signal, wlen, inc, NIS, alpha, beta);

output = output / max(abs(output));
len = min(length(output), length(x));

% Three speech signals
x = x(1:len);
signal = signal(1:len);
output = output(1:len);
output2 = output2(1:len);
output3 = output3(1:len);
output4 = output4(1:len);
output5 = output5(1:len);

output = real(output);
output2 = real(output2);
output3 = real(output3);
output4 = real(output4);
output5 = real(output5);

%---------------------------------------------
audiowrite("pesq_nolimiterout.wav", output, fs)

% % Create a compressor object and set parameters
% compressorObj = compressor('Threshold', -10, 'Ratio', 4, 'KneeWidth', 1);
% 
% % Create a limiter object and set parameters (suitable for Lookahead Limiter)
% limiterObj = limiter('Threshold', -1, 'AttackTime', 0.01, 'ReleaseTime', 0.1);
% 
% % Apply compressor and limiter
% audioCompressed = compressorObj(output);
% audioLimited = limiterObj(audioCompressed);
% 
% % Write audio file
% audiowrite('output_toolbox_limiter.wav', audioLimited, fs);

%------------------------------------------------------
% Calculate max and min levels
%maxLevel = max(abs(audioLimited));
%minLevel = min(abs(audioLimited));

% Dynamic range
%dynamicRange = 20 * log10(maxLevel / minLevel);
%disp(['Dynamic Range: ', num2str(dynamicRange), ' dB']);

% Calculate max and min levels
%maxLevel = max(abs(output));
%minLevel = min(abs(output));

% Dynamic range
%dynamicRange2 = 20 * log10(maxLevel / minLevel);
%disp(['Dynamic Range 2: ', num2str(dynamicRange2), ' dB']);

%----------------------------------------------------------------
% STOI
noisySpeechSTOI = stoi(signal, x, fs);
enhanceSpeech = stoi(output, x, fs);
enhanceSpeech2 = stoi(output2, x, fs);
enhanceSpeech3 = stoi(output3, x, fs);
enhanceSpeech4 = stoi(output4, x, fs);
enhanceSpeech5 = stoi(output5, x, fs);

%audioout = stoi(audioLimited, x, fs);

sub = enhanceSpeech - noisySpeechSTOI;
disp(['noise short-time objective intelligibility index:', num2str(noisySpeechSTOI)]);
disp(['enhance short-time objective intelligibility index:', num2str(enhanceSpeech)]);
disp(['sub', num2str(sub)]);
%disp(['audioout', num2str(audioout)]);
disp(['enhance short-time objective intelligibility index2:', num2str(enhanceSpeech2)]);
disp(['enhance short-time objective intelligibility index3:', num2str(enhanceSpeech3)]);
disp(['enhance short-time objective intelligibility index4:', num2str(enhanceSpeech4)]);
disp(['enhance short-time objective intelligibility index5:', num2str(enhanceSpeech5)]);

% PESQ
noisySpeechpesq = pesq('clean_pesq.wav', 'noise_pesq.wav');
enhanceSpeechpesq = pesq('clean_pesq.wav', 'pesq_nolimiterout.wav');

%audiooutpesq = stoi(audioLimited, x, fs);

%sub = enhanceSpeechpesq - noisySpeechpesq;
disp(['noise PESQ:', num2str(noisySpeechpesq)]);
disp(['enhance PESQ:', num2str(enhanceSpeechpesq)]);
%disp(['sub', num2str(sub)]);
%disp(['audioout', num2str(audiooutpesq)]);

% Signal-to-noise ratio
snr1 = SNR_Calc(x, signal);           % Calculate initial SNR
snr2 = SNR_Calc(x, output);           % Calculate SNR after noise reduction
snr = snr2 - snr1;
fprintf('snr1=%5.4f   snr2=%5.4f   snr=%5.4f\n', snr1, snr2, snr);

% Segmental SNR
segSnr1 = seg_SNR(x, signal, 200);    % Calculate initial segmental SNR
segSnr2 = seg_SNR(x, output, 200);    % Calculate segmental SNR after noise reduction
segSnr = segSnr2 - segSnr1;
fprintf('segSnr1=%5.4f   segSnr2=%5.4f   segSnr=%5.4f\n', segSnr1, segSnr2, segSnr);

% Plotting
figure;
time = (0:len-1) / fs;                % Set time
% subplot 311; plot(time,x,'k'); grid; axis tight;
% title('Pure Speech Waveform'); ylabel('Amplitude')
% subplot 312; plot(time,signal,'k'); grid; axis tight;
% title(['Noisy Speech SNR=' num2str(SNR) 'dB']); ylabel('Amplitude')
% subplot 313; plot(time,output,'k'); grid;
% title('Filtered Waveform'); ylabel('Amplitude'); xlabel('Time/s');

subplot 311; plot(time, x, 'k'); grid; axis tight;
title('Original Speech Signal'); ylabel('Amplitude')
subplot 312; plot(time, signal, 'k'); grid; axis tight;
title(['Noisy Speech Signal-to-Noise Ratio =' num2str(SNR) 'dB']); ylabel('Amplitude')
subplot 313; plot(time, output, 'k'); grid;
title('Filtered Signal'); ylabel('Amplitude'); xlabel('Time/s');

% x, signal, and output are your three audio signals, Fs is the sampling rate

figure;

% Spectrum of the first signal
subplot(3, 1, 1);
%spectrogram(x, 256, 250, 256, fs, 'yaxis');
spectrogram(x, 200, 80, 256, fs, 'yaxis');
xlabel('Time');
ylabel('Magnitude');  % Set Y-axis label to English
colorbar('off');
title('Original Speech Spectrogram');

% Spectrum of the second signal
subplot(3, 1, 2);
spectrogram(signal, 200, 80, 256, fs, 'yaxis');
xlabel('Time');
ylabel('Magnitude');  % Set Y-axis label to English
colorbar('off');
title('Noisy Speech Spectrogram');

% Spectrum of the third signal
subplot(3, 1, 3);
spectrogram(output, 200, 80, 256, fs, 'yaxis');
xlabel('Time');
ylabel('Magnitude');  % Set Y-axis label to English
colorbar('off');
title('Enhanced Speech Spectrogram');

% Set color limits for the entire plot to maintain consistency
clim = get(gca, 'CLim');  % Get current axis color limits
subplot(3, 1, 1);
set(gca, 'CLim', clim);
subplot(3, 1, 2);
set(gca, 'CLim', clim);

% Set colormap
colormap jet;

if length(output) > length(x)
    output = output(1:length(x)); % Trim output to the same length as x
end
%output = real(output);
sound(output, fs);

%audiowrite("est_clean.wav", output, fs)
