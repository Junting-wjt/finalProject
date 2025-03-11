% -------------------------------------------------------------------- 
% Mainly used to generate sample sound examples
% --------------------------------------------------------------------  
clc;
clear;
close all;
warning off;
addpath(genpath(pwd));

%[xx, fs] = audioread('lathe.wav');       % Read data file
[xx, fs] = audioread("audio\clean\clean\sp03.wav");

xx = xx - mean(xx);                      % Remove DC component
x = xx / max(abs(xx));                   % Amplitude normalization

%audiowrite("clean_pesq.wav", x, fs);
audiowrite("sample_clean.mp3", x, fs);
IS = 0.25;                               % Set leading silence length
wlen = 200;                              % Set frame length to 25ms
inc = 80;                                % Set frame shift to 10ms

%SNR = 5;                                 % Set signal-to-noise ratio (SNR)
NIS = fix((IS * fs - wlen) / inc + 1);   % Calculate the number of silent frames
alpha = 3;
beta = 0.01;

%%
%signal = awgn(x, SNR, 'measured', 'db');  % Add noise

%%
[signal, fs] = audioread("audio\car_0dB\0dB\sp03_car_sn0.wav");

%%
% pink_noise = pink(length(x));  % Generate pink noise
% 
% % Calculate the power of the original signal to set an appropriate noise level
% signalPower = sum(x.^2) / length(x);
% %noisePower = signalPower / 10^(0/10);  % Assume desired SNR is 0dB
% noisePower = signalPower / 10^(5/10);  % Assume desired SNR is 5dB
% 
% % Adjust noise level
% adjustedPinkNoise = pink_noise * sqrt(noisePower / (sum(pink_noise.^2) / length(pink_noise)));
% 
% % Add noise to signal
% signal = x + adjustedPinkNoise;

%%
%audiowrite("noise_pesq.wav", signal, fs)
audiowrite("sample_noisy.mp3", signal, fs)

%%
output = Weina_Norm_EN(signal, wlen, inc, NIS, alpha, beta);

%-------------------- Improvement
output2 = Weina_Norm(signal, wlen, inc, NIS, alpha, beta);
output3 = Weina_Norm(signal, wlen, inc, NIS, alpha, beta);
output4 = Weina_Norm(signal, wlen, inc, NIS, alpha, beta);

output = output / max(abs(output));
len = min(length(output), length(x));

% Three speech signals
x = x(1:len);
signal = signal(1:len);
output = output(1:len);
output2 = output2(1:len);
output3 = output3(1:len);
output4 = output4(1:len);

output = real(output);
output2 = real(output2);
output3 = real(output3);
output4 = real(output4);

%---------------------------------------------
%audiowrite("pesq_nolimiterout.wav", output, fs)

%----------------------------------------------------------------
% STOI
noisySpeechSTOI = stoi(signal, x, fs);
enhanceSpeech = stoi(output, x, fs);
enhanceSpeech2 = stoi(output2, x, fs);
enhanceSpeech3 = stoi(output3, x, fs);
enhanceSpeech4 = stoi(output4, x, fs);

sub = enhanceSpeech - noisySpeechSTOI;
disp(['noise STOI index:', num2str(noisySpeechSTOI)]);
disp(['enhance STOI index:', num2str(enhanceSpeech)]);
disp(['sub', num2str(sub)]);
disp(['enhance STOI index2:', num2str(enhanceSpeech2)]);
disp(['enhance STOI index3:', num2str(enhanceSpeech3)]);
disp(['enhance STOI index4:', num2str(enhanceSpeech4)]);

%% PESQ
% noisySpeechpesq = pesq('clean_pesq.wav', 'noise_pesq.wav');
% 
% enhanceSpeechpesq = pesq('clean_pesq.wav', 'pesq_nolimiterout.wav');
% 
% %sub = enhanceSpeechpesq - noisySpeechpesq;
% disp(['noise PESQ:', num2str(noisySpeechpesq)]);
% disp(['enhance PESQ:', num2str(enhanceSpeechpesq)]);
% %disp(['sub', num2str(sub)]);

%%
% Signal-to-noise ratio
snr1 = SNR_Calc(x, signal);             % Calculate initial SNR
snr2 = SNR_Calc(x, output);             % Calculate SNR after noise reduction
snr = snr2 - snr1;
fprintf('snr1=%5.4f   snr2=%5.4f   snr=%5.4f\n', snr1, snr2, snr);

%%
% Segmental SNR
segSnr1 = seg_SNR(x, signal, 200);      % Calculate initial segmental SNR
segSnr2 = seg_SNR(x, output, 200);      % Calculate segmental SNR after noise reduction
segSnr = segSnr2 - segSnr1;
fprintf('segSnr1=%5.4f   segSnr2=%5.4f   segSnr=%5.4f\n', segSnr1, segSnr2, segSnr);

%%
% Plotting
time = (0:len-1) / fs;                  % Set time

% subplot 311; plot(time, x, 'k'); grid; axis tight;
% title('original speech signal'); ylabel('Amplitude')
% subplot 312; plot(time, signal, 'k'); grid; axis tight;
% title(['noisy speech  signal-to-noise ratio =' num2str(SNR) 'dB']); ylabel('Amplitude')
% subplot 313; plot(time, output, 'k'); grid;
% title('filtered signal'); ylabel('Amplitude'); xlabel('time/s');

%%
%x, signal, and output are your three audio signals, fs is the sampling rate

figure;

% Spectrum of the first signal
subplot(3, 1, 1);
spectrogram(x, 200, 80, 256, fs, 'yaxis');
colorbar('off');
title('Original Speech Spectrogram');

% Spectrum of the second signal
subplot(3, 1, 2);
spectrogram(signal, 200, 80, 256, fs, 'yaxis');
colorbar('off');
title('Noisy Speech Spectrogram');

% Spectrum of the third signal
subplot(3, 1, 3);
spectrogram(output, 200, 80, 256, fs, 'yaxis');
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

%%
if length(output) > length(x)
    output = output(1:length(x));       % Trim output to the same length as x
end
sound(output2, fs);

audiowrite("sample_output.mp3", output, fs);
