%% Neurophotonics HW1 - Run Script
% Executes CalcNIRS on both provided recordings, generates hemoglobin
% concentration plots for channels 1-2, and performs spectral analysis
% (FFT + SNR) on channel 1 of the first file.

close all; clc;

% Resolve paths relative to this script's location
scriptFolder = fileparts(mfilename('fullpath'));
cd(scriptFolder);

%% Measurement parameters
sourceDetectorSep = 3;             % [cm]
tissue = 'adult_head';
channelsToPlot = [1 2];

%% Process File 1
[dHbR1, dHbO1, fig1] = CalcNIRS('FN_031_V2_Postdose2_Nback.mat', ...
    sourceDetectorSep, tissue, channelsToPlot);
saveas(fig1, 'FN_031_V2_Postdose2_Nback_HbR_HbO.png');
fprintf('File 1: dHbR [%d x %d], dHbO [%d x %d]\n', ...
    size(dHbR1,1), size(dHbR1,2), size(dHbO1,1), size(dHbO1,2));

%% Process File 2
[dHbR2, dHbO2, fig2] = CalcNIRS('FN_032_V1_Postdose1_Nback.mat', ...
    sourceDetectorSep, tissue, channelsToPlot);
saveas(fig2, 'FN_032_V1_Postdose1_Nback_HbR_HbO.png');
fprintf('File 2: dHbR [%d x %d], dHbO [%d x %d]\n', ...
    size(dHbR2,1), size(dHbR2,2), size(dHbO2,1), size(dHbO2,2));

%% FFT analysis - Channel 1, File 1
% The cardiac pulse modulates arterial blood volume, producing a periodic
% component in the fNIRS signal. I use the FFT to identify this frequency
% and evaluate signal quality via the SNR.

fileData = load('FN_031_V2_Postdose2_Nback.mat');
timeVec = double(fileData.t(:));
dt = timeVec(2) - timeVec(1);
fs = 1 / dt;                          % ~7.81 Hz
nSamples = length(timeVec);

ch1Signal = dHbO1(:, 1);
fftVals = fft(ch1Signal);

% Single-sided magnitude spectrum (normalized)
nHalf = floor(nSamples/2) + 1;
freqAxis = (0:nHalf-1) * fs / nSamples;
magSpectrum = abs(fftVals(1:nHalf)) / nSamples;

%% SNR computation
% "Noise" = mean spectral magnitude above 2.5 Hz (beyond hemodynamic and
% cardiac bands; dominated by instrumentation noise at these frequencies).
% "Signal" = peak magnitude in the cardiac frequency range.

NOISE_CUTOFF = 2.5;                   % [Hz]
noiseRegion = freqAxis > NOISE_CUTOFF;
avgNoiseMag = mean(magSpectrum(noiseRegion));

HR_BAND = [0.5  2.5];                 % expected cardiac range [Hz]
hrMask = freqAxis >= HR_BAND(1) & freqAxis <= HR_BAND(2);
hrMagnitudes = magSpectrum(hrMask);
hrFreqs = freqAxis(hrMask);
[peakMag, peakIdx] = max(hrMagnitudes);
heartbeatFreq = hrFreqs(peakIdx);

snrValue = peakMag / avgNoiseMag;

fprintf('\n--- Spectral Analysis (Channel 1, File 1) ---\n');
fprintf('Heartbeat frequency : %.3f Hz  (%.1f BPM)\n', ...
    heartbeatFreq, heartbeatFreq * 60);
fprintf('SNR                 : %.2f\n', snrValue);

%% FFT plot
fftFig = figure;
plot(freqAxis, magSpectrum, 'k', 'LineWidth', 0.8);
hold on;
xline(heartbeatFreq, 'r--', ...
    sprintf('Heartbeat = %.2f Hz', heartbeatFreq), 'LineWidth', 1.2);
xline(NOISE_CUTOFF, ':', 'Noise cutoff (2.5 Hz)', ...
    'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title(sprintf('FFT of \\DeltaHbO - Ch 1 - FN\\_031\nSNR = %.2f', snrValue));
legend('FFT magnitude', sprintf('Heartbeat = %.2f Hz', heartbeatFreq), ...
    'Noise cutoff');
xlim([0 5]);
grid on;
saveas(fftFig, 'FFT_Channel1_File1.png');
fprintf('All plots saved.\n');
