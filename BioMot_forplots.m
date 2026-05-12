%%% Script to create amplitude spectra plots - run separately on 3M/6M data

close all
clear all

eeglab

%path = 'Y:\hoehl\projects\biomot\EEG_analysis\fixed_look_markers\3MO_clean\rereferenced'; % 3 months
path = 'path\to\EEGfiles'; % enter path to EEG files (3M or 6M)

files = dir(fullfile(path, '*.set'));

all_spectra = zeros(length(files), 626, 2); %last dimension: 1 - up, 2 -inv

for g = 1:length(files)
    file = files(g).name;
    partip_ID = file(1:6);
    EEG = pop_loadset(fullfile(path, file));
    EEG_roi = EEG.data(22:31, :, :); %22:31    

    % average signal depending on condition

    roi_channels = size(EEG_roi, 1);
    no_trials = size(EEG_roi, 3);

    signal_up = zeros(roi_channels, 834); %834 -- 1.668s in samples
    signal_inv = zeros(roi_channels, 834);


    up_count = 0;
    inv_count = 0;

    for k = 1:no_trials
        trial_type = EEG.epoch(k).eventtype{1};
        if strcmp(trial_type, 'S 15') || strcmp(trial_type, 'S 17')
            up_count = up_count + 1;
            for j = 1:roi_channels
                signal_up(j, :) = signal_up(j, :) + EEG_roi(j, :, k);
            end
        elseif strcmp(trial_type, 'S 16') || strcmp(trial_type, 'S 18')
            inv_count = inv_count + 1;
            for j = 1:roi_channels
                signal_inv(j, :) = signal_inv(j, :) + EEG_roi(j, :, k);
            end

        end
    end

    signal_up = signal_up./up_count;
    signal_inv = signal_inv./inv_count;

    signal_up_chanavg = mean(signal_up, 1);
    signal_inv_chanavg = mean(signal_inv,1);


    % compute FFT of the signal

    fs = 500;
    sig_length = 834;
    target_resolution = 0.4;

    % calculate required FFT length
    N = ceil(fs / target_resolution);

    zero_padding = N - sig_length;

    % apply zero-padding
    padded_signal_up = [signal_up_chanavg, zeros(1, zero_padding)];
    padded_signal_inv = [signal_inv_chanavg, zeros(1, zero_padding)];



    N = length(padded_signal_up);          % number of samples
    X = fft(padded_signal_up);             % FFT of the signal
    X2 = fft(padded_signal_inv); 
    % 
    % compute the power spectrum
    P2 = abs(X/N);          % two-sided spectrum
    P1 = P2(1:N/2+1);       % one-sided spectrum %UPRIGHT!
    P1(2:end-1) = 2*P1(2:end-1); % adjust amplitude for single-sided spectrum

    P4 = abs(X2/N);          % two-sided spectrum
    P3 = P4(1:N/2+1);       % one-sided spectrum %INVERTED!
    P3(2:end-1) = 2*P3(2:end-1); % adjust amplitude for single-sided spectrum

     all_spectra(g, :, 1) = P1;
     all_spectra(g, :, 2) = P3;
end

all_spect_mean = squeeze(mean(all_spectra, 1));

spect_mean_up = squeeze(all_spect_mean(:,1))';
spect_mean_inv = squeeze(all_spect_mean(:,2))';

fs = EEG.srate;
f = (0:(N/2))*fs/N;     % Frequency vector

% Plot the power spectrum
figure;
plot(f,spect_mean_up, 'g');
hold on
plot(f,spect_mean_inv, 'r');
title('Single-Sided Amplitude Spectrum of the Signal');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
grid on;
xlim([0 20])

% 
% calculate SNRs
% 
SNR_up = zeros(length(files), 626);
SNR_inv = zeros(length(files), 626);

for l = 1:length(files)
    for h = 4:624
        SNR_up(l,h) = all_spectra(l, h, 1) - (all_spectra(l, h-2, 1)+all_spectra(l, h+2, 1))./2;
        SNR_inv(l,h) = all_spectra(l,h, 2) - (all_spectra(l, h-2, 2)+all_spectra(l, h+2, 2))./2;
    end
end

avg_SNR_up = mean(SNR_up, 1);
avg_SNR_inv = mean(SNR_inv, 1);

figure;
plot(f,avg_SNR_up, 'Color', '#1A85FF', 'LineWidth', 3);
hold on
plot(f,avg_SNR_inv, 'Color', '#D41159', 'LineWidth', 3);
ax = gca;
ax.FontSize = 16; 
title(' ');
xlabel('Frequency [Hz]', 'fontsize',18);
ylabel('Amplitude [\muV]', 'fontsize',18);
grid on;
xlim([1.0 11.2])
ylim([-0.4 0.4])
lgd = legend('upright','inverted');
fontsize(lgd,18,'points')
set(gcf, 'PaperUnits', 'inches');
x_width=9.125 ;y_width=7.25;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,'6MO_paper.png')


% average across all conditions

all_SNR = (avg_SNR_up + avg_SNR_inv)./2;

figure;
plot(f, all_SNR, 'b')
title(' ');
xlabel('Frequency (Hz)');
ylabel('Baseline-Subtracted Amplitude');
grid on;
xlim([1 11])
%ylim([-0.25 0.3])
set(gca,'fontsize',14)