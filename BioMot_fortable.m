% Script to create data structures - run separately on 3M/6M data
close all
clear all

eeglab

path = 'path\to\EEGfiles'; % enter path to EEG files (3M or 6M)

files = dir(fullfile(path, '*.set')); %access files

results = table( ... %create results table
    strings(0,1), ...  % participant_ID (string)
    strings(0,1), ...  % orientation (string)
    strings(0,1), ...  % direction (string)
    strings(0,1), ...  % harmonic (string)
    strings(0,1), ...  % hemisphere (string)
    zeros(0,1),   ...  % no_epochs (numeric)
    zeros(0,1),   ...  % SNR (numeric)
    'VariableNames', {'participant_ID', 'orientation', 'direction', ...
                      'harmonic', 'hemisphere', 'no_epochs', 'SNR'});

row_counter = 2;
for g = 1:length(files)
    file = files(g).name; 
    partip_ID = file(1:6);
    EEG = pop_loadset(fullfile(path, file)); % load the preprocessed EEG set
    EEG_roi_right = EEG.data([25, 26, 30, 31], :, :); % right side parietal-occipital
    EEG_roi_left = EEG.data([22, 23, 27, 28], :, :); % left side parietal-occipital
    
    roi_channels = size(EEG_roi_right, 1); % number of channels in ROI, is the same for left
    no_trials = size(EEG_roi_right, 3); % number of epochs

    signal_up_left_lhemi = zeros(roi_channels, 834); % create a zeros matrix for signals in the upright/left condition, left hemisphere ROI; below analogous
    signal_up_left_rhemi = zeros(roi_channels, 834);
    signal_up_right_lhemi = zeros(roi_channels, 834);
    signal_up_right_rhemi = zeros(roi_channels, 834);
    signal_inv_left_lhemi = zeros(roi_channels, 834);
    signal_inv_left_rhemi = zeros(roi_channels, 834);
    signal_inv_right_lhemi = zeros(roi_channels, 834);
    signal_inv_right_rhemi = zeros(roi_channels, 834);

    up_left_count = 0; % count how many epochs in the upright-left condition
    up_right_count = 0; % analogous
    inv_left_count = 0;
    inv_right_count = 0;

    for k = 1:no_trials
        trial_type = EEG.epoch(k).eventtype{1}; % extract event marker
        if strcmp(trial_type, 'S 15') % upright & right condition
            up_right_count = up_right_count + 1;
            for j = 1:roi_channels
                signal_up_right_lhemi(j, :) = signal_up_right_lhemi(j, :) + EEG_roi_left(j, :, k); %add signal to the matrix
                signal_up_right_rhemi(j, :) = signal_up_right_rhemi(j, :) + EEG_roi_right(j, :, k);
            end 
        elseif strcmp(trial_type, 'S 16') % inverted & right condition
            inv_right_count = inv_right_count + 1;
            for j = 1:roi_channels
                signal_inv_right_lhemi(j, :) = signal_inv_right_lhemi(j, :) + EEG_roi_left(j, :, k); %add signal to the matrix
                signal_inv_right_rhemi(j, :) = signal_inv_right_rhemi(j, :) + EEG_roi_right(j, :, k);
            end
        elseif strcmp(trial_type, 'S 17') % upright & left condition
            up_left_count = up_left_count + 1;
            for j = 1:roi_channels
                signal_up_left_lhemi(j, :) = signal_up_left_lhemi(j, :) + EEG_roi_left(j, :, k); %add signal to the matrix
                signal_up_left_rhemi(j, :) = signal_up_left_rhemi(j, :) + EEG_roi_right(j, :, k);
            end
        elseif strcmp(trial_type, 'S 18') % inverted & left condition 
            inv_left_count = inv_left_count + 1;
            for j = 1:roi_channels
                signal_inv_left_lhemi(j, :) = signal_inv_left_lhemi(j, :) + EEG_roi_left(j, :, k);
                signal_inv_left_rhemi(j, :) = signal_inv_left_rhemi(j, :) + EEG_roi_right(j, :, k);
            end

        end
    end
    
    % average the signals across epochs
    signal_up_left_lhemi = signal_up_left_lhemi./up_left_count;
    signal_up_left_rhemi = signal_up_left_rhemi./up_left_count;
    signal_up_right_lhemi = signal_up_right_lhemi./up_right_count;
    signal_up_right_rhemi = signal_up_right_rhemi./up_right_count;
    signal_inv_left_lhemi = signal_inv_left_lhemi./inv_left_count;
    signal_inv_left_rhemi = signal_inv_left_rhemi./inv_left_count;
    signal_inv_right_lhemi = signal_inv_right_lhemi./inv_right_count;
    signal_inv_right_rhemi = signal_inv_right_rhemi./inv_right_count;
    
    % average the signals across channels
    signal_up_left_lhemi_chanavg = mean(signal_up_left_lhemi, 1);
    signal_up_left_rhemi_chanavg = mean(signal_up_left_rhemi, 1);
    signal_up_right_lhemi_chanavg = mean(signal_up_right_lhemi, 1);
    signal_up_right_rhemi_chanavg = mean(signal_up_right_rhemi, 1);
    signal_inv_left_lhemi_chanavg = mean(signal_inv_left_lhemi,1);
    signal_inv_left_rhemi_chanavg = mean(signal_inv_left_rhemi,1);
    signal_inv_right_lhemi_chanavg = mean(signal_inv_right_lhemi,1);
    signal_inv_right_rhemi_chanavg = mean(signal_inv_right_rhemi,1);

    % zero-pad the signal
    fs = 500.;
    sig_length = 834;
    target_resolution = 0.4;

    N = ceil(fs/target_resolution); % required FFT length
    zero_padding = N - sig_length;

    padded_up_left_lhemi = [signal_up_left_lhemi_chanavg, zeros(1, zero_padding)];
    padded_up_left_rhemi = [signal_up_left_rhemi_chanavg, zeros(1, zero_padding)]; 
    padded_up_right_lhemi = [signal_up_right_lhemi_chanavg, zeros(1, zero_padding)]; 
    padded_up_right_rhemi = [signal_up_right_rhemi_chanavg, zeros(1, zero_padding)]; 
    padded_inv_left_lhemi = [signal_inv_left_lhemi_chanavg, zeros(1, zero_padding)]; 
    padded_inv_left_rhemi = [signal_inv_left_rhemi_chanavg, zeros(1, zero_padding)];
    padded_inv_right_lhemi = [signal_inv_right_lhemi_chanavg, zeros(1, zero_padding)];
    padded_inv_right_rhemi = [signal_inv_right_rhemi_chanavg, zeros(1, zero_padding)];



    % compute FFTs of the averaged signals
    N = length(padded_up_left_lhemi); % will be the same for all conditions
    Xul_left = fft(padded_up_left_lhemi); % upright left, left hemisphere, etc.
    Xul_right = fft(padded_up_left_rhemi);
    Xur_left = fft(padded_up_right_lhemi);
    Xur_right = fft(padded_up_right_rhemi);
    Xil_left = fft(padded_inv_left_lhemi);
    Xil_right = fft(padded_inv_left_rhemi);
    Xir_left = fft(padded_inv_right_lhemi);
    Xir_right = fft(padded_inv_right_rhemi);

    % compute the power spectrum
    Pp_ul_left = abs(Xul_left/N);          % two-sided spectrum %upright left, left hemisphere, etc.
    P_ul_left = Pp_ul_left(1:N/2+1);       % one-sided spectrum
    P_ul_left(2:end-1) = 2*P_ul_left(2:end-1); % adjust amplitude for single-sided spectrum

    Pp_ul_right = abs(Xul_right/N);          % two-sided spectrum
    P_ul_right = Pp_ul_right(1:N/2+1);       % one-sided spectrum
    P_ul_right(2:end-1) = 2*P_ul_right(2:end-1); % adjust amplitude for single-sided spectrum

    Pp_ur_left = abs(Xur_left/N);          % two-sided spectrum
    P_ur_left = Pp_ur_left(1:N/2+1);       % one-sided spectrum
    P_ur_left(2:end-1) = 2*P_ur_left(2:end-1); % adjust amplitude for single-sided spectrum

    Pp_ur_right = abs(Xur_right/N);          % two-sided spectrum
    P_ur_right = Pp_ur_right(1:N/2+1);       % one-sided spectrum
    P_ur_right(2:end-1) = 2*P_ur_right(2:end-1); % adjust amplitude for single-sided spectrum

    Pp_il_left = abs(Xil_left/N);          % two-sided spectrum
    P_il_left = Pp_il_left(1:N/2+1);       % one-sided spectrum
    P_il_left(2:end-1) = 2*P_il_left(2:end-1); % adjust amplitude for single-sided spectrum

    Pp_il_right = abs(Xil_right/N);          % two-sided spectrum
    P_il_right = Pp_il_right(1:N/2+1);       % one-sided spectrum
    P_il_right(2:end-1) = 2*P_il_right(2:end-1); % adjust amplitude for single-sided spectrum

    Pp_ir_left = abs(Xir_left/N);          % two-sided spectrum
    P_ir_left = Pp_ir_left(1:N/2+1);       % one-sided spectrum
    P_ir_left(2:end-1) = 2*P_ir_left(2:end-1); % adjust amplitude for single-sided spectrum

    Pp_ir_right = abs(Xir_right/N);          % two-sided spectrum
    P_ir_right = Pp_ir_right(1:N/2+1);       % one-sided spectrum
    P_ir_right(2:end-1) = 2*P_ir_right(2:end-1); % adjust amplitude for single-sided spectrum
    spectra{g, 1} = partip_ID; % create a cell array where you store all 8 spectra per participant
    spectra{g, 2} = P_ur_right;
    spectra{g, 3} = P_ur_left;
    spectra{g, 4} = P_ul_right;
    spectra{g, 5} = P_ul_left;
    spectra{g, 6} = P_ir_right;
    spectra{g, 7} = P_ir_left;
    spectra{g, 8} = P_il_right;
    spectra{g, 9} = P_il_left;

    % compute SNRs (baseline-subtracted amplitudes)

    SNR_ur_right = zeros(1, 626);
    SNR_ur_left = zeros(1, 626);
    SNR_ul_right = zeros(1, 626);
    SNR_ul_left = zeros(1, 626);
    SNR_ir_right = zeros(1, 626);
    SNR_ir_left = zeros(1, 626);
    SNR_il_right = zeros(1, 626);
    SNR_il_left = zeros(1, 626);

    for h = 4:624
        SNR_ur_right(h) = P_ur_right(h) - (P_ur_right(h-2) + P_ur_right(h+2))./2;
        SNR_ur_left(h) = P_ur_left(h) - (P_ur_left(h-2) + P_ur_left(h+2))./2;
        SNR_ul_right(h) = P_ul_right(h) - (P_ul_right(h-2) + P_ul_right(h+2))./2;
        SNR_ul_left(h) = P_ul_left(h) - (P_ul_left(h-2) + P_ul_left(h+2))./2;
        SNR_ir_right(h) = P_ir_right(h) - (P_ir_right(h-2) + P_ir_right(h+2))./2;
        SNR_ir_left(h) = P_ir_left(h) - (P_ir_left(h-2) + P_ir_left(h+2))./2;
        SNR_il_right(h) = P_il_right(h) - (P_il_right(h-2) + P_il_right(h+2))./2;
        SNR_il_left(h) = P_il_left(h) - (P_il_left(h-2) + P_il_left(h+2))./2;
    end

    SNRS{g, 1} = partip_ID;
    SNRS{g, 2} = SNR_ur_right;
    SNRS{g, 3} = SNR_ur_left;
    SNRS{g, 4} = SNR_ul_right;
    SNRS{g, 5} = SNR_ul_left;
    SNRS{g, 6} = SNR_ir_right;
    SNRS{g, 7} = SNR_ir_left;
    SNRS{g, 8} = SNR_il_right;
    SNRS{g, 9} = SNR_il_left;

    % up & right condition
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "2_4Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = up_right_count;
    results{row_counter, 7} = SNR_ur_right(7);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "4_8Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = up_right_count;
    results{row_counter, 7} = SNR_ur_right(13);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "7_2Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = up_right_count;
    results{row_counter, 7} = SNR_ur_right(19);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "2_4Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = up_right_count;
    results{row_counter, 7} = SNR_ur_left(7);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "4_8Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = up_right_count;
    results{row_counter, 7} = SNR_ur_left(13);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "7_2Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = up_right_count;
    results{row_counter, 7} = SNR_ur_left(19);
    row_counter = row_counter + 1;

    %up & left condition
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "2_4Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = up_left_count;
    results{row_counter, 7} = SNR_ul_right(7);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "4_8Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = up_left_count;
    results{row_counter, 7} = SNR_ul_right(13);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "7_2Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = up_left_count;
    results{row_counter, 7} = SNR_ul_right(19);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "2_4Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = up_left_count;
    results{row_counter, 7} = SNR_ul_left(7);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "4_8Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = up_left_count;
    results{row_counter, 7} = SNR_ul_left(13);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "up";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "7_2Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = up_left_count;
    results{row_counter, 7} = SNR_ul_left(19);
    row_counter = row_counter + 1;

    % inv&right condition
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "2_4Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = inv_right_count;
    results{row_counter, 7} = SNR_ir_right(7);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "4_8Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = inv_right_count;
    results{row_counter, 7} = SNR_ir_right(13);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "7_2Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = inv_right_count;
    results{row_counter, 7} = SNR_ir_right(19);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "2_4Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = inv_right_count;
    results{row_counter, 7} = SNR_ir_left(7);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "4_8Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = inv_right_count;
    results{row_counter, 7} = SNR_ir_left(13);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "right";
    results{row_counter, 4} = "7_2Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = inv_right_count;
    results{row_counter, 7} = SNR_ir_left(19);
    row_counter = row_counter + 1;

    % inverted & left
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "2_4Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = inv_left_count;
    results{row_counter, 7} = SNR_il_right(7);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "4_8Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = inv_left_count;
    results{row_counter, 7} = SNR_il_right(13);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "7_2Hz";
    results{row_counter, 5} = "right";
    results{row_counter, 6} = inv_left_count;
    results{row_counter, 7} = SNR_il_right(19);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "2_4Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = inv_left_count;
    results{row_counter, 7} = SNR_il_left(7);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "4_8Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = inv_left_count;
    results{row_counter, 7} = SNR_il_left(13);
    row_counter = row_counter + 1;
    results{row_counter, 1} = string(partip_ID);
    results{row_counter, 2} = "inverted";
    results{row_counter, 3} = "left";
    results{row_counter, 4} = "7_2Hz";
    results{row_counter, 5} = "left";
    results{row_counter, 6} = inv_left_count;
    results{row_counter, 7} = SNR_il_left(19);
    row_counter = row_counter + 1;

end


writetable(results, '6M_results.csv'); %change file name depending on age group

%save('SNRs_6MO.csv','SNRS' ); %change file name depending on age group





      


