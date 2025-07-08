% script to calculate the power modulation indices, absolute, and relative
% power for a given age group and frequency band

clear

% load segment list: table with indices of clean useable segments, different for T1
% (6 months) and T2 (12 months)
load('REPLACE\WITH\FILEPATH\kept_indices_TX.mat');

%load the existing results.csv table
results = readtable('REPLACE\WITH\FILEPATH\\results.csv'); 


%% EEG

% Specify the folder path containing the datasets: separate for T1 and T2
% data
folderPath = 'REPLACE\WITH\FILEPATH\';

% Get a list of all .set files in the folder
fileList = dir(fullfile(folderPath, '*.set'));

% variable vor rowIndex
allRows = [];

%define ROI
fronto_central_channels = [1, 2, 6, 7, 27, 28]; % F3, Fz, F4, FC3, FCz, FC4

eeglab; % Start EEGLAB

% Loop through each dataset
for i = 1:length(fileList)
    % Load dataset in EEGLAB
    
    EEG = pop_loadset('filename', fileList(i).name, 'filepath', folderPath);
    eeglab redraw;

    % Extract the dataset name
    datasetName = EEG.filename;

    % Find the specific row in the CSV dataset
    % Convert the first column to a cell array of strings
    csvDataCol1 = table2cell(results(:, 1));
    csvDataCol1 = cellstr(csvDataCol1);
    % Find row
    rowIndex = find(strncmpi(csvDataCol1, datasetName, 5), 1);
    allRows = [allRows, rowIndex];

    % show number of epochs in the dataset
    num_epochs = EEG.trials;
    % save to dataset
    colIndex_seg = 5;
    % Assign the value to cell
    results{rowIndex, colIndex_seg} = num_epochs;

    disp(['Number of epochs in the dataset: ' num2str(num_epochs)]);

    % Find the coresponding row in the kept_indices cell
    % Convert the first column to a cell array of string
    kept_indicesCol1 = cellstr(kept_indices2(: , 1)); 
    % Find row
    indicesRow = find(strncmpi(kept_indicesCol1, datasetName, 5), 1);
    % get segment numbers from kept_indices
    segment_num = kept_indices2{indicesRow, 2}; 

    EEG_cont = eeg_epoch2continuous(EEG); %concatenate the 1s non-overlapping segments
    %removing boundaries to disable discontinuity rejection
    EEG_cont.event(strcmp({EEG_cont.event.type}, 'boundary')) = [];
    EEG_cont = eeg_checkset(EEG_cont, 'eventconsistency'); 

    %insert events for epoching every 0.5 sec (for creating 1-second
    %segments with 50% overlap)
    EEG_cont = eeg_regepochs(EEG_cont, 'recurrence', 0.5, 'eventtype','3','extractepochs','off');


    %epoch according to the new events
    EEG_over = pop_epoch(EEG_cont, {'3'}, [0, 1]);
    eeglab redraw

    % Select fronto-central electrode channels
    EEG = pop_select(EEG, 'channel', fronto_central_channels);
    EEG_over = pop_select(EEG_over, 'channel', fronto_central_channels);


    %% FFT
    npts = EEG.pnts;

    fs = EEG.srate;
    segment_length = fs;

    %Define parameters for pwelch

    noverlap = 0; %0% overlap
    nfft = segment_length; % FFT length

    num_segments = size(EEG.data, 3);
    num_segmentsover = size(EEG_over.data, 3);

    % Define the frequency range of interest
    freq_range = [6 9];  % Frequency range in Hz, change for different frequency bands! [3 6] for theta, [6 9] for alpha
    % Find the corresponding frequency indices
    freq_idx = find(EEG.srate/npts * (0:npts/2) >= freq_range(1) & EEG.srate/npts * (0:npts/2) <= freq_range(2));

    % Initialize average power per segment matrix
    avg_power = zeros(length(EEG.epoch), 1);

    %Initialize matrix to store pwelch psd for each channel
    power_pwelch = zeros(6, nfft/2+1, num_segmentsover); %channels x psd points x epochs


    % Loop through each trial
    for trial = 1:length(EEG.epoch)
        % Extract data for the current trial and electrode channels
        data = EEG.data(:, :, trial);

        % Compute FFT for each channel
        fft_data = fft(data, [], 2);

        % Compute amplitude power spectrum
        power_spectrum = (abs((fft_data)./npts)).^2;

        % Average power across electrode channels for frequencies 3, 4, 5 and 6 Hz
        avg_powerele = mean(power_spectrum(:, freq_idx), 1);

        % compute average fft_power for frequency range per trial
        avg_power(trial) = mean(avg_powerele);
    end

    %Loop through each channel and epoch and compute pwelch psd
    for chan = 1:6
        for k = 1:num_segmentsover
            [Pxx, f] = pwelch(EEG_over.data(chan, :,k), hamming(segment_length), noverlap, nfft, fs);
            power_pwelch(chan,:,k) = Pxx; % Store PSD for averaging
        end
    end

    %% correlation average power and segment number

    % correlate average theta/alpha power per segment with the segment number
    corr_coeff = corrcoef(avg_power(:,1), segment_num);
    power_change = corr_coeff(1, 2); %pearson

    % for reliability testing: taking every other segment

    seg_odd = segment_num(1:2:end); %taking only odd segments
    seg_even = segment_num(2:2:end); %taking only even segments

    corr_coeffr1 = corrcoef(avg_power(1:2:end, 1), seg_odd);
    corr_coeffr2 = corrcoef(avg_power(2:2:end, 1), seg_even);

    corr1 = corr_coeffr1(1, 2);
    corr2 = corr_coeffr2(1, 2);

    %calculate absolute and relative power

    frequency_band = freq_range;
    total_range = [1 45]; %for relative power
    psd_mean = mean(power_pwelch, 3); %psd averaged across trials
    psd_averaged_channels = log(mean(psd_mean, 1)); %psd averaged across channels
    absolute_power = trapz(f(frequency_band), psd_averaged_channels(frequency_band));

    psd_averaged_channels_nolog = mean(psd_mean, 1);
    total_power_nolog = trapz(f(f >= total_range(1) & f <= total_range(2)), ...
                    psd_averaged_channels_nolog(f >= total_range(1) & f <= total_range(2)));
    power_nolog = trapz(f(frequency_band), psd_averaged_channels_nolog(frequency_band));

    relative_power = power_nolog / total_power_nolog;


    % save new correlation coefficient for each dataset
    % fill cell with correlation
    colmodulation = 2;
    colr1modulation = 3;
    colr2modulation = 4;
    colabspow = 6;
    colrelpow = 7;
    v_modulation = power_change(:,1);
    v_r1modulation = corr1(:, 1);
    v_r2modulation = corr2(:, 1);
    % Assign the value to cell
    results{rowIndex, colmodulation} = v_modulation;
    results{rowIndex, colr1modulation} = v_r1modulation;
    results{rowIndex, colr2modulation} = v_r2modulation;
    results{rowIndex, colabspow} = absolute_power;
    results{rowIndex, colrelpow} = relative_power;

    % Save the updated table
    writetable(results, 'alphaT2_R.csv'); %change the name depending on frequency band/timepoint

    % Display the correlation coefficient
    disp(['Modulation index: ', num2str(power_change)]);


end
