%% Preprocessing Pipeline
%   1.  Import raw csv files -> band-pass -> save to 2_filter path 
%   2.  ASR bad-channel & segment removal -> save to 3_asr path
%   3.  Interpolate bad channels + re-reference to common average -> save
%   to 4_interp
%   4.  ICA (Picard) on retained channels -> save to 5_ica
%   5.   ICLabel components rejection -> save to 6_clean

%% Step 0: House-keeping
clear; 
clc; 
eeglab;

%% Step 0a: Setup paths
data_path   = '/path/to/EEG_2';
raw_path    = fullfile(data_path,'1_segmented');
filter_path = fullfile(data_path,'2_filter');
asr_path    = fullfile(data_path,'3_asr');
interp_path = fullfile(data_path,'4.1_interp');
ica_path    = fullfile(data_path,'5.1_ica');
clean_path  = fullfile(data_path,'6.1_clean');
chanloc_file = '/path/to/Standard-10-5-Cap385_witheog.elp';


required_channels = {'Fp1','F7','F8','T4','T6','T5','T3','Fp2','O1','P3','Pz', ...
                    'F3','Fz','F4','C4','P4','POz','C3','Cz','O2'};
%% Step 1: Import raw eeg csv files -> filter -> save to 2_filter
raw_files = dir(fullfile(raw_path,'*.csv'));
for i = 1:numel(raw_files)
    infile = fullfile(raw_path, raw_files(i).name);
    fprintf('Step1 [%d/%d]: Import %s\n', i, numel(raw_files), raw_files(i).name);

    % Read channel names (row 1) and numeric data (rows 2:end)
    opts = detectImportOptions(infile);
    opts.VariableNamesLine = 1;
    chNames = readtable(infile, opts, 'ReadVariableNames', true).Properties.VariableNames;
    opts.DataLines = [2 Inf];
    mat = readmatrix(infile, opts)';  % transpose to [chan x time]

    % Create EEGLAB dataset
    EEG = eeg_emptyset();
    EEG.data = mat;
    EEG.srate = 256;
    EEG.nbchan = numel(chNames);
    EEG.pnts = size(mat,2);
    EEG.setname = erase(raw_files(i).name, '.csv');
    EEG.filename = [EEG.setname '.set'];
    EEG.filepath = filter_path;

    % Assign channel labels and locations
    EEG.chanlocs = struct('labels', chNames);
    EEG = pop_chanedit(EEG, 'lookup', chanloc_file);

    % Store original channel structure
    EEG.etc.orig_chanlocs = EEG.chanlocs;
    EEG.etc.orig_nbchan = EEG.nbchan;

    % Band-pass and notch filters
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.5); % highpass filter 
    EEG = pop_eegfiltnew(EEG, 'hicutoff',40); % lowpass filter
    EEG = pop_eegfiltnew(EEG, 'locutoff',58,'hicutoff',62,'revfilt',1); % notch filter

    % Save filtered dataset
    EEG = pop_saveset(EEG, 'filename', EEG.filename, 'filepath', filter_path);
end

%% Step 2: ASR bad-channel & segment removal -> save to 3_asr
filter_sets = dir(fullfile(filter_path,'*.set'));
for i = 1:numel(filter_sets)
    fname = filter_sets(i).name;
    EEG = pop_loadset('filename',fname,'filepath',filter_path);

    % Restore original chanlocs if missing
    if ~isfield(EEG.etc,'orig_chanlocs')
        EEG = pop_chanedit(EEG,'lookup',chanloc_file);
        EEG.etc.orig_chanlocs = EEG.chanlocs;
        EEG.etc.orig_nbchan = EEG.nbchan;
    end

    % Run ASR
    pre_labels = {EEG.chanlocs.labels};
    EEG = pop_clean_rawdata(EEG, ...
        'FlatlineCriterion', 4, ...
        'ChannelCriterion', 0.9, ...
        'LineNoiseCriterion', 4, ...
        'Highpass', 'off', ...
        'BurstCriterion', 20, ...
        'WindowCriterion', 0.25, ...
        'BurstRejection', 'on', ...
        'WindowCriterionTolerances', [-Inf 7], ...
        'Distance', 'Euclidian');

    post_labels = {EEG.chanlocs.labels};
    EEG.etc.removed_channels = setdiff(pre_labels,post_labels);

    % Save ASR-cleaned set
    EEG = pop_saveset(EEG,'filename',fname,'filepath',asr_path);
end

%% Step 3: Interpolate bad channels + re-reference to common average -> 4_interp/
asr_sets = dir(fullfile(asr_path,'*.set'));
for i = 1:numel(asr_sets)
    fname = asr_sets(i).name;
    fprintf('Step3 [%d/%d]: %s\n', i, numel(asr_sets), fname);
    EEG = pop_loadset('filename',fname,'filepath',asr_path);

    orig = EEG.etc.orig_chanlocs;
    if EEG.nbchan < numel(orig)
        EEG = pop_interp(EEG,orig,'spherical');
    end

    curr = {EEG.chanlocs.labels};
    miss = setdiff(required_channels,curr);
    if ~isempty(miss), warning('Missing: %s',strjoin(miss,', ')); end

    EEG = pop_reref(EEG, []); % re-reference

    pop_saveset(EEG,'filename',fname,'filepath',interp_path);
end

%% Step 4: Run Independent Component Analysis → 5_ica/
interp_sets = dir(fullfile(interp_path,'*.set'));
for i = 1:numel(interp_sets)
    fname = interp_sets(i).name;
    fprintf('Step4 [%d/%d]: %s\n', i, numel(interp_sets), fname);
    EEG = pop_loadset('filename',fname,'filepath',interp_path);

    % Check Data Rank after interpolation and re-reference
    % Count interpolated channels 
    if isfield(EEG.etc, 'removed_channels') && ~isempty(EEG.etc.removed_channels)
        num_interpolated = length(EEG.etc.removed_channels);
        theoretical_rank = EEG.nbchan - num_interpolated;
        fprintf('Interpolated channels: %d\n', num_interpolated);
    else
        theoretical_rank = EEG.nbchan;
        num_interpolated = 0;
    end
    
    % Empirical rank from eigenvalues
    eigenvalues = eig(cov(double(EEG.data'))); % check number of non-zero eigenvalues
    tolerance = max(size(EEG.data)) * eps(max(eigenvalues));
    empirical_rank = sum(eigenvalues > tolerance);
    
    % Use the more conservative estimate of rank 
    dataRank = min(empirical_rank, theoretical_rank);
    
    % Check valid PCA dimension
    pcaDim = min(dataRank, EEG.nbchan - 1);
   
    % Run runica with PCA
    EEG = pop_runica(EEG,'icatype','picard','pca',pcaDim,'maxiter',500);

    EEG = pop_saveset(EEG,'filename',fname,'filepath',ica_path);
end

%% Step 5: ICLabel & components rejection ->  6_clean/
ica_sets = dir(fullfile(ica_path,'*.set'));
flagMat = [ NaN NaN; 0.9 1.0; 0.9 1.0; 0.9 1.0; NaN NaN; NaN NaN; NaN NaN ]; % eye, heart, muscle
for i = 1:numel(ica_sets)
    fname = ica_sets(i).name;
    fprintf('Step5 [%d/%d]: %s\n', i, numel(ica_sets), fname);
    EEG = pop_loadset('filename',fname,'filepath',ica_path);

    if ~isempty(EEG.icaweights)
        EEG = pop_iclabel(EEG,'default');
        EEG = pop_icflag(EEG,flagMat);
        bad = find(EEG.reject.gcompreject);
        if ~isempty(bad), EEG = pop_subcomp(EEG,bad,0); end
    end

    EEG = pop_rejcont(EEG, 'elecrange', 1:EEG.nbchan, ...
                             'threshold', 100, ...
                             'verbose','off');   % legacy continuous spectral rejection (threshold is in dB)

    % Final verify channels
    curr = {EEG.chanlocs.labels};
    miss = setdiff(required_channels,curr);
    if ~isempty(miss), warning('Missing: %s',strjoin(miss,', ')); end

    EEG = pop_saveset(EEG,'filename',fname,'filepath',clean_path);
end
