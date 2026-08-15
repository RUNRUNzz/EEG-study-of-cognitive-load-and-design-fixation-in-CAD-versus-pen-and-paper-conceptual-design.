%% EEG feature extraction used for the current CHB analyses
% Produces:
%   1. Frequency_bands_power.csv  (channel-level band power for TRP/MCA)
%   2. eeg_features.csv           (session-level PSD summaries and CLI)
%
% IMPORTANT: This file preserves the numerical settings in the uploaded
% analysis source that generated the current derived results. See
% VALIDATION_NOTES.md before public release because several settings differ
% from the wording in the current manuscript.

clear;
clc;
eeglab;

%% Paths
project_root = '/path/to/EEG_2';
data_path = fullfile(project_root, '6.1_clean');
output_path = fullfile(project_root, '7.1_initial_analysis');
if ~exist(output_path, 'dir'), mkdir(output_path); end

%% Spectral settings used in the original analysis
Fs = 256;
window_sec = 1;
window_len = Fs * window_sec;
nfft = 2^nextpow2(window_len);  % 256 points in the uploaded analysis

bands = struct( ...
    'delta', [1, 4], ...
    'theta', [4, 8], ...
    'alpha', [8, 13], ...
    'beta',  [13, 30], ...
    'gamma', [30, 40]);
band_names = fieldnames(bands);

% Channel sets used in the uploaded CLI calculation
frontal_ch = {'Fp1','Fp2','F3','Fz','F4','F7','F8'};
parietal_ch = {'P3','Pz','P4'};

%% 1. Band power for each channel
files = dir(fullfile(data_path, '*.set'));
results = cell(0, 9);

for i = 1:numel(files)
    filename = files(i).name;
    EEG = pop_loadset('filename', filename, 'filepath', data_path);

    [~, basename] = fileparts(filename);
    parts = strsplit(basename, '_');
    participant = parts{1};
    condition = parts{2};
    phase = strjoin(parts(3:end), '_');

    for ch = 1:EEG.nbchan
        channel = strtrim(EEG.chanlocs(ch).labels);
        n_windows = floor(EEG.pnts / window_len);
        all_psd = zeros(n_windows, nfft / 2 + 1);
        frequencies = [];

        for w = 1:n_windows
            idx = (w - 1) * window_len + 1 : w * window_len;
            segment = double(EEG.data(ch, idx));
            segment = segment - mean(segment);

            % The uploaded analysis used zero overlap within each 1-s segment.
            [pxx, f] = pwelch(segment, hamming(window_len), 0, nfft, Fs, 'psd');
            all_psd(w, :) = pxx';
            if isempty(frequencies), frequencies = f; end
        end

        mean_psd = mean(all_psd, 1);
        power = struct();

        for b = 1:numel(band_names)
            range = bands.(band_names{b});
            use = frequencies >= range(1) & frequencies <= range(2);
            power.(band_names{b}) = trapz(frequencies(use), mean_psd(use));
        end

        results(end + 1, :) = {participant, condition, phase, channel, ... %#ok<SAGROW>
            power.delta, power.theta, power.alpha, power.beta, power.gamma};
    end
end

if isempty(results), error('No EEG feature results were generated.'); end

band_power = cell2table(results, 'VariableNames', ...
    {'Participant','Condition','Phase','Channel','Delta','Theta','Alpha','Beta','Gamma'});
writetable(band_power, fullfile(output_path, 'Frequency_bands_power.csv'));

%% 2. Session-level PSD summaries
band_power.PSD_Aggregated = band_power.Delta + band_power.Theta + ...
    band_power.Alpha + band_power.Beta + band_power.Gamma;

features = groupsummary(
    band_power,
    {'Participant','Condition','Phase'},
    'mean',
    {'PSD_Aggregated','Delta','Theta','Alpha','Beta','Gamma'});

features.GroupCount = [];
features = renamevars(features, ...
    {'mean_PSD_Aggregated','mean_Delta','mean_Theta','mean_Alpha','mean_Beta','mean_Gamma'}, ...
    {'PSD_Aggregated','PSD_Delta','PSD_Theta','PSD_Alpha','PSD_Beta','PSD_Gamma'});

%% 3. Cognitive Load Index (CLI)
frontal = band_power(ismember(band_power.Channel, frontal_ch), :);
frontal = groupsummary(frontal, {'Participant','Condition','Phase'}, 'mean', 'Theta');
frontal = renamevars(frontal, 'mean_Theta', 'ThetaFrontal');
frontal.GroupCount = [];

parietal = band_power(ismember(band_power.Channel, parietal_ch), :);
parietal = groupsummary(parietal, {'Participant','Condition','Phase'}, 'mean', 'Alpha');
parietal = renamevars(parietal, 'mean_Alpha', 'AlphaParietal');
parietal.GroupCount = [];

cli = outerjoin(frontal, parietal, ...
    'Keys', {'Participant','Condition','Phase'}, 'MergeKeys', true);
cli.CLI = cli.ThetaFrontal ./ cli.AlphaParietal;

features = outerjoin(features, ...
    cli(:, {'Participant','Condition','Phase','CLI'}), ...
    'Keys', {'Participant','Condition','Phase'}, 'MergeKeys', true);

writetable(features, fullfile(output_path, 'eeg_features.csv'));
