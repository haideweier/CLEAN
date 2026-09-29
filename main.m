% CLEAN denoising example.
% Input:  data/data.mat (one anonymized subject, muscle, and nerve)
%         data/background_replacement.mat
% Output: result/cleaned_data.mat

clear; clc;
project_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(project_dir, 'preprocess'));

input_file = fullfile(project_dir, 'data', 'data.mat');
background_file = fullfile(project_dir, 'data', 'background_replacement.mat');
output_dir = fullfile(project_dir, 'result');
output_file = fullfile(output_dir, 'cleaned_data.mat');

% Sampling frequency and CLEAN parameters from the original analysis.
Fs = 1 / 0.000156;
late_window_ms = 55;
early_window_ms = 45;
similarity_threshold = 0.75;
neighbor_search_range = 2;
attenuation_factor = 10;
override_similarity_threshold = 0.90;
replacement_window_multiplier = 15.0;
prominence_method = 'std';
std_mad = 1;
window_multiplier = 5.0;
fwhm_width_threshold_ms = 0.5;

% Load the anonymized trials and background replacement signal.
loaded = load(input_file, 'data');
if ~isfield(loaded, 'data') || ...
        ~isfield(loaded.data, 'Level_01') || ...
        ~isfield(loaded.data.Level_01, 'Subject_01') || ...
        ~isfield(loaded.data.Level_01.Subject_01, 'Muscle_01') || ...
        ~isfield(loaded.data.Level_01.Subject_01.Muscle_01, 'Nerve_01')
    error('Input file does not contain the expected anonymized data fields.');
end
data = loaded.data;
trials = data.Level_01.Subject_01.Muscle_01.Nerve_01;

loaded_bg = load(background_file, 'background_data');
if ~isfield(loaded_bg, 'background_data') || ...
        ~isfield(loaded_bg.background_data, 'Background_01')
    error('Background file must contain background_data.Background_01.');
end
background_signal = loaded_bg.background_data.Background_01(:)';
if isempty(background_signal) || ~isnumeric(background_signal)
    error('The background signal must be a nonempty numeric vector.');
end

% Step 1: suppress narrow peaks. Step 2: apply CLEAN.
[preprocessed, ~] = removeNarrowPeaksByFWHM( ...
    trials, Fs, background_signal, replacement_window_multiplier, ...
    fwhm_width_threshold_ms, std_mad);
[cleaned, ~, ~, ~] = cleanEarlySignal_HATA_threshold( ...
    preprocessed, Fs, late_window_ms, early_window_ms, ...
    similarity_threshold, prominence_method, window_multiplier, ...
    replacement_window_multiplier, neighbor_search_range, ...
    override_similarity_threshold, std_mad, attenuation_factor, ...
    background_signal);

% Preserve trial names and Stimulus_Amplitude; replace only EMG values.
trial_names = fieldnames(trials);
for i = 1:numel(trial_names)
    trial_name = trial_names{i};
    if ~isfield(trials.(trial_name), 'EMG')
        error('Trial %s does not contain EMG.', trial_name);
    end
    if ~isfield(cleaned, trial_name) || ...
            ~isfield(cleaned.(trial_name), 'cleaned_raw')
        error('CLEAN did not return a signal for trial %s.', trial_name);
    end
    data.Level_01.Subject_01.Muscle_01.Nerve_01.(trial_name).EMG = ...
        cleaned.(trial_name).cleaned_raw;
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(output_file, 'data');
fprintf('Saved %d cleaned trials to %s\n', numel(trial_names), output_file);
