% Suppress the single tallest narrow peak based on its FWHM.
function [processed_data, fwhm_removed_trials] = removeNarrowPeaksByFWHM(nerve_data, Fs, bg_signal, replacement_multiplier, width_thresh_ms, std_mad)
    fwhm_removed_trials = {};
    % Start with a copy of the input trials.
    processed_data = nerve_data; 
    
    trials = fieldnames(nerve_data);
    dt = 1 / Fs;

    for i = 1:length(trials)
        t_name = trials{i};
        
        % Only process trials with an EMG field.
        if ~isstruct(nerve_data.(t_name)) || ~isfield(nerve_data.(t_name), 'EMG')
            continue; 
        end
        
        % Read the original EMG signal.
        raw_sig = nerve_data.(t_name).EMG;
        N = length(raw_sig);

        % Find the tallest peak in the absolute signal.
        abs_sig = abs(raw_sig);
        [max_val, max_idx] = max(abs_sig);
        
        % =======================================================
        % Compute an adaptive threshold, as in the CLEAN step.
        % =======================================================
        % HATA threshold retained for reference.
        sigma_noise = median(abs_sig) / 0.6745; 
        log2_n = log2(N);
        term_inside_sqrt = (2 * (log2_n^4)) / N;
        heuristic_factor = sqrt(term_inside_sqrt);
        dynamic_thresh = sigma_noise * heuristic_factor * std_mad;
        
        % Alternative mean-plus-standard-deviation threshold.
        % sig_mean = mean(abs_sig);
        % sig_std  = std(abs_sig);
        % dynamic_thresh = sig_mean + std_mad * sig_std;
        % =======================================================
        
        % Criterion 1: the tallest peak must exceed the amplitude threshold.
        if max_val > 30 % dynamic_thresh
            half_max = max_val / 2;

            % Find the half-maximum crossings.
            left_idx = find(abs_sig(1:max_idx) <= half_max, 1, 'last');
            right_idx = find(abs_sig(max_idx:end) <= half_max, 1, 'first') + max_idx - 1;

            if isempty(left_idx), left_idx = 1; end
            if isempty(right_idx), right_idx = N; end

            fwhm_pts = right_idx - left_idx;
            fwhm_ms = fwhm_pts * dt * 1000;

            % Criterion 2: replace peaks narrower than the FWHM limit.
            if fwhm_ms < width_thresh_ms
                % Determine the replacement interval.
                radius = round((fwhm_pts / 2) * replacement_multiplier);
                start_idx = max(1, max_idx - radius);
                end_idx = min(N, max_idx + radius);
                seg_len = end_idx - start_idx + 1;
                indices = start_idx:end_idx;

                % Replace the interval with external background noise.
                if ~isempty(bg_signal) && length(bg_signal) >= seg_len
                    max_rand = length(bg_signal) - seg_len + 1;
                    r_start = randi(max_rand);
                    fill_raw = bg_signal(r_start : r_start + seg_len - 1);

                    raw_sig(indices) = fill_raw;
                    
                    % Store the updated EMG signal.
                    processed_data.(t_name).EMG = raw_sig; 
                    
                    % Record the processed trial.
                    fwhm_removed_trials{end+1} = t_name; 
                end
            end
        end
    end
end
