function [cleaned_data, mean_noise_template, noise_waveforms_cell, fallback_trials] = cleanEarlySignal_HATA_threshold(nerve_data, Fs, late_window_ms, early_window_ms, similarity_threshold, prominence_method, window_multiplier, replacement_window_multiplier, neighbor_search_range, override_similarity_threshold, std_mad, attenuation_factor, external_bg_signal)
% CLEAN processing with HATA detection, external background replacement, and late-window cleaning.
% 1. Preserve the peak detection and template matching logic.
% 2. Return the preprocessed signals if too few noise templates are found.
% 3. Replace detected artifacts with the external background signal.
% 4. Apply the global noise check with boundary protection.
% 5. Detect and replace artifacts in the late window.

% --- Initialization ---
amplititude_range = 1.5; % Default amplitude ratio
cleaned_data = struct();
mean_noise_template = [];
noise_waveforms_cell = {};
fallback_trials = {}; 
% Store late-window noise locations for possible later replacement.
all_late_info = struct();
all_fields = fieldnames(nerve_data);
x_fields = all_fields(strncmp(all_fields, 'X_', 2));
if isempty(x_fields), warning('No trial fields starting with "X_" were found.'); return; end

% The caller loads and validates the anonymized background signal.

% --- Sort trials by their numeric suffix ---
trial_nums = zeros(length(x_fields), 1);
for i = 1:length(x_fields)
    num = sscanf(x_fields{i}, 'X_%d');
    if ~isempty(num)
        trial_nums(i) = num;
    else
        trial_nums(i) = NaN;
    end
end
[~, sorted_indices] = sort(trial_nums);
sorted_x_fields = x_fields(sorted_indices);
num_trials = length(sorted_x_fields);

%% =================== Step 1: Learn HATA noise templates ===================
fprintf('>>> Step 1: Learning noise templates from the final %d ms of %d trials (HATA threshold)...\n', late_window_ms, num_trials);
max_noise_len = 0;
for i = 1:num_trials
    current_x = sorted_x_fields{i}; 
    emg_raw = nerve_data.(current_x).EMG;
    emg_rectified = abs(emg_raw);
    N = length(emg_rectified);    
    
    % Prepare the late-window samples.
    late_samples = round(late_window_ms / 1000 * Fs);
    region_start_idx = N - late_samples + 1;
    if region_start_idx < 1, region_start_idx = 1; end
    late_region = emg_rectified(region_start_idx:end);
    
    % Compute the HATA threshold.
    sigma_noise = median(emg_rectified) / 0.6745; 
    log2_n = log2(N);
    term_inside_sqrt = (2 * (log2_n^4)) / N;
    heuristic_factor = sqrt(term_inside_sqrt);
    q = std_mad; 
    min_prominence = sigma_noise * heuristic_factor * q;
    
    [~, locs_local, w] = findpeaks(late_region, 'MinPeakProminence', min_prominence);
    locs_global = locs_local + region_start_idx - 1;

    current_late_peaks = []; % Late-window noise peaks in this trial
    
    for j = 1:length(locs_global)
        peak_loc = locs_global(j); 
        radius = round(w(j) / 2 * window_multiplier);
        start_idx = peak_loc - radius; end_idx = peak_loc + radius;
        if start_idx >= 1 && end_idx <= N
            segment = emg_rectified(start_idx:end_idx);
            if max(segment) > 10
                noise_waveforms_cell{end+1} = segment;
                if length(segment) > max_noise_len, max_noise_len = length(segment); end

                if length(segment) > max_noise_len, max_noise_len = length(segment); end
                % Record the peak location for later replacement.
                l_pk.location = peak_loc;
                l_pk.width = w(j);
                current_late_peaks = [current_late_peaks, l_pk];
            end
        end
    end
    all_late_info.(current_x) = current_late_peaks; % Save this trial's late-window peaks
end

num_noise_templates = length(noise_waveforms_cell);
min_required_templates = 0.10 * num_trials;

% If no templates or too few templates are found, pass the preprocessed signals through.
if isempty(noise_waveforms_cell) || num_noise_templates < min_required_templates
    
    % Report why template matching is skipped.
    if isempty(noise_waveforms_cell)
        warning('No noise peaks were found in the late window; skipping HATA template matching.');
    else
        warning('Only %d noise templates were found, below 10%% of %d trials; skipping template matching.', num_noise_templates, num_trials);
    end
    
    % Return the FWHM-preprocessed signal for each trial.
    mean_noise_template = []; noise_waveforms_cell = {};
    for i = 1:num_trials
        current_x = sorted_x_fields{i}; 
        emg_raw = nerve_data.(current_x).EMG; 
        N_current = length(emg_raw);
        emg_rectified = abs(emg_raw);
        
        cleaned_data.(current_x).original_rectified = emg_rectified;
        cleaned_data.(current_x).cleaned_rectified = emg_rectified;
        cleaned_data.(current_x).original_raw = emg_raw;
        cleaned_data.(current_x).cleaned_raw = emg_raw; % Preserve the FWHM preprocessing result
        cleaned_data.(current_x).t_ms = (0:N_current-1) / Fs * 1000;
        cleaned_data.(current_x).peaks_info = [];
    end
    return;
end
fprintf('Template learning complete: %d noise templates above the amplitude limit of 10.\n\n', length(noise_waveforms_cell));

% Compute the mean template for output; matching uses the template cell array.
noise_matrix = NaN(length(noise_waveforms_cell), max_noise_len);
for i = 1:length(noise_waveforms_cell)
    waveform = noise_waveforms_cell{i}; len = length(waveform);
    start_col = floor((max_noise_len - len) / 2) + 1;
    noise_matrix(i, start_col:(start_col+len-1)) = waveform;
end
mean_noise_template = nanmean(noise_matrix, 1);
mean_noise_template(isnan(mean_noise_template)) = [];


%% =================== Step 2: Detect candidate noise peaks ===================
fprintf('>>> Step 2: Detecting candidate noise peaks in the early and late windows...\n');
all_trials_info = struct();
trial_is_noisy_neighbor = struct();

for i = 1:num_trials
    current_x = sorted_x_fields{i};
    emg_raw = nerve_data.(current_x).EMG;
    emg_rectified = abs(emg_raw); 
    N = length(emg_rectified);
    
    % HATA threshold for the full segment.
    sigma_noise = median(emg_rectified) / 0.6745; 
    log2_n = log2(N);
    term_inside_sqrt = (2 * (log2_n^4)) / N;
    heuristic_factor = sqrt(term_inside_sqrt);
    q = std_mad; 
    min_prominence = sigma_noise * heuristic_factor * q;
    
    % Store all detected peaks for this trial.
    peaks_info = [];
    
    % -------------------------------------------------------------
    % Region A: search the early window.
    % -------------------------------------------------------------
    early_samples = round(early_window_ms / 1000 * Fs);
    buffer_ms = 5; 
    buffer_samples = round(buffer_ms / 1000 * Fs);
    search_end_idx = min(early_samples + buffer_samples, N);
    
    early_region_buffered = emg_rectified(1:search_end_idx);
    
    % Detect early peaks, then apply shadowing and filtering.
    padded_signal = [0; early_region_buffered(:)];
    [pks_raw, locs_raw, w_raw] = findpeaks(padded_signal, 'MinPeakProminence', min_prominence);
    locs_raw = locs_raw - 1;
    valid_mask = locs_raw > 0;
    
    locs_candidates = locs_raw(valid_mask);
    pks_candidates = pks_raw(valid_mask);
    w_candidates = w_raw(valid_mask);
    
    % Shadowing: a larger peak suppresses nearby smaller peaks.
    keep_mask = true(size(locs_candidates));
    check_radius = 25; amplitude_ratio = 1.2; 
    for k = 1:length(locs_candidates)
        my_amp = pks_candidates(k); my_loc = locs_candidates(k);
        dist_vec = abs(locs_candidates - my_loc);
        neighbors_idx = find(dist_vec > 0 & dist_vec <= check_radius);
        if ~isempty(neighbors_idx)
            max_neighbor_amp = max(pks_candidates(neighbors_idx));
            if max_neighbor_amp > (my_amp * amplitude_ratio), keep_mask(k) = false; end
        end
    end
    
    locs_survivors = locs_candidates(keep_mask);
    pks_survivors = pks_candidates(keep_mask);
    w_survivors = w_candidates(keep_mask);
    
    % Min Distance Logic
    final_keep_mask = true(size(locs_survivors));
    [~, sort_idx] = sort(pks_survivors, 'descend');
    target_min_dist = 25; 
    for k = 1:length(sort_idx)
        idx_curr = sort_idx(k);
        if ~final_keep_mask(idx_curr), continue; end
        loc_curr = locs_survivors(idx_curr);
        dist_vec = abs(locs_survivors - loc_curr);
        conflict_idx = find(dist_vec > 0 & dist_vec < target_min_dist);
        final_keep_mask(conflict_idx) = false;
    end
    
    locs_early = locs_survivors(final_keep_mask);
    w_early = w_survivors(final_keep_mask);
    
    % Add the early peaks.
    for j = 1:length(locs_early)
        peak_loc = locs_early(j);
        if peak_loc > early_samples, continue; end % Ignore peaks outside the early window
        
        % Match the peak against the learned templates.
        % Use the shared template-matching function.
        radius = round(w_early(j) / 2 * window_multiplier);
        start_idx = peak_loc - radius; end_idx = peak_loc + radius;
        if start_idx < 1, start_idx = 1; end
        if start_idx >= 1 && end_idx <= N
             % Match each candidate against the learned noise templates.
             [final_idx, final_sim, is_noise, best_tmpl_amp] = matchTemplate(emg_rectified(start_idx:end_idx), noise_waveforms_cell, similarity_threshold, override_similarity_threshold);
             
             peak_struct = struct();
             peak_struct.location = peak_loc;
             peak_struct.width = w_early(j);
             peak_struct.similarity = round(final_sim * 100) / 100;
             peak_struct.best_template_idx = final_idx;
             peak_struct.amplitude = pks_survivors(find(locs_survivors==peak_loc,1)); % Approximate amplitude
             if isempty(peak_struct.amplitude), peak_struct.amplitude = max(emg_rectified(start_idx:end_idx)); end
             peak_struct.best_template_amplitude = best_tmpl_amp;
             peak_struct.is_potential_noise = is_noise;
             
             peaks_info = [peaks_info, peak_struct];
        end
    end
    
    % -------------------------------------------------------------
    % Region B: search the late window.
    % -------------------------------------------------------------
    late_samples = round(late_window_ms / 1000 * Fs);
    late_start_global_idx = N - late_samples + 1;
    if late_start_global_idx < 1, late_start_global_idx = 1; end
    
    late_region = emg_rectified(late_start_global_idx:end);
    
    % Apply the same peak detection to the late region.
    % Zero-pad the region to retain peaks near its boundaries.
    padded_late = [0; late_region(:)]; 
    [pks_late_raw, locs_late_raw, w_late_raw] = findpeaks(padded_late, 'MinPeakProminence', min_prominence);
    locs_late_raw = locs_late_raw - 1;
    valid_late = locs_late_raw > 0;
    
    locs_late_c = locs_late_raw(valid_late);
    pks_late_c = pks_late_raw(valid_late);
    w_late_c = w_late_raw(valid_late);
    
    % Shadowing Logic (Late)
    keep_mask_l = true(size(locs_late_c));
    for k = 1:length(locs_late_c)
        my_amp = pks_late_c(k); my_loc = locs_late_c(k);
        dist_vec = abs(locs_late_c - my_loc);
        neighbors_idx = find(dist_vec > 0 & dist_vec <= check_radius);
        if ~isempty(neighbors_idx)
            max_neighbor_amp = max(pks_late_c(neighbors_idx));
            if max_neighbor_amp > (my_amp * amplitude_ratio), keep_mask_l(k) = false; end
        end
    end
    
    locs_late_surv = locs_late_c(keep_mask_l);
    pks_late_surv = pks_late_c(keep_mask_l);
    w_late_surv = w_late_c(keep_mask_l);
    
    % Add the late peaks.
    for j = 1:length(locs_late_surv)
        % Convert the location to a global sample index.
        local_loc = locs_late_surv(j);
        global_loc = local_loc + late_start_global_idx - 1;
        
        radius = round(w_late_surv(j) / 2 * window_multiplier);
        start_idx = global_loc - radius; end_idx = global_loc + radius;
        if start_idx < 1, start_idx = 1; end
        if end_idx > N, end_idx = N; end
        
        if start_idx >= 1 && end_idx <= N
             % Late peaks often match well because the templates come from this region.
             [final_idx, final_sim, is_noise, best_tmpl_amp] = matchTemplate(emg_rectified(start_idx:end_idx), noise_waveforms_cell, similarity_threshold, override_similarity_threshold);
             
             peak_struct = struct();
             peak_struct.location = global_loc; % Store the global sample index
             peak_struct.width = w_late_surv(j);
             peak_struct.similarity = round(final_sim * 100) / 100;
             peak_struct.best_template_idx = final_idx;
             peak_struct.amplitude = pks_late_surv(j);
             peak_struct.best_template_amplitude = best_tmpl_amp;
             peak_struct.is_potential_noise = is_noise;
             
             peaks_info = [peaks_info, peak_struct];
        end
    end

    all_trials_info.(current_x) = peaks_info;

    % Identify noisy neighboring trials.
    is_noisy_flag = false;
    % Recompute early_samples for the neighbor check.
    early_limit_samples = round(early_window_ms / 1000 * Fs);
    
    if ~isempty(peaks_info)
        for k = 1:length(peaks_info)
            pk = peaks_info(k);
            
            % Late-window noise must not affect the noisy-neighbor decision.
            % Mark a trial as a noisy neighbor only for early-window artifacts.
            if pk.location > early_limit_samples
                continue; 
            end
            
            if pk.is_potential_noise
                is_amplitude_ok = false;
                if ~isnan(pk.best_template_amplitude) && pk.best_template_amplitude > 0
                    if pk.amplitude < amplititude_range * pk.best_template_amplitude && amplititude_range * pk.amplitude > pk.best_template_amplitude
                        is_amplitude_ok = true;
                    end
                end
                if is_amplitude_ok
                    is_noisy_flag = true; 
                    break;
                end
            end
        end
    end
    trial_is_noisy_neighbor.(current_x) = is_noisy_flag;
end
fprintf('Peak detection complete.\n\n');

%% =================== Step 3: Replace artifacts with background ===================
% The replacement logic uses global peak locations.
% Late-window peaks added in Step 2 are handled here.

fprintf('>>> Step 3: Cleaning with external background replacement, including the late window...\n');
for i = 1:num_trials
    current_x = sorted_x_fields{i};
    target_emg_raw = nerve_data.(current_x).EMG;
    target_emg_rect = abs(target_emg_raw);
    N_target = length(target_emg_raw);
    
    emg_cleaned_rectified = target_emg_rect;
    emg_cleaned_raw = target_emg_raw;            
    
    is_replaced_mask = false(1, N_target); 
    
    current_peaks_info = all_trials_info.(current_x);
    
    % Determine which peaks to replace.
    for j = 1:length(current_peaks_info)
        peak = current_peaks_info(j);
        
        % 1. Amplitude check.
        is_amplitude_valid = false;
        if isfield(peak, 'best_template_amplitude') && ~isnan(peak.best_template_amplitude) && peak.best_template_amplitude > 0
            is_amplitude_valid = (peak.amplitude < 1.5 * peak.best_template_amplitude && 1.5 * peak.amplitude > peak.best_template_amplitude);
        end
        % is_amplitude_valid = true; % Ablation: omit the Phase 2 amplitude check
        % 2. Similarity check.
        is_high_sim = (peak.similarity >= override_similarity_threshold) && is_amplitude_valid;
        % is_high_sim = (peak.similarity >= similarity_threshold) && is_amplitude_valid; % Ablation: use a single Phase 2 threshold
        % 3. Context check.
        is_context_noise = false;
        if peak.is_potential_noise
            has_noisy_neighbor = false;
            for k = 1:neighbor_search_range
                if (i-k >= 1) && trial_is_noisy_neighbor.(sorted_x_fields{i-k}), has_noisy_neighbor = true; break; end
            end
            if ~has_noisy_neighbor
                for k = 1:neighbor_search_range
                    if (i+k <= num_trials) && trial_is_noisy_neighbor.(sorted_x_fields{i+k}), has_noisy_neighbor = true; break; end
                end
            end
            if has_noisy_neighbor && is_amplitude_valid, is_context_noise = true; end
        end
        % is_context_noise = false; % Ablation: omit the Phase 2 context check

        % 4. Location check: before 2 ms or in the late window.
        peak_time_ms = peak.location / Fs * 1000; 
        is_too_early = peak_time_ms < 2.0;
        
        % Identify late-window peaks, for which replacement criteria are relaxed.
        % Replace late-window peaks that match a noise template.
        total_ms = N_target / Fs * 1000;
        is_in_late_window = peak_time_ms > (total_ms - late_window_ms);
        
        % if peak.best_template_idx == -1  % Newly added condition
        %     current_peaks_info(j).to_be_replaced = true;
            
        if is_high_sim || is_context_noise || is_too_early || (is_in_late_window && peak.is_potential_noise)
            current_peaks_info(j).to_be_replaced = true;
        else
            current_peaks_info(j).to_be_replaced = false;
        end
    end
    
    % =====================================================================
    % Step 3a: replace only early-window peaks.
    % Keep late-window noise until the global low-SNR check is complete.
    % =====================================================================
    final_peaks_info = current_peaks_info;
    replaced_count_early = 0;
    late_peaks_indices = []; % Save late-window peak indices for later
    
    late_start_global_idx = N_target - round(late_window_ms / 1000 * Fs) + 1;
    
    for j = 1:length(final_peaks_info)
        if final_peaks_info(j).to_be_replaced
            
            % Defer late-window peaks and save their indices.
            if final_peaks_info(j).location >= late_start_global_idx
                late_peaks_indices = [late_peaks_indices, j];
                continue; 
            end
            
            % Replace the early-window peak.
            [emg_cleaned_raw, emg_cleaned_rectified, is_replaced_mask, success] = ...
                performReplacement(emg_cleaned_raw, emg_cleaned_rectified, is_replaced_mask, final_peaks_info(j), N_target, replacement_window_multiplier, external_bg_signal);
            
            final_peaks_info(j).replaced = success;
            if success, final_peaks_info(j).replacement_source = 'External BG'; replaced_count_early = replaced_count_early + 1; end
        else
            final_peaks_info(j).replaced = false;
        end
    end
    
    if replaced_count_early > 0
        fprintf('  - %s (Step 3a): Replaced %d early-window peak(s).\n', current_x, replaced_count_early);
    end
    
    % =====================================================================
    % Step 3b: global noise check with boundary protection.
    % =====================================================================
    early_samples_count = round(early_window_ms / 1000 * Fs);
    late_samples_count = round(late_window_ms / 1000 * Fs);
    
    should_fallback = false;
    fallback_reason = '';
    global_check_thr = 10;
    
    if length(emg_cleaned_rectified) > max(early_samples_count, late_samples_count)
        
        % Maximum early-window amplitude after cleaning, excluding replaced samples.
        early_region_post = emg_cleaned_rectified(1:early_samples_count);
        early_mask_post = is_replaced_mask(1:early_samples_count);
        valid_indices = find(~early_mask_post);
        if isempty(valid_indices), max_early_post = 0; else, max_early_post = max(early_region_post(valid_indices)); end
        
        % Maximum amplitude in the original late region.
        late_region_dirty = emg_cleaned_rectified(late_start_global_idx:end);
        [max_late_dirty, idx_late_local] = max(late_region_dirty);
            
        % Branch 2: low SNR (late peak > 10 and 1.5 x late peak > early peak).
        if (max_early_post > global_check_thr) && (max_late_dirty > global_check_thr) && (1.5 * max_late_dirty > max_early_post)
            should_fallback = true; 
            fallback_reason = sprintf('low SNR (Late=%.1f vs Early=%.1f)', max_late_dirty, max_early_post);
            
            % Protect waves that cross from the early into the late window.
            % Convert the late peak location to a global sample index.
            idx_late_global = late_start_global_idx + idx_late_local - 1;
            
            % Search backward for the wave start, using amplitude < 5 as baseline.
            wave_start_idx = idx_late_global;
            for s = idx_late_global:-1:1
                if emg_cleaned_rectified(s) < 5
                    wave_start_idx = s;
                    break;
                end
                wave_start_idx = s; 
            end
            
            % A wave starting in the early window is not independent late-window noise.
            if wave_start_idx <= early_samples_count
                should_fallback = false; % Cancel whole-window replacement
                fprintf('>> Preserving %s: late peak at sample %d starts in the early window at sample %d.\n', current_x, idx_late_global, wave_start_idx);
            end
            % ---------------------------------------------
        end
        
        if should_fallback
            fprintf('!! Warning: %s has %s. Replacing the entire early window...\n', current_x, fallback_reason);
            segment_len_needed = early_samples_count;
            indices = 1:early_samples_count;
            if ~isempty(external_bg_signal) && length(external_bg_signal) >= segment_len_needed
                r_start = randi(length(external_bg_signal) - segment_len_needed + 1);
                best_raw = external_bg_signal(r_start : r_start + segment_len_needed - 1);
                emg_cleaned_raw(indices) = best_raw;
                emg_cleaned_rectified(indices) = abs(best_raw);
                fallback_trials{end+1} = current_x; 
                is_replaced_mask(indices) = true;
                fprintf('>> %s: Replaced the entire early window with external background.\n', current_x);
            end
        end
    end

    % =====================================================================
    % Step 3c: replace late-window peaks.
    % The global check is complete, so late-window peaks can now be replaced.
    % =====================================================================
    replaced_count_late = 0;
    if ~isempty(late_peaks_indices)
        for idx = late_peaks_indices
            % Replace the late-window peak.
            [emg_cleaned_raw, emg_cleaned_rectified, is_replaced_mask, success] = ...
                performReplacement(emg_cleaned_raw, emg_cleaned_rectified, is_replaced_mask, final_peaks_info(idx), N_target, replacement_window_multiplier, external_bg_signal);

            final_peaks_info(idx).replaced = success;
            if success, final_peaks_info(idx).replacement_source = 'External BG (Late)'; replaced_count_late = replaced_count_late + 1; end
        end
    end

    if replaced_count_late > 0
        fprintf('  - %s (Step 3c): Replaced %d late-window peak(s).\n', current_x, replaced_count_late);
    end
    
    cleaned_data.(current_x).original_rectified = target_emg_rect;
    cleaned_data.(current_x).cleaned_rectified = emg_cleaned_rectified;
    cleaned_data.(current_x).original_raw = target_emg_raw;
    cleaned_data.(current_x).cleaned_raw = emg_cleaned_raw; 
    cleaned_data.(current_x).t_ms = (0:N_target-1) / Fs * 1000;
    cleaned_data.(current_x).peaks_info = final_peaks_info;
end
fprintf('Cleaning complete.\n');
end

% Helper: replace one peak interval.
function [raw_out, rect_out, mask_out, success] = performReplacement(raw_in, rect_in, mask_in, peak, N, multiplier, bg_signal)
    raw_out = raw_in; rect_out = rect_in; mask_out = mask_in; success = false;
    
    radius = round(peak.width / 2 * multiplier);
    start_idx = max(1, peak.location - radius);
    end_idx = min(N, peak.location + radius);
    seg_len = end_idx - start_idx + 1;
    indices = start_idx:end_idx;
    
    if ~isempty(bg_signal) && length(bg_signal) >= seg_len
        max_rand = length(bg_signal) - seg_len + 1;
        r_start = randi(max_rand);
        fill_raw = bg_signal(r_start : r_start + seg_len - 1);
        
        raw_out(indices) = fill_raw;
        rect_out(indices) = abs(fill_raw);
        mask_out(indices) = true;
        success = true;
    end
end

% Helper: match a waveform against noise templates.
function [final_idx, final_sim, is_noise, best_tmpl_amp] = matchTemplate(waveform, noise_templates, sim_thresh, override_thresh)
    global_max_sim = -1; 
    global_best_idx = 0;
    
    candidate_indices = [];
    candidate_sims = [];
    candidate_amps = [];
    
    wave_amp = max(waveform);
    
    for i = 1:length(noise_templates)
        tmpl = noise_templates{i};
        
        % Align lengths and smooth the waveforms.
        len_w = length(waveform); len_t = length(tmpl);
        if len_w <= len_t
             center_t = floor(len_t/2)+1; half_w = floor(len_w/2);
             s_t = center_t - half_w; e_t = s_t + len_w - 1;
             if s_t < 1, s_t=1; e_t=s_t+len_w-1; end
             if e_t > len_t, e_t=len_t; s_t=e_t-len_w+1; end
             tmpl_seg = tmpl(s_t:e_t);
        else
             pad = len_w - len_t;
             tmpl_seg = [zeros(1, floor(pad/2)), tmpl, zeros(1, ceil(pad/2))];
        end
        
        a = smoothdata(waveform, 'gaussian', 5);
        b = smoothdata(tmpl_seg, 'gaussian', 5);
        a_n = a - mean(a); b_n = b - mean(b);
        
        if any(isnan(a_n)) || any(isnan(b_n)) || all(a_n==0) || all(b_n==0)
            sim = 0;
        else
            sim = sum(a_n.*b_n) / (sqrt(sum(a_n.^2)) * sqrt(sum(b_n.^2)));
        end
        
        if sim > global_max_sim, global_max_sim = sim; global_best_idx = i; end
        if sim >= sim_thresh
            candidate_indices(end+1) = i;
            candidate_sims(end+1) = sim;
            candidate_amps(end+1) = max(tmpl);
        end
    end
    
    % Select the best template match.
    if ~isempty(candidate_indices)
        if max(candidate_sims) >= override_thresh, v_idx = find(candidate_sims >= override_thresh);
        else, v_idx = find(candidate_sims >= sim_thresh); end
        
        if ~isempty(v_idx)
             [~, m_idx] = min(abs(candidate_amps(v_idx) - wave_amp));
             real_ptr = v_idx(m_idx);
             final_idx = candidate_indices(real_ptr); final_sim = candidate_sims(real_ptr);
        else, final_idx = global_best_idx; final_sim = global_max_sim; 
        end
    else, final_idx = global_best_idx; final_sim = global_max_sim; 
    end
    
    % Decide whether the candidate is noise.
    if final_idx > 0
        best_tmpl_amp = max(noise_templates{final_idx});
        
        % Require similarity to reach sim_thresh (for example, 0.75).
        is_noise = (wave_amp > 10) && (final_sim >= sim_thresh);
        
    else
        best_tmpl_amp = NaN; 
        is_noise = false; 
    end
end
