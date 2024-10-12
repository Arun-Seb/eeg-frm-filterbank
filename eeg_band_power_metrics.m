function eeg_band_power_metrics(band_signals, true_components, params, band_map)
%EEG_BAND_POWER_METRICS  Compute clinically relevant EEG power metrics.
%
%  Computes widely used EEG indices:
%
%  1. ABSOLUTE BAND POWER (µV²)
%     Total signal power in each frequency band.
%     Useful for comparing across conditions.
%
%  2. RELATIVE BAND POWER (%)
%     Each band's power as % of total power.
%     More robust to inter-subject amplitude differences.
%
%  3. CLINICALLY RELEVANT RATIOS:
%
%     a) Alpha/Beta Ratio
%        High ratio => relaxed state, low arousal
%        Low ratio  => stressed, cognitively active
%        Reference: Ratio > 1 suggests relaxed, < 0.5 suggests stress
%
%     b) Theta/Alpha Ratio
%        High ratio => drowsiness, early sleep onset
%        Used in fatigue monitoring (driver drowsiness detection)
%
%     c) Delta/Alpha Ratio
%        High ratio => deep sleep or pathological slowing
%        Used in sleep staging, encephalopathy detection
%
%     d) (Theta + Alpha) / (Alpha + Beta) — Engagement Index
%        Low => high cognitive engagement
%        High => low engagement / drowsiness
%
%     e) Gamma / (Alpha + Beta) — Gamma relative power
%        Elevated in certain cognitive tasks and pathologies
%
%  4. INDIVIDUAL ALPHA FREQUENCY (IAF)
%     Frequency of peak alpha power (usually 8-13 Hz).
%     Marker of cognitive status — slows with fatigue/aging.
%
%  5. SPECTRAL EDGE FREQUENCY (SEF95)
%     Frequency below which 95% of total spectral power lies.
%     Used in anesthesia monitoring.

    fprintf('--- EEG Band Power Metrics ---\n\n');

    Fs    = params.Fs;
    NFFT  = params.NFFT;

    band_names_ext = {'delta','theta','alpha','beta','gamma'};
    band_ranges    = {[0.5,4],[4,8],[8,13],[13,30],[30,60]};
    band_labels    = {'Delta','Theta','Alpha','Beta','Gamma'};
    band_colors    = {[0.2,0.4,0.9],[0.1,0.7,0.3],[0.9,0.7,0.1], ...
                      [0.9,0.3,0.2],[0.7,0.1,0.7]};

    %% ---- Get extracted band signals from band_signals (FRM subbands) ----
    % Use true_components as proxy for clean extracted bands
    % (In real use, replace with extracted.(band) from eeg_extract_bands)
    sigs = struct();
    sigs.delta = true_components.delta;
    sigs.theta = true_components.theta;
    sigs.alpha = true_components.alpha;
    sigs.beta  = true_components.beta;
    sigs.gamma = true_components.gamma;

    %% ---- Compute PSD for each band ----
    win_len = min(round(Fs*2), length(sigs.delta));
    [~, f]  = pwelch(sigs.delta, hann(win_len), round(win_len/2), NFFT, Fs);

    power   = zeros(1, length(band_names_ext));
    psd_all = zeros(length(band_names_ext), length(f));

    for b = 1 : length(band_names_ext)
        sig = sigs.(band_names_ext{b});
        [pxx, ~] = pwelch(sig, hann(win_len), round(win_len/2), NFFT, Fs);
        psd_all(b,:) = pxx;

        % Integrate PSD over band range to get power
        idx = (f >= band_ranges{b}(1)) & (f <= band_ranges{b}(2));
        if any(idx)
            power(b) = trapz(f(idx), pxx(idx));
        end
    end

    total_power = sum(power);
    rel_power   = 100 * power / total_power;

    %% ---- Print absolute and relative power ----
    fprintf('  Absolute Band Power:\n');
    for b = 1 : length(band_names_ext)
        fprintf('    %-8s : %10.4f µV²\n', band_labels{b}, power(b));
    end
    fprintf('    %-8s : %10.4f µV²\n', 'TOTAL', total_power);

    fprintf('\n  Relative Band Power:\n');
    for b = 1 : length(band_names_ext)
        fprintf('    %-8s : %6.2f %%\n', band_labels{b}, rel_power(b));
    end

    %% ---- Compute clinical ratios ----
    delta_p = power(1);  theta_p = power(2);
    alpha_p = power(3);  beta_p  = power(4);  gamma_p = power(5);

    ratios = struct();
    ratios.alpha_beta      = alpha_p / max(beta_p, 1e-10);
    ratios.theta_alpha     = theta_p / max(alpha_p, 1e-10);
    ratios.delta_alpha     = delta_p / max(alpha_p, 1e-10);
    ratios.engagement      = (theta_p + alpha_p) / max(alpha_p + beta_p, 1e-10);
    ratios.gamma_rel       = gamma_p / max(alpha_p + beta_p, 1e-10);

    fprintf('\n  Clinical EEG Ratios:\n');
    fprintf('    Alpha/Beta ratio          : %.4f\n', ratios.alpha_beta);
    fprintf('      > 1 => relaxed; < 0.5 => stressed/active\n');
    fprintf('    Theta/Alpha ratio         : %.4f\n', ratios.theta_alpha);
    fprintf('      High => drowsiness/fatigue\n');
    fprintf('    Delta/Alpha ratio         : %.4f\n', ratios.delta_alpha);
    fprintf('      High => sleep or pathological slowing\n');
    fprintf('    Engagement index          : %.4f\n', ratios.engagement);
    fprintf('      (Theta+Alpha)/(Alpha+Beta); Low => engaged\n');
    fprintf('    Gamma relative power      : %.4f\n', ratios.gamma_rel);

    %% ---- Individual Alpha Frequency (IAF) ----
    [pxx_full, f_full] = pwelch(sigs.alpha, hann(win_len), round(win_len/2), NFFT, Fs);
    alpha_idx = (f_full >= 8) & (f_full <= 13);
    [~, iaf_idx] = max(pxx_full(alpha_idx));
    f_alpha = f_full(alpha_idx);
    iaf = f_alpha(iaf_idx);
    fprintf('\n  Individual Alpha Frequency (IAF): %.2f Hz\n', iaf);
    fprintf('    Normal IAF: 9-11 Hz; slows with fatigue/aging\n');

    %% ---- Spectral Edge Frequency (SEF95) ----
    [pxx_raw, f_raw] = pwelch(sigs.delta + sigs.theta + sigs.alpha + sigs.beta + sigs.gamma, ...
                              hann(win_len), round(win_len/2), NFFT, Fs);
    cum_power = cumsum(pxx_raw);
    cum_power = cum_power / cum_power(end);
    sef95_idx = find(cum_power >= 0.95, 1, 'first');
    sef95     = f_raw(sef95_idx);
    fprintf('  Spectral Edge Frequency (SEF95)  : %.2f Hz\n', sef95);
    fprintf('    Used in anesthesia depth monitoring\n\n');

    %% ---- Plot 1: Absolute and relative power bar charts ----
    fig1 = figure('Name', 'EEG Band Power Analysis', 'Position', [50, 100, 1100, 500]);

    subplot(1,3,1);
    b_colors_mat = cell2mat(band_colors');
    bar_h = bar(power, 'FaceColor', 'flat');
    bar_h.CData = b_colors_mat;
    set(gca, 'XTickLabel', band_labels, 'XTickLabelRotation', 20);
    ylabel('Power (µV²)'); title('Absolute Band Power');
    grid on;

    subplot(1,3,2);
    pie(rel_power, band_labels);
    title('Relative Band Power (%)');
    colormap(b_colors_mat);

    subplot(1,3,3);
    ratio_names = {'α/β','θ/α','δ/α','Engage','γ-rel'};
    ratio_vals  = [ratios.alpha_beta, ratios.theta_alpha, ratios.delta_alpha, ...
                   ratios.engagement, ratios.gamma_rel];
    bar_r = bar(ratio_vals, 'FaceColor', [0.4, 0.6, 0.8]);
    set(gca, 'XTickLabel', ratio_names);
    ylabel('Ratio Value'); title('Clinical EEG Ratios');
    grid on;
    yline(1, 'r--', 'LineWidth', 1.2);

    sgtitle('EEG Band Power Metrics and Clinical Ratios', 'FontSize', 13);
    saveas(fig1, 'EEG_Band_Power_Metrics.png');

    %% ---- Plot 2: PSD overlay with band shading ----
    fig2 = figure('Name', 'EEG PSD with Band Shading', 'Position', [100, 100, 1000, 500]);

    hold on;
    for b = 1 : length(band_names_ext)
        idx = (f >= band_ranges{b}(1)) & (f <= band_ranges{b}(2));
        if any(idx)
            % Fill band area
            fill([f(idx); flipud(f(idx))], ...
                 [psd_all(b, idx)'; zeros(sum(idx),1)], ...
                 band_colors{b}, 'FaceAlpha', 0.5, 'EdgeColor', band_colors{b}, ...
                 'LineWidth', 1.5, 'DisplayName', band_labels{b});
        end
    end

    % IAF marker
    xline(iaf, 'k--', sprintf('IAF=%.1f Hz', iaf), 'LineWidth', 1.5, 'FontSize', 10);
    xline(sef95, 'm--', sprintf('SEF95=%.1f Hz', sef95), 'LineWidth', 1.5, 'FontSize', 10);
    hold off;

    xlabel('Frequency (Hz)'); ylabel('PSD (µV²/Hz)');
    title('EEG Power Spectral Density — Band Contributions', 'FontSize', 12);
    xlim([0, 70]); grid on;
    legend('Location', 'northeast', 'FontSize', 9);

    saveas(fig2, 'EEG_PSD_BandShaded.png');

    fprintf('  Band power metric plots saved.\n\n');
end
