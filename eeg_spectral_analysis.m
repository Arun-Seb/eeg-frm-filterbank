function eeg_spectral_analysis(eeg_signal, extracted, band_signals, t, params, band_map)
%EEG_SPECTRAL_ANALYSIS  Power spectral density and spectrogram analysis.
%
%  Produces:
%    1. PSD of raw EEG vs extracted bands (Welch's method)
%    2. Spectrogram of raw EEG (time-frequency representation)
%    3. Band power over time (sliding window)
%    4. Topographic-style power bar chart

    fprintf('--- Spectral Analysis ---\n');

    Fs        = params.Fs;
    N         = length(eeg_signal);
    NFFT      = params.NFFT;

    band_names  = {'delta','theta','alpha','beta','gamma'};
    band_labels = {'Delta\newline0.5-4 Hz','Theta\newline4-8 Hz', ...
                   'Alpha\newline8-13 Hz','Beta\newline13-30 Hz','Gamma\newline30-60 Hz'};
    band_ranges = {[0.5,4],[4,8],[8,13],[13,30],[30,60]};
    band_colors = {[0.2,0.4,0.9],[0.1,0.7,0.3],[0.9,0.7,0.1], ...
                   [0.9,0.3,0.2],[0.7,0.1,0.7]};

    %% ---- Plot 1: Power Spectral Density ----
    fig1 = figure('Name', 'EEG Power Spectral Density', 'Position', [50, 100, 1100, 600]);

    % Welch PSD of raw EEG
    win_len    = round(Fs * 2);   % 2-second window
    win_len    = min(win_len, N);
    [pxx, f]  = pwelch(eeg_signal, hann(win_len), round(win_len/2), NFFT, Fs);

    subplot(2, 1, 1);
    semilogy(f, pxx, 'k', 'LineWidth', 1.5);
    hold on;
    % Shade each EEG band
    for b = 1 : length(band_names)
        idx = (f >= band_ranges{b}(1)) & (f <= band_ranges{b}(2));
        if any(idx)
            fill([f(idx); flipud(f(idx))], ...
                 [pxx(idx); ones(sum(idx),1)*min(pxx)*0.5], ...
                 band_colors{b}, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
                 'DisplayName', band_names{b});
        end
    end
    hold off;
    xlabel('Frequency (Hz)'); ylabel('PSD (µV²/Hz)');
    title('Raw EEG Power Spectral Density (Welch)', 'FontSize', 11);
    xlim([0, min(80, Fs/2)]); grid on;
    legend(band_names, 'Location', 'northeast', 'FontSize', 8);

    % PSD comparison: raw vs extracted bands
    subplot(2, 1, 2);
    semilogy(f, pxx, 'k-', 'LineWidth', 2, 'DisplayName', 'Raw EEG');
    hold on;
    for b = 1 : length(band_names)
        sig = extracted.(band_names{b});
        [pxx_b, ~] = pwelch(sig, hann(win_len), round(win_len/2), NFFT, Fs);
        semilogy(f, pxx_b, 'Color', band_colors{b}, 'LineWidth', 1.5, ...
                 'DisplayName', band_names{b});
    end
    hold off;
    xlabel('Frequency (Hz)'); ylabel('PSD (µV²/Hz)');
    title('PSD: Raw EEG vs Extracted Bands', 'FontSize', 11);
    xlim([0, min(80, Fs/2)]); grid on;
    legend('Location', 'northeast', 'FontSize', 8);

    sgtitle('EEG Power Spectral Density Analysis', 'FontSize', 13);
    saveas(fig1, 'EEG_PSD_Analysis.png');

    %% ---- Plot 2: Spectrogram ----
    fig2 = figure('Name', 'EEG Spectrogram', 'Position', [100, 100, 1100, 500]);

    win_spec = hann(round(Fs * 1));    % 1-second window
    noverlap = round(length(win_spec) * 0.75);

    spectrogram(eeg_signal, win_spec, noverlap, NFFT, Fs, 'yaxis');
    ylim([0, min(80, Fs/2)]);
    colormap(jet);
    colorbar;
    title('EEG Spectrogram (Short-Time Fourier Transform)', 'FontSize', 12);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    % Add horizontal lines for band boundaries
    hold on;
    boundaries = [4, 8, 13, 30];
    b_labels   = {'  δ/θ','  θ/α','  α/β','  β/γ'};
    for b = 1 : length(boundaries)
        yline(boundaries(b), 'w--', b_labels{b}, 'LineWidth', 1.5, 'FontSize', 10);
    end
    hold off;

    saveas(fig2, 'EEG_Spectrogram.png');

    %% ---- Plot 3: Band power over time (sliding window) ----
    fig3 = figure('Name', 'EEG Band Power Over Time', 'Position', [150, 100, 1100, 600]);

    win_s   = round(Fs * 1.0);    % 1-second window
    step_s  = round(Fs * 0.25);   % 250 ms step (75% overlap)
    starts  = 1 : step_s : N - win_s + 1;
    t_win   = t(starts + round(win_s/2));   % center time of each window

    band_power = zeros(length(band_names), length(starts));

    for b = 1 : length(band_names)
        sig = extracted.(band_names{b});
        for w = 1 : length(starts)
            idx = starts(w) : starts(w) + win_s - 1;
            idx = min(idx, length(sig));
            band_power(b, w) = sum(sig(idx).^2) / win_s;
        end
    end

    % Normalize each band's power to 0-1 for comparison
    for b = 1 : length(band_names)
        bp = band_power(b,:);
        band_power(b,:) = (bp - min(bp)) / max(max(bp) - min(bp), 1e-10);
    end

    subplot(2,1,1);
    hold on;
    for b = 1 : length(band_names)
        plot(t_win, band_power(b,:), 'Color', band_colors{b}, ...
             'LineWidth', 1.5, 'DisplayName', band_names{b});
    end
    hold off;
    xlabel('Time (s)'); ylabel('Normalized Power');
    title('EEG Band Power Over Time (1-second sliding window, 75% overlap)');
    legend('Location', 'eastoutside', 'FontSize', 9); grid on;
    xlim([0, params.duration]);

    subplot(2,1,2);
    area(t_win, band_power', 'LineWidth', 0.5);
    colororder(cell2mat(band_colors'));
    xlabel('Time (s)'); ylabel('Relative Power (stacked)');
    title('Stacked Band Power Over Time');
    legend(band_names, 'Location', 'eastoutside', 'FontSize', 9); grid on;
    xlim([0, params.duration]);

    sgtitle('EEG Band Power Dynamics', 'FontSize', 13);
    saveas(fig3, 'EEG_Band_Power_Time.png');

    fprintf('  Spectral analysis plots saved.\n\n');
end
