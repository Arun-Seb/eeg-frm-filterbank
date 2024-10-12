function eeg_analyse_results(eeg_signal, t, extracted, band_signals, true_components, params, band_map)
%EEG_ANALYSE_RESULTS  Plot and compare extracted EEG bands vs ground truth.
%
%  Produces:
%    1. Multi-panel plot: raw EEG + all extracted bands (stacked)
%    2. Band-by-band comparison: extracted vs ground truth
%    3. Correlation and SNR metrics table

    fprintf('--- Analysing Extraction Results ---\n');

    Fs       = params.Fs;
    duration = params.duration;
    N        = length(eeg_signal);

    %% ---- Plot 1: Stacked EEG decomposition ----
    band_names  = {'delta','theta','alpha','beta','gamma'};
    band_labels = {'Delta (0.5-4 Hz)','Theta (4-8 Hz)', ...
                   'Alpha (8-13 Hz)','Beta (13-30 Hz)','Gamma (30-60 Hz)'};
    band_colors = {[0.2,0.4,0.9],[0.1,0.7,0.3],[0.9,0.7,0.1], ...
                   [0.9,0.3,0.2],[0.7,0.1,0.7]};
    truth_fields = {'delta','theta','alpha','beta','gamma'};

    fig1 = figure('Name', 'EEG Band Decomposition', 'Position', [50, 50, 1200, 900]);

    n_rows = length(band_names) + 1;

    % Raw EEG
    subplot(n_rows, 1, 1);
    plot(t, eeg_signal, 'k', 'LineWidth', 0.7);
    ylabel('µV'); title('Raw Synthetic EEG', 'FontSize', 10);
    grid on; xlim([0, duration]);
    set(gca, 'XTickLabel', []);

    % Each extracted band
    for b = 1 : length(band_names)
        subplot(n_rows, 1, b+1);
        sig = extracted.(band_names{b});
        plot(t, sig, 'Color', band_colors{b}, 'LineWidth', 1.2);
        ylabel('µV');
        title(sprintf('Extracted: %s', band_labels{b}), 'FontSize', 9);
        grid on; xlim([0, duration]);
        if b < length(band_names)
            set(gca, 'XTickLabel', []);
        else
            xlabel('Time (s)');
        end
    end

    sgtitle('EEG Signal Decomposition — FRM Filter Bank', 'FontSize', 13);
    saveas(fig1, 'EEG_Band_Decomposition.png');

    %% ---- Plot 2: Extracted vs Ground Truth (first 5 seconds) ----
    fig2 = figure('Name', 'Extracted vs Ground Truth', 'Position', [100, 50, 1300, 900]);
    t_show = t(t <= 5);   % show first 5 seconds
    idx5   = t <= 5;

    for b = 1 : length(band_names)
        subplot(length(band_names), 1, b);
        truth = true_components.(truth_fields{b});

        % Normalize both to same scale for visual comparison
        ext_sig = extracted.(band_names{b});

        hold on;
        plot(t_show, truth(idx5),    'b-',  'LineWidth', 1.8, 'DisplayName', 'Ground Truth');
        plot(t_show, ext_sig(idx5),  'r--', 'LineWidth', 1.4, 'DisplayName', 'FRM Extracted');
        hold off;

        ylabel('µV');
        title(band_labels{b}, 'FontSize', 9);
        grid on;
        if b == 1
            legend('Location', 'northeast', 'FontSize', 8);
        end
        if b < length(band_names)
            set(gca, 'XTickLabel', []);
        else
            xlabel('Time (s)');
        end
    end

    sgtitle('Extracted Bands vs Ground Truth (first 5 seconds)', 'FontSize', 13);
    saveas(fig2, 'EEG_Extracted_vs_Truth.png');

    %% ---- Compute and report quality metrics ----
    fprintf('\n  %-12s | %-10s | %-10s | %-10s | %-10s\n', ...
            'Band', 'Corr (r)', 'SNR (dB)', 'RMSE (µV)', 'RMS ratio');
    fprintf('  %s\n', repmat('-', 1, 62));

    metrics = struct();
    for b = 1 : length(band_names)
        bname = band_names{b};
        truth = true_components.(truth_fields{b});
        ext   = extracted.(bname);

        % Trim to same length
        L   = min(length(truth), length(ext));
        tr  = truth(1:L);
        ex  = ext(1:L);

        % Pearson correlation
        cc  = corrcoef(tr, ex);
        r   = cc(1,2);

        % SNR: signal power / noise power
        noise    = ex - tr;
        snr_val  = 10 * log10(sum(tr.^2) / max(sum(noise.^2), 1e-10));

        % RMSE
        rmse_val = sqrt(mean((tr - ex).^2));

        % RMS ratio (should be ~1 for good extraction)
        rms_ratio = rms(ex) / max(rms(tr), 1e-10);

        metrics.(bname).r         = r;
        metrics.(bname).snr       = snr_val;
        metrics.(bname).rmse      = rmse_val;
        metrics.(bname).rms_ratio = rms_ratio;

        fprintf('  %-12s | %10.4f | %10.2f | %10.4f | %10.4f\n', ...
                bname, r, snr_val, rmse_val, rms_ratio);
    end

    %% ---- Plot metrics summary ----
    fig3 = figure('Name', 'Extraction Quality Metrics', 'Position', [200, 100, 1000, 450]);

    r_vals   = cellfun(@(b) metrics.(b).r,    band_names);
    snr_vals = cellfun(@(b) metrics.(b).snr,  band_names);

    subplot(1,2,1);
    bar(r_vals, 'FaceColor', [0.2, 0.5, 0.8]);
    set(gca, 'XTickLabel', band_labels, 'XTickLabelRotation', 25);
    ylabel('Pearson Correlation (r)');
    title('Extraction Correlation vs Ground Truth');
    ylim([0, 1.1]); grid on;
    yline(0.9, 'r--', 'r=0.9 threshold', 'LineWidth', 1.2);

    subplot(1,2,2);
    bar(snr_vals, 'FaceColor', [0.2, 0.7, 0.4]);
    set(gca, 'XTickLabel', band_labels, 'XTickLabelRotation', 25);
    ylabel('SNR (dB)');
    title('Signal-to-Noise Ratio of Extracted Bands');
    grid on;
    yline(10, 'r--', '10 dB threshold', 'LineWidth', 1.2);

    sgtitle('EEG Band Extraction Quality Metrics', 'FontSize', 13);
    saveas(fig3, 'EEG_Extraction_Metrics.png');

    fprintf('\n  Plots saved: EEG_Band_Decomposition.png, EEG_Extracted_vs_Truth.png\n\n');
end
