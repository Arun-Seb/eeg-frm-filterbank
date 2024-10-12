function eeg_plot_filterbank(H_bank, freq_axis, params, band_map)
%EEG_PLOT_FILTERBANK  Visualise the EEG-tuned FRM filter bank response.
%
%  Shows:
%    1. All subband magnitude responses overlaid (linear + dB)
%    2. EEG band boundaries overlaid on filter bank
%    3. Individual EEG band filter responses

    fprintf('--- Plotting EEG Filter Bank ---\n');

    Fs        = params.Fs;
    NUM_BANDS = params.NUM_BANDS;

    % Limit display to EEG-relevant range
    f_max_disp = min(80, Fs/2);
    f_idx      = freq_axis <= f_max_disp;
    f_disp     = freq_axis(f_idx);

    colors     = lines(NUM_BANDS);

    %% ---- Plot 1: All subbands — linear and dB ----
    fig1 = figure('Name', 'EEG FRM Filter Bank Response', 'Position', [50, 100, 1200, 600]);

    subplot(1,2,1);
    hold on;
    for i = 1 : NUM_BANDS
        plot(f_disp, H_bank(i, f_idx), 'Color', colors(i,:), 'LineWidth', 1.8);
    end
    % Add EEG band boundary lines
    boundaries = [4, 8, 13, 30];
    b_labels   = {'δ/θ  ','θ/α  ','α/β  ','β/γ  '};
    for b = 1 : length(boundaries)
        if boundaries(b) <= f_max_disp
            xline(boundaries(b), 'k--', b_labels{b}, 'LineWidth', 1.2, 'FontSize', 9);
        end
    end
    hold off;
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 11);
    ylabel('Magnitude', 'FontSize', 11);
    title('FRM Filter Bank — Linear Magnitude', 'FontSize', 11);
    xlim([0, f_max_disp]); ylim([0, 1.3]);
    legend(arrayfun(@(i) sprintf('Sub-%d',i), 1:NUM_BANDS, 'UniformOutput',false), ...
           'Location','eastoutside','FontSize',7);

    subplot(1,2,2);
    hold on;
    for i = 1 : NUM_BANDS
        mag_dB = 20*log10(max(H_bank(i, f_idx), 1e-10));
        plot(f_disp, mag_dB, 'Color', colors(i,:), 'LineWidth', 1.8);
    end
    for b = 1 : length(boundaries)
        if boundaries(b) <= f_max_disp
            xline(boundaries(b), 'k--', b_labels{b}, 'LineWidth', 1.2, 'FontSize', 9);
        end
    end
    yline(-80, 'r:', '−80 dB', 'LineWidth', 1);
    hold off;
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 11);
    ylabel('Magnitude (dB)', 'FontSize', 11);
    title('FRM Filter Bank — dB Magnitude', 'FontSize', 11);
    xlim([0, f_max_disp]); ylim([-100, 5]);

    sgtitle(sprintf('EEG-Tuned %d-Band FRM Non-Uniform Filter Bank (Fs=%d Hz)', ...
            NUM_BANDS, Fs), 'FontSize', 13);
    saveas(fig1, 'EEG_FilterBank_Response.png');

    %% ---- Plot 2: EEG band overlay ----
    fig2 = figure('Name', 'EEG Bands on Filter Bank', 'Position', [100, 100, 1100, 550]);

    eeg_bands  = {'Delta','Theta','Alpha','Beta','Gamma'};
    eeg_ranges = {[0.5,4],[4,8],[8,13],[13,30],[30,60]};
    eeg_colors = {[0.2,0.4,0.9,0.15],[0.1,0.7,0.3,0.15],[0.9,0.7,0.1,0.15], ...
                  [0.9,0.3,0.2,0.15],[0.7,0.1,0.7,0.15]};

    subplot(2,1,1);
    hold on;
    % Shade EEG bands first
    for b = 1 : length(eeg_bands)
        if eeg_ranges{b}(1) <= f_max_disp
            c = eeg_colors{b};
            patch([eeg_ranges{b}(1), eeg_ranges{b}(2), eeg_ranges{b}(2), eeg_ranges{b}(1)], ...
                  [0, 0, 1.3, 1.3], c(1:3), 'FaceAlpha', 0.2, ...
                  'EdgeColor', c(1:3), 'LineWidth', 1, 'DisplayName', eeg_bands{b});
        end
    end
    % Then plot subbands
    for i = 1 : NUM_BANDS
        plot(f_disp, H_bank(i, f_idx), 'Color', [colors(i,:), 0.8], 'LineWidth', 1.5, ...
             'HandleVisibility', 'off');
    end
    hold off;
    xlabel('Frequency (Hz)'); ylabel('Magnitude');
    title('FRM Filter Subbands with EEG Band Overlay', 'FontSize', 11);
    xlim([0, f_max_disp]); ylim([0, 1.3]); grid on;
    legend('Location','eastoutside','FontSize',9);

    subplot(2,1,2);
    % Show which subbands contribute to each EEG band
    contribution = zeros(length(eeg_bands), NUM_BANDS);
    for b = 1 : length(eeg_bands)
        f_idx_band = (freq_axis >= eeg_ranges{b}(1)) & (freq_axis <= eeg_ranges{b}(2));
        for i = 1 : NUM_BANDS
            if any(f_idx_band)
                contribution(b,i) = sum(H_bank(i, f_idx_band).^2);
            end
        end
        % Normalize per band
        contribution(b,:) = contribution(b,:) / max(sum(contribution(b,:)), 1e-10);
    end

    imagesc(contribution);
    colormap(hot);
    colorbar;
    xlabel('Subband Index'); ylabel('EEG Band');
    set(gca, 'YTickLabel', eeg_bands, 'YTick', 1:length(eeg_bands));
    title('Subband Contribution to Each EEG Band (normalized energy)');
    xticks(1:NUM_BANDS);

    sgtitle('EEG Band — FRM Subband Mapping', 'FontSize', 13);
    saveas(fig2, 'EEG_Band_Subband_Mapping.png');

    fprintf('  Filter bank plots saved.\n\n');
end
