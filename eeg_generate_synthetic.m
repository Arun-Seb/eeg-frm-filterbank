function [eeg_signal, t, true_components] = eeg_generate_synthetic(params)
%EEG_GENERATE_SYNTHETIC  Generate a realistic synthetic EEG signal.
%
%  The synthetic signal models realistic EEG by combining:
%
%  1. DELTA (0.5-4 Hz): Large amplitude, slow waves
%     - Simulates deep sleep slow-wave activity
%     - Amplitude: 50-100 µV (dominant in sleep)
%
%  2. THETA (4-8 Hz): Medium amplitude
%     - Simulates hippocampal theta during drowsiness/memory
%     - Amplitude: 20-40 µV
%
%  3. ALPHA (8-13 Hz): Prominent at rest (eyes closed)
%     - The classic "alpha rhythm" — 10 Hz dominant
%     - Amplitude: 30-50 µV (strong when relaxed)
%
%  4. BETA (13-30 Hz): Low amplitude, fast activity
%     - Active thinking, concentration
%     - Amplitude: 5-20 µV
%
%  5. GAMMA (30-60 Hz): Very low amplitude
%     - High-level cognitive processing
%     - Amplitude: 2-8 µV
%
%  6. 1/f BACKGROUND NOISE: Pink noise (realistic EEG background)
%     - Real EEG has a 1/f power spectrum
%     - Simulated by filtering white noise
%
%  7. LINE NOISE (50 Hz): Power line interference
%     - Common artifact in EEG recordings
%     - Included to test filter performance
%
%  AMPLITUDE NOTE:
%  EEG amplitudes are in microvolts (µV). The signal is normalized
%  to unit variance for filter bank processing, then scaled back.
%
%  Inputs:
%    params - parameter struct
%
%  Outputs:
%    eeg_signal      - combined synthetic EEG [1 x N] (µV)
%    t               - time vector (seconds)
%    true_components - struct with individual band signals (ground truth)

    fprintf('--- Generating Synthetic EEG Signal ---\n');

    Fs       = params.Fs;
    duration = params.duration;
    N        = Fs * duration;
    t        = (0 : N-1) / Fs;

    rng(42);   % fixed seed for reproducibility

    %% ---- Delta (0.5 - 4 Hz) ----
    % Multiple harmonics for realism
    delta = 80 * sin(2*pi*1.0*t) + ...
            50 * sin(2*pi*2.0*t + 0.5) + ...
            30 * sin(2*pi*3.5*t + 1.2);

    %% ---- Theta (4 - 8 Hz) ----
    theta = 35 * sin(2*pi*5.0*t) + ...
            20 * sin(2*pi*6.5*t + 0.8) + ...
            15 * sin(2*pi*7.5*t + 0.3);

    %% ---- Alpha (8 - 13 Hz) ----
    % Alpha is typically amplitude-modulated (waxes and wanes)
    alpha_carrier = sin(2*pi*10.0*t) + ...
                    0.6*sin(2*pi*9.5*t + 0.4) + ...
                    0.4*sin(2*pi*11.0*t + 0.2);
    % Amplitude modulation: slow variation
    alpha_env   = 1 + 0.4 * sin(2*pi*0.3*t);
    alpha       = 45 * alpha_env .* alpha_carrier;

    %% ---- Beta (13 - 30 Hz) ----
    beta = 15 * sin(2*pi*18.0*t) + ...
           10 * sin(2*pi*22.0*t + 0.6) + ...
            8 * sin(2*pi*26.0*t + 1.0) + ...
            5 * sin(2*pi*28.0*t + 0.2);

    %% ---- Gamma (30 - 60 Hz) ----
    gamma_comp = 6 * sin(2*pi*35.0*t) + ...
                 4 * sin(2*pi*42.0*t + 0.7) + ...
                 3 * sin(2*pi*55.0*t + 0.4);

    %% ---- 1/f Pink Background Noise ----
    % Real EEG has power spectrum ~1/f^alpha (alpha ~1-2)
    % Generate by filtering white noise
    white_noise = randn(1, N);
    % Simple 1/f approximation via integration (cumsum)
    pink_noise  = cumsum(white_noise);
    pink_noise  = pink_noise - mean(pink_noise);
    pink_noise  = pink_noise / std(pink_noise) * 15;   % scale to ~15 µV std

    %% ---- 50 Hz Line Noise (power line artifact) ----
    line_noise = 8 * sin(2*pi*50.0*t);

    %% ---- Combine all components ----
    eeg_signal = delta + theta + alpha + beta + gamma_comp + pink_noise + line_noise;

    %% ---- Store true components for validation ----
    true_components.delta       = delta;
    true_components.theta       = theta;
    true_components.alpha       = alpha;
    true_components.beta        = beta;
    true_components.gamma       = gamma_comp;
    true_components.pink_noise  = pink_noise;
    true_components.line_noise  = line_noise;

    %% ---- Report signal statistics ----
    fprintf('  Duration          : %d seconds\n', duration);
    fprintf('  Samples           : %d\n', N);
    fprintf('  Total RMS (µV)    : %.2f\n', rms(eeg_signal));
    fprintf('  Component RMS (µV):\n');
    fprintf('    Delta           : %.2f\n', rms(delta));
    fprintf('    Theta           : %.2f\n', rms(theta));
    fprintf('    Alpha           : %.2f\n', rms(alpha));
    fprintf('    Beta            : %.2f\n', rms(beta));
    fprintf('    Gamma           : %.2f\n', rms(gamma_comp));
    fprintf('    Pink noise      : %.2f\n', rms(pink_noise));
    fprintf('    Line noise      : %.2f\n\n', rms(line_noise));

    %% ---- Plot synthetic EEG ----
    fig = figure('Name', 'Synthetic EEG Signal', 'Position', [50, 100, 1200, 700]);

    % Full signal
    subplot(4, 2, [1 2]);
    plot(t, eeg_signal, 'k', 'LineWidth', 0.8);
    xlabel('Time (s)'); ylabel('Amplitude (µV)');
    title('Synthetic EEG Signal (all components combined)', 'FontSize', 11);
    grid on; xlim([0, duration]);

    % Individual components
    comp_names  = {'Delta (1-3.5 Hz)', 'Theta (5-7.5 Hz)', ...
                   'Alpha (9.5-11 Hz)', 'Beta (18-28 Hz)', ...
                   'Gamma (35-55 Hz)', 'Pink + Line Noise'};
    comp_data   = {delta, theta, alpha, beta, gamma_comp, pink_noise + line_noise};
    comp_colors = {[0.2,0.4,0.9],[0.1,0.7,0.3],[0.9,0.7,0.1], ...
                   [0.9,0.3,0.2],[0.7,0.1,0.7],[0.5,0.5,0.5]};

    for k = 1 : 6
        subplot(4, 2, k+2);
        plot(t, comp_data{k}, 'Color', comp_colors{k}, 'LineWidth', 1);
        xlabel('Time (s)'); ylabel('µV');
        title(comp_names{k}, 'FontSize', 9);
        grid on; xlim([0, min(3, duration)]);  % show first 3s for clarity
    end

    sgtitle('Synthetic EEG — Individual Components (Ground Truth)', 'FontSize', 13);
    saveas(fig, 'EEG_Synthetic_Signal.png');
    fprintf('  Synthetic EEG plot saved.\n\n');
end
