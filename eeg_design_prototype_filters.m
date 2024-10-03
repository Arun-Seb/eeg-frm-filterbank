function [Hz, MFz, filter_info] = eeg_design_prototype_filters(params)
%EEG_DESIGN_PROTOTYPE_FILTERS  Design H(z) and MF(z) for EEG sampling rates.
%
%  KEY DIFFERENCE FROM HEARING AID VERSION:
%  ----------------------------------------
%  Hearing aid: Fs = 16000 Hz, bands span 0-8000 Hz
%  EEG:         Fs = 256 Hz,  bands span 0.5-60 Hz
%
%  The normalized transition bandwidth (norm_tbw) is the same concept,
%  but because EEG bands are much narrower relative to Fs, we need
%  tighter filters (smaller norm_tbw) => longer filter lengths.
%
%  Example:
%    Delta/Theta boundary: 4 Hz out of 128 Hz Nyquist
%    Normalized: 4/128 = 0.031 => very narrow => needs sharp filter
%
%  Half-band filter H(z):
%    - Cutoff at Fs/4 (Nyquist/2), i.e., 64 Hz for Fs=256 Hz
%    - All odd-indexed coefficients (except center) are zero
%    - Symmetric => linear phase
%    - Very sparse coefficient structure => low multiplier count
%
%  Masking filter MF(z):
%    - Relaxed transition band (does not need sharp cutoff)
%    - Short length (15 coefficients) keeps complexity low
%    - Cutoff chosen to shape the cascaded response correctly
%
%  FRM equation (Lim 1986):
%    Ha(z) = H(z^M) * MF(z) + Hc(z^M) * MFc(z)
%
%  Inputs:
%    params - parameter struct from main script
%
%  Outputs:
%    Hz          - H(z) filter coefficients
%    MFz         - MF(z) filter coefficients
%    filter_info - struct with design details

    fprintf('--- Designing EEG Prototype Filters ---\n');

    Fs       = params.Fs;
    norm_tbw = params.norm_tbw;
    MF_len   = params.MF_len;

    % ---- Compute H(z) length from transition bandwidth ----
    % Parks-McClellan estimate: N ≈ -10*log10(delta_p * delta_s) / (14.6 * delta_f)
    % Simplified: N ≈ 3.3 / norm_tbw (good approximation for equiripple)
    H_len = round(3.3 / norm_tbw / 2) * 2 + 1;   % force odd length
    H_len = max(H_len, 21);                         % minimum practical length

    fprintf('  Fs              : %d Hz\n', Fs);
    fprintf('  Nyquist         : %d Hz\n', Fs/2);
    fprintf('  Norm. TBW       : %.4f\n', norm_tbw);
    fprintf('  H(z) length     : %d taps\n', H_len);
    fprintf('  MF(z) length    : %d taps\n', MF_len);

    %% ---- Design H(z): Half-band lowpass filter ----
    % Half-band: passband edge at Wn=0.5*(1 - norm_tbw)
    %            stopband edge at Wn=0.5*(1 + norm_tbw)
    % Normalized to [0,1] where 1 = Nyquist
    Wc    = 0.5;                          % half-band center
    f_pb  = Wc - norm_tbw / 2;           % passband edge (normalized)
    f_sb  = Wc + norm_tbw / 2;           % stopband edge (normalized)

    f_pb  = max(f_pb, 0.01);             % numerical safety
    f_sb  = min(f_sb, 0.99);

    % Ripple specifications (80 dB stopband attenuation)
    delta_p = 0.01;     % passband ripple  (~0.09 dB)
    delta_s = 0.0001;   % stopband ripple  (~80 dB)

    try
        % firpm: optimal equiripple (Parks-McClellan)
        Hz = firpm(H_len - 1, ...
                   [0, f_pb, f_sb, 1], ...
                   [1, 1, 0, 0], ...
                   [1/delta_p, 1/delta_s]);
    catch ME
        fprintf('  firpm failed (%s), using Kaiser window design.\n', ME.message);
        % Kaiser window fallback
        beta  = 0.1102 * (80 - 8.7);   % Kaiser beta for 80 dB
        Hz    = fir1(H_len - 1, Wc, 'low', kaiser(H_len, beta));
    end

    % Enforce half-band symmetry: zero all odd-index coefficients except center
    Hz = enforce_halfband_symmetry(Hz);

    %% ---- Design MF(z): Masking filter ----
    % MF(z) shapes the final passband. Its transition band can be relaxed
    % because the interpolated H(z^M) provides the sharp band edges.
    % Cutoff at ~0.25 normalized frequency (quarter-band)
    f_pb_mf = 0.18;
    f_sb_mf = 0.32;

    try
        MFz = firpm(MF_len - 1, ...
                    [0, f_pb_mf, f_sb_mf, 1], ...
                    [1, 1, 0, 0]);
    catch
        MFz = fir1(MF_len - 1, 0.25, 'low', hamming(MF_len));
    end

    %% ---- Compute and report filter properties ----
    NFFT      = params.NFFT;
    H_resp    = freqz(Hz,  1, NFFT/2, 'whole');
    MF_resp   = freqz(MFz, 1, NFFT/2, 'whole');

    % Stopband attenuation
    freq_n    = (0 : NFFT/4 - 1) / (NFFT/2);
    sb_idx    = freq_n > (Wc + norm_tbw);
    if any(sb_idx)
        sb_atten = -20 * log10(max(abs(H_resp(sb_idx))));
    else
        sb_atten  = 0;
    end

    % Non-zero multiplier count (after half-band zeroing)
    Hz_nz  = sum(abs(Hz)  > 1e-8);
    MFz_nz = sum(abs(MFz) > 1e-8);
    n_mult = ceil(Hz_nz/2) + ceil(MFz_nz/2);   % symmetry halves count

    fprintf('  H(z) non-zero coeffs  : %d / %d\n', Hz_nz, H_len);
    fprintf('  MF(z) non-zero coeffs : %d / %d\n', MFz_nz, MF_len);
    fprintf('  Stopband attenuation  : %.1f dB\n', sb_atten);
    fprintf('  Estimated multipliers : ~%d (shared across bank)\n', n_mult);

    %% ---- Store filter info ----
    filter_info.H_len    = H_len;
    filter_info.MF_len   = MF_len;
    filter_info.Hz_nz    = Hz_nz;
    filter_info.MFz_nz   = MFz_nz;
    filter_info.sb_atten = sb_atten;
    filter_info.n_mult   = n_mult;
    filter_info.norm_tbw = norm_tbw;

    %% ---- Plot prototype filter responses ----
    fig = figure('Name', 'EEG Prototype Filters H(z) and MF(z)', ...
                 'Position', [50, 100, 1000, 500]);

    freq_Hz = (0 : NFFT/4 - 1) / (NFFT/2) * (Fs/2);

    subplot(2, 2, 1);
    plot(freq_Hz, abs(H_resp(1:NFFT/4)), 'b-', 'LineWidth', 2);
    xlabel('Frequency (Hz)'); ylabel('|H(z)|');
    title('H(z) — Magnitude (Linear)'); grid on;
    xlim([0, Fs/2]);

    subplot(2, 2, 2);
    plot(freq_Hz, 20*log10(max(abs(H_resp(1:NFFT/4)), 1e-10)), 'b-', 'LineWidth', 2);
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title('H(z) — Magnitude (dB)'); grid on;
    xlim([0, Fs/2]); ylim([-100, 5]);
    yline(-80, 'r--', '−80 dB', 'LineWidth', 1.2);

    subplot(2, 2, 3);
    plot(freq_Hz, abs(MF_resp(1:NFFT/4)), 'r-', 'LineWidth', 2);
    xlabel('Frequency (Hz)'); ylabel('|MF(z)|');
    title('MF(z) — Magnitude (Linear)'); grid on;
    xlim([0, Fs/2]);

    subplot(2, 2, 4);
    plot(freq_Hz, 20*log10(max(abs(MF_resp(1:NFFT/4)), 1e-10)), 'r-', 'LineWidth', 2);
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title('MF(z) — Magnitude (dB)'); grid on;
    xlim([0, Fs/2]); ylim([-100, 5]);

    sgtitle('EEG FRM Prototype Filters', 'FontSize', 13);
    saveas(fig, 'EEG_Prototype_Filters.png');
    fprintf('  Prototype filter plot saved.\n\n');
end


function h = enforce_halfband_symmetry(h)
%ENFORCE_HALFBAND_SYMMETRY  Zero odd-indexed coefficients except center.
%
%  In a true half-band FIR filter (cutoff at exactly pi/2):
%    h(n) = 0 for all n where |n - center| is odd  (0-based from center)
%  This halves the number of non-zero coefficients.
%
%  This property is exact for ideal half-band filters and approximately
%  holds for designed filters. We enforce it explicitly to reduce complexity.

    N   = length(h);
    ctr = ceil(N / 2);    % 1-based center index
    for k = 1 : N
        dist = abs(k - ctr);
        if mod(dist, 2) == 1
            h(k) = 0;
        end
    end
    % Re-normalize: ensure unity DC gain
    dc_gain = sum(h);
    if abs(dc_gain) > 1e-6
        h = h / dc_gain;
    end
end
