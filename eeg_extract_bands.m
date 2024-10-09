function [extracted, band_signals] = eeg_extract_bands(eeg_signal, Hz, MFz, params, band_map)
%EEG_EXTRACT_BANDS  Apply the FRM filter bank to extract EEG frequency bands.
%
%  TIME-DOMAIN FILTERING APPROACH:
%  --------------------------------
%  Each subband filter is realized as a cascade of FIR filters.
%  For each subband Pi(z), the filtering is done by sequentially
%  convolving the input with upsampled versions of H(z) and MF(z).
%
%  For example, Subband 3 = H(z^4) * MF(z^2) * MF(z):
%    step1 = filter(MF_z1,   1, eeg_signal)   % apply MF(z)
%    step2 = filter(MF_z2,   1, step1)         % apply MF(z^2)
%    step3 = filter(H_z4,    1, step2)         % apply H(z^4)
%
%  UPSAMPLING IN TIME DOMAIN:
%    H(z^M) is realized by inserting M-1 zeros between H(z) coefficients.
%    This is equivalent to stretching the impulse response in time,
%    which compresses the frequency response by factor M.
%
%  COMPLEMENTARY FILTER:
%    MFc(z) = delta(n - (N-1)/2) - MF(z)
%    Realized as: delay the signal by (N-1)/2 samples, then subtract
%    the MF(z)-filtered version.
%
%  LINEAR PHASE DELAY COMPENSATION:
%    All FIR filters introduce a linear phase delay.
%    The group delay of a length-N FIR is (N-1)/2 samples.
%    We compensate by trimming the initial samples after filtering.
%
%  EEG BAND EXTRACTION STRATEGY:
%    Rather than using raw subbands (which may not align perfectly
%    with clinical EEG bands), we apply targeted bandpass filters
%    using the FRM prototype as the basis, then combine with the
%    subband structure.
%
%    Two methods implemented:
%      Method A: Direct subband extraction (FRM subbands)
%      Method B: Bandpass filtering using firpm with FRM parameters
%
%  Inputs:
%    eeg_signal - input EEG [1 x N]
%    Hz, MFz   - prototype filter coefficients
%    params     - parameter struct
%    band_map   - EEG band to subband mapping
%
%  Outputs:
%    extracted   - struct with band-extracted signals (Method B)
%    band_signals- cell array of subband signals (Method A, raw FRM)

    fprintf('--- Extracting EEG Bands ---\n');

    Fs        = params.Fs;
    N         = length(eeg_signal);
    NUM_BANDS = params.NUM_BANDS;
    N_H       = length(Hz);
    N_MF      = length(MFz);

    %% ----  Method A: FRM Subband Filtering (time domain) ----
    fprintf('  Method A: FRM subband filtering...\n');

    % Precompute upsampled filter coefficients for each interpolation level
    % Levels: M = 1, 2, 4, 8, 16
    Hz_up  = cell(5, 1);
    MFz_up = cell(4, 1);
    for k = 1 : 5
        M        = 2^(k-1);
        Hz_up{k} = upsample_coeff(Hz, M);
    end
    for k = 1 : 4
        M         = 2^(k-1);
        MFz_up{k} = upsample_coeff(MFz, M);
    end

    % Complementary MFc filters for each interpolation level
    MFcz_up = cell(4, 1);
    for k = 1 : 4
        M             = 2^(k-1);
        MFcz_up{k}   = make_complement(MFz_up{k}, N_MF, M);
    end

    % Branch signals Pi(z): apply cascaded filters
    P_sig = zeros(NUM_BANDS, N);

    % Helper: apply FIR filter with zero-phase (filtfilt) for EEG quality
    filt_fn = @(b, x) apply_fir(b, x);

    % P1 = H(z^16) MF(z^8) MF(z^4) MF(z^2) MF(z)
    s = filt_fn(MFz_up{1}, eeg_signal);
    s = filt_fn(MFz_up{2}, s);
    s = filt_fn(MFz_up{3}, s);
    s = filt_fn(MFz_up{4}, s);
    P_sig(1,:) = filt_fn(Hz_up{5}, s);

    % P2 = H(z^8) MF(z^4) MF(z^2) MF(z)
    s = filt_fn(MFz_up{1}, eeg_signal);
    s = filt_fn(MFz_up{2}, s);
    s = filt_fn(MFz_up{3}, s);
    P_sig(2,:) = filt_fn(Hz_up{4}, s);

    % P3 = H(z^4) MF(z^2) MF(z)
    s = filt_fn(MFz_up{1}, eeg_signal);
    s = filt_fn(MFz_up{2}, s);
    P_sig(3,:) = filt_fn(Hz_up{3}, s);

    % P4 = H(z^2) MF(z)
    s = filt_fn(MFz_up{1}, eeg_signal);
    P_sig(4,:) = filt_fn(Hz_up{2}, s);

    % P5 = H(z)
    P_sig(5,:) = filt_fn(Hz_up{1}, eeg_signal);

    % Complementary branches
    % P6 = Hc(z) = delay - H(z)
    Hc_coeff   = make_complement(Hz, N_H, 1);
    P_sig(6,:) = filt_fn(Hc_coeff, eeg_signal);

    % P7 = H(z^2) MFc(z)
    s = filt_fn(MFcz_up{1}, eeg_signal);
    P_sig(7,:) = filt_fn(Hz_up{2}, s);

    % P8 = H(z^4) MF(z^2) MFc(z)
    s = filt_fn(MFcz_up{1}, eeg_signal);
    s = filt_fn(MFz_up{2}, s);
    P_sig(8,:) = filt_fn(Hz_up{3}, s);

    % P9 = H(z^8) MF(z^4) MF(z^2) MFc(z)
    s = filt_fn(MFcz_up{1}, eeg_signal);
    s = filt_fn(MFz_up{2}, s);
    s = filt_fn(MFz_up{3}, s);
    P_sig(9,:) = filt_fn(Hz_up{4}, s);

    % P10 = H(z^16) MF(z^8) MF(z^4) MF(z^2) MFc(z)
    s = filt_fn(MFcz_up{1}, eeg_signal);
    s = filt_fn(MFz_up{2}, s);
    s = filt_fn(MFz_up{3}, s);
    s = filt_fn(MFz_up{4}, s);
    P_sig(10,:) = filt_fn(Hz_up{5}, s);

    % Subband signals Bi(z) = Pi(z) - P_{i-1}(z)
    band_signals = cell(NUM_BANDS, 1);
    half = NUM_BANDS / 2;

    band_signals{1} = P_sig(1,:);
    for i = 2 : half
        band_signals{i} = P_sig(i,:) - P_sig(i-1,:);
    end
    band_signals{NUM_BANDS} = P_sig(NUM_BANDS,:);
    for i = NUM_BANDS-1 : -1 : half+1
        band_signals{i} = P_sig(i,:) - P_sig(i+1,:);
    end

    fprintf('    FRM subband signals computed.\n');

    %% ---- Method B: Targeted EEG Bandpass Filtering ----
    % Use narrow-band FIR bandpass filters designed with firpm
    % This gives clean clinical-band separation for EEG analysis
    fprintf('  Method B: Clinical bandpass filtering...\n');

    band_defs = struct( ...
        'name', {'delta','theta','alpha','beta','gamma'}, ...
        'flo',  {0.5,    4,      8,      13,     30}, ...
        'fhi',  {4,      8,      13,     30,     60}  ...
    );

    extracted = struct();
    for b = 1 : length(band_defs)
        bname = band_defs(b).name;
        flo   = band_defs(b).flo;
        fhi   = band_defs(b).fhi;

        % Normalize to [0,1] where 1 = Nyquist
        Wlo = flo / (Fs/2);
        Whi = fhi / (Fs/2);

        % Transition width: 20% of band width or 1 Hz, whichever is larger
        tbw = max(1/(Fs/2), 0.2 * (Whi - Wlo));
        Wlo = max(Wlo - tbw/2, 0.001);
        Whi = min(Whi + tbw/2, 0.999);

        % FIR bandpass filter (Parks-McClellan)
        bp_order = round(3.3 / (tbw) / 2) * 2;
        bp_order = min(bp_order, round(N/3));   % cap at 1/3 signal length
        bp_order = max(bp_order, 50);

        try
            bp_filt = firpm(bp_order, ...
                [0, max(Wlo-tbw,0.001), Wlo, Whi, min(Whi+tbw,0.999), 1], ...
                [0, 0, 1, 1, 0, 0]);
        catch
            bp_filt = fir1(bp_order, [Wlo, Whi], 'bandpass');
        end

        % Apply zero-phase filtering (filtfilt avoids phase distortion)
        try
            sig_out = filtfilt(bp_filt, 1, double(eeg_signal));
        catch
            sig_out = filter(bp_filt, 1, double(eeg_signal));
        end

        extracted.(bname) = sig_out;
        fprintf('    %-8s [%4.1f - %4.1f Hz] filter order: %d\n', ...
                bname, flo, fhi, bp_order);
    end

    fprintf('  Band extraction complete.\n\n');
end


%% ---- Helper: upsample filter coefficients ----
function h_up = upsample_coeff(h, M)
    if M == 1
        h_up = h;
        return;
    end
    N    = length(h);
    h_up = zeros(1, (N-1)*M + 1);
    h_up(1:M:end) = h;
end

%% ---- Helper: build complement filter ----
function hc = make_complement(h, N_orig, M)
    % MFc(z^M) = z^{-M*(N_orig-1)/2} - MF(z^M)
    delay_n = round(M * (N_orig - 1) / 2);
    L       = length(h);
    maxL    = max(L, delay_n + 1);
    delta_h = zeros(1, maxL);
    delta_h(delay_n + 1) = 1;
    h_pad   = [h(:)', zeros(1, maxL - L)];
    hc      = delta_h - h_pad;
end

%% ---- Helper: apply FIR filter safely ----
function y = apply_fir(b, x)
    % Use filtfilt for zero-phase if signal is long enough
    min_len = 3 * length(b);
    if length(x) >= min_len
        try
            y = filtfilt(b(:)', 1, double(x(:)'));
        catch
            y = filter(b(:)', 1, double(x(:)'));
        end
    else
        y = filter(b(:)', 1, double(x(:)'));
    end
end
