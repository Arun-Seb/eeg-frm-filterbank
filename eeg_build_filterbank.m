function [H_bank, freq_axis, band_map] = eeg_build_filterbank(Hz, MFz, params)
%EEG_BUILD_FILTERBANK  Build the FRM filter bank tuned for EEG bands.
%
%  FILTER BANK ARCHITECTURE:
%  --------------------------
%  The bank has NUM_BANDS subbands arranged symmetrically around
%  the half-band frequency (Fs/4 = 64 Hz for Fs=256).
%
%  Lower subbands (1..NUM_BANDS/2): original MF(z) path
%  Upper subbands (NUM_BANDS/2+1..NUM_BANDS): complementary MFc(z) path
%
%  Each subband Pi(z) is formed by cascading interpolated versions:
%    P_k = H(z^{2^k}) * MF(z^{2^{k-1}}) * ... * MF(z)
%
%  The complementary filter:
%    MFc(z) = z^{-(N-1)/2} - MF(z)    [Eq. 2 in Sebastian & James 2014]
%
%  BAND MAPPING FOR EEG:
%  ---------------------
%  The 10 non-uniform subbands map to EEG bands as follows
%  (for Fs=256 Hz, Nyquist=128 Hz):
%
%  Subband | Approx range (Hz) | EEG Band
%  --------|-------------------|----------
%    1     |  0 –  4           | Delta
%    2     |  4 –  8           | Theta
%    3     |  8 – 16           | Alpha + low Beta
%    4     | 16 – 32           | Beta
%    5     | 32 – 64           | Low Gamma
%    6     | 64 – 96           | High Gamma / noise
%    7     | 96 – 112          | (above EEG interest)
%    8     |112 – 120          | (above EEG interest)
%    9     |120 – 124          | (above EEG interest)
%   10     |124 – 128          | (above EEG interest)
%
%  For EEG we typically use subbands 1-5 directly.
%
%  Inputs:
%    Hz, MFz - prototype filter coefficients
%    params  - parameter struct
%
%  Outputs:
%    H_bank   - [NUM_BANDS x NFFT/2+1] magnitude response matrix
%    freq_axis- frequency vector (Hz)
%    band_map - struct mapping EEG bands to subbands

    fprintf('--- Building EEG FRM Filter Bank ---\n');

    Fs        = params.Fs;
    NUM_BANDS = params.NUM_BANDS;
    NFFT      = params.NFFT;

    freq_axis = (0 : NFFT/2) / NFFT * Fs;   % frequency in Hz

    %% ---- Compute DTFT of prototype filters ----
    H_H  = freqz(Hz,  1, NFFT, 'whole');
    H_MF = freqz(MFz, 1, NFFT, 'whole');

    %% ---- Complementary filters ----
    N_H   = length(Hz);
    N_MF  = length(MFz);
    w     = (0 : NFFT-1)' * 2*pi / NFFT;

    % Hc(z)  = z^{-(N_H-1)/2}  - H(z)
    % MFc(z) = z^{-(N_MF-1)/2} - MF(z)
    H_Hc  = exp(-1j * w * (N_H  - 1)/2) - H_H;
    H_MFc = exp(-1j * w * (N_MF - 1)/2) - H_MF;

    %% ---- Build interpolated filter responses ----
    % For 10-band: max interpolation = 2^4 = 16
    % H(z^M) in time domain = upsample H by M (insert M-1 zeros)
    max_exp = 4;   % 2^4 = 16

    H_Hm  = cell(max_exp + 1, 1);
    H_MFm = cell(max_exp,     1);

    H_Hm{1}  = H_H;
    H_MFm{1} = H_MF;

    for k = 2 : max_exp + 1
        M     = 2^(k-1);
        hz_up = upsample_coeff(Hz,  M);
        H_Hm{k} = freqz(hz_up, 1, NFFT, 'whole');
    end
    for k = 2 : max_exp
        M      = 2^(k-1);
        mf_up  = upsample_coeff(MFz, M);
        H_MFm{k} = freqz(mf_up, 1, NFFT, 'whole');
    end

    % Complementary masking filter interpolated versions
    H_MFcm    = cell(max_exp, 1);
    H_MFcm{1} = H_MFc;
    for k = 2 : max_exp
        M     = 2^(k-1);
        mfc_up = upsample_complement(MFz, N_MF, M, NFFT);
        H_MFcm{k} = mfc_up;
    end

    %% ---- Build branch transfer functions Pi(z) ----
    P = zeros(NFFT, NUM_BANDS);

    % Lower branches (original MF path)
    % P1 = H(z^16) MF(z^8) MF(z^4) MF(z^2) MF(z)
    P(:,1)  = H_Hm{5} .* H_MFm{4} .* H_MFm{3} .* H_MFm{2} .* H_MFm{1};
    % P2 = H(z^8)  MF(z^4) MF(z^2) MF(z)
    P(:,2)  = H_Hm{4} .* H_MFm{3} .* H_MFm{2} .* H_MFm{1};
    % P3 = H(z^4)  MF(z^2) MF(z)
    P(:,3)  = H_Hm{3} .* H_MFm{2} .* H_MFm{1};
    % P4 = H(z^2)  MF(z)
    P(:,4)  = H_Hm{2} .* H_MFm{1};
    % P5 = H(z)
    P(:,5)  = H_Hm{1};

    % Upper branches (complementary MFc path)
    % P6 = Hc(z)
    P(:,6)  = H_Hc;
    % P7 = H(z^2) MFc(z)
    P(:,7)  = H_Hm{2} .* H_MFcm{1};
    % P8 = H(z^4) MF(z^2) MFc(z)
    P(:,8)  = H_Hm{3} .* H_MFm{2} .* H_MFcm{1};
    % P9 = H(z^8) MF(z^4) MF(z^2) MFc(z)
    P(:,9)  = H_Hm{4} .* H_MFm{3} .* H_MFm{2} .* H_MFcm{1};
    % P10= H(z^16) MF(z^8) MF(z^4) MF(z^2) MFc(z)
    P(:,10) = H_Hm{5} .* H_MFm{4} .* H_MFm{3} .* H_MFm{2} .* H_MFcm{1};

    %% ---- Compute subband outputs Bi(z) ----
    half = NUM_BANDS / 2;
    B    = zeros(NFFT, NUM_BANDS);

    % Lower half: B1=P1, Bi = Pi - P_{i-1}
    B(:,1) = P(:,1);
    for i = 2 : half
        B(:,i) = P(:,i) - P(:,i-1);
    end
    % Upper half: B_last=P_last, Bi = Pi - P_{i+1}
    B(:, NUM_BANDS) = P(:, NUM_BANDS);
    for i = NUM_BANDS-1 : -1 : half+1
        B(:,i) = P(:,i) - P(:,i+1);
    end

    %% ---- Extract one-sided magnitude response ----
    H_bank = zeros(NUM_BANDS, NFFT/2 + 1);
    for i = 1 : NUM_BANDS
        H_bank(i, :) = abs(B(1:NFFT/2+1, i));
    end

    %% ---- Map EEG bands to subbands ----
    % Find which subband has maximum energy in each EEG band
    band_map = map_eeg_bands(H_bank, freq_axis, params);

    %% ---- Print subband center frequencies ----
    fprintf('  Subband center frequencies (Hz):\n');
    for i = 1 : NUM_BANDS
        [~, pk] = max(H_bank(i,:));
        fc = freq_axis(pk);
        fprintf('    Subband %2d: peak at %6.2f Hz\n', i, fc);
    end
    fprintf('\n');
end


%% ---- Helper: upsample by M (insert M-1 zeros between samples) ----
function h_up = upsample_coeff(h, M)
    N    = length(h);
    h_up = zeros(1, (N-1)*M + 1);
    h_up(1 : M : end) = h;
end

%% ---- Helper: compute DTFT of MFc(z^M) ----
function H_out = upsample_complement(mf, N_MF, M, NFFT)
    mf_up    = upsample_coeff(mf, M);
    delay_n  = round(M * (N_MF - 1) / 2);
    L        = length(mf_up);
    delay_imp = zeros(1, max(L, delay_n + 1));
    delay_imp(delay_n + 1) = 1;
    mf_padded = [mf_up, zeros(1, length(delay_imp) - L)];
    mfc       = delay_imp - mf_padded;
    H_out     = freqz(mfc, 1, NFFT, 'whole');
end

%% ---- Map EEG clinical bands to subbands ----
function band_map = map_eeg_bands(H_bank, freq_axis, params)
    bands     = params.bands;
    NUM_BANDS = params.NUM_BANDS;
    n_bands   = length(bands);
    band_map  = struct();

    fprintf('  EEG Band → Subband mapping:\n');
    for b = 1 : n_bands
        flo   = bands(b).flo;
        fhi   = bands(b).fhi;
        f_idx = (freq_axis >= flo) & (freq_axis <= fhi);

        if ~any(f_idx)
            band_map(b).subband_idx = 1;
            band_map(b).name        = bands(b).name;
            continue;
        end

        % Find subband with most energy in this band
        energy = zeros(NUM_BANDS, 1);
        for i = 1 : NUM_BANDS
            energy(i) = sum(H_bank(i, f_idx).^2);
        end
        [~, best_sb] = max(energy);

        band_map(b).name        = bands(b).name;
        band_map(b).flo         = flo;
        band_map(b).fhi         = fhi;
        band_map(b).subband_idx = best_sb;
        band_map(b).color       = bands(b).color;

        fprintf('    %-12s [%3d-%3d Hz] => Subband %2d\n', ...
                bands(b).name, flo, fhi, best_sb);
    end
    fprintf('\n');
end
