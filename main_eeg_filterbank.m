%% =========================================================
%  EEG Band Extraction using FRM Non-Uniform Filter Bank
%  =========================================================
%  Extracts EEG frequency bands (Delta, Theta, Alpha, Beta,
%  Gamma) using the Frequency Response Masking (FRM) technique
%  adapted from the hearing aid filter bank design.
%
%  Based on:
%    Lim, Y.C. (1986). Frequency-response masking approach for
%    the synthesis of sharp linear phase digital filters.
%    IEEE Trans. Circuits Systems, 33(4), 357-364.
%
%    Sebastian & James (2014). A Low Complex 12-Band Non-uniform
%    FIR Digital Filter Bank Using FRM Technique For Hearing Aid.
%    IJRITCC, 2(9), 2786-2790.
%
%  Key adaptations for EEG:
%    - Sampling rate: 256 Hz (standard EEG, vs 16 kHz in hearing aid)
%    - Band boundaries aligned to clinical EEG bands
%    - Non-uniform spacing: fine resolution at low freqs (delta/theta)
%    - Synthetic EEG signal generation for demonstration
%
%  Usage:
%    Simply run this script. All parameters are set below.
%    Figures and results are saved automatically.
%
%  Requirements:
%    - MATLAB R2016b+
%    - Signal Processing Toolbox
%    - Optimization Toolbox (optional, for lsqnonneg)
% =========================================================

clear; clc; close all;

fprintf('=====================================================\n');
fprintf('  EEG Band Extraction — FRM Non-Uniform Filter Bank\n');
fprintf('=====================================================\n\n');

%% ============================================================
%  SECTION 1: PARAMETERS
%  All key parameters in one place for easy modification
%  ============================================================

params.Fs           = 256;      % EEG sampling frequency (Hz)
                                % Common values: 128, 256, 512, 1024 Hz
params.duration     = 10;       % Synthetic EEG duration (seconds)
params.norm_tbw     = 0.04;     % Normalized transition bandwidth
                                % Smaller = sharper filters, longer length
                                % Range tested: 0.02 to 0.08
params.MF_len       = 15;       % Masking filter length (odd integer)
                                % Increase for sharper masking
params.NFFT         = 4096;     % FFT size for frequency analysis
params.NUM_BANDS    = 10;       % Filter bank size (10 recommended for EEG)

% EEG Band definitions (Hz) — standard clinical boundaries
params.bands = struct( ...
    'name',  {'Delta','Theta','Alpha','Beta_Low','Beta_High','Gamma_Low','Gamma_Mid'}, ...
    'flo',   {0.5,    4,      8,      13,        20,         30,         45}, ...
    'fhi',   {4,      8,      13,     20,         30,         45,         60}, ...
    'color', {[0.2,0.4,0.9],[0.1,0.7,0.3],[0.9,0.7,0.1],[0.9,0.3,0.2],[0.7,0.1,0.7],[0.1,0.8,0.8],[0.5,0.5,0.5]} ...
);

fprintf('Parameters:\n');
fprintf('  Sampling frequency : %d Hz\n', params.Fs);
fprintf('  Signal duration    : %d seconds\n', params.duration);
fprintf('  Norm. TBW          : %.3f\n', params.norm_tbw);
fprintf('  Filter bank bands  : %d\n\n', params.NUM_BANDS);

%% ============================================================
%  SECTION 2: DESIGN PROTOTYPE FILTERS
%  ============================================================
[Hz, MFz, filter_info] = eeg_design_prototype_filters(params);

%% ============================================================
%  SECTION 3: BUILD EEG FILTER BANK
%  ============================================================
[H_bank, freq_axis, band_map] = eeg_build_filterbank(Hz, MFz, params);

%% ============================================================
%  SECTION 4: PLOT FILTER BANK RESPONSE
%  ============================================================
eeg_plot_filterbank(H_bank, freq_axis, params, band_map);

%% ============================================================
%  SECTION 5: GENERATE SYNTHETIC EEG SIGNAL
%  ============================================================
[eeg_signal, t, true_components] = eeg_generate_synthetic(params);

%% ============================================================
%  SECTION 6: EXTRACT EEG BANDS
%  ============================================================
[extracted, band_signals] = eeg_extract_bands(eeg_signal, Hz, MFz, params, band_map);

%% ============================================================
%  SECTION 7: ANALYSE AND PLOT RESULTS
%  ============================================================
eeg_analyse_results(eeg_signal, t, extracted, band_signals, true_components, params, band_map);

%% ============================================================
%  SECTION 8: POWER SPECTRAL ANALYSIS
%  ============================================================
eeg_spectral_analysis(eeg_signal, extracted, band_signals, t, params, band_map);

%% ============================================================
%  SECTION 9: BAND POWER METRICS
%  ============================================================
eeg_band_power_metrics(band_signals, true_components, params, band_map);

fprintf('\n=====================================================\n');
fprintf('  All done! Check generated figures and output files.\n');
fprintf('=====================================================\n');
