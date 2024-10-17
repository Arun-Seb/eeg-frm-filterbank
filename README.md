# EEG Band Extraction using FRM Non-Uniform Digital Filter Bank

A MATLAB implementation of **EEG frequency band extraction** using the **Frequency Response Masking (FRM)** technique, adapted from the hearing aid non-uniform digital filter bank design by Sebastian & James (2014).

---

## Table of Contents

1. [Background](#background)
2. [EEG Frequency Bands](#eeg-frequency-bands)
3. [Why FRM for EEG?](#why-frm-for-eeg)
4. [How the FRM Filter Bank Works](#how-the-frm-filter-bank-works)
5. [Key Adaptations from Hearing Aid to EEG](#key-adaptations-from-hearing-aid-to-eeg)
6. [Signal Flow](#signal-flow)
7. [File Structure](#file-structure)
8. [Usage](#usage)
9. [Parameters Guide](#parameters-guide)
10. [Output Files](#output-files)
11. [Clinical EEG Metrics Computed](#clinical-eeg-metrics-computed)
12. [Using Your Own EEG Data](#using-your-own-eeg-data)
13. [Limitations and Future Work](#limitations-and-future-work)
14. [References](#references)

---

## Background

### What is EEG?

**Electroencephalography (EEG)** records electrical activity of the brain through electrodes placed on the scalp. The raw EEG signal is a superposition of neural oscillations occurring at different frequencies, each associated with different brain states and cognitive processes.

### The Core Problem

The raw EEG signal looks like noise — it contains all frequency bands mixed together. To analyse brain state, diagnose conditions, or build brain-computer interfaces (BCIs), we need to **separate the EEG into its constituent frequency bands**. This is the band extraction (or band decomposition) problem.

---

## EEG Frequency Bands

| Band | Range (Hz) | Brain State / Function |
|------|-----------|----------------------|
| **Delta** (δ) | 0.5 – 4 | Deep sleep, unconsciousness, brain injury. Largest amplitude (~50-100 µV). |
| **Theta** (θ) | 4 – 8 | Drowsiness, memory encoding, REM sleep, meditation. Prominent in hippocampus. |
| **Alpha** (α) | 8 – 13 | Relaxed wakefulness, eyes closed. The dominant "resting" rhythm (~10 Hz). |
| **Beta** (β) | 13 – 30 | Active thinking, concentration, anxiety. Low amplitude, fast. |
| **Gamma** (γ) | 30 – 100 | High-level cognition, feature binding, perception. Very low amplitude (~2-8 µV). |

> **Audiogram analogy:** Just as a hearing aid audiogram shows hearing loss at different frequencies and requires a filter bank to apply selective amplification, an EEG recording shows brain activity at different frequencies and requires a filter bank to separate those activities.

---

## Why FRM for EEG?

### Traditional EEG filtering approaches

| Method | Pros | Cons |
|--------|------|------|
| Butterworth IIR | Simple, fast | Non-linear phase, instability risk |
| FIR bandpass | Linear phase, stable | High order needed for sharp cutoff |
| Wavelet decomposition | Multi-resolution | Not aligned to clinical EEG bands |
| FFT-based | Exact frequency resolution | No time-domain output |

### Why FRM is a good fit

**1. Non-uniform band spacing is natural for EEG**

EEG bands are not equally spaced. The clinically important transitions are:
- Delta/Theta boundary: 4 Hz (very narrow, requires sharp filter)
- Theta/Alpha boundary: 8 Hz (narrow)
- Beta spans 17 Hz (wider)
- Gamma can span 70+ Hz (very wide)

This mirrors the hearing aid problem where more resolution was needed at low frequencies — and is exactly the motivation for the non-uniform FRM approach.

**2. Linear phase = no phase distortion**

FIR filters have exact linear phase. For EEG signals where the timing relationships between different frequency components carry information (e.g., cross-frequency coupling, phase-amplitude coupling), this is essential. IIR filters distort phase differently at each frequency.

**3. Low hardware complexity — critical for wearable EEG**

The FRM bank needs only ~10 multipliers for the entire filter bank. This is directly applicable to:
- **Wearable EEG headsets** (limited battery)
- **Implantable BCIs** (very constrained power budget)
- **Real-time seizure detection devices**

**4. Guaranteed stability**

FIR filters are unconditionally stable. IIR filters (like Butterworth or Chebyshev) can become unstable, which is dangerous in medical devices.

---

## How the FRM Filter Bank Works

### The FRM Principle (Lim, 1986)

A sharp linear-phase FIR filter is built by combining a **prototype filter** with **masking filters**:

```
Ha(z) = H(z^M) · MF(z) + Hc(z^M) · MFc(z)
```

Where:
- `H(z)` — prototype **half-band** filter (the "band-edge shaping" filter)
- `MF(z)` — **masking filter** (shapes the passband/stopband boundary)
- `MFc(z) = z^{-(N-1)/2} - MF(z)` — **complementary masking filter**
- `H(z^M)` — **interpolated** H(z): insert M-1 zeros between each coefficient
- `M` — interpolation factor (2, 4, 8, 16, 32...)

### Half-band Filter Properties (Why H(z) is efficient)

A half-band FIR filter has two special properties that make it very efficient:

1. **Zero odd coefficients**: All coefficients at odd distance from center are exactly zero. Only about half the coefficients are non-zero.
2. **Symmetry**: The non-zero coefficients are symmetric, so each unique value is used twice.

Combined: a length-19 half-band filter needs only ~5 multiplications instead of 19.

### Building Subbands

The 10-band filter bank builds each subband by cascading interpolated filters:

```
Subband 1 (lowest freq):  H(z^16) → MF(z^8) → MF(z^4) → MF(z^2) → MF(z)
Subband 2:                H(z^8)  → MF(z^4) → MF(z^2) → MF(z)
Subband 3:                H(z^4)  → MF(z^2) → MF(z)
Subband 4:                H(z^2)  → MF(z)
Subband 5:                H(z)
──────────────── mid-band ────────────────
Subband 6:                Hc(z)
Subband 7:                H(z^2)  → MFc(z)
Subband 8:                H(z^4)  → MF(z^2) → MFc(z)
Subband 9:                H(z^8)  → MF(z^4) → MF(z^2) → MFc(z)
Subband 10 (highest):     H(z^16) → MF(z^8) → MF(z^4) → MF(z^2) → MFc(z)
```

The actual **subband outputs** are differences between adjacent branches:

```
B₁(z) = P₁(z)
Bᵢ(z) = Pᵢ(z) − Pᵢ₋₁(z),    i = 2, 3, 4, 5
B₆(z) = P₆(z) − P₇(z)
Bᵢ(z) = Pᵢ(z) − Pᵢ₊₁(z),    i = 7, 8, 9
B₁₀(z) = P₁₀(z)
```

This subtraction creates the non-overlapping subband structure.

### Why Interpolation Compresses Frequency?

When you replace `z` with `z^M` (insert M-1 zeros between coefficients), the frequency response **compresses by factor M**:

- `H(z)` has transition at normalized frequency 0.5 (Nyquist/2)
- `H(z^2)` has transitions at 0.25 and 0.75
- `H(z^4)` has transitions at 0.125, 0.375, 0.625, 0.875

This is how a single prototype filter creates multiple narrow bands at different frequencies.

---

## Key Adaptations from Hearing Aid to EEG

| Aspect | Hearing Aid | EEG |
|--------|-------------|-----|
| Sampling rate | 16,000 Hz | 256 Hz (typical) |
| Frequency range of interest | 0 – 8,000 Hz | 0.5 – 60 Hz |
| Band boundaries | Logarithmically spaced (audiogram) | Delta/Theta/Alpha/Beta/Gamma |
| Transition bandwidth | norm_tbw = 0.15 | norm_tbw = 0.04 (much tighter) |
| Filter length H(z) | ~19 taps | ~83 taps |
| Filter length MF(z) | 11 taps | 15 taps |
| Stopband attenuation | 80 dB | 80 dB (same) |
| Signal amplitude | Acoustic (dB SPL) | EEG microvolts (µV) |
| Gain adjustment | Per audiogram | Per band of interest |

### The Critical Difference: Normalized Transition Bandwidth

In hearing aids, the delta-theta boundary is a large fraction of the sampling rate:
- 4 Hz / 8000 Hz = 0.0005 → but multiple bands over entire range, so TBW ≈ 0.15

In EEG, with Fs = 256 Hz:
- Delta/Theta transition: 4 Hz / 128 Hz = 0.031 → requires norm_tbw ≈ 0.04
- This means **much longer filters** are needed for sharp EEG band separation

This is why EEG filter design is more challenging computationally — the bands are much narrower relative to the sampling rate.

---

## Signal Flow

```
EEG Signal x(t)  [256 Hz, µV]
        │
        ▼
┌───────────────────────────────────────────────────────────────┐
│            FRM Filter Bank (10 subbands)                      │
│                                                               │
│  ┌─────────────────────────────────────────────────────┐     │
│  │  Prototype Filters: H(z) [half-band], MF(z)         │     │
│  │  Complementary:     Hc(z), MFc(z)                   │     │
│  └─────────────────────────────────────────────────────┘     │
│                                                               │
│  Subband 1  →  0 – ~4 Hz    (Delta band)                     │
│  Subband 2  →  ~4 – 8 Hz    (Theta band)                     │
│  Subband 3  →  ~8 – 16 Hz   (Alpha + low Beta)               │
│  Subband 4  →  ~16 – 32 Hz  (Beta band)                      │
│  Subband 5  →  ~32 – 64 Hz  (Low Gamma)                      │
│  Subbands 6-10 → upper half (above EEG interest)             │
└───────────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────────┐
│          Clinical Bandpass Filtering (Method B)               │
│  firpm FIR bandpass applied for precise clinical boundaries:  │
│    Delta:  0.5–4 Hz  │  Theta: 4–8 Hz  │  Alpha: 8–13 Hz    │
│    Beta: 13–30 Hz    │  Gamma: 30–60 Hz                       │
└───────────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────────┐
│          Analysis & Metrics                                    │
│  • Time-domain waveforms per band                             │
│  • Power Spectral Density (Welch)                             │
│  • Spectrogram (STFT)                                         │
│  • Band power over time (sliding window)                      │
│  • Alpha/Beta ratio, Theta/Alpha ratio, Engagement index      │
│  • Individual Alpha Frequency (IAF)                           │
│  • Spectral Edge Frequency (SEF95)                            │
└───────────────────────────────────────────────────────────────┘
```

---

## File Structure

```
EEG_FRM_FilterBank/
├── main_eeg_filterbank.m           % Entry point — runs complete pipeline
├── eeg_design_prototype_filters.m  % Design H(z) and MF(z) for EEG Fs
├── eeg_build_filterbank.m          % Build 10-band FRM bank, map to EEG bands
├── eeg_plot_filterbank.m           % Plot subband responses + EEG band overlay
├── eeg_generate_synthetic.m        % Generate synthetic EEG with known components
├── eeg_extract_bands.m             % Apply filter bank to EEG signal (2 methods)
├── eeg_analyse_results.m           % Compare extracted vs truth, compute metrics
├── eeg_spectral_analysis.m         % PSD, spectrogram, band power over time
├── eeg_band_power_metrics.m        % Clinical ratios, IAF, SEF95
└── README.md
```

---

## Usage

### Run the full pipeline

```matlab
run('main_eeg_filterbank.m')
```

This runs the complete pipeline:
1. Designs prototype filters for Fs = 256 Hz
2. Builds the 10-band FRM filter bank
3. Generates 10-second synthetic EEG with known Delta/Theta/Alpha/Beta/Gamma components
4. Extracts EEG bands using the filter bank
5. Computes all metrics and saves all figures

### Run individual components

```matlab
% Just design filters and view their responses
params.Fs = 256; params.norm_tbw = 0.04; params.MF_len = 15; params.NFFT = 4096;
[Hz, MFz, info] = eeg_design_prototype_filters(params);

% Generate a synthetic EEG signal only
params.duration = 10;
[eeg, t, truth] = eeg_generate_synthetic(params);
```

---

## Parameters Guide

All parameters are set at the top of `main_eeg_filterbank.m`:

| Parameter | Default | Effect |
|-----------|---------|--------|
| `params.Fs` | `256` | EEG sampling frequency. Common: 128, 256, 512, 1024 Hz. |
| `params.duration` | `10` | Synthetic EEG duration in seconds. |
| `params.norm_tbw` | `0.04` | Normalized transition bandwidth. **Smaller = sharper filter = longer filter**. Recommended: 0.03–0.06 for EEG. |
| `params.MF_len` | `15` | Masking filter length (odd integer). Longer = slightly sharper masking. |
| `params.NFFT` | `4096` | FFT size for frequency analysis. Power of 2. |
| `params.NUM_BANDS` | `10` | Number of subbands. Fixed at 10 for this implementation. |

### Choosing `norm_tbw` for EEG

The normalized transition bandwidth determines filter sharpness:

| norm_tbw | H(z) length | Sharpness | Use case |
|----------|-------------|-----------|----------|
| 0.06 | ~55 taps | Moderate | Low-power wearable, coarse separation |
| 0.04 | ~83 taps | Good | Standard EEG analysis (recommended) |
| 0.03 | ~111 taps | Very good | Research grade, plenty of compute |
| 0.02 | ~165 taps | Excellent | High-resolution, offline analysis |

---

## Output Files

| File | Description |
|------|-------------|
| `EEG_Prototype_Filters.png` | H(z) and MF(z) magnitude responses (linear + dB) |
| `EEG_FilterBank_Response.png` | All 10 subbands overlaid with EEG band boundaries |
| `EEG_Band_Subband_Mapping.png` | Which subbands contribute to each EEG band (heat map) |
| `EEG_Synthetic_Signal.png` | Raw synthetic EEG + all 5 known components |
| `EEG_Band_Decomposition.png` | Raw EEG stacked above each extracted band |
| `EEG_Extracted_vs_Truth.png` | Extracted signal vs ground truth (overlay, first 5 s) |
| `EEG_Extraction_Metrics.png` | Correlation and SNR bar charts |
| `EEG_PSD_Analysis.png` | Welch PSD of raw EEG and each extracted band |
| `EEG_Spectrogram.png` | STFT spectrogram with band boundary lines |
| `EEG_Band_Power_Time.png` | Band power over time + stacked area chart |
| `EEG_Band_Power_Metrics.png` | Absolute power, relative power pie, clinical ratios |
| `EEG_PSD_BandShaded.png` | PSD with each band shaded, IAF and SEF95 markers |

---

## Clinical EEG Metrics Computed

### Band Power

- **Absolute power (µV²)**: Integral of PSD over each band
- **Relative power (%)**: Each band's share of total power — more robust across subjects

### Clinical Ratios

| Ratio | Formula | Clinical Meaning |
|-------|---------|-----------------|
| Alpha/Beta | α / β | > 1: relaxed; < 0.5: stressed/active |
| Theta/Alpha | θ / α | High: drowsiness, fatigue |
| Delta/Alpha | δ / α | High: deep sleep, pathological slowing |
| Engagement Index | (θ+α) / (α+β) | Low: cognitively engaged |
| Gamma Relative | γ / (α+β) | Elevated in certain cognitive tasks |

### Individual Alpha Frequency (IAF)

The frequency of peak alpha power (typically 9–11 Hz). IAF slows with:
- Fatigue and sleep deprivation
- Aging
- Some neurological conditions

### Spectral Edge Frequency (SEF95)

The frequency below which 95% of total EEG power lies. Used in:
- Anesthesia depth monitoring
- Sleep staging

---

## Using Your Own EEG Data

Replace the synthetic EEG generation step with your recorded data:

```matlab
% Option 1: Load from .mat file
load('my_eeg_data.mat');           % assumes variable 'eeg' exists
eeg_signal = eeg(1, :);            % use first channel
params.Fs  = 256;                  % set your actual sampling rate

% Option 2: Load from EDF file (requires EEGLAB or fieldtrip)
% eeg_signal = edfread('recording.edf');

% Option 3: Load from CSV
data       = csvread('eeg_data.csv');
eeg_signal = data(1, :);

% Then run extraction:
[Hz, MFz, ~]           = eeg_design_prototype_filters(params);
[H_bank, freq_axis, band_map] = eeg_build_filterbank(Hz, MFz, params);
[extracted, band_sigs] = eeg_extract_bands(eeg_signal, Hz, MFz, params, band_map);
eeg_spectral_analysis(eeg_signal, extracted, band_sigs, t, params, band_map);
eeg_band_power_metrics(band_sigs, [], params, band_map);
```

### Important preprocessing for real EEG

Before running the filter bank on real EEG data:

1. **Remove DC offset**: `eeg_signal = eeg_signal - mean(eeg_signal);`
2. **Artifact rejection**: Remove epochs with eye blinks, muscle artifacts (amplitude > 100 µV)
3. **Notch filter**: Remove power line noise at 50 Hz (EU) or 60 Hz (US)
   ```matlab
   % 50 Hz notch filter
   notch = designfilt('bandstopfir','FilterOrder',200, ...
                      'CutoffFrequency1',49,'CutoffFrequency2',51, ...
                      'SampleRate',Fs);
   eeg_signal = filtfilt(notch, eeg_signal);
   ```
4. **Re-referencing**: Common average reference or linked mastoids

---

## Limitations and Future Work

### Current Limitations

1. **Single channel**: This implementation processes one EEG channel at a time. Multi-channel processing (for source localization) would need to loop over channels.

2. **Fixed band boundaries**: The standard clinical band definitions (delta 0.5–4 Hz, etc.) are used. Real research sometimes uses subject-specific boundaries.

3. **Subband-to-EEG-band alignment**: The FRM subband boundaries don't perfectly align with clinical EEG band edges. Method B (targeted bandpass) gives more precise clinical boundaries.

4. **Real-time processing**: The current implementation uses `filtfilt` (zero-phase, non-causal). For real-time BCI applications, replace with causal `filter()` and account for group delay.

### Potential Extensions

- **Multi-channel EEG**: Apply filter bank to all channels, compute connectivity metrics (coherence, phase synchrony) between bands
- **Event-related desynchronization (ERD/ERS)**: Detect motor imagery for BCI by tracking alpha/beta power changes
- **Seizure detection**: Monitor delta power increase and alpha/beta decrease — hallmarks of seizure onset
- **Sleep staging**: Use delta/theta/alpha ratios to classify sleep stages automatically
- **Neurofeedback**: Real-time alpha/beta feedback for attention training
- **Adaptive filter bank**: Adjust band boundaries to the subject's Individual Alpha Frequency

---

## Requirements

- MATLAB R2016b or later
- Signal Processing Toolbox (`firpm`, `freqz`, `pwelch`, `spectrogram`, `filtfilt`)
- Optimization Toolbox (optional — `lsqnonneg`; falls back to `\` operator)

---

## References

1. **Lim, Y.C.** (1986). Frequency-response masking approach for the synthesis of sharp linear phase digital filters. *IEEE Transactions on Circuits and Systems*, 33(4), 357–364.

2. **Sebastian, A. & James, T.G.** (2014). A Low Complex 12-Band Non-uniform FIR Digital Filter Bank Using Frequency Response Masking Technique For Hearing Aid. *IJRITCC*, 2(9), 2786–2790.

3. **Sebastian, A. & James, T.G.** (2014). A Low Complex 10-Band Non-uniform FIR Digital Filter Bank Using Frequency Response Masking Technique For Hearing Aid. *ICCSC 2014*.

4. **Lian, Y. & Wei, Y.** (2005). A Computationally Efficient Non-Uniform FIR Digital Filter Bank for Hearing Aid. *IEEE Transactions on Circuits and Systems I*, 52, 2754–2762.

5. **Niedermeyer, E. & da Silva, F.L.** (2004). *Electroencephalography: Basic Principles, Clinical Applications, and Related Fields* (5th ed.). Lippincott Williams & Wilkins.

6. **Sanei, S. & Chambers, J.A.** (2007). *EEG Signal Processing*. Wiley-Interscience.

7. **Cohen, M.X.** (2014). *Analyzing Neural Time Series Data: Theory and Practice*. MIT Press.

---

## License

This code is provided for educational and research purposes. The FRM filter bank design follows the methodology in the referenced papers. If you use this in published research, please cite the original FRM paper (Lim 1986) and the hearing aid adaptation papers (Sebastian & James 2014).
