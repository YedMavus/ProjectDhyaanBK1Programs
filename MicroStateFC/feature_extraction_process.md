🧠 Microstate + Functional Connectivity Feature Extraction Pipeline
📌 Overview

This pipeline performs:

Full-signal spectral analysis (0–100 Hz) using multitaper FFT
Microstate backfitting on band-limited signals (alpha, gamma)
Microstate temporal feature extraction
Hilbert-based functional connectivity (PLI, WPLI)
Organized saving of all outputs for downstream analysis
📥 Input Requirements
1. Raw EEG / LFP Data

Expected directory structure:

rootPath/
  subject/
    EEG/
      date/
        protocol (M1/G2)/
          segmentedData/
            LFP/
              lfpInfo.mat
              elec1.mat
              elec2.mat
              ...
              elec70.mat
File Details
elecX.mat
Contains:
analogData → [numTrials × numSamples]
Example: 360 × 2500
lfpInfo.mat
Contains:
electrode indices
optional channel labels
optional sampling rate (Fs)
2. Microstate Files

Required per stage and band:

013AR_M1_Microstates_alpha.mat
013AR_M1_Microstates_gamma.mat
013AR_G2_Microstates_alpha.mat
013AR_G2_Microstates_gamma.mat
Contents:
microstateMaps → [M × C]
Example: 5 × 70
⚙️ Pipeline Steps
🧠 A) Full-Signal Spectral Features
Method
FieldTrip: ft_freqanalysis
Configuration:
method = 'mtmfft'
output = 'fourier'
taper = 'dpss'
tapsmofrq = 4
foilim = [0 100]
pad = nextpow2 → nfft = 4096
Outputs
Power Spectrum
powspctrm_full → [C × F]
Connectivity (per frequency)
COH_full  → [C × C × F]
PPC_full  → [C × C × F]
PLI_full  → [C × C × F]
WPLI_full → [C × C × F]
Node Strength
nodeStrength_* → [C × F]
🎯 B) Band-Specific Processing (Alpha / Gamma)
Step 1: Bandpass Filtering
FFT-based filtering (same as microstate extraction)
Applied on continuous signal
Step 2: Microstate Backfitting

For each time point:

Normalize EEG topography
Normalize microstate templates
Compute correlation
Assign label with highest similarity
Step 3: Microstate Smoothing
Problem:
Backfitting causes rapid switching
Solution:
Minimum segment duration = 30 ms
Algorithm:
Run-length encoding
Replace short segments with neighbors:
Prefer longer neighbor
Iterate until stable (max 20 iterations)
Step 4: Microstate Features

For each microstate m:

Coverage
Fraction of time spent in state
Duration
Average segment length (seconds)
Occurrence
Number of segments per second
Dwell Time Std
Standard deviation of segment lengths
Transition Matrix
TransitionCount → [M × M]
TransitionProb  → [M × M]
Entropy
Label entropy
Transition entropy
GEV (Global Explained Variance)
GEV_m → contribution of each microstate
TotalGEV → sum over all states
🔁 C) Hilbert-Based Functional Connectivity
Step 1: Analytic Signal

Using Hilbert transform:

U = x + j * H(x)
Step 2: Downsampling
Default: every 10 samples (reduces computation)
Step 3: Connectivity Computation
PLI
PLI = |mean(sign(imag(Sxy)))|
WPLI
WPLI = |sum(imag(Sxy))| / sum(|imag(Sxy)|)
Outputs
Full Signal
PLI_full  → [C × C]
WPLI_full → [C × C]
Microstate-Specific
PLI_ms  → [C × C × M]
WPLI_ms → [C × C × M]
Node Strength
Sum of connections per node
💾 Output Structure
extracted_features/
  subject/
    stage/
      full_signal_mtmfft_spectral_fc.mat

      alpha/
        microstate_backfit_and_features.mat
        fc_hilbert_bandlimited.mat

      gamma/
        microstate_backfit_and_features.mat
        fc_hilbert_bandlimited.mat
📦 Saved Variables
1. Full Signal File
spectral_full:
  freqVec                → [1 × F]
  powspctrm_full         → [C × F]
  COH_full               → [C × C × F]
  PPC_full               → [C × C × F]
  PLI_full               → [C × C × F]
  WPLI_full              → [C × C × F]
  nodeStrength_*         → [C × F]
2. Microstate File
microstateMaps           → [M × C]
z_trial_backfit_smooth   → [T × S]
z_cont_backfit_smooth    → [1 × (T×S)]

ms_features:
  Duration_m
  Coverage_m
  Occurrence_m
  GEV_m
  TransitionProb
  Entropy metrics
3. Hilbert FC File
hilbert_fc:
  PLI_full   → [C × C]
  WPLI_full  → [C × C]
  PLI_ms     → [C × C × M]
  WPLI_ms    → [C × C × M]
  nodeStrength_*
🧩 Design Choices
Multitaper FFT (mtmfft)
Reduces spectral noise
More stable than plain FFT
Hilbert for Microstates
Time-domain aligned with microstate segmentation
Avoids heavy memory cost of TF analysis
Smoothing
Enforces physiological plausibility (~30 ms minimum state duration)
⚠️ Limitations
No time-frequency FC (memory constraints)
FFT-based band filtering (not FIR/IIR)
Hilbert assumes narrowband signals
✅ Summary

This pipeline produces:

Full-spectrum connectivity (0–100 Hz)
Microstate segmentation (alpha & gamma)
Temporal microstate dynamics
State-specific connectivity (PLI/WPLI)
Structured outputs for further analysis
