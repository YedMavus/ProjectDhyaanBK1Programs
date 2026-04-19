# Microstate and Functional Connectivity Feature Extraction Pipeline

## Overview

This pipeline performs the following steps:

1. Computes full-signal spectral features from raw EEG using multitaper FFT.
2. Recreates alpha- and gamma-band signals using the same FFT-based filtering used during microstate extraction.
3. Backfits saved microstate templates to every sample of the band-limited signal.
4. Smooths rapid microstate switching using a minimum segment duration rule.
5. Computes temporal microstate features such as duration, coverage, occurrence, GEV, transition probabilities, and entropy.
6. Computes Hilbert-based PLI and WPLI for:
   - the full band-limited signal
   - each microstate separately
7. Saves all outputs in a structured folder layout for later analysis and visualization.

This version of the pipeline does **not** compute time-frequency FC and does **not** compute microstate-specific coherence or PPC.

---

## Input Data

The code expects the following directory structure for each subject and stage:

```text
<rootPath>/<subjectName>/EEG/<expDate>/<protocolName>/segmentedData/LFP/
```

For example:

```text
/Users/anatapmitra/Desktop/NSP project/data/Segmented Data/013AR/EEG/280122/M1/segmentedData/LFP/
```

Inside that folder, the following files are required.

### 1. `lfpInfo.mat`

This file is used to obtain:

- the electrode indices
- optional channel labels
- optional sampling frequency

The code checks for these possible field names:

- `analogChannelsStored`
- `Fs`, `fs`, or `samplingRate`
- `analogChannelLabels` or `channelLabels`

If no channel labels are present, labels are assigned as:

```text
elec1, elec2, elec3, ...
```

If sampling frequency is not found, the code assumes:

```text
Fs = 1000 Hz
```

### 2. `elec*.mat`

For each electrode listed in `lfpInfo.analogChannelsStored`, one file must exist:

```text
elec1.mat
elec2.mat
...
elec70.mat
```

Each file must contain:

- `analogData` of shape:

```text
[numTrials x numSamples]
```

Example:

```text
360 x 2500
```

After loading all channels, the code builds:

```text
eegDataAll: [C x T x S]
```

where:

- `C` = number of channels
- `T` = number of trials
- `S` = number of samples per trial

### 3. Saved microstate files

For each stage and band, the code expects a microstate file inside the same `LFP` folder:

```text
<subject>_<protocol>_Microstates_alpha.mat
<subject>_<protocol>_Microstates_gamma.mat
```

Examples:

```text
013AR_M1_Microstates_alpha.mat
013AR_M1_Microstates_gamma.mat
013AR_G2_Microstates_alpha.mat
013AR_G2_Microstates_gamma.mat
```

Each microstate file must contain at least:

- `microstateMaps` of shape:

```text
[M x C]
```

Example:

```text
5 x 70
```

---

## Processing Steps

## A. Load raw stage data

For each protocol such as `M1` or `G2`, the code:

1. loads `lfpInfo.mat`
2. loads all `elec*.mat` files
3. stacks them into:

```text
eegDataAll: [C x T x S]
```

It also reshapes the same data into a continuous signal:

```text
eegContinuousRaw: [C x (T*S)]
```

This continuous version is later used for FFT-based band filtering and backfitting.

---

## B. Full-signal spectral features using multitaper FFT

Before band-specific processing, the code computes full-signal spectral features from the raw EEG.

### Method

FieldTrip function used:

```matlab
ft_freqanalysis
```

Configuration:

```matlab
cfg.method      = 'mtmfft';
cfg.output      = 'fourier';
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 4;
cfg.foilim      = [0 100];
cfg.keeptrials  = 'yes';
cfg.keeptapers  = 'yes';
cfg.pad         = 'nextpow2';
```

### What this produces

The code obtains Fourier coefficients trial-by-trial and taper-by-taper, then manually computes:

- PSD
- Coherence
- PPC
- PLI
- WPLI

### Saved spectral outputs

These are stored in `spectral_full`.

#### Frequency vector

```text
freqVec: [1 x F]
```

Because `pad = nextpow2`, the FFT length becomes 4096 for 2500-sample trials at 1000 Hz, which leads to approximately 0.244 Hz frequency spacing. Therefore, in the 0-100 Hz range, `freqVec` typically has 411 bins.

#### PSD

```text
powspctrm_full: [C x F]
```

This is computed as:

```text
mean(abs(Fourier).^2, across observations)
```

#### Connectivity matrices

```text
COH_full:  [C x C x F]
PPC_full:  [C x C x F]
PLI_full:  [C x C x F]
WPLI_full: [C x C x F]
```

#### Node strength

For each metric, node strength is computed by summing all off-diagonal connections for a channel:

```text
nodeStrength_COH_full:  [C x F]
nodeStrength_PPC_full:  [C x F]
nodeStrength_PLI_full:  [C x F]
nodeStrength_WPLI_full: [C x F]
```

### Mathematical definitions used

For a pair of channels `i,j` at frequency `f`:

- Cross-spectrum:

  ```text
  Sxy = Xi * conj(Xj)
  ```

- Coherence:

  ```text
  COH = |mean(Sxy)| / sqrt(mean(Sxx) * mean(Syy))
  ```

- PLI:

  ```text
  PLI = |mean(sign(imag(Sxy)))|
  ```

- WPLI:

  ```text
  WPLI = |sum(imag(Sxy))| / sum(|imag(Sxy)|)
  ```

- PPC:

  ```text
  PPC = (|sum(u)|^2 - N) / (N*(N-1))
  ```

  where `u = Sxy / |Sxy|` and `N` is the number of valid observations.

---

## C. Band-limited signal reconstruction

For each band (`alpha`, `gamma`), the code reconstructs a band-limited continuous signal using the same FFT masking procedure as the original microstate extraction script.

### Bands

- Alpha: 8-13 Hz
- Gamma: 30-80 Hz

### Method

For each channel:

1. take FFT of the continuous signal
2. keep only the selected frequency band and its mirrored negative-frequency portion
3. inverse FFT back to time domain

This gives:

```text
eegContinuousBand: [C x (T*S)]
```

and trialwise:

```text
eegBandTrials: [C x T x S]
```

---

## D. Microstate backfitting

The saved microstate templates are backfitted sample-by-sample to the band-limited signal.

### Inputs

- `eegContinuousBand: [C x N]`
- `microstateMaps: [M x C]`

### Algorithm

1. Normalize each template to unit norm.
2. For each time sample:
   - take the scalp topography vector `x`
   - normalize `x`
   - compute similarity with each template using polarity-invariant absolute dot product
   - assign the label of the best matching template

Formally:

```text
z(n) = argmax_m |A(:,m)' * x(:,n)|
```

### Output

Raw backfitted labels:

```text
z_cont_raw:  [1 x (T*S)]
z_trial_raw: [T x S]
```

---

## E. Microstate smoothing

Backfitting often produces unrealistically fast switching. The code smooths labels using a minimum segment duration rule.

### Threshold

```text
cfgMain.minSegmentMs = 30
```

At `Fs = 1000 Hz`, this becomes:

```text
minSegSamples = 30
```

### Smoothing algorithm

The function `smooth_microstate_labels` works separately for each trial:

1. Run-length encode the microstate label sequence.
2. Find all segments shorter than `minSegSamples`.
3. Replace each short segment by a neighboring state:
   - if only the right neighbor exists, assign the right label
   - if only the left neighbor exists, assign the left label
   - if both neighbors exist, compare the lengths of the neighboring segments
   - merge into the longer neighboring segment
   - if lengths are equal, merge into the left neighbor
4. Repeat until no changes are made or until 20 iterations are reached.

### Output

Smoothed labels:

```text
z_trial_smooth: [T x S]
z_cont_smooth:  [1 x (T*S)]
```

---

## F. Microstate temporal features

After smoothing, the code computes temporal microstate statistics.

These are stored in `ms_features`.

### Coverage

Fraction of total time spent in each microstate:

```text
Coverage_m: [M x 1]
```

### Duration

Mean segment length in seconds:

```text
Duration_m: [M x 1]
```

### Occurrence

Number of segments per second:

```text
Occurrence_m: [M x 1]
```

### Dwell time variability

Standard deviation of segment lengths:

```text
DwellStd_m: [M x 1]
```

### Number of segments

```text
NumSegments_m: [M x 1]
```

### Transition counts and probabilities

```text
TransitionCount: [M x M]
TransitionProb:  [M x M]
```

### Label entropy

Computed from the state occupancy probabilities.

### Transition entropy

Computed from transition probabilities, weighted by row mass.

### GEV

Global explained variance is computed from the GFP-weighted squared correlation between each sample topography and its assigned template.

Outputs:

```text
GEV_m:    [M x 1]
TotalGEV: scalar
```

---

## G. Hilbert-based PLI and WPLI

For each band-limited signal, the code computes connectivity using the analytic signal.

### Step 1: Analytic signal

For each trial:

```text
U = hilbert(X')'
```

This gives the complex analytic signal for all channels.

### Step 2: Downsampling

To reduce computational cost, the code uses:

```text
cfgMain.hilbertDownsampleFactor = 10
```

So only every 10th sample is used.

### Step 3: Full-signal PLI/WPLI

For all retained samples across all trials:

- Cross-product:

  ```text
  crossObs = Ui * conj(Uj)
  ```

- PLI:

  ```text
  PLI = |mean(sign(imag(crossObs)))|
  ```

- WPLI:

  ```text
  WPLI = |sum(imag(crossObs))| / sum(|imag(crossObs)|)
  ```

Outputs:

```text
PLI_full:  [C x C]
WPLI_full: [C x C]
```

### Step 4: Microstate-specific PLI/WPLI

Using the smoothed microstate labels, the code groups downsampled observations by microstate and recomputes PLI/WPLI separately for each state.

Outputs:

```text
PLI_ms:  [C x C x M]
WPLI_ms: [C x C x M]
```

### Step 5: Node strength

For each connectivity matrix, node strength is computed by summing off-diagonal weights.

Outputs:

```text
nodeStrength_PLI_full:  [C x 1]
nodeStrength_WPLI_full: [C x 1]

nodeStrength_PLI_ms:  [C x M]
nodeStrength_WPLI_ms: [C x M]
```

---

## Saved Output Structure

Outputs are saved under:

```text
<outputBase>/<subjectName>/<protocolName>/
```

Example:

```text
/Users/anatapmitra/Desktop/NSP project/extracted_features/013AR/M1/
```

### Stage-level file

Saved once per stage:

```text
full_signal_mtmfft_spectral_fc.mat
```

Contains:

- `metadata_stage`
- `spectral_full`

### Band-level folders

For each band:

```text
alpha/
gamma/
```

Each band folder contains:

#### 1. `microstate_backfit_and_features.mat`

Contains:

- `metadata_band`
- `microstateMaps`
- `ms_features`
- `z_trial_backfit_smooth`
- `z_cont_backfit_smooth`
- optionally `z_trial_backfit_raw`

#### 2. `fc_hilbert_bandlimited.mat`

Contains:

- `metadata_band`
- `microstateMaps`
- `hilbert_fc`

---

## Metadata

Each saved file includes metadata describing:

- subject name
- experiment date
- protocol name
- band name and band range where applicable
- number of electrodes
- number of trials
- number of samples per trial
- sampling frequency
- electrode indices
- channel labels
- channel pair definitions

This ensures that saved matrices can be correctly interpreted later.

---

## Design Choices

### Why `mtmfft` for full-signal spectral features

A plain FFT PSD estimate is noisy. The code uses multitaper FFT with DPSS tapers to obtain more stable spectral estimates for PSD, coherence, PPC, PLI, and WPLI.

### Why Hilbert for microstate-specific FC

Microstates are defined in time. A Hilbert-based approach is well matched to time-domain state masking and avoids the memory cost of time-frequency decomposition.

### Why smoothing is necessary

Samplewise backfitting often creates extremely short, unrealistic microstate bursts. The minimum-segment smoothing step makes the label sequence more physiologically plausible.

---

## Limitations

This reduced version of the pipeline does **not** include:

- time-frequency FC
- microstate-specific coherence
- microstate-specific PPC

Band filtering is performed via FFT masking rather than FIR/IIR filtering.

Hilbert-based connectivity assumes narrowband signals.

---

## Summary

This pipeline produces:

- full-signal PSD from 0-100 Hz
- full-signal coherence, PPC, PLI, and WPLI from 0-100 Hz
- alpha- and gamma-specific microstate backfitting
- smoothed microstate label sequences
- temporal microstate features
- Hilbert-based full-signal and microstate-specific PLI/WPLI
- node strength summaries
- structured `.mat` files for later visualization and analysis
