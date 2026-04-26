# Stage-wise Gamma-band Microstate Analysis  
## State-specific microstate signatures in meditators vs controls

This study is done on **gamma-band (typically ~30–45 Hz) microstate templates**.

---

## Motivation

This study investigates whether meditators and controls show different EEG microstate signatures across experimental stages, focusing specifically on the **gamma band**, which is often associated with:

- high-level cognitive processing,
- attention and awareness,
- meditation-related synchrony and integration.

We aim to determine whether meditation modulates:

1. microstate occupancy (coverage),
2. temporal stability (duration),
3. functional connectivity within microstates,

in a **stage-dependent manner**.

The analysis is performed across:

`EO1, EC1, G1, M1, G2, EO2, EC2, M2`

with comparisons between controls and meditators at each stage, along with modeling of overall Stage, Group, and Stage × Group effects.

---

## Dataset

- 5 control subjects  
- 5 meditator subjects  

Subject selection is loaded from:

```matlab
selected_subjects.mat
````

The analysis is performed on **gamma-band filtered EEG data**.

---

## Main analysis pipeline

### 1. Load subject-level microstate templates

Each subject has `K = 5` **gamma-band microstate maps**.

All subject-level templates are concatenated and clustered to derive a **common gamma-band template set**.

Output:

```text
commonTemplates_gamma_study2.mat
P01_common_templates_gamma_study2.png
```

---

### 2. Backfitting

For each subject and stage:

* EEG data is reshaped to `channels × time`
* Preprocessing includes:

  * removal of invalid samples,
  * channel mean removal,
  * **gamma-band filtering (~30–45 Hz)**,
  * average referencing.

Common templates are backfitted in 30-second windows using spatial correlation.

Each time point is assigned a microstate label → produces a label sequence per window.

---

### 3. Window-level feature extraction

For each window and microstate:

#### Coverage

Fraction of time spent in a microstate.

#### Duration

Average duration of continuous microstate segments (in seconds).

#### Global WPLI

* computed only on time points belonging to a microstate,
* phase-based connectivity in the **gamma band**,
* averaged across all channel pairs.

Output:

```text
window_features_gamma_study2.csv
window_features_gamma_study2.mat
```

---

### 4. Subject-level aggregation

Window-level features are averaged within:

```text
Subject × Group × Band × Stage × Microstate
```

Output:

```text
subject_level_features_gamma_study2.csv
subject_level_features_gamma_study2.mat
```

---

## Statistical analysis

### 1. Stage-wise permutation tests

For each feature, stage, and microstate:

* control mean vs meditator mean,
* difference (meditator − control),
* Cohen’s d,
* permutation p-value (`nPerm = 5000`).

Output:

```text
stagewise_group_tests_gamma_study2.csv
```

---

### 2. Linear mixed-effects model

For each feature and microstate:

```text
Response ~ Stage * Group + (1 | Subject)
```

Tests:

* Stage effect
* Group effect
* Stage × Group interaction

Output:

```text
LMM_summary_gamma_study2.csv
```

---

## Dominance analysis

Determines dominant microstates based on **gamma-band coverage** per group and stage.

Output:

```text
coverage_dominance_gamma_study2.csv
```

---

## Plots generated

### P01: Common gamma microstate templates

```text
P01_common_templates_gamma_study2.png
```

---

### P02–P04: Group-difference heatmaps

For:

* Coverage
* Duration
* GlobalWPLI

```text
P02_heatmap_delta_Coverage_gamma.png
P03_heatmap_delta_Duration_gamma.png
P04_heatmap_delta_GlobalWPLI_gamma.png
```

---

### P05: Dominance pattern

```text
P05_dominance_pattern_gamma.png
```

---

### P06–P08: Key stage profiles

* M1 coverage profile
* EC average coverage
* G-stage GlobalWPLI

---

## Interpretation goal

This analysis aims to determine whether meditators exhibit **gamma-specific microstate dynamics**, particularly:

* increased/decreased occupancy of specific microstates,
* altered temporal stability,
* enhanced or reduced functional connectivity within microstates.

Given the link between gamma activity and higher-order processing, the focus is on whether meditation induces **stage-dependent changes in network synchronization and microstate organization**.

Key outputs:

1. Coverage heatmaps
2. GlobalWPLI heatmaps
3. Dominance patterns
4. M1 and EC profiles
5. LMM interaction effects

These collectively reveal whether meditation modulates **gamma-band microstate dynamics and connectivity across stages**.
