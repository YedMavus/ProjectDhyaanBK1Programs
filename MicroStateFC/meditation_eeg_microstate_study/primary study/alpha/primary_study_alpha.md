# Stage-wise Alpha-band Microstate Analysis  
## State-specific microstate signatures in meditators vs controls

This study is conducted on **alpha-band (8–12 Hz) microstate templates**.

---

## Motivation

This study investigates whether meditators and controls exhibit different EEG microstate signatures across experimental stages. The central idea is that meditation may not only influence brain activity during meditation itself, but may also induce **stage-specific and carryover effects** before and after meditation.

We focus specifically on the **alpha-band EEG signal**, which is strongly associated with:

- relaxed wakefulness,  
- internal attention,  
- meditation-related neural processes.  

We examine whether the two groups differ in:

1. **Microstate coverage** (time spent in each microstate),  
2. **Microstate duration** (temporal stability),  
3. **Functional connectivity within microstates**.  

The analysis is performed stage-wise across:

```

EO1, EC1, G1, M1, G2, EO2, EC2, M2

````

The goal is to:

- compare controls and meditators at each stage,  
- model overall **Stage**, **Group**, and **Stage × Group interaction effects**.  

---

## Dataset

The study uses the `full10` subject set:

- **5 control subjects**  
- **5 meditator subjects**  

Subject metadata is loaded from:

```matlab
selected_subjects.mat
````

The analysis is performed on **alpha-band filtered EEG data**.

---

## Main Analysis Pipeline

### 1. Subject-level Microstate Templates

* Each subject has **K = 5 alpha-band microstate maps**.
* All subject templates are concatenated and reclustered to obtain a **common alpha-band template set** shared across subjects and groups.

**Outputs:**

```
commonTemplates_alpha_study2.mat  
P01_common_templates_alpha_study2.png  
```

---

### 2. Backfitting to EEG Data

For each subject and stage:

* EEG is reshaped into:

```
channels × time
```

* Preprocessing steps:

  * removal of invalid time points
  * channel-wise mean removal
  * **alpha-band filtering (8–12 Hz)**
  * average reference subtraction

* Data is divided into **30-second windows**

* Microstate templates are backfitted using spatial correlation

* Each time point is assigned a microstate label

**Output:** label sequence per window

---

### 3. Window-level Feature Extraction

For each window and each microstate:

#### Coverage

Fraction of time spent in a microstate:

```
Coverage = samples assigned / total samples
```

#### Duration

Average continuous run length (in seconds):

```
Duration = average run length / Fs
```

#### Global WPLI

* Computed using only time points belonging to the microstate
* Phase-based connectivity in the **alpha band**
* Averaged across all channel pairs

**Outputs:**

```
window_features_alpha_study2.csv  
window_features_alpha_study2.mat  
```

---

### 4. Subject-level Aggregation

Window-level features are averaged within:

```
Subject × Group × Band × Stage × Microstate
```

This produces one value per subject, stage, and microstate.

**Outputs:**

```
subject_level_features_alpha_study2.csv  
subject_level_features_alpha_study2.mat  
```

---

## Statistical Analysis

### 1. Stage-wise Group Comparison (Permutation Test)

For each feature, stage, and microstate:

* control mean
* meditator mean
* difference (meditator − control)
* Cohen’s d
* permutation p-value

**Parameters:**

```matlab
nPerm = 5000
```

**Outputs:**

```
stagewise_group_tests_alpha_study2.csv  
stagewise_group_tests_alpha_study2.mat  
```

---

### 2. Linear Mixed-Effects Model

Model:

```
Response ~ Stage * Group + (1 | Subject)
```

Tests:

* main effect of Stage
* main effect of Group
* Stage × Group interaction

This captures **stage-dependent modulation of alpha microstates**.

**Outputs:**

```
LMM_summary_alpha_study2.csv  
LMM_ANOVA_<Feature>_MS<Microstate>_alpha_study2.csv  
LMM_<Feature>_MS<Microstate>_alpha_study2.mat  
```

---

## Dominance Analysis

Determines dominant microstates based on **alpha-band coverage**.

For each group and stage:

* most dominant microstate
* second most dominant microstate
* corresponding coverage values

**Outputs:**

```
coverage_dominance_alpha_study2.csv  
coverage_dominance_alpha_study2.mat  
```

---

## Plots Generated

### P01: Common Alpha Microstate Templates

```
P01_common_templates_alpha_study2.png
```

---

### P02–P04: Group Difference Heatmaps

Heatmaps of:

```
Meditator mean − Control mean
```

Generated for:

* Coverage
* Duration
* GlobalWPLI

```
P02_heatmap_delta_Coverage_alpha.png  
P03_heatmap_delta_Duration_alpha.png  
P04_heatmap_delta_GlobalWPLI_alpha.png  
```

---

### P05: Dominance Pattern

```
P05_dominance_pattern_alpha.png
```

---

### P06: M1 Coverage Profile

```
P06_M1_coverage_profile_alpha.png
```

---

### P07: EC Average Coverage Profile

```
P07_EC_avg_coverage_profile_alpha.png
```

---

### P08: G-stage GlobalWPLI Profile

```
P08_G_avg_GlobalWPLI_profile_alpha.png
```

---

## Interpretation Goal

This analysis aims to determine whether meditators exhibit **alpha-specific microstate signatures** compared to controls.

Key questions:

* Do meditators show different microstate dominance patterns?
* Are alpha microstates more or less temporally stable?
* Does functional connectivity differ within microstates across stages?

### Key outputs for interpretation:

1. Coverage difference heatmaps
2. GlobalWPLI difference heatmaps
3. Dominance patterns
4. M1 coverage profile
5. EC1/EC2 coverage profile
6. LMM Stage × Group interaction

Together, these results help determine whether meditation modulates
**alpha-band microstate dynamics and connectivity across the experimental trajectory**.
