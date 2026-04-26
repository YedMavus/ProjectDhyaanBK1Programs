# Stage-wise Broadband Microstate Analysis  
## State-specific microstate signatures in meditators vs controls
This study is done on **broadband microstate templates.**

## Motivation

This study investigates whether meditators and controls show different EEG microstate signatures across different experimental stages. The main idea is that meditation may not only affect EEG activity during meditation itself, but may also produce stage-specific or carryover differences before and after meditation-related states.

Here, we focus on the broadband EEG signal and ask whether the two groups differ in:

1. how much time they spend in each microstate,
2. how long each microstate remains stable,
3. how strongly functional connectivity is expressed within each microstate.

The analysis is performed stage-wise across:

`EO1, EC1, G1, M1, G2, EO2, EC2, M2`

where the goal is to compare controls and meditators at each stage separately, and also model the overall effect of stage, group, and their interaction.

---

## Dataset

This script uses the `full10` subject set:

- 5 control subjects
- 5 meditator subjects

The selected subject list and group labels are loaded from:

```matlab
selected_subjects.mat
````

The analysis is run for the broadband case.

---

## Main analysis pipeline

### 1. Load subject-level microstate templates

For each subject, the script loads the saved broadband subject-level microstate templates.

Each subject has `K = 5` microstate maps.

The templates from all subjects are concatenated and clustered again to derive one common set of broadband templates shared across all subjects and groups.

Output:

```text
commonTemplates_broadband_study2.mat
P01_common_templates_broadband_study2.png
```

---

### 2. Backfit common templates to EEG data

For every subject and every stage, the EEG data is loaded from the segmented LFP files.

The data is reshaped into a continuous 2D matrix:

```text
channels × time
```

The broadband preprocessing includes:

* removing invalid time points,
* channel-wise mean removal,
* 0.5–45 Hz bandpass filtering,
* average reference subtraction.

Then, for every 30-second window, the common microstate templates are backfitted to the EEG signal.

At each time point, the EEG topography is assigned to the microstate template with the highest absolute spatial correlation.

This gives one microstate label sequence per window.

---

### 3. Extract window-level features

For each 30-second window and each microstate, the script computes:

### Coverage

Fraction of time spent in a microstate.

```text
Coverage = samples assigned to microstate / total samples
```

### Duration

Average continuous run length of a microstate, converted to seconds.

```text
Duration = average run length / Fs
```

### Global WPLI

For each microstate, the script takes only the time points assigned to that microstate and computes phase-based connectivity using WPLI.

The WPLI matrix is computed across channels, and then averaged over all unique channel pairs to obtain one global WPLI value per microstate.

Output:

```text
window_features_broadband_study2.csv
window_features_broadband_study2.mat
```

---

### 4. Aggregate to subject-level features

The window-level values are averaged within each:

```text
Subject × Group × Band × Stage × Microstate
```

This produces one feature value per subject, stage, and microstate.

Output:

```text
subject_level_features_broadband_study2.csv
subject_level_features_broadband_study2.mat
```

---

## Statistical analysis

### 1. Stage-wise group permutation tests

For each feature, stage, and microstate, the script compares controls and meditators.

The tested features are:

```text
Coverage
Duration
GlobalWPLI
```

For each comparison, it computes:

* control mean,
* meditator mean,
* meditator minus control difference,
* Cohen’s d,
* permutation p-value.

The permutation test uses:

```matlab
nPerm = 5000
```

Output:

```text
stagewise_group_tests_broadband_study2.csv
stagewise_group_tests_broadband_study2.mat
```

---

### 2. Linear mixed-effects model

For each feature and each microstate, the script fits:

```text
Response ~ Stage * Group + (1 | Subject)
```

This tests:

* main effect of Stage,
* main effect of Group,
* Stage × Group interaction.

This is useful for checking whether group differences depend on experimental stage.

Output:

```text
LMM_summary_broadband_study2.csv
LMM_ANOVA_<Feature>_MS<Microstate>_broadband_study2.csv
LMM_<Feature>_MS<Microstate>_broadband_study2.mat
```

---

## Dominance analysis

The script also computes the dominant microstate for each group and stage based on coverage.

For every group-stage pair, it identifies:

* the most dominant microstate,
* the second most dominant microstate,
* their coverage values.

Output:

```text
coverage_dominance_broadband_study2.csv
coverage_dominance_broadband_study2.mat
```

---

## Plots generated

The script generates a minimal figure set for reporting.

### P01: Common broadband microstate templates

Shows the final common K=5 broadband templates.

```text
P01_common_templates_broadband_study2.png
```

### P02–P04: Group-difference heatmaps

Heatmaps showing:

```text
Meditator mean − Control mean
```

for each stage and microstate.

Generated separately for:

* Coverage,
* Duration,
* GlobalWPLI.

```text
P02_heatmap_delta_Coverage_broadband.png
P03_heatmap_delta_Duration_broadband.png
P04_heatmap_delta_GlobalWPLI_broadband.png
```

### P05: Dominant microstate pattern

Shows which microstate is dominant in each group and stage.

```text
P05_dominance_pattern_broadband.png
```

### P06: M1 coverage profile

Compares control and meditator coverage profiles during the M1 stage.

```text
P06_M1_coverage_profile_broadband.png
```

### P07: EC average coverage profile

Compares average coverage across EC1 and EC2.

```text
P07_EC_avg_coverage_profile_broadband.png
```

### P08: Gamma-stage GlobalWPLI profile

Compares average GlobalWPLI across G1 and G2.

```text
P08_G_avg_GlobalWPLI_profile_broadband.png
```

Note: although this plot is named “gamma average” in the script, the current analysis is broadband. Here, G1 and G2 refer to experimental stages, not gamma-band filtering.

---

## Interpretation goal

This analysis is designed to identify whether meditators show stage-specific microstate signatures compared with controls.

The most important outputs for interpretation are:

1. Coverage group-difference heatmap,
2. GlobalWPLI group-difference heatmap,
3. dominance pattern plot,
4. M1 coverage profile,
5. EC1/EC2 average coverage profile,
6. LMM Stage × Group interaction table.

Together, these results help determine whether meditators and controls differ in the temporal dominance and connectivity structure of broadband EEG microstates across the experimental trajectory.
