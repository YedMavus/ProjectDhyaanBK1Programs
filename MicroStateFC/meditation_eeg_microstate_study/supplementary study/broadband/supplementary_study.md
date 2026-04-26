#  EEG Microstate Trajectory Analysis (Broadband)

##  Motivation

This project investigates how **brain-state dynamics differ between controls and long-term meditators** using EEG microstate analysis.

Rather than focusing on *within-session meditation effects*, we analyze how **state transitions across stages** differ between groups. This is important because:

* There may be **unknown time gaps between stages**, making causal claims about meditation unreliable.
* Instead, we aim to identify **stable differences in brain dynamics** that may reflect **long-term meditation practice**.

---

##  Studies Overview

We analyze two trajectory-based studies:

### **Study 1: EC1 → M1 → EC2**

* Eyes-closed baseline → Meditation → Post-meditation rest
* Goal:
  Understand how **resting-state brain dynamics evolve across stages**

---

### **Study 2: G1 → M1 → G2**

* Gamma (task/active state) → Meditation → Post-meditation gamma
* Goal:
  Analyze how **non-resting brain dynamics transition across stages**

---

##  What the Code Does

The pipeline is identical across both studies:

### 1. **Template Construction**

* Load subject-level microstate templates
* Derive **common group-level templates (K = 5)** using clustering

---

### 2. **Backfitting**

* For each subject and stage:

  * EEG signals are segmented into windows
  * Each time point is assigned a microstate label via template matching

---

### 3. **Feature Extraction**

For each microstate (A–E), we compute:

* **Coverage** → fraction of time spent in a state
* **Duration** → average temporal stability of the state
* **Global WPLI** → functional connectivity within that state

---

### 4. **Statistical Analysis**

* Linear Mixed Models (LMM):

  * Test **Stage × Group interactions**
* Compare trajectories across:

  * Controls vs Meditators

---

### 5. **Change Score Analysis**

Compute key transitions:

* Study 1:

  * M1 − EC1
  * EC2 − M1
  * EC2 − EC1
* Study 2:

  * M1 − G1
  * G2 − M1
  * G2 − G1

---

### 6. **Visualization**

* Trajectory plots (stage-wise evolution)
* Change-score scatter plots
* Window-level distributions

---

##  Key Idea

We do **not claim**:

>  “Meditation directly causes these changes”

Instead, we analyze:

>  How **brain-state trajectories differ between groups**

---

##  Main Findings (High-Level)

* **Controls:**

  * Stronger **microstate redistribution**
  * More **reactive / rebound-like transitions**

* **Meditators:**

  * More **selective microstate modulation**
  * Clearer **connectivity regulation (WPLI suppression → recovery)**
  * Consistent **MS-C enhancement** across conditions

---

##  Files

* `study1_broadband.m` → EC1 → M1 → EC2 analysis
* `study2_broadband.m` → G1 → M1 → G2 analysis

---

##  Summary

This supplementary study provides a **trajectory-based framework** to study EEG microstates, showing that:

> Brain dynamics differ systematically between controls and meditators, reflecting differences in large-scale neural organization rather than just transient state effects.


