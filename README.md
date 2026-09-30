# CLEAN EMG Denoising Example

This MATLAB example loads EMG trials, removes artifacts, and saves the cleaned signals. It does not perform boundary-continuity evaluation or generate Excel reports.

## Requirements

- MATLAB (tested with R2024b)
- `findpeaks` and `nanmean` available in MATLAB

## Repository structure

```text
.
├── main.m
├── preprocess/
│   ├── removeNarrowPeaksByFWHM.m
│   └── cleanEarlySignal_HATA_threshold.m
└── data/
    ├── data.mat
    └── background_replacement.mat
```

## Input data

`data/data.mat` contains a variable named `data`. Each trial has an EMG signal and a stimulus amplitude:

```matlab
data.Level_01.Subject_01.Muscle_01.Nerve_01.X_0.EMG
data.Level_01.Subject_01.Muscle_01.Nerve_01.X_0.Stimulus_Amplitude
```

The example retains trial names and `Stimulus_Amplitude` values. The level, subject, muscle, and nerve field names use generic labels.

`data/background_replacement.mat` contains the background signal used for artifact replacement:

```matlab
background_data.Background_01
```

## Run

In MATLAB, open the repository directory and run:

```matlab
main
```

The pipeline first suppresses narrow, high-amplitude peaks. It then learns noise templates, detects candidate artifacts, and replaces them with segments from the background signal.

The cleaned trials are saved to `result/cleaned_data.mat`. The output variable is named `data` and retains the input hierarchy, trial names, and `Stimulus_Amplitude` values. Each trial’s `EMG` field contains its cleaned signal.
