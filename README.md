# Mental_Imagery_Pathway_Specific_Automaticity

Code repository for my honors thesis: *Pathway-Specific Signatures in Visual Mental Imagery: Pupil and Gaze Measures During Occlusion*.

This project tests whether internally generated visual representations during temporary occlusion differ in automaticity across two dimensions:

- a ventral / surface-feature signature indexed by pupil-based luminance responses
- a dorsal / spatiotemporal signature indexed by gaze-based trajectory tracking

In the task, participants viewed a rotating coin that briefly disappeared behind a background-matched occluder. The main comparison is between `implicit` trials (passive viewing, minimal intentional demand) and `explicit` trials (instructed imagery), with `perception` used as a reference phase.

The central result of the thesis is that gaze-based trajectory tracking was reduced during passive viewing relative to instructed imagery, whereas pupil-based luminance responses were comparatively stable across phases. This supports pathway-specific differences in the automaticity of internal visual representations.

## Repository Structure

- `DataPreprocessing.R`: preprocess raw EyeLink/DataViewer `.txt` exports into participant-level `.rds` files
- `raw_graph/`: time-course and summary plotting scripts for pupil and gaze data
- `analysis/Main_LMM/`: main ventral and dorsal summary scripts plus mixed-effects models
- `analysis/TrialLevel_Coupling/`: trial-level coupling analyses between dorsal and ventral signatures
- `analysis/ParticipantLevel_Coupling/`: participant-level coupling summaries and visualization
- `analysis/Dorsal_Analysis/`: alternate dorsal metrics and derived plots (`corridor`, `cosine`, `RMSE`)
- `analysis/Dorsal_Subgroup_analysis/`: group-specific dorsal analyses and the corridor hit-rate pipeline
- `analysis/Questionnaire/`: questionnaire-performance association analyses (`VVIQ`, `VOSI`)
- `StimulusPhaseSummary.R`: participant-by-phase HTML/CSV summaries for ventral and dorsal measures
- `VentralTrialResponsesMerged.R`: merged HTML report for trial-level pupil responses

## Expected Local Data Layout

This repository currently contains code only. Raw eye-tracking data, processed participant files, questionnaire spreadsheets, reference trajectory files, and generated outputs are not included.

Most scripts assume a local directory layout like:

```text
data/
  raw/Raw_Data_Final/
  processed/RDS_Data/
  processed/*.csv
  metadata/Participants_info.xlsx
  reference/coin_path_white.xlsx
  reference/coin_path_black.xlsx
outputs/
```

These defaults can usually be overridden with command-line arguments.

## Typical Workflow

1. Preprocess raw eye-tracking exports into participant-level `.rds` files:

```bash
Rscript DataPreprocessing.R data/raw/Raw_Data_Final data/processed/RDS_Data
```

2. Compute ventral trial-level pupil responses:

```bash
Rscript analysis/Main_LMM/VentrualAutomaticityIndex.R data/processed/RDS_Data outputs/PupilSize_Raw
```

3. Derive the main dorsal trajectory metric (corridor hit rate) and fit its main model:

```bash
Rscript analysis/Dorsal_Subgroup_analysis/GroupMixedModel_CorridorHitRate.R data/processed/RDS_Data data/metadata/Participants_info.xlsx outputs/Gaze_corrodor
```

4. Fit the main automaticity models combining ventral and dorsal signatures:

```bash
Rscript analysis/Main_LMM/AutomaticityMixedModels_Corridor.R outputs/PupilSize_Raw/ventral_trial_responses.csv outputs/Gaze_corrodor/gaze_corridor_hit_rate_trials.csv outputs/Gaze_corrodor
```

5. Optional downstream analyses include:

- cosine-based gaze alignment: `analysis/Dorsal_Analysis/GazeTarget_Cosine.R`
- trial-level coupling: `analysis/TrialLevel_Coupling/`
- participant-level coupling: `analysis/ParticipantLevel_Coupling/`
- questionnaire-performance correlations: `analysis/Questionnaire/`
- figure generation: `raw_graph/`

## R Dependencies

Core packages used across the pipeline include:

- `dplyr`
- `tidyr`
- `ggplot2`
- `readxl`
- `lme4`
- `lmerTest`
- `emmeans`
- `data.table`
- `PupilPre`
- `eyetrackingR`
- `stringr`
- `purrr`
- `rlang`
- `scales`
- `tibble`

## Notes

- In this codebase, `implicit` corresponds to passive viewing and `explicit` corresponds to instructed imagery.
- Some filenames preserve legacy spelling or naming conventions, such as `VentrualAutomaticityIndex.R` and `Gaze_corrodor`.
- The repository is intended for code sharing and reproducibility support; raw participant data may be restricted from public release.
