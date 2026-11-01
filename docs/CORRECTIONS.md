# Corrections to the original coursework

Two errors were found by rechecking figures and file metadata against the raw model output. Both are fixed in this repository. The originally submitted coursework was not changed.

## 1. Duplicate hours in the MYNN2.5 run (Project 01)

**Symptom:** in the original time series, the MYNN2.5 line ran straight backwards across the plot from Jan 19 17:00 to Jan 19 01:00.

**Cause:** the run had two output files, and the data loader joined them without checking for overlap.

| File | Records | Valid times |
|------|---------|-------------|
| `wrfout_d02_2024-01-16_00_00_00` | 90 | Jan 16 00:00 to Jan 19 17:00 |
| `wrfout_d02_2024-01-19_01_00_00` | 17 | Jan 19 01:00 to Jan 19 17:00 (repeats) |

Those 17 hours were counted twice in the MYNN2.5 statistics. A spot check of U at one grid point matched in all 17 hours, so the repeats look like restart copies rather than different data.

**Fix:** `add_time_coord()` in [analysis.py](../projects/01_dunkelflaute_pbl_sensitivity/src/analysis.py) now drops repeated valid times and raises an error if times are out of order.

**Effect:**

| MYNN2.5 vs ERA5 at 100 m | Original | Corrected |
|--------------------------|----------|-----------|
| MAE (m/s) | 1.88 | 1.82 |
| Mean bias (m/s) | +1.52 | +1.51 |
| RMSE (m/s) | 2.37 | 2.24 |

At the evaluation point, the conclusion did not change: MYNN2.5 is still closer to ERA5 than YSU. The YSU statistics and the point difference between schemes (1.01 m/s) were not affected.

The map of time-averaged differences changed more. The original showed MYNN2.5 slower almost everywhere, by up to about 1.5 m/s. The corrected map shows differences within about ±0.8 m/s, with MYNN2.5 faster over the southeast ocean. The double-counted late hours had pulled down the MYNN2.5 average.

## 2. Mismatched baseline in the wake analysis (Project 02)

**Symptom:** the 18 UTC "wake" map showed strong banded differences of several m/s far from any turbine.

**Cause:** wakes are measured by subtracting a no-turbine run from a turbine run, so the two runs should differ only by the turbines. File metadata showed that they did not:

| Run | Original location | `BL_PBL_PHYSICS` | Scheme |
|-----|-------------------|------------------|--------|
| With turbines | `ICA11/wrf_files` | 5 | MYNN |
| No turbines, original baseline | `ICA11/NWF` | 1 | YSU |
| No turbines, corrected baseline | `ICA07/Prof` | 5 | MYNN |

**Fix:** the wake maps were regenerated with the matching MYNN no-turbine run. It has the same grid, and the other physics settings compared also match.

**Effect** (113 m wind-speed difference over the mapped area, turbines minus no turbines):

| Time | Original RMS (m/s) | Corrected RMS (m/s) | Original mean (m/s) | Corrected mean (m/s) |
|------|--------------------|---------------------|---------------------|----------------------|
| 06 UTC | 1.15 | 1.08 | -0.65 | -0.64 |
| 18 UTC | 1.53 | 0.41 | -0.38 | -0.02 |

Most of the original 18 UTC signal came from the difference between schemes. The 06 UTC wake pattern was broadly right, but it is now much cleaner.

## Figure and wording changes

- Project 01: "ERA5 ground truth" became "ERA5 reference". ERA5 is a reanalysis that also supplied the model's boundary conditions, so it is not an independent observation.
- Project 01: moved labels that overlapped data.
- Project 02: values beyond the color scale are now shaded instead of appearing as white holes.
- Scripts and output files were renamed from assignment numbers to descriptive names.
