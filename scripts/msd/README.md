# MSD Analysis

This directory contains scripts for calculating the mean squared displacement (MSD) of bulk water trajectories using MDAnalysis.

## Script

### `msd_parallel.py`

Computes the mean squared displacement (MSD) of water oxygen atoms using multiple time-origin averaging with parallel processing.

### Method

The script:

* Loads an **unwrapped** molecular dynamics trajectory.
* Selects all water oxygen atoms (`OW`).
* Computes the displacement of each oxygen atom relative to its position at a chosen time origin.
* Calculates the mean squared displacement (MSD) for every frame within a specified lag time.
* Repeats the calculation for multiple independent time origins separated by a user-defined frame interval.
* Averages the MSD over all time origins to improve statistical sampling.
* Reports the standard deviation of the MSD across all time origins.

The final output is written as a CSV file containing:

* Time (ps)
* Mean MSD (Å²)
* Standard deviation (Å²)

## Input

* `md.gro`
* `unwrapped.xtc`

The trajectory **must be unwrapped** before running the analysis so that molecular displacements remain continuous across periodic boundaries.

## User Parameters

| Parameter        | Description                                                                |
| ---------------- | -------------------------------------------------------------------------- |
| `MSD_TIME`       | Maximum lag time (ps) over which the MSD is computed for each time origin. |
| `FRAME_INTERVAL` | Spacing (in frames) between successive time origins used for averaging.    |
| `N_CORES`        | Number of CPU cores used by Joblib (`-1` uses all available cores).        |

## Output

`msd.csv`

| Column    | Description                                              |
| --------- | -------------------------------------------------------- |
| Time (ps) | Lag time                                                 |
| MSD (Å²)  | Mean squared displacement averaged over all time origins |
| STD (Å²)  | Standard deviation across all time origins               |

## Notes

* The script assumes Cartesian coordinates are expressed in Å.
* Multiple time-origin averaging significantly reduces statistical noise compared with using a single reference frame.
* Longer values of `MSD_TIME` increase computational cost because each time origin is propagated over more trajectory frames.
* Smaller values of `FRAME_INTERVAL` produce more time origins and better statistics, at the expense of increased runtime.
