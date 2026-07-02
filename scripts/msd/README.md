# Mean Squared Displacement (MSD) Analysis

This directory contains scripts for calculating the mean squared displacement (MSD) of water molecules from molecular dynamics trajectories.

## Available Scripts

- `msd_parallel.py` – Computes the MSD of water oxygen atoms using multiple time-origin averaging with parallel processing.

## Method

The MSD is averaged over multiple independent time origins to improve statistical sampling.

## Output

The script generates a CSV file containing:

- Time (ps)
- Mean squared displacement (Å²)
- Standard deviation (Å²)
