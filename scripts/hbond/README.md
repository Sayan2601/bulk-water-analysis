# Average Hydrogen Bonds

This script computes the average number of hydrogen bonds per water molecule throughout a molecular dynamics trajectory using MDAnalysis.

Hydrogen bonds are identified using geometric criteria based on donor–acceptor distance and hydrogen-bond angle. For each trajectory frame, the total number of hydrogen bonds involving every water oxygen atom is counted, and the average number of hydrogen bonds per water molecule is calculated.

The script outputs the frame-wise average hydrogen-bond count, allowing analysis of the temporal evolution and overall average hydrogen-bonding behavior in bulk water.

## Features

- Calculates hydrogen bonds using geometric criteria.
- Supports periodic boundary conditions.
- Computes the average hydrogen bonds per water molecule for every frame.
- Saves the results as a CSV file.

## Input

- `md.gro`
- `md.xtc`

| Column | Description |
|--------|-------------|
| Frame | Trajectory frame number |
| Average H-bonds | Average number of hydrogen bonds per water molecule |

## Hydrogen-bond criteria

The default geometric criteria are

- Donor–acceptor distance ≤ 3.5 Å
- Hydrogen-bond angle ≤ 30°

These values can be modified in the `ini_hbond_params` dictionary.
