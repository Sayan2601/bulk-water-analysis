# Bulk Water Analysis

A collection of Python scripts for post-processing and analyzing molecular dynamics (MD) simulations of bulk water.

The repository contains analysis tools developed for studying the structural and dynamical properties of liquid water using MD trajectories. The scripts are primarily written using MDAnalysis and are designed to be easily adapted to different simulation systems.

## Features

Current analyses include:

- **Radial Distribution Function (RDF)**
  - Oxygen–oxygen (O–O) RDF
  - Intermolecular oxygen–hydrogen (O–H) RDF
  - Hydrogen-centered oxygen (H–O) RDF
  - First-shell coordination number

- **Mean Squared Displacement (MSD)**
  - Multiple time-origin averaging
  - Parallel processing using Joblib
  - Statistical averaging over independent time origins

## Repository Structure

```
bulk-water-analysis/
├── scripts/
│   ├── rdf/
│   └── msd/
├── LICENSE
└── README.md
```

## Requirements

- Python 3.10+
- MDAnalysis
- NumPy
- Pandas
- Matplotlib
- Joblib
- tqdm

Install the required packages with:

```bash
pip install MDAnalysis numpy pandas matplotlib joblib tqdm
```

## Input Files

The scripts are designed for molecular dynamics trajectories and typically require:

- Structure file (`.gro`)
- Trajectory file (`.xtc`)

Modify the input filenames in each script before running.

## Future Development

Planned additions include:

- Hydrogen-bond analysis
- Rotational time correlation functions
- Density profile calculations
- Diffusion coefficient analysis
- Structural order parameter analysis

## License

This project is released under the MIT License.
