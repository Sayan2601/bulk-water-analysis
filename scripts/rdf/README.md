# Radial Distribution Function (RDF) Analysis

This directory contains Python scripts for calculating radial distribution functions (RDFs) from molecular dynamics trajectories of bulk water using MDAnalysis.

The RDF characterizes the local structure of a liquid by describing the probability of finding neighboring atoms at a given distance from a reference atom relative to an ideal homogeneous fluid.

## Available Scripts

### `rdf_oo.py`

Computes the oxygen–oxygen (O–O) radial distribution function of bulk water.

The script:

- calculates the O–O RDF using periodic boundary conditions,
- writes the RDF to a CSV file,
- estimates the first-shell coordination number by integrating the RDF up to a user-defined cutoff corresponding to the first minimum.

**Output**

- `rdf_oo.csv`
- First-shell coordination number (printed to the terminal)

---

### `rdf_oh.py`

Computes the intermolecular oxygen–hydrogen (O–H) or hydrogen–oxygen (H–O) radial distribution function.

The reference and neighbor atom selections can be modified at the beginning of the script to calculate either:

- O(reference) – H(neighbor)
- H(reference) – O(neighbor)

Intramolecular O–H pairs are excluded so that only intermolecular structural correlations are included.

**Output**

- `rdf_oh.csv` or `rdf_ho.csv`

## Requirements

- Python 3
- MDAnalysis
- NumPy
- Pandas
- tqdm

## Input

The scripts require:

- GROMACS structure file (`md.gro`)
- GROMACS trajectory file (`md.xtc`)

Update the filenames in each script if different input files are used.
