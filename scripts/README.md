# Scripts

This directory contains Python scripts for post-processing and analyzing bulk water molecular dynamics (MD) simulations.

## Available Analysis

### RDF Analysis (`rdf/`)

- `rdf_oo.py` — Calculates the oxygen–oxygen (O–O) radial distribution function (RDF) and computes the first-shell coordination number.
- `rdf_oh.py` — Calculates the intermolecular oxygen–hydrogen (O–H) and hydrogen–oxygen (H–O) radial distribution functions.

### MSD Analysis (`msd/`)

- `msd_parallel.py` — Computes the mean squared displacement (MSD) of water oxygen atoms using multiple time-origin averaging with parallel processing.
