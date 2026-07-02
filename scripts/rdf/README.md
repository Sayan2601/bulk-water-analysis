# Radial Distribution Function (RDF) Analysis

This directory contains Python scripts for calculating radial distribution functions (RDFs) from molecular dynamics trajectories of bulk water.

## Available Scripts

- `rdf_oo.py` – Calculates the oxygen–oxygen RDF and estimates the first-shell coordination number.
- `rdf_oh.py` – Calculates the intermolecular oxygen–hydrogen and hydrogen-centered oxygen RDFs.

## Output

The scripts generate CSV files containing:

- Distance (Å)
- RDF, g(r)

The oxygen–oxygen RDF script also computes the coordination number by integrating g(r) up to the first minimum.
