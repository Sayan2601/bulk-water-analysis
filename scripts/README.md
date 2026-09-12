# Scripts

This directory contains Python scripts for post-processing and analysis of bulk water molecular dynamics (MD) simulations.

## Available Analysis

### RDF Analysis (`rdf/`)

* `rdf_oo.py` — Calculates the oxygen–oxygen (O–O) radial distribution function and first-shell coordination number.
* `rdf_oh.py` — Calculates the intermolecular oxygen–hydrogen (O–H) and hydrogen–oxygen (H–O) radial distribution functions.

### MSD Analysis (`msd/`)

* `msd.py` — Computes the mean squared displacement (MSD) of water oxygen atoms using multiple time-origin averaging.

### Hydrogen-Bond Analysis (`hbond/`)

* `avg_hbond.py` — Computes the average number of hydrogen bonds per water molecule using geometric hydrogen-bond criteria.

Further analysis workflows will be added progressively.
