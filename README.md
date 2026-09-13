# `bulk-water-analysis`

> Molecular dynamics analysis workflows for structural and dynamical characterization of liquid water.

[![Python](https://img.shields.io/badge/Python-3.x-3776AB?logo=python\&logoColor=white)](https://www.python.org/)
[![MDAnalysis](https://img.shields.io/badge/MDAnalysis-trajectory%20analysis-222222)](https://www.mdanalysis.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## `scripts/`

### `rdf/`

Structural analysis through pair correlation functions.

* **O–O RDF**
* **O–H / H–O RDF**
* **First-shell coordination number**

### `hbond/`

Geometric hydrogen-bond analysis and hydrogen-bond dynamics.

* **Hydrogen-bond population**
* **Average H-bonds per water molecule**
* **Continuous H-bond lifetime correlation**
* **Intermittent H-bond lifetime correlation**

### `msd/`

Translational dynamics from trajectory coordinates.

* **Oxygen MSD**
* **Multiple time-origin averaging**
* **Statistical averaging across time origins**

---

## Example results

Representative figures produced by the analysis workflows are shown below.

<p align="center">
  <img src="examples/figures/rdf_example.png" width="45%">
  <img src="examples/figures/hbond_example.png" width="45%">
</p>

<p align="center">
  <img src="examples/figures/hb_lifetime_example.png" width="45%">
  <img src="examples/figures/msd_example.png" width="45%">
</p>

---

## Repository structure

```text
bulk-water-analysis/
│
├── scripts/
│   ├── rdf/
│   │   ├── rdf_oo.py
│   │   └── rdf_oh.py
│   │
│   ├── hbond/
│   │   ├── hbond_average.py
│   │   └── hbond_lifetime.py
│   │
│   └── msd/
│       └── msd.py
│
├── input/
│   └── README.md
│
├── examples/
│   └── figures/
│       ├── rdf_example.png
│       ├── hbond_example.png
│       ├── hb_lifetime_example.png
│       └── msd_example.png
│
├── requirements.txt
├── LICENSE
└── README.md
```

---

## Stack

```text
Python
├── MDAnalysis
├── NumPy
├── Pandas
├── Matplotlib
├── Joblib
└── tqdm
```

The analysis workflows are developed for atomistic molecular dynamics trajectories and are primarily used with GROMACS-generated data.

---

## Status

```text
[active development]

rdf        ██████████  complete
hbond      ██████████  complete
msd        ██████████  complete
additional analyses   → incoming
```

Further analysis workflows will be added as they are developed and finalized.

---

## License

MIT
