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
* Example RDF plot

### `hbond/`

Geometric hydrogen-bond analysis and hydrogen-bond dynamics.

* **Hydrogen-bond population**
* **Average H-bonds per water molecule**
* **Continuous H-bond lifetime correlation**
* **Intermittent H-bond lifetime correlation**
* Example H-bond lifetime plot

### `msd/`

Translational dynamics from trajectory coordinates.

* **Oxygen MSD**
* **Multiple time-origin averaging**
* **Statistical averaging across time origins**
* Example MSD plot

---

## `repository/`

```text
bulk-water-analysis/
│
├── scripts/
│   ├── rdf/
│   │   ├── rdf_oo.py
│   │   ├── rdf_oh.py
│   │   └── rdf_example.png
│   │
│   ├── hbond/
│   │   ├── hbond_average.py
│   │   ├── hbond_lifetime.py
│   │   └── hb_lifetime_example.png
│   │
│   └── msd/
│       ├── msd_parallel.py
│       └── msd_example.png
│
├── LICENSE
└── README.md
```

---

## `stack`

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

## `status`

```text
[active development]

rdf        ██████████  complete
hbond      ██████████  complete
msd        ██████████  complete
additional analyses   → incoming
```

Further analysis workflows will be added as they are developed and finalized.

---

## `license`

MIT
