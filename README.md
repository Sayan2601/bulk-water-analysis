# `bulk-water-analysis`

> Molecular dynamics analysis workflows for structural and dynamical characterization of liquid water.

[![Python](https://img.shields.io/badge/Python-3.x-3776AB?logo=python\&logoColor=white)](https://www.python.org/)
[![MDAnalysis](https://img.shields.io/badge/MDAnalysis-trajectory%20analysis-222222)](https://www.mdanalysis.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## `analysis/`

### `rdf/`

Structural analysis through pair correlation functions.

* **O–O RDF**
* **O–H / H–O RDF**
* **First-shell coordination number**

### `hbond/`

Geometric hydrogen-bond analysis.

* **Hydrogen-bond population**
* **Average H-bonds per water molecule**

### `msd/`

Translational dynamics from trajectory coordinates.

* **Oxygen MSD**
* **Multiple time-origin averaging**
* **Statistical averaging across time origins**

---

## `repository/`

```text
bulk-water-analysis/
│
├── scripts/
│   ├── rdf/
│   ├── hbond/
│   └── msd/
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
hbond      ██████████  current
msd        ██████████  current
additional analyses   → incoming
```

Further analysis workflows will be added as they are developed and finalized.

---

## `license`

MIT
