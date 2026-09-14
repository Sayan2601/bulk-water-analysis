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

### `correlation/`

Rotational dynamics from molecular orientation.

* **O–H rotational correlation functions**
* **First-, second-, and third-rank correlation functions (P1, P2, P3)**
* **O–H bond-vector self correlation**
* **The same correlation procedure can be adapted to normalized molecular dipole vectors**

---

## Installation

Clone the repository and install the Python dependencies:

```bash
git clone https://github.com/Sayan2601/bulk-water-analysis.git
cd bulk-water-analysis
python -m pip install -r requirements.txt
```

The workflows are intended for Python 3.x and GROMACS-generated molecular dynamics trajectories.

---

## Input files

The scripts currently use fixed input filenames. Place the required files in the directory from which the script is run.

| Analysis | Required input | Output |
|---|---|---|
| O–O RDF | `md.gro`, `md.xtc` | `rdf_oo.csv` |
| O–H / H–O RDF | `md.gro`, `md.xtc` | `rdf_oh.csv` or `rdf_ho.csv` |
| H-bond average | `md.gro`, `md.xtc` | printed average |
| H-bond lifetime | `dynamic.gro`, `dynamic.xtc` | continuous/intermittent CSV files |
| MSD | `md.gro`, `unwrapped.xtc` | `msd.csv` |
| O–H rotational correlation | `dynamic.gro`, `dynamic.xtc` | `tcf-15ps-p1-p2-p3.csv` |

Large production trajectories are intentionally not stored in the repository. See [`input/`](input/) for details.

---

## Usage

Run the scripts from a directory containing the required input files.

### O–O RDF

```bash
python scripts/rdf/rdf_oo.py
```

Calculates the oxygen–oxygen RDF and prints the coordination number using the current first-shell cutoff.

### O–H RDF

```bash
python scripts/rdf/rdf_oh.py
```

The default calculation is intermolecular O–H. The atom selections in the script can be swapped to calculate H–O instead.

### Hydrogen-bond average

```bash
python scripts/hbond/hbond_average.py
```

Calculates the geometric hydrogen-bond population and reports the average number of hydrogen bonds per water molecule.

### Hydrogen-bond lifetime

```bash
python scripts/hbond/hbond_lifetime.py
```

Calculates continuous and intermittent hydrogen-bond lifetime correlation functions using multiple analysis windows.

### Mean-squared displacement

```bash
python scripts/msd/msd.py
```

Calculates oxygen MSD using multiple time origins. The trajectory must contain unwrapped coordinates.

### O–H rotational correlation

```bash
python scripts/correlation/orientational-correlation.py
```

Calculates first-, second-, and third-rank rotational correlation functions (P1, P2, and P3) of O–H bond vectors. The calculation correlates the same O–H vector across time origins and averages over the O–H vectors.

The same correlation procedure can be applied to normalized molecular dipole vectors by replacing the O–H vector construction while keeping the correlation calculation unchanged.

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
│   ├── msd/
│   │   └── msd.py
│   │
│   └── correlation/
│       └── orientational-correlation.py
│
├── input/
│   └── README.md
│
├── examples/
│   └── figures/
│       ├── README.md
│       ├── rdf_example.png
│       ├── hbond_example.png
│       ├── hb_lifetime_example.png
│       ├── msd_example.png
│       └── p1_p2_p3_correlation.png
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
├── SciPy
├── Seaborn
├── Joblib
└── tqdm
```

---

## Example figures

Representative outputs are collected in [`examples/figures/`](examples/figures/). The figure README displays the plots directly.

### Orientational correlation

![O–H rotational correlation example](examples/figures/p1_p2_p3_correlation.png)

Example first-, second-, and third-rank rotational correlation functions (P1, P2, and P3) of O–H bond vectors.

---

## Status

```text
[active development]

rdf          ██████████  complete
hbond        ██████████  complete
msd          ██████████  complete
correlation  ██████████  complete
additional analyses     → incoming
```

Further analysis workflows will be added as they are developed and finalized.

---

## License

MIT
