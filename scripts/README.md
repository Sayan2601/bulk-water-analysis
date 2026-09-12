# `scripts`

Python-based trajectory analysis workflows for bulk-water molecular dynamics simulations.

## `modules`

```text id="v3k8qm"
scripts/
│
├── rdf/
│   ├── rdf_oo.py
│   ├── rdf_oh.py
│   └── rdf_example.png
│
├── hbond/
│   └── avg_hbond.py
│
└── msd/
    └── msd.py
```

### `rdf/`

Pair-distribution analysis of liquid water.

* `rdf_oo.py` — O–O radial distribution function and first-shell coordination number.
* `rdf_oh.py` — intermolecular O–H / H–O radial distribution functions.

### `hbond/`

Geometric hydrogen-bond analysis.

* `avg_hbond.py` — frame-resolved average hydrogen-bond population per water molecule.

### `msd/`

Translational dynamics.

* `msd.py` — oxygen-atom mean squared displacement using multiple time origins.

## `dependencies`

```text
MDAnalysis
NumPy
Pandas
Matplotlib
Joblib
tqdm
```

The scripts are intended to remain modular and adaptable to different molecular dynamics trajectories and simulation conditions.

---

`status: active development`
