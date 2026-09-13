# Input files

The analysis scripts expect molecular dynamics input files such as:

- `md.gro` — GROMACS structure/topology coordinate file
- `md.xtc` — GROMACS trajectory for structural analyses
- `unwrapped.xtc` — unwrapped trajectory for MSD analysis
- `dynamic.gro` / `dynamic.xtc` — trajectory used by the hydrogen-bond lifetime analysis

Large production trajectories are **not stored in this repository**. They can be several GB in size and are only needed locally when running the scripts.

For example, the current working data include trajectories such as `dynamic.xtc` and `unwrapped.xtc`, which are several GB and should remain outside GitHub.

Small example input files may be added here later if they are useful for reproducing a demonstration without making the repository unnecessarily large.
