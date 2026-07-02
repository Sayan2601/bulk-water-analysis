"""
rdf_oo.py

Compute the oxygen–oxygen radial distribution function (RDF) and
estimate the first-shell coordination number from bulk water
molecular dynamics trajectories using MDAnalysis.
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm

# Input files
GRO_FILE = "md.gro"
XTC_FILE = "md.xtc"
OUTPUT_FILE = "rdf_oo.csv"

# Analysis parameters
R_MAX = 10.0          # Å
N_BINS = 200
FRAME_STEP = 1
R_CUT = 3.5           # Å

u = mda.Universe(GRO_FILE, XTC_FILE)
oxygen = u.select_atoms("resname SOL and name OW")
dr = R_MAX / N_BINS
bins = np.linspace(0.0, R_MAX, N_BINS + 1)
rdf_hist = np.zeros(N_BINS)
n_atoms = len(oxygen)
n_frames = 0
volumes = []

for ts in tqdm(u.trajectory[::FRAME_STEP], desc="Calculating O-O RDF"):
    n_frames += 1
    positions = oxygen.positions
    box = u.dimensions
    distances = distance_array(positions, positions, box=box)
    distances = distances[np.triu_indices(n_atoms, k=1)]
    hist, _ = np.histogram(distances, bins=bins)
    rdf_hist += hist
    volumes.append(np.prod(box[:3]))

r = 0.5 * (bins[:-1] + bins[1:])
shell_volume = 4.0 * np.pi * r**2 * dr
average_volume = np.mean(volumes)
rho = n_atoms / average_volume
g_r = rdf_hist / (shell_volume * n_frames * n_atoms * rho)
g_r *= 2.0
rdf_df = pd.DataFrame({"r (Å)": r,"g(r)": g_r})
rdf_df.to_csv(OUTPUT_FILE, index=False)

idx = np.where(r <= R_CUT)[0]
coordination_number = (4.0* np.pi* rho* np.trapz(g_r[idx] * r[idx]**2, r[idx]))
print(f"First-shell cutoff   : {R_CUT:.2f} Å")
print(f"Coordination number  : {coordination_number:.3f}")
