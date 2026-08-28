"""
Calculate the mean squared displacement (MSD) of water oxygen atoms
using multiple time origins.

The trajectory must contain unwrapped coordinates so that molecular
displacements remain continuous across periodic boundaries.
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm

# Input/output files
GRO_FILE = "md.gro"
XTC_FILE = "unwrapped.xtc"
OUTPUT_FILE = "msd.csv"

# Analysis parameters
MSD_TIME = 50.0       # ps
FRAME_INTERVAL = 1000

u = mda.Universe(GRO_FILE, XTC_FILE)
oxygen = u.select_atoms("name OW")
dt = u.trajectory.dt
n_frames = len(u.trajectory)
frame_lag = int(MSD_TIME / dt)
if frame_lag >= n_frames:
    raise ValueError("MSD_TIME is longer than the trajectory duration.")

positions = np.empty((n_frames, len(oxygen), 3),dtype=np.float32)
for i, ts in enumerate(tqdm(u.trajectory, desc="Reading trajectory")):
    positions[i] = oxygen.positions

origins = range(0,n_frames - frame_lag,FRAME_INTERVAL)
msd_results = []
for start in tqdm(origins, desc="Calculating MSD"):
    trajectory = positions[start:start + frame_lag]
    initial = positions[start]
    displacement = trajectory - initial
    squared_displacement = np.sum(displacement**2, axis=2)
    msd = np.mean(squared_displacement, axis=1)
    msd_results.append(msd)

msd_results = np.asarray(msd_results)
mean_msd = np.mean(msd_results, axis=0)
std_msd = np.std(msd_results, axis=0)
time = np.arange(frame_lag) * dt

results = pd.DataFrame({
    "Time (ps)": time,
    "MSD (Å²)": mean_msd,
    "STD (Å²)": std_msd
})
results.to_csv(OUTPUT_FILE,index=False)

print(f"Average MSD calculated using {len(msd_results)} time origins.")
