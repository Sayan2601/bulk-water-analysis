"""
Calculate continuous or intermittent hydrogen-bond lifetime
correlation functions for bulk water.

Hydrogen bonds are identified using geometric distance and angle
criteria. Multiple time windows are used for statistical averaging.
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, calc_angles, calc_bonds, distance_array
from tqdm import tqdm

def calculate_hydrogen_bonds(u, donors_sel, hydrogens_sel, acceptors_sel, d_h_cutoff=1.2, d_a_cutoff=3.5, angle_cutoff=30,oh_cutoff=2.45):
    # donor-hydrogen pairs within the cutoff distance
    d_h_indices, _ = capped_distance(donors_sel.positions, hydrogens_sel.positions, max_cutoff=d_h_cutoff, box=u.dimensions)
    donors = donors_sel[d_h_indices.T[0]]
    hydrogens = hydrogens_sel[d_h_indices.T[1]]
    # donor-acceptor pairs within the cutoff distance
    d_a_indices, d_a_distances = capped_distance(donors.positions, acceptors_sel.positions, max_cutoff=d_a_cutoff, min_cutoff=1.0, box=u.dimensions)
    tmp_donors = donors[d_a_indices.T[0]]
    tmp_hydrogens = hydrogens[d_a_indices.T[0]]
    tmp_acceptors = acceptors_sel[d_a_indices.T[1]]
    # Calculate angles between hydrogen, donor, and acceptor
    h_d_a_angles = np.rad2deg(calc_angles(tmp_hydrogens.positions, tmp_donors.positions, tmp_acceptors.positions, box=u.dimensions))
    # Filter bonds based on angle cutoff
    hbond_indices = np.where(h_d_a_angles < angle_cutoff)[0]
    hbond_donors = tmp_donors[hbond_indices]
    hbond_hydrogens = tmp_hydrogens[hbond_indices]
    hbond_acceptors = tmp_acceptors[hbond_indices]
    hbond_distances = d_a_distances[hbond_indices]
    hbond_angles = h_d_a_angles[hbond_indices]
    # hydrogen bond tuples
    dist_arr = distance_array(hbond_hydrogens.positions, hbond_acceptors.positions, box=u.dimensions)
    hbond_tuples = list(zip(hbond_donors.resids, hbond_hydrogens.resids, hbond_acceptors.resids, dist_arr.diagonal()))
    # Filter tuples based on final distance cutoff
    filtered_hbond_tuples = [tup for tup in hbond_tuples if tup[3] <= oh_cutoff]
    return [(tup[0], tup[1], tup[2]) for tup in filtered_hbond_tuples]

gro = f'dynamic.gro'
xtc = f'dynamic.xtc' # 500ps xtc file with 2fs saving freq
u = mda.Universe(gro, xtc)
analysis_time=10 #ps
skipping_time=5  #ps
total_frame=len(u.trajectory)  
saving_frequency=u.trajectory.dt
total_simulation_time = int(total_frame * saving_frequency)
pieces_frame = int(analysis_time) * int(total_frame / total_simulation_time)
skipping_frame = int(skipping_time) * int(total_frame / total_simulation_time)
frame_pairs = []
start = 0
while start < total_frame-pieces_frame:
    end = start + pieces_frame
    frame_pairs.append((start, min(end, total_frame))) 
    start = end + skipping_frame
print(pieces_frame,skipping_frame)
print(frame_pairs)
print(len(frame_pairs))

def hb_lifetime(frame_pair,gro,xtc,cal_type):
    # Load the Universe
    u = mda.Universe(gro, xtc)
    water = u.select_atoms("resname SOL")
    ow = water.select_atoms("name OW")
    #~~~~~~~~~~~~~~~~~~~~initial-Frame~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    u.trajectory[frame_pair[0]]
    #~~~~~~~~~WATER SELECTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ini_resid_oxygen = ow.resids
    hydrogens_sel = u.select_atoms(f"resname SOL and name H* and resid {' '.join(map(str, ini_resid_oxygen))}")
    donors_sel = u.select_atoms(f"resname SOL and name OW and resid {' '.join(map(str, ini_resid_oxygen))}")
    acceptors_sel = u.select_atoms(f"resname SOL and name OW and resid {' '.join(map(str, ini_resid_oxygen))}")
    hbond_tuples_initial = calculate_hydrogen_bonds(u, donors_sel,hydrogens_sel, acceptors_sel,
        d_h_cutoff=1.2, d_a_cutoff=3.5,angle_cutoff=30,oh_cutoff=2.45)
    pair_presence = {pair: [] for pair in hbond_tuples_initial}
    #~~~~~~~~~~~~~~~~~~~~~~~~~~CURRENT_FRAME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for ts in tqdm(u.trajectory[frame_pair[0]:frame_pair[1] + 1]):
        curr_resid_oxygen = ow.resids
        curr_hydrogens_sel = u.select_atoms(f"resname SOL and name H* and resid {' '.join(map(str, curr_resid_oxygen))}")
        curr_donors_sel = u.select_atoms(f"resname SOL and name OW and resid {' '.join(map(str, curr_resid_oxygen))}")
        curr_acceptors_sel = u.select_atoms(f"resname SOL and name OW and resid {' '.join(map(str, curr_resid_oxygen))}")
        hbond_tuples_current = calculate_hydrogen_bonds(u, curr_donors_sel, curr_hydrogens_sel, curr_acceptors_sel,
                d_h_cutoff=1.2, d_a_cutoff=3.5, angle_cutoff=30, oh_cutoff=2.45)
        for pair in hbond_tuples_initial:
            pair_presence[pair].append(1 if pair in hbond_tuples_current else 0)
    pair_presence_df = pd.DataFrame(pair_presence)
    # Continuous or intermittent calculation
    if cal_type == 'continuous':
        cumulative_product = pair_presence_df.cumprod(axis=0)
        cumulative_product = cumulative_product.transpose()
        mean_correlation = cumulative_product.mean(axis=0).values.tolist()
    elif cal_type == 'intermittent':
        transpose_df=pair_presence_df.transpose()
        mean_correlation = transpose_df.mean(axis=0).values.tolist()
    else:
        raise ValueError("cal_type must be either 'continuous' or 'intermittent'")
    return mean_correlation

from joblib import Parallel, delayed
from tqdm import tqdm
cal_type='continuous'
c_hb_lifetime = Parallel(n_jobs=-1)(delayed(hb_lifetime)(frame_pair,gro,xtc,cal_type) for frame_pair in tqdm(frame_pairs))
df_c_hb_lifetime=pd.DataFrame(c_hb_lifetime)
c_hb_final_df=pd.DataFrame()
c_hb_final_df['time']=[i*(u.trajectory.dt) for i in range(len(df_c_hb_lifetime.columns))]
c_hb_final_df['c_hb']=df_c_hb_lifetime.mean(axis=0).values.tolist()
c_hb_final_df.to_csv(f'spce-300k-{cal_type}-HB-dyn-{analysis_time}ps.csv',index=False)
cal_type='intermittent'
i_hb_lifetime = Parallel(n_jobs=-1)(delayed(hb_lifetime)(frame_pair,gro,xtc,cal_type) for frame_pair in tqdm(frame_pairs))
df_i_hb_lifetime=pd.DataFrame(i_hb_lifetime)
i_hb_final_df=pd.DataFrame()
i_hb_final_df['time']=[i*(u.trajectory.dt) for i in range(len(df_i_hb_lifetime.columns))]
i_hb_final_df['c_hb']=df_i_hb_lifetime.mean(axis=0).values.tolist()
i_hb_final_df.to_csv(f'spce-300k_{cal_type}_HB_dyn_{analysis_time}ps.csv',index=False)
