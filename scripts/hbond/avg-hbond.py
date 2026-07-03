import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, calc_angles, calc_bonds, distance_array
from tqdm import tqdm

def calculate_hydrogen_bonds(u,donors_sel,hydrogens_sel,acceptors_sel,params,return_tuple=False):
    d_h_cutoff = params["d_h_cutoff"]
    d_a_cutoff = params["d_a_cutoff"]
    angle_cutoff = params["angle_cutoff"]
    oh_cutoff = params["oh_cutoff"]

    lower_angle = params["lower_angle"]
    higher_angle = params["higher_angle"]

    lower_dist = params["lower_dist"]
    higher_dist = params["higher_dist"]

    # Donor-Hydrogen pairs
    d_h_indices, _ = capped_distance(donors_sel.positions,hydrogens_sel.positions,max_cutoff=d_h_cutoff,box=u.dimensions)
    donors = donors_sel[d_h_indices.T[0]]
    hydrogens = hydrogens_sel[d_h_indices.T[1]]

    # Donor-Acceptor pairs
    d_a_indices, d_a_distances = capped_distance(
        donors.positions,acceptors_sel.positions,max_cutoff=d_a_cutoff,min_cutoff=1.0,box=u.dimensions)
    if len(d_a_indices) == 0:
        if return_tuple:
            return [(donor.index, hydrogen.index, None) for donor, hydrogen in zip(donors, hydrogens)]

        return [(donor.index, hydrogen.index, None, None, None, None) for donor, hydrogen in zip(donors, hydrogens)]

    tmp_donors = donors[d_a_indices.T[0]]
    tmp_hydrogens = hydrogens[d_a_indices.T[0]]
    tmp_acceptors = acceptors_sel[d_a_indices.T[1]]

    h_d_a_angles = np.rad2deg(
        calc_angles(tmp_hydrogens.positions,tmp_donors.positions,tmp_acceptors.positions,box=u.dimensions))
    hbond_indices = np.where(h_d_a_angles < angle_cutoff)[0]

    if len(hbond_indices) == 0:
        if return_tuple:
            return [(donor.index, hydrogen.index, None)
                for donor, hydrogen in zip(donors, hydrogens)]
        return [
            (donor.index, hydrogen.index, None, None, None, None)
            for donor, hydrogen in zip(donors, hydrogens)]

    hbond_donors = tmp_donors[hbond_indices]
    hbond_hydrogens = tmp_hydrogens[hbond_indices]
    hbond_acceptors = tmp_acceptors[hbond_indices]
    hbond_distances = d_a_distances[hbond_indices]
    hbond_angles = h_d_a_angles[hbond_indices]

    oh_distances = distance_array(
        hbond_hydrogens.positions,
        hbond_acceptors.positions,
        box=u.dimensions).diagonal()

    hbond_tuples = list(zip(
        hbond_donors.indices,
        hbond_hydrogens.indices,
        hbond_acceptors.indices,
        hbond_distances,
        hbond_angles,
        oh_distances))

    filtered_hbond_tuples = [
        tup for tup in hbond_tuples
        if (lower_dist < tup[3] < higher_dist)
        and (lower_angle < tup[4] < higher_angle)]

    if len(filtered_hbond_tuples) == 0:
        if return_tuple:
            return [(donor.index, hydrogen.index, None) for donor, hydrogen in zip(donors, hydrogens)]
        return [(donor.index, hydrogen.index, None, None, None, None) for donor, hydrogen in zip(donors, hydrogens)]

    if return_tuple:
        return [(tup[0],tup[1],tup[2],tup[3],tup[4]) for tup in filtered_hbond_tuples]
    return filtered_hbond_tuples

ini_hbond_params = {"d_h_cutoff": 1.2,"d_a_cutoff": 3.5,"angle_cutoff": 30,"oh_cutoff": 2.45,"lower_angle": 0,"higher_angle": 30,"lower_dist": 0,"higher_dist": 3.5}

gro = f'md.gro'
xtc = f'md.xtc'
u = mda.Universe(gro, xtc)  
avg_hbonds = []
for ts in tqdm(u.trajectory):
    donors = u.select_atoms(f"resname SOL and name OW")
    hydrogens = u.select_atoms(f"resname SOL and name HW1 HW2")
    acceptors = donors
    hbond_tuples = calculate_hydrogen_bonds(u,donors,hydrogens,acceptors,ini_hbond_params,return_tuple=True)
    df = pd.DataFrame(hbond_tuples,columns=["donor", "hydrogen", "acceptor", "dist", "angle"])
    counts = []
    for idx in oxygen_indices:
        k = df[(df["donor"] == idx) |(df["acceptor"] == idx)]
        counts.append(len(k))
    avg_hbonds.append([ts.frame, np.mean(counts)])
dfs = pd.DataFrame(avg_hbonds,columns=["Frame", "Average H-bonds"])
print(df["Average H-bonds"].mean())        
