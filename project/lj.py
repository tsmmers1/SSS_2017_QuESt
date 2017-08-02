import numpy as np
import psi4

# .mol returns a psi4 geometry
# setmol can be given a string and it sets up a molecule

# atom - Input molecule (should be one atom)
# returnAll - Returns coefficients and energies/distances in addition to sigma if True
def lj_fit(molecule, returnAll=False):
    # stuff to return (change later)
    sigma = 0.5
    A = 2
    B = 3
    # set up distances for PES
    start = 2.0
    stop = 10.0
    step = 0.25
    distances = np.arange(start, stop, step)
    energies = np.zeros(distances.size)
    # set up geometry stuff -- make new geom object to avoid overwriting the old one
    #atom_str = molecule.mol().create_psi4_string_from_molecule()
    atom_str = molecule.create_psi4_string_from_molecule().splitlines()[2]
    print(atom_str)
    new_mol = molecule
    geom_string = '{:s}\n{:s} 1   {:2.5f}'
    # do MP2 on each distance in the array
    for i, distance in enumerate(distances):
        # construct molecule w/ correct distance
        mol_geom = geom_string.format(atom_str, atom_str, distance)
        #new_mol.set_geometry(mol_geom)
        # call MP2 on molecule and get energy
        # scf_wfn = scf( some arguments )
        # mp2_wfn = mp2(scf_wfn)
        energy = distance
        # add MP2 energy to energies list
        energies[i] = energy
    # doing the fit
    if(returnAll):
        return sigma, A, B, energies, distances
    else:
        return sigma



