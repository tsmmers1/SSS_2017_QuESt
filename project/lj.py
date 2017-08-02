import numpy as np

# .mol returns a psi4 geometry
# setmol can be given a string and it sets up a molecule

# atom - Input molecule (should be one atom)
# returnAll - Returns coefficients and energies/distances in addition to sigma if True
def lj(molecule, returnAll=False):
    # stuff to return (change later)
    sigma = 0
    A = 0
    B = 0
    # set up distances for PES
    start = 2.0
    stop = 10.0
    step = 0.25
    distances = np.arange(start, stop, step)
    energies = np.zeros(distances.size)
    # do MP2 on each distance in the array
    for i, distance in enumerate(distances):
        # construct molecule w/ correct distance
        #molecule_updated = .format()
        #molecule.set_geometry( # psi4 
        # call MP2 on molecule
        energy = distance
        # add MP2 energy to energies list
        energies[i] = energy
    if(returnAll):
        return sigma, A, B, energies, distances
    else:
        return sigma


def get_atom(molecule):
    # call
    # molecule.mol()
    # figure out how to get identity of molecule here!!
    return "He"
