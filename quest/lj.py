import numpy as np
import molecule
import psi4

# .mol returns a psi4 geometry
# setmol can be given a string and it sets up a molecule

# atom - Input molecule (should be one atom)
# returnAll - Returns coefficients and energies/distances in addition to sigma if True
def lj_fit(mol, returnAll=False):
    # set up distances for PES
    start = 2.0
    stop = 10.0
    step = 0.25
    distances = np.arange(start, stop, step)
    energies = np.zeros(distances.size)
    # set up geometry stuff -- make new geom object to avoid overwriting the old one
    atom_str = mol.mol().create_psi4_string_from_molecule()
    #atom_str = mol.create_psi4_string_from_molecule().splitlines()[2]
    print(atom_str)
    new_mol = molecule.Molecule(mol=None, bas=mol.bas())
    geom_string = '{:s}\n{:s} 1   {:2.5f}'
    # do MP2 on each distance in the array
    for i, distance in enumerate(distances):
        # construct molecule w/ correct distance
        mol_geom = geom_string.format(atom_str, atom_str, distance)
        new_mol.set_geometry(mol_geom)
        # call MP2 on molecule and get energy
        #rhf_object = RHF(new_mol, new_mol.basisset())
        #scf_wfn = rhf.compute_energy()
        # scf_wfn = scf( some arguments )
        # mp2_wfn = mp2(scf_wfn)
        energy = distance
        # add MP2 energy to energies list
        energies[i] = energy
    # doing the fit
    A,B = get_coeffs(distances, energies)
    # calculate sigma
    sigma = np.power(A/(-B), 1.0/6.0)
    # returning the correct values
    if(returnAll):
        return sigma, A, B, energies, distances
    else:
        return sigma

def get_coeffs(distances, energies):
    """
    Takes a list of distances and corresponding energies 
    perform a list squares fitting and 
    returns coefficients A for the r^12 term and B for the r^6 term
    """

    powers = [-12.0, -6.0]
    r_power = np.power(np.array(distances).reshape(-1,1), powers) 
    A,B = np.linalg.lstsq(r_power, energies)[0]
    return A,B

if __name__=='__main__':
    distances = [2.0,3.0,4.0,5.0,6.0,7.0]
    energies = [2.0,0.0,-1.0,-1.0,2.0,4.0]
    A,B = get_coeffs(distances, energies)
    print(A, B)
    mol_psi = psi4.geometry("He")
    test = molecule.Molecule(mol_psi, "sto-3g")
    s, a, b, e, d = lj_fit(test, True)
    print("sigma", s)
    print("A", a)
    print("B", b)


