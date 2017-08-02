import psi4

class Molecule(object):
    '''
    A molecule class that is used to store basic parameters including
    psi4 geometry, basis set, number of electrons, and number of doubly
    occupied orbitals.
    '''

    def __init__(self, mol=None, bas=None):
        '''
        Initialize molecule with geometry mol and basis
        set bas. Calculate number of electrons nel and
        number of doubly occupied orbitals ndocc.
        '''
        self.mol = mol
        self.bas = bas

        # Calculate the number of doubly occupied orbitals and the number of electrons
        self.nel = sum(self.molecule.Z(n) for n in range(self.molecule.natom()))
        self.nel -= self.molecule.molecular_charge()
        
        if not (self.nel / 2.0).is_integer():
            raise ValueError("Molecule must have an even number of electrons to perform RHF.")
        self.ndocc = int(self.nel / 2.0)

    def set_geometry(self, geom_str):
        '''
        Set the geometry of the molecule by providing
        a psi4 style string or psi4 geometry.
        '''
        if isinstance(mol, geom_str):
            self.mol = psi4.geometry(mol)
        elif isinstance(mol, psi4.core.Molecule):
            self.mol = mol
        else:
            raise TypeError("Input molecule of type %s is not understood" % type(mol))

        self.mol.update_geometry()

    def set_basis(bas_str):
        '''
        set the basis set that is used in any further calculations
        '''
        if(self.mol is None):
            raise Exception('Error: Geometry not defined')
        else:
            self.bas = psi4.core.BasisSet.build(self.mol, target=bas_str)

