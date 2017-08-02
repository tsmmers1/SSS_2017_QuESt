"""
A basic wave function class
"""

import psi4


class WaveFunction(object):
    """
    Basic wave function class

    Members:
        options: dictionary with various options, including the basis name.
        mints: contains a psi4 mints object built from the basis set.
        energies: a dictionary of various energies that have been
                  calculated (including scf, nuclear, mp2, etc.).
        arrays: a dictionary of various arrays (including a coefficient
                matrix, a fock matrix, etc.)
    """

    def __init__(self, mol):
        """
        Initialize the WaveFunction class.

        Args:
            mol: a Molecule class from the QuESt repository. Needs to
                 have a basis set at least.
        """

        # Set whichever options you would like. The only default is the
        # basis set name.
        self.options = {'basis_name': mol.bas.name}

        # Build the mints object
        self.mints = psi4.core.MintsHelper(mol.bas)

        # Energies: scf_energy, mp2_os_corr, mp2_ss_corr, E_nuc
        self.energies = {}

        # coefficients, fock_matrix, density,
        # and whatever else you wanna put in it... communicate about it
        self.arrays = {}
