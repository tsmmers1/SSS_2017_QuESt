"""
A basic wave function class
"""

import numpy as np
import psi4

class WaveFunction(object):
    """Basic wave function class"""

    def __init__(self, psi4mol, basisname="sto-3g", scf_type="df",
                 e_conv=1.e-8, d_conv=1.e-8):
        psi4.set_options({'basis': basisname,
                          'scf_type': scf_type,
                          'e_convergence': e_conv,
                          'd_convergence': d_conv})

        bas = psi4.core.BasisSet.build(psi4mol, target=basis)

        # Build the mints object
        self.mints = psi4.core.MintsHelper(bas)

        self.basisname = basisname
        self.scf_type = scf_type
        self.e_conv = e_conv
        self.d_conv = d_conv

        # Energies: scf_energy, mp2_os_corr, mp2_ss_corr, E_nuc
        self.energies = {}

        # coefficients, fock_matrix, density,
        # and whatever else you wanna put in it... communicate about it
        self.arrays = {}
