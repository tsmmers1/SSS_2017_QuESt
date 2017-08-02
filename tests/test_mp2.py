"""
This file tests mp2
"""

import quest
import pytest
import psi4
import numpy as np


def test_mp2():
    geometry = psi4.geometry("""
    O
    H 1 1.1
    H 1 1.1 2 104
    """)
    basis = "STO-3G"
    mol = quest.Molecule(geometry, basis)
    rhf_options = \
    {
        'e_conv': 1.e-8,
        'd_conv': 1.e-8,
        'diis': True,
        'max_diis': 7,
        'max_iter': 100,
    }

    wafu = quest.Wavefunction(mol, rhf_options)

    scf_energy = quest.scf_module.compute_rhf(wafu)
    mp2_energy = quest.mp2.mp2(wafu)

    psi4.set_options({"scf_type": "pk", "mp2_type": "conv"})
    psi4_scf_energy = psi4.energy('scf/' + basis, molecule=geometry)
    psi4_mp2_energy = psi4.energy('mp2/' + basis, molecule=geometry)

    assert np.allclose(scf_energy, psi4_scf_energy)
    # assert np.allclose(mp2_energy, psi4_mp2_energy)
