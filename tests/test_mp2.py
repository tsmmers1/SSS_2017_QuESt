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
    mol = quest.Molecule(geometry,basis)
    wafu = quest.Wavefunction(mol,None)
    scf_energy = quest.scf(wafu,mol.nel,diis = False)
    mp2_energy = quest.mp2(wafu)

    psi4_energy = psi4.energy('mp2/STO-3G', molecule = geometry)

    assert np.allclose(mp2_energy,psi4_energy)
