import pytest
import project
import psi4
import numpy as np
def test_get_coeffs():
    A = 148.0
    B = 1.0
    distances = np.arange(2.0, 10.0, 0.25)
    energies = A/(distances ** 12.0) - B/(distances ** 6.0)
    A, B = project.lj.get_coeffs(distances, energies)
    assert np.isclose(A, 148.0) and np.isclose(B, -1.0)

def test_lj_fit():
    mol = psi4.geometry("He")
    assert project.lj.lj_fit(mol) == 0.5
