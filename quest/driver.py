"""
Drives modules
"""

from . import mp2
from . import scf


def compute_mp2(molecule, basis="aug-cc-pvdz"):
	mp2.mp2(molecule, basis)


def compute_rhf(molecule, basis_name="aug-cc-pvdz", numpy_memory=1.e9, maxiter=12, E_conv=1.e-6, D_conv=1.e-4):
	RHF_object = scf.RHF(molecule, basis_name, numpy_memory)
	RHF_object.compute_energy(maxiter, E_conv, D_conv)


def compute_mc(sigma, epsilon, coord_file_path, params):
	pass
