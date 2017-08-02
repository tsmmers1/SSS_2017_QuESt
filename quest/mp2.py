import numpy as np



def mp2(wavefunction):
    """MP2 energy using conventional ERIs

    Parameters
    ----------
    wavefunction: Wavefunction object
        Input `wavefunction` should have MO coefficients and orbital energies in
        `wavefunction.arrays` dictionary.

    Returns
    --------
    mp2_energy: float
        The mp2 energy, scf_energy + mp2 correlation energy
    wavefunction: Wavefunction object
        On return this will be the input `wavefunction` with MP2 energy quantities added to the
        `wavefunction.energies` dictionary.

    """
    pass


def df_mp2(wavefunction):
    """MP2 energy using Density fitted ERIs

    Parameters
    ----------
    wavefunction: Wavefunction object
        Input `wavefunction` should have MO coefficients and orbital energies in
        `wavefunction.arrays` dictionary.

    Returns
    --------
    mp2_energy: float
        The mp2 energy, scf_energy + mp2 correlation energy
    wavefunction: Wavefunction object
        On return this will be the input `wavefunction` with MP2 energy quantities added to the
        `wavefunction.energies` dictionary.
    """
    pass


def _mo_transform(g, C, nocc):
    """
    Transforms the MP2 relevant subset of the ERI tensor to MO basis
    In the future, wavefunction class will have attributes g, C, and nocc
    so the only argument will need to be a wfn class instance
    C = wfn.C
    g = wfn.g
    nocc = wfn.nocc
    """
    O = slice(None, nocc)
    V = slice(nocc, None)
    g_iajb = np.einsum('pQRS, pP -> PQRS',
        np.einsum('pqRS, qQ -> pQRS',
        np.einsum('pqrS, rR -> pqRS',
        np.einsum('pqrs, sS -> pqrS', g, C[:,V]), C[:,O]), C[:,V]), C[:,O])
    return g_iajb


def _denom():
    pass


def _compute_conv_e(I_iajb, D_ijab):
    # OS = (ia|jb)(ia|jb)/(ei+ej-ea-eb)
    OS = np.einsum("iajb,iajb,ijab->", I_iajb, I_iajb, D_ijab)
    # build (ib|ja)
    I_ibja = I_iajb.swapaxes(1, 3)
    # SS = [(ia|jb)-(ib|ja)](ia|jb)/(ei+ej-ea-eb)
    SS = np.einsum("iajb,iajb,ijab->", (I_iajb - I_ibja), I_iajb, D_ijab)
    # opposite spin, same spin
    E_mp2 = OS - SS
    return E_mp2


def _compute_DF_e():
    pass
