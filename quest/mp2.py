import numpy as np



def mp2(wavefunction):
    """
    Conventional MP2 energy

    Takes a quest.wavefunction.Wavefunction object, returns a wavefunction object
    modifies the energy dict with the MP2_corr_E, MP2_E variables.
    """
    pass


def df_mp2(wavefunction):
    """
    Density-Fitted MP2 energy

    Takes a quest.wavefunction.Wavefunction object, returns a wavefunction object
    modifies the energy dict with the MP2_corr_E, MP2_E variables
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


def _denom(wfn):
    #get energies from fock matrix)
    #multiply C and F to get 
    eps = wfn.arrays.get("EPSILON", None )
    if eps != None:
        #must pull nocc from wavefunction when implemented
        nocc = 5
        eocc = eps[:nocc]
        evir = eps[nocc:]
        e_denom = 1 / (eocc.reshape(-1, 1, 1, 1) - evir.reshape(-1, 1, 1) + eocc.reshape( -1, 1) - evir)
   else:
        raise Exception ("orbital energy array  'EPSILON' doesn't exist")
   return e_denom



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
