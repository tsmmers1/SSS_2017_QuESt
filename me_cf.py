import numpy as np


def make_f_mo(g, kappa, F, nocc, nvirt):
    """
    Make F_ia^kappa in MO basis
    """

    kappa.reshape((nocc, nvirt))
    # Einsum pieces of Eqn 14 in response handout
    one = np.einsum('ab,jb->ja', F[:nocc, :nocc], kappa)
    two = np.einsum('ij,jb->ib', F[nocc:, nocc:], kappa)
    thr = np.einsum('jb,iajb->ia', kappa, g) 
    fou = np.einsum('jb,ijab->ia', kappa, g)
    fiv = np.einsum('jb,bjai->ia', kappa, g)
    
    Fiakappa = - one + two + 4 * thr - fou - fiv

    Fiakappa.ravel()

    return Fiakappa




