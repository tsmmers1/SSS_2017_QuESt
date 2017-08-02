import numpy as np
import scipy.sparse.linalg as spla

# modules that need interaction w/ other grps
from scf_utils import xform_2, xform_4

def get_E_bf(g, C, F, nbas, nocc):
    """
    Get four-tensor E in a brute-force way, reshape it in to
    nocc * nvirt-by-nocc * nvirt
    """
    nvirt = nbas - nocc
    F = xform_2(F, C)
    g = xform_4(g, C)
    print(nbas, nocc, F.shape, g.shape)

    # original version
    E0 = np.zeros([nbas, nbas, nbas, nbas])
    E0 += np.einsum('ij,ab->iajb', np.eye(nbas), F)
    E0 -= np.einsum('ab,ij->iajb', np.eye(nbas), F)
    E0 += 4 * np.einsum('aijb->iajb', g) - np.einsum("jiab->iajb", g) - \
        np.einsum("ajbi->iajb", g)
    E0 = E0[:nocc, nocc:, :nocc, nocc:]
    E0 = E0.reshape(nocc * nvirt, nocc * nvirt)

    # concise version
    """
    E = -F.reshape(1, nbas, 1, nbas)[:, :nvirt, :, :nvirt] + \
        F.reshape(nbas, 1, nbas, 1)[:nocc, :, :nocc, :] + \
        4. * g[:nocc, :nvirt, :nocc, :nvirt] - \
        g[:nocc, :nocc, :nvirt, :nvirt].swapaxes(1, 2) -\
        g[:nvirt, :nocc, :nvirt, :nocc].swapaxes(0, 3).swapaxes(1, 2)
    """

    E = E.reshape(nocc * nvirt, nocc * nvirt)

    return E0
