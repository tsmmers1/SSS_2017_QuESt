import numpy as np



def mp2(wfn):
    pass


def df_mp2(wfn):
    pass


def _mo_transform():
    pass


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
