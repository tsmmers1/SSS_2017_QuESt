import numpy as np
import scipy.sparse.linalg as spla
import scipy as sp
from hessian_builders import *



def ao_to_mo(tensor, transform):
    if len(g.shape) == 2:
        return transform.T @ tensor @ transform
    if len(g.shape) == 4:
        tensor = np.einsum('pqrs,pt->tqrs', tensor, transform)
        tensor = np.einsum('pqrs,qt->ptrs', tensor, transform)
        tensor = np.einsum('pqrs,rt->pqts', tensor, transform)
        tensor = np.einsum('pqrs,st->pqrt', tensor, transform)
        return tensor
    raise Exception("Not a 2 or 4 tensor")

def response(g, F, C, L, R, nocc):
    '''
    Calculates the CPHF response.
    Expects the eri 4-tensor, g,  the Fock matrix, F, the MO coeefcient matrix, C, the left and right response tensors, L and R, and the occupation number, nocc.
    '''

    F = ao_to_mo(F, C)
    g = ao_to_mo(g, C)

    nbas = F.shape[0]
    nvirt = nbas - nocc 
    opeartor_vs = 1
    if operator_vs == 0:
        E = get_E(F, g, nocc, nbas)
        E_inv_R = np.linalg.solve(E,R.T)
    elif operator_vs == 1:
        E_kappa = E_kappa_MO(F, g, nocc, nbas)
        E_inv_R = spla.cg(E_kappa, R)
    elif operator_vs == 2:
        E_kappa = E_kappa_AO(F, g, C, get_JK, nocc, nbas)
        E_inv_R = spla.cg(E_kappa, R)
    return L @ E_inv_R
