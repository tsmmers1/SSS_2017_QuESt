"""
An SCF module
"""

import numpy as np
import psi4


def compute_JK(g, D, version):
    if version == "conv":
        J = np.einsum("pqrs,rs->pq", g, D)
        K = np.einsum("prqs,rs->pq", g, D)
    else:
        raise KeyError("Key %s not recognized" % version)

    return J, K


def compute_rhf(wfn, df=True, diis=True, maxiter=25, e_conv=1.e-6, d_conv=1.e-6):
    """
    Add docs here!
    """
    nbf = wfn.mints.nbf()

    if (nbf > 100):
        raise Exception("More than 100 basis functions!")

    # Build core hamiltonian and overlap matrix
    V = np.array(wfn.mints.ao_potential())
    T = np.array(wfn.mints.ao_kinetic())
    H = V + T
    S = np.array(wfn.mints.ao_overlap())

    # Build ERI matrix
    g = np.array(wfn.mints.ao_eri())

    # Build orthogonalization matrix
    A = wfn.mints.ao_overlap()
    A.power(-0.5, 1.e-14)
    A = np.array(A)

    # Grab the number of doubly occupied orbs
    ndocc = wfn.options["ndocc"]

    # An internal diaognalize function
    def diag(F, A):
        Fp = A.T @ F @ A
        eps, Cp = np.linalg.eigh(Fp)
        C = A @ Cp
        return eps, C

    # Form a desnity
    eps, C = diag(H, A)
    Cocc = C[:, :ndocc]
    D = Cocc @ Cocc.T

    # Start with a few defaults
    F_list = []
    DIIS_grad = []
    grad_rms_list = []
    E_old = 0.0
    F_old = None

    # Roothan iterations
    for iteration in range(maxiter):
        J, K = compute_JK(g, D, "conv")

        # Fock Matrix and gradient
        F = H + 2.0 * J - K
        grad = F @ D @ S - S @ D @ F
        grad_rms = np.mean(grad**2)**0.5

        if diis:
            # Do DIIS
            F_list.append(F)
            diis_grad = A.T @ grad @ A
            grad_rms = np.mean(diis_grad**2)**0.5
            grad_rms_list.append(grad_rms)
            DIIS_grad.append(diis_grad)

            # Manipulate the DIIS graidents
            if len(diis_grad) > wfn.options['max_diis']:
                index = grad_rms_list.index(max(grad_rms_list[:-1]))
                F_list.pop(index)
                DIIS_grad.pop(index)
                grad_rms_list.pop(index)

            F = diis(F_list, DIIS_grad)

        else:
            # Use damping of the Fock matrix
            damp_value = 0.20
            damp_start = 5
            if iteration >= damp_start:
                F_new = F
                F = damp_value * F_old + (1.0 - damp_value) * F_new


        E = np.sum((F + H) * D)

        eps, C = diag(F, A)
        Cocc = C[:, :ndocc]
        D = Cocc @ Cocc.T

        if (E - E_old < e_conv) and (grad_rms < d_conv):
            break

        E_old = E
        F_old = F

    wfn.energies['scf_energy'] = E
    wfn.arrays['fock_matrix'] = F
    wfn.arrays['density'] = D
    wfn.arrays['coefficients'] = C

    return wfn.energies["scf_energy"]
    # wfn.arrays['grad_rms'] = grad_rms
