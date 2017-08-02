"""
This is a file where DIIS and CG would go.
"""

def DIIS_step(fock_list, r_list, diis_vectors = 6):
    B = -1*np.ones((len(r_list)+1, len(r_list)+1))
    B[-1,-1] = 0
    for i in range(len(r_list)):
        for j in range(i, len(r_list)):
            ri_dot_rj = np.vdot(r_list[i],r_list[j])
            B[i,j] = B[j,i] = ri_dot_rj
    vec = np.zeros((len(r_list)+1))
    vec[-1] = -1
    coeff =  np.linalg.solve(B, vec)
    F = np.zeros_like(fock_list[-1])
    for i, fock in enumerate(fock_list):
        F += coeff[i]*fock
    return F
