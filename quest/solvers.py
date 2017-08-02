"""
This is a file where DIIS and CG would go.
"""


def DIIS_step(state_list, error_list):
    """Routine that returns a guess for a new state vector
    using DIIS

    Parameters
    ----------
    state_list : list of numpy arrays
                 n previous states to form the new state from
    error_list : list of numpy arrays
                 n previous errors corresponding to `state_list`

    Returns
    -------
    new_state : numpy array
                new state array from DIIS update

    Notes
    -----
    Whichever method calls this routine must keep track of its `state_list` and
    `error_list` and update them.

    Examples
    --------
    new_state = solvers.DIIS_step(state_list, error_list)

    """

    if len(state_list) != len(error_list):
        raise Exception("State and error list sizes don't match")

    if len(state_list) < 2:
        raise Exception("Must have at least two previous state vectors for DIIS")

    B = -1 * np.ones((len(error_list) + 1, len(error_list) + 1))
    B[-1, -1] = 0
    for i in range(len(error_list)):
        for j in range(i, len(error_list)):
            ri_dot_rj = np.vdot(error_list[i], error_list[j])
            B[i, j] = B[j, i] = ri_dot_rj

    vec = np.zeros((len(error_list) + 1))
    vec[-1] = -1
    coeff = np.linalg.solve(B, vec)
    new_state = np.zeros_like(state_list[-1])

    for i, prev_state in enumerate(state_list):
        new_state += coeff[i] * prev_state

    return new_state
