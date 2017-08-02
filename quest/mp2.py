import numpy as np



def mp2(wfn):
    pass


def df_mp2(wfn):
    pass


def _mo_transform(orbital_energies):
    
    pass


def _denom(fock_matrix):
    #get energies from fock matrix)
    #multiply C and F to get 
    eps = wfn.arrays.get("EPSILON", None )
    if eps != None:
        nocc = 5 #must pull from wavefunction when implemented
        eocc = eps[:nocc]
        evir = eps[nocc:]
        e_denom = 1 / (eocc.reshape(-1, 1, 1, 1) - evir.reshape(-1, 1, 1) + eocc.reshape( -1, 1) - evir)
   else:
        raise Exception ("orbital energy array  'EPSILON' doesn't exist")
   return e_denom



def _compute_conv_e():
    pass


def _compute_DF_e():
    pass
