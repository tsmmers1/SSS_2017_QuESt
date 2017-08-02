"""
This is to calculate the F matrix give kappa in AO basis

"""

from . import core

import numpy as np
import psi4  # not complete 

def E_kappa(kappa,C,F,nocc):
	Cocc = C[:nocc,:]
	kocc = kapp[:nocc,:]
	Focc  = F[:nocc,:nocc]
	V = np.einsum('ij,ik,kl->jl',Cocc,kocc,C,optimize=True) #Eq. (14) in the handout
	J,K = get_JK(V);

	JK = 4 * J - K - K.T; # denote the centor part of Eq.(16)
	FF = np.einsum('ij,jk->ik',Focc,kocc) - np.einsum('ij,jk->ik',kocc,F)

	FF +=np.einsum('ij,jk,lk->il',Cocc,JK,C)
	
	return FF




