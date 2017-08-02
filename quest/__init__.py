"""
The primary init for the project.
"""

from . import molecule
from . import scf
from . import jk
from . import solvers
from . import core
from . import mp2
from . import mollib
from . import cli


# Make sure Psi4 respects the global OMP_NUM_THREADS
import psi4
import os
if "OMP_NUM_THREADS" in list(os.environ):
    psi4.set_num_threads(int(os.environ["OMP_NUM_THREADS"]))
