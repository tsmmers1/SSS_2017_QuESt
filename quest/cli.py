"""
QuESt command line interface.
"""
import os
import argparse
import yaml
from molecule import Molecule
from mollib import mollib
# from . import driver


parser = argparse.ArgumentParser(
    description="""
-------------------------------------------

  ██████╗         ███████╗███████╗
 ██╔═══██╗        ██╔════╝██╔════╝  ██║
 ██║   ██║██║  ██║█████╗  ███████╗██████╗
 ██║▄▄ ██║██║  ██║██╔══╝  ╚════██║╚═██╔═╝
 ╚██████╔╝╚█████╔╝███████╗███████║  ███║
  ╚══▀▀═╝  ╚════╝ ╚══════╝╚══════╝  ╚══╝

QuESt: Quantum Energy and Stuff
-------------------------------------------
    """,
    epilog="""
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter)

default_params = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'default_params.yml')

parser.add_argument('--molecule', '-m', type=str, default='h2o', metavar='',
                    help='Molecule file name (default: h2o)')
parser.add_argument('--parameters', '-p', type=str, default=default_params, metavar='',
                    help='Parameters file (default: parameters.yml)')
parser.add_argument('-qm', action='store_true', default=True,
                    help='Quantum mechanics (default: True)')
parser.add_argument('-mm', action='store_true', default=True,
                    help='Molecular mechanics (default: True)')

# Parse arguments
args = parser.parse_args()

# Read parameters
with open(args.parameters, 'r') as inp_file:
    params = yaml.load(inp_file)

# Create molecule object, assign basis set and name
mol = Molecule(mol=mollib[args.molecule], bas=params['qm']['basis_set'])
mol.name = os.path.splitext(os.path.basename(args.molecule))[0]
mol.mol.print_out()

# driver.compute_mp2(molecule, "aug-cc-pvdz")
