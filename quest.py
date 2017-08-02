import argparse


parser = argparse.ArgumentParser(
    description="""
-------------------------------------------

  ██████╗         ███████╗███████╗
 ██╔═══██╗        ██╔════╝██╔════╝  ██║
 ██║   ██║██║  ██║█████╗  ███████╗██████╗
 ██║▄▄ ██║██║  ██║██╔══╝  ╚════██║╚═██╔═╝
 ╚██████╔╝╚█████╔╝███████╗███████║  ███║
  ╚══▀▀═╝  ╚════╝ ╚══════╝╚══════╝  ╚══╝

-------------------------------------------
    """,
    epilog="""
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--molecule', '-m', type=str, default='tests/water', metavar='',
                    help='Molecule file name (default: water)')
parser.add_argument('-qm', action='store_true', default=True,
                    help='Quantum mechanics (default: True)')
parser.add_argument('-mm', action='store_true', default=True,
                    help='Molecular mechanics (default: True)')


args = parser.parse_args()
print(args)
