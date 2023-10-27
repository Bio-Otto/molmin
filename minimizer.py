import sys
from distutils.util import strtobool
from Classic_MD import *
from openmm import unit
from mut_and_fix import *

from openmm import *
import time

if __name__ == '__main__':
    import argparse
    import multiprocessing as mp

    ### PARSER
    parser = argparse.ArgumentParser(description='The Program fully applying molecule aa-uto fixing and minimizing. '
                                                 'OpenMM Molecular Dynamic Toolkit, which also supports the Cuda '
                                                 'platform.')

    parser.add_argument('-p', '--topology', type=str, help='Need *.pdb file for loading trajectory file',
                        required=True)

    parser.add_argument('-pff', '--protein_ff', choices=['amber03', 'amber10', 'amber96', 'amber99sb', 'amber99sbildn',
                                                         'charmm36'], default='amber96', nargs='?', type=str,
                        help='Protein Forcefield (The program defaultly will use "amber96" forcefield)', required=False)

    parser.add_argument('-wff', '--water_ff', choices=['tip3p', 'tip5p', 'spce', 'tip4pew'], default='tip3p', nargs='?',
                        type=str, help='Water Forcefield (The program defaultly will use "tip3p" forcefield)',
                        required=False)

    parser.add_argument('-ts', '--long_md_total_step', default=300000, nargs='?', type=int,
                        help='It is the total number of steps the simulation wants to run. (The program defaultly will '
                             'use "300000")', required=False)

    parser.add_argument('-temp', '--temperature', nargs='?', default=310, type=float, required=False,
                        help='The temperature unit is kelvin. (The program defaultly will use "310 Kelvin")')

    parser.add_argument('-minim_step', '--minimize_step', default=500, nargs='?', type=int, required=False,
                        help='The program defaultly will minimize system for 500 steps if mimimize option is not '
                             '"False"')

    parser.add_argument('-mut', '--mut_region', nargs='+', type=str, default=None, required=False,
                        help='The program defaultly will mutate selected residue')

    parser.add_argument('-mut_ch', '--mut_chain', nargs='?', type=str, default=None, required=False,
                        help='The program defaultly will mutate selected residue according to indicated chain')

    parser.add_argument("--mutate_only", action="store_true",
                        help="The program will perform mutation only")

    parsed = parser.parse_args()
    start_time = time.time()

    mut_file_path = None
    if parsed.mut_region and parsed.mut_chain is not None:
        if parsed.mutate_only:
            print("ONlY MUTATING")
            mut_file_path = mutate(pdb_path=parsed.topology, mut_region=parsed.mut_region, chain_id=parsed.mut_chain)

        elif not parsed.mutate_only:
            mut_file_path = mutate(pdb_path=parsed.topology, mut_region=parsed.mut_region, chain_id=parsed.mut_chain)
            Classic_MD_Engine(pdb_path=mut_file_path, protein_ff=parsed.protein_ff, water_ff=parsed.water_ff,
                              total_Steps=parsed.long_md_total_step, temp=parsed.temperature,
                              minimize_steps=parsed.minimize_step)

    elif parsed.mut_region is not None and parsed.mut_chain is None:

        if parsed.mutate_only:
            print("ONlY MUTATING -- > to ALLL")
            mut_file_path = mutate(pdb_path=parsed.topology, mut_region=parsed.mut_region, chain_id=None)

        elif not parsed.mutate_only:
            mut_file_path = mutate(pdb_path=parsed.topology, mut_region=parsed.mut_region, chain_id=None)
            Classic_MD_Engine(pdb_path=mut_file_path, protein_ff=parsed.protein_ff, water_ff=parsed.water_ff,
                              total_Steps=parsed.long_md_total_step, temp=parsed.temperature,
                              minimize_steps=parsed.minimize_step)

    elif parsed.mut_region is None and parsed.mut_chain is None:
        print("No Mutation, just minimization applying")
        Classic_MD_Engine(pdb_path=parsed.topology, protein_ff=parsed.protein_ff, water_ff=parsed.water_ff,
                          total_Steps=parsed.long_md_total_step, temp=parsed.temperature, minimize=False,
                          minimize_steps=parsed.minimize_step)

    else:
        print("You must specify mutation region and chain together !")
        sys.exit()
    # -->

    print("\n--- %s seconds ---" % (time.time() - start_time))
    # ------------------------------------------------------------------------------------------------------------ #
