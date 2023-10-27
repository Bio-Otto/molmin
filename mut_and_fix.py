"""
    Author: Halil ibrahim özdemir
    Loc: Marmara University / Bioengineering
"""

import pdbfixer
from openmm import app
import os


def mutate(pdb_path, mut_region=None, chain_id=None):
    """
    Make a mutant protein easy.

        Parameters
        ----------

        pdb_path: Give your pdb whole path to this parameter

        mut_region : list of strings
            Each string must include the resName (original), index,
            and resName (target).  For example, ALA-133-GLY will mutate
            alanine 133 to glycine.

        chain_id : str
            Chain ID to apply mutation.

        Example
        ----------

        mutate('C:/Users/HIbrahim/Desktop/MolDynAnalyze/test/last.pdb', mut_region=['ASP-306-ARG'], chain_id='A')

    """

    try:
        print(pdb_path, mut_region, chain_id)
        pdb_name = os.path.basename(pdb_path).split('.')[0]
        pdb_directory = os.path.dirname(pdb_path)

        if chain_id is None:
            mut_file_name = pdb_name + '_chain' + '_all' + '_' + str(mut_region[0]) + '.pdb'
            mut_file_path = os.path.join(pdb_directory, mut_file_name)
        else:
            mut_file_name = pdb_name + '_chain' + chain_id + '_' + str(mut_region[0]) + '.pdb'
            mut_file_path = os.path.join(pdb_directory, mut_file_name)

        fixer = pdbfixer.PDBFixer(pdb_path+'.pdb')

        if chain_id is None:
            chains = [chain.id for chain in fixer.topology.chains()]
            for _chain in chains:
                fixer.applyMutations(mut_region, _chain)

        else:
            fixer.applyMutations(mut_region, chain_id)

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        with open(mut_file_path, 'w') as w_file:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, w_file, keepIds=True)

        return mut_file_path

    except ValueError as error:
        print(error)
        print('Please Check Your Input Parameters !!')


