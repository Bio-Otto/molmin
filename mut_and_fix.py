"""
    Author: Halil ibrahim özdemir
    Loc: Marmara University / Bioengineering
"""

import pdbfixer
from openmm import app
import os, sys

Res_Syn_ID = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


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
        pdb_name = os.path.basename(pdb_path).split('.')[0]
        pdb_directory = os.path.dirname(pdb_path)
        mut_file_name = pdb_name + '_' + Res_Syn_ID[mut_region[0][0:3]] + mut_region[0].split('-')[1] + Res_Syn_ID[
            mut_region[0][-3:]] + '.pdb'

        if chain_id is None:
            mut_file_path = os.path.join(pdb_directory, mut_file_name)
        else:
            mut_file_path = os.path.join(pdb_directory, mut_file_name)

        try:
            fixer = pdbfixer.PDBFixer(pdb_path + '.pdb')
        except FileNotFoundError:
            fixer = pdbfixer.PDBFixer(pdbid=os.path.splitext(os.path.basename(pdb_path))[0])

        if chain_id is None:
            chains = list(set([chain.id for chain in fixer.topology.chains()]))
            print("While mutating the PDB ID {}, the program is searching for chain IDs {}".format(
                os.path.splitext(os.path.basename(pdb_path))[0].upper(), chains))
            for _chain in chains:
                try:
                    fixer.applyMutations(mut_region, _chain)
                except Exception as Chain_Error:
                    pass
                    #print(Chain_Error, "So, passing!")

        else:
            print("While mutating the PDB ID {}, the program is searching for chain IDs {}".format(
                os.path.splitext(os.path.basename(pdb_path))[0].upper(), chain_id))
            fixer.applyMutations(mut_region, chain_id)

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        with open(mut_file_path, 'w') as w_file:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, w_file, keepIds=True)

        return mut_file_path

    except Exception as error:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        print('Please Check Your Input Parameters !!')
