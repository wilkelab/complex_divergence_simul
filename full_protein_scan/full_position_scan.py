import re, shutil, random, os, math, sys, glob, subprocess
from Bio import *
import Bio.PDB as PDB

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

rev_resdict = { 'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', \
                'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', \
                'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', \
                'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}


def main():
  parser = PDB.PDBParser()
  structure = parser.get_structure('working_pdb', '2eke_optimized.bak')

  total_length = 0
  total_sequence = ''

  count = 0
  ppb = PDB.PPBuilder()
  for pp in ppb.build_peptides(structure):
    total_length += len(pp.get_sequence())
    total_sequence += pp.get_sequence()
    if count == 0:
      first_chain_length = total_length
    count += 1
 
  chain_letters = ''
  residue_numbers = []
  for chains in structure.get_chains():
    chain_letters += chains.get_id()
  for chains in structure.get_residues():
    residue_numbers.append(str(chains.get_id()[1]))

  for i in range(0, len(residue_numbers)):
   if int(residue_numbers[i]) > 156:
     mutant = total_sequence[i] + chain_letters[1] + residue_numbers[i] + 'a'
   else:
     continue

   runFoldxSimpleMutator(mutant, ['2eke_optimized.bak'])
    

def runFoldxSimpleMutator(mutant, pdbs):
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()

  name = 'run_' + str(mutant) + '.foldx'
  makeFoldxPositionScan(mutant, name, 'false')
  command = '/home/austin/local/bin/foldx64Linux -runfile ' + name
  status = -1
  while status < 0:
    status = subprocess.call(command, shell=True)
  print('\n\nMutation exit status = ' + str(status))

def makeFoldxPositionScan(mutant, name, output_pdb):
  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<PositionScan>#,' + mutant + ';\n'+\
            '<END>#;\n'+\
            '<OPTIONS>FOLDX_optionfile;\n'+\
            '<Temperature>298;\n'+\
            '<R>#;\n'+\
            '<pH>7;\n'+\
            '<IonStrength>0.050;\n'+\
            '<water>-CRYSTAL;\n'+\
            '<metal>-CRYSTAL;\n'+\
            '<VdWDesign>2;\n'+\
            '<OutPDB>' + output_pdb + ';\n'+\
            '<pdb_hydrogens>;\n'+\
            '<END>#;\n'+\
            '<JOBEND>#;\n'+\
            '<ENDFILE>#;\n'

  output = open(name, 'w')
  output.write(to_file)
  output.close()

#Run main program
if __name__ == '__main__':
   main()


