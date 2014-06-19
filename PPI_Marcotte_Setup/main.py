#!/usr/bin/python

import foldx, re, shutil, random, os, math, sys
import numpy as np
from Bio import *
import Bio.PDB as PDB

def main():
  prefix = '2eke'
      
  output = open('data.txt', 'w')
  to_file = 'mutant\tbinding\tstability1\tstability2\tprobability'
  output.write(to_file)
  output.close()
  
  for i in range(0, 1000):
    sys.stdout.flush()
    (mutation_code, site) = generate_mutation_code(prefix)
    foldx.runFoldxSimpleMutator(mutation_code, [prefix + '.pdb'])
    (new_mutant_name, old_mutant_name) = recode_mutant_pdb(mutation_code, site, prefix)

    #Make sure both pdbs exist
    if not os.path.isfile(new_mutant_name) or not os.path.isfile(old_mutant_name):
      score_ob.cleanUp(['*' + new_mutant_name[0:-4] + '*'])
      continue

    score_ob = foldx.Scores()
    foldx.runFoldxAnalyzeComplex(new_mutant_name[0:-4] + '_complex', [old_mutant_name, new_mutant_name])
    print('\n\nThese are the names: ' + old_mutant_name + ' ' + new_mutant_name)
    score_ob.parseAnalyzeComplex()

    ids = score_ob.getIds()
    binding = score_ob.getInteractionEnergies()
    stab1 = [score_ob.getStability1()[0], score_ob.getStability2()[0]]
    stab2 = [score_ob.getStability1()[1], score_ob.getStability2()[1]]
    
    probability = (binding_probability(binding) * stability_probability1(stab1) * stability_probability2(stab2))

    if random.random() < probability:
      print('\n\nPassing to the next round...\n')
      score_ob.cleanUp([])
      output = open('data.txt', 'a')
      to_file = '\n' + str(ids[1]) + '\t' + str(binding[1]) + '\t' + str(stab1[1]) + '\t' + str(stab2[1]) + '\t' + str(probability)
      output.write(to_file)
      output.close()
      prefix = new_mutant_name[0:-4]
    else:
      print('\n\nMutation is being reverted...\n')
      score_ob.cleanUp(['*' + new_mutant_name[0:-4] + '*'])
  
def generate_mutation_code(prefix):
  start_name = prefix + '.pdb'

  parser = PDB.PDBParser()
  structure = parser.get_structure('working_pdb', start_name)
  
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

  site = random.randint(0, total_length - 1)
  
  chain = 0
  if site > first_chain_length - 1:
    chain = 1
    
  chain_letters = ''
  residue_numbers = []
  for chains in structure.get_chains():
    chain_letters += chains.get_id()
  for chains in structure.get_residues():
    residue_numbers.append(str(chains.get_id()[1]))

  mutation = total_sequence[site]
  
  while( mutation == total_sequence[site] ):
    mutation = random.choice(foldx.rev_resdict.keys())

  mutation_code = total_sequence[site] + chain_letters[chain] + residue_numbers[site] + mutation
  
  return(mutation_code, residue_numbers[site])
  
def binding_probability(binding):
    binding_ddG = float(binding[1]) - float(binding[0])
    if binding_ddG <= 0:
      return(1.0)
    elif binding_ddG > 0:
      return(math.exp(-binding_ddG) * (binding_ddG))
      
def stability_probability1(stab1):
    stab1_ddG = float(stab1[1]) - float(stab1[0])
    if stab1_ddG <= 0:
      return(1.0)
    elif stab1_ddG > 0:
      return(math.exp(-stab1_ddG) * (stab1_ddG))
      
def stability_probability2(stab2):
    stab2_ddG = float(stab2[1]) - float(stab2[0])
    if stab2_ddG <= 0:
      return(1.0)
    elif stab2_ddG > 0:
      return(math.exp(-stab2_ddG) * (stab2_ddG))
  
def recode_mutant_pdb(mutation_code, site, prefix):
  mutant = foldx.rev_resdict[mutation_code[-1]]
  recoded_mutant = mutation_code[0] + site + mutation_code[-1] + '.pdb'
  pdb_code = mutant + site + '_' + prefix + '.pdb'
  
  initial_mutant_name = prefix + '_1.pdb'
  wt_comparison_name = 'WT_' + prefix + '_1.pdb'
  if os.path.isfile(initial_mutant_name) and os.path.isfile(wt_comparison_name):
    shutil.move(initial_mutant_name, recoded_mutant)
    shutil.move(wt_comparison_name, recoded_mutant[0:-4] + '.wt.pdb')
    print('\n\nMoved both files.\n')
  else:
    raise Exception('\n\nOne of the two required files do not exist.\n')
  
  return(recoded_mutant, recoded_mutant[0:-4] + '.wt.pdb')

def capture_mutant_pdb(out_name, mutant, chain_letter):
  parser = PDB.PDBParser()
  structure = parser.get_structure('working_pdb', mutant)
  
  writer = PDB.PDBIO()
  writer.set_structure(structure)
  writer.save(out_name, select=SelectChains(chain_letter))
    
class SelectChains(PDB.Select):
  """ Only accept the specified chains when saving. """
  def __init__(self, chain_letters):
    self.chain_letters = chain_letters

  def accept_chain(self, chain):
    return (chain.get_id() in self.chain_letters)
        
#Run main program
if __name__ == '__main__':
   main()

