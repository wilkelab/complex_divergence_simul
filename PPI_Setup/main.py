#!/usr/bin/python

import foldx, re, shutil, random, os, math, sys, glob
import numpy as np
from Bio import *
import Bio.PDB as PDB

def main():
  prefix = '2eke'
      
  output = open('data.txt', 'w')
  to_file = 'mutant\tbinding\tstability1\tstability2\tprobability'
  output.write(to_file)
  output.close()

  all_kept_mutants = []
  all_mutants_tried = []

  output = open('all_mutants_tried.txt', 'w')
  to_file = 'count\tmutant\tstability1\tstability2\tbinding\n'
  output.write(to_file)
  output.close()

  count = 0

  foldx.runFoldxRepair(prefix, [prefix + '.bak'])
  score_ob = foldx.Scores()
  score_ob.cleanUp([])
  repair_file = glob.glob('RepairPDB_*pdb')
  if len(repair_file) == 1:
    shutil.move(repair_file[0], prefix + '.pdb')
  else:
    raise Exception('No output from RepairPDB.')

  for i in range(0, 10000):
    sys.stdout.flush()
    
    #Make sure the pdb exists
    if not os.path.isfile(prefix + '.pdb') and i > 0:
      all_kept_mutants = all_kept_mutants[0:-1]
      prefix = all_kept_mutants[-1]
      all_mutants_tried = all_mutants_tried[0:-1]
      count -= 1
      continue

    (mutation_code, site) = generate_mutation_code(prefix)
    foldx.runFoldxSimpleMutator(mutation_code, [prefix + '.pdb'])
    proceed = foldx.checkOutputMutator(prefix)

    #See if we got the files we needed from the mutator
    if not proceed:
      score_ob = foldx.Scores()
      score_ob.cleanUp(['*_*.pdb', '*energies*'])
      continue
      
    (new_mutant_name, old_mutant_name) = recode_mutant_pdb(mutation_code, site, prefix)
    foldx.runFoldxAnalyzeComplex(new_mutant_name[0:-4] + '_complex', [old_mutant_name, new_mutant_name])
    proceed = foldx.checkOutputAnalyzeComplex(new_mutant_name[0:-4])

    #See if we got the files we needed from Analyze Complex
    if not proceed:
      score_ob = foldx.Scores()
      score_ob.cleanUp(['*' + new_mutant_name[0:-4] + '*', '*energies*'])
      continue

    print('\n\nThe problem came after foldx.\n')

    score_ob = foldx.Scores()
    score_ob.parseAnalyzeComplex()

    ids = score_ob.getIds()
    binding = score_ob.getInteractionEnergies()
    stab1 = [score_ob.getStability1()[0], score_ob.getStability2()[0]]
    stab2 = [score_ob.getStability1()[1], score_ob.getStability2()[1]]
    
    probability = (calc_prob(binding)[1] * calc_prob(stab1)[1] * calc_prob(stab2)[1])

    print('\n\nThe problem came after probability calculation\n')
    
    all_mutants_tried.append(new_mutant_name[0:-4])
    count += 1

    output = open('all_mutants_tried.txt', 'a')
    to_file = str(count) + '\t' + str(ids[1]) + '\t' + str(stab1[1]) + '\t' + str(stab2[1]) + str(binding[1]) + '\t' + str(probability) + '\n'
    output.write(to_file)
    output.close()

    if random.random() < probability:
      print('\n\nPassing to the next round...\n')
      score_ob.cleanUp(['*energies*'])
      output = open('data.txt', 'a')
      to_file = '\n' + str(ids[1]) + '\t' + str(binding[1]) + '\t' + str(stab1[1]) + '\t' + str(stab2[1]) + '\t' + str(probability)
      output.write(to_file)
      output.close()
      prefix = new_mutant_name[0:-4]
      all_kept_mutants.append(prefix)
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
  
def calc_prob(data):
    ddG = float(data[1]) - float(data[0])
    if ddG <= 0.0:
      return((ddG, 1.0))
    else:
      return((ddG, math.exp(-ddG) * (ddG)))
      
def recode_mutant_pdb(mutation_code, site, prefix):
  recoded_mutant = mutation_code[0] + site + mutation_code[-1]
  
  files = glob.glob('*_' + prefix + '.pdb')

  for a_file in files:
    if foldx.rev_resdict[recoded_mutant[0]] in a_file:
      shutil.move(a_file, recoded_mutant + '.wt.pdb')

      #This line makes it so the energy is evaluated relative to the last best mutant
      #Comment the line out if you want to evaluate the energy versus the repositioned WT output from positionscan
      shutil.copy(prefix + '.pdb', recoded_mutant + '.wt.pdb')
    else:
      shutil.move(a_file, recoded_mutant + '.pdb')

  return(recoded_mutant + '.pdb', recoded_mutant + '.wt.pdb')

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

