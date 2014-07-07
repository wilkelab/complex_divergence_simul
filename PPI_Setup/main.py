#!/share/apps/python-2.7.2/bin/python

import foldx, re, shutil, random, os, math, sys, glob
from Bio import *
import Bio.PDB as PDB

def main():
  prefix = '2eke'
      
  output = open('data.txt', 'w')
  to_file = 'mutant\tcount\tbinding\tstability1\tstability2\tprobability'
  output.write(to_file)
  output.close()

  all_kept_mutants = []
  all_mutants_tried = []

  output = open('all_mutants_tried.txt', 'w')
  to_file = 'count\tmutant\tstability1\tstability2\tbinding\tprobability\n'
  output.write(to_file)
  output.close()

  #Should we select for binding on both, on one, or on none
  count = 0
  both = True
  one = True

  foldx.runFoldxRepair(prefix, [prefix + '.bak'])
  score_ob = foldx.Scores()
  score_ob.cleanUp([])
  repair_file = glob.glob('RepairPDB_*pdb')
  if len(repair_file) == 1:
    shutil.move(repair_file[0], prefix + '.pdb')
  else:
    raise Exception('No output from RepairPDB.')

  for i in range(0, 1000):
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

    #Declare the score parsing object
    score_ob = foldx.Scores()
    score_ob.parseAnalyzeComplex()

    #Grab the scores to be used in the probability calculations
    ids = score_ob.getIds()
    stab1 = [score_ob.getStability1()[0], score_ob.getStability2()[0]]
    stab2 = [score_ob.getStability1()[1], score_ob.getStability2()[1]]
    binding = score_ob.getInteractionEnergies()
    thresholds = [-19.7, -4.2, -10]
    
    if both:
      #To this function you need 6 variables: stab1, stab2, binding, N, beta, and threshold
      probability = calc_prob(stab1, stab2, binding, 10000, 1, thresholds)
    else:
      raise Exception("We're not doing both?")
    
    all_mutants_tried.append(new_mutant_name[0:-4])
    count += 1

    output = open('all_mutants_tried.txt', 'a')
    to_file = str(count) + '\t' + str(ids[1]) + '\t' + str(stab1[1]) + '\t' + str(stab2[1]) + '\t' + str(binding[1]) + '\t' + str(probability) + '\n'
    output.write(to_file)
    output.close()

    if random.random() < probability:
      print('\n\nPassing to the next round...\n')
      score_ob.cleanUp(['*energies*'])
      output = open('data.txt', 'a')
      to_file = '\n' + str(ids[1]) + '\t' + str(count) + '\t' + str(binding[1]) + '\t' + str(stab1[1]) + '\t' + str(stab2[1]) + '\t' + str(probability)
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
  
def calc_prob(stab1, stab2, binding, N, beta, thresholds):
  '''In other to use this function, you need to provide a number of parameters.
     The stab1, stab2, and binding data should be coming from the foldx values
     and they need to be ABSOLUTE energy not differences.  The N, beta and 
     threshold numbers need to specified for the theoretical population size,
     the beta distribution constant, and the soft threshold for survival of 
     each protein in the complex.
  
     At this point, the function cannot be used if binding on both chains is
     not desired.'''

  mutant = [stab1[1], stab2[1], binding[1]]
  origin = [stab1[0], stab2[0], binding[0]]

  xi = calc_x(origin, beta, thresholds)
  xj = calc_x(mutant, beta, thresholds)

  if xj > xi:
    return((1.0))
  else:
    #Need to make sure you check numbers that are too big for the math library
    exponent = -2 * float(N) * (xi - xj)
    
    return(safe_calc(exponent))

def calc_x(data, beta, thresholds):
  total = 0
  for i in range(0, len(data)):
    #Need to make sure you check numbers that are too big for the math library
    exponent = float(beta) * (float(data[i]) - float(thresholds[i]))

    total += -math.log(safe_calc(exponent) + 1)

  return(total)

def safe_calc(exponent):
  if exponent > 500:
    return(1.0)
  else:
    return(math.exp(exponent))
      
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

