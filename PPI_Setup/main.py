#!/share/apps/python-2.7.2/bin/python

import foldx, re, shutil, random, os, math, sys, glob
from Bio import *
import Bio.PDB as PDB

def main():

  args =  sys.argv
  prefix          = args[1]
  population_size = float(args[2])
  beta            = float(args[3])
  dGt1            = float(args[4])
  dGt2            = float(args[5])
  dGt3            = float(args[6])
  out_file        = args[7]
      
  output = open(out_file, 'w')
  to_file = 'file\tmutant\tcount\tbinding\tstability1\tstability2\tprobability\n'
  output.write(to_file)
  output.close()

  all_kept_mutants = []
  all_mutants_tried = []

  output = open('all_mutants_tried.txt', 'w')
  to_file = 'file\tmutant\tcount\tstability1\tstability2\tbinding\tprobability\n'
  output.write(to_file)
  output.close()

  #Should we select for binding on both, on one, or on none
  count = 0
  both = True
  one = True

  #foldx.runFoldxRepair(prefix, [prefix + '.bak'])
  #score_ob = foldx.Scores()
  #score_ob.cleanUp([])
  #repair_file = glob.glob('RepairPDB_' + prefix + '*.pdb')
  #if len(repair_file) == 1:
  #  shutil.move(repair_file[0], prefix + '.pdb')
  #else:
  #  raise Exception('No output from RepairPDB.')

  output_dict = {}

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

    (new_mutant_name, old_mutant_name) = recode_mutant_pdb(mutation_code, site, prefix)
 
    foldx.runFoldxRepair(new_mutant_name[0:-4], [new_mutant_name])
    repair_file = glob.glob('RepairPDB_' + new_mutant_name[0:-4] + '*.pdb')
    
    shutil.move(repair_file[0], new_mutant_name)
      
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
    thresholds = [dGt1, dGt2, dGt3]
    
    if both:
      #To this function you need 6 variables: stab1, stab2, binding, N, beta, and threshold
      probability = calc_prob(stab1, stab2, binding, population_size, beta, thresholds)
    else:
      raise Exception("We're not doing both?")
    
    all_mutants_tried.append(new_mutant_name[0:-4])

    count += 1

    output = open('all_mutants_tried.txt', 'a')
    to_file = str(count) + '.pdb' + '\t' + str(ids[1][0:-4]) + '\t' + str(count) + '\t' + str(stab1[1]) + '\t' + str(stab2[1]) + '\t' + str(binding[1]) + '\t' + str(probability) + '\n'
    output.write(to_file)
    output.close()

    if random.random() < probability:
      print('\n\nPassing to the next round...\n')
      score_ob.cleanUp(['*energies*', 'WT_*'])
      
      output = open(out_file, 'a')
      to_file = str(count) + '.pdb' + '\t' + str(ids[1][0:-4]) + '\t' + str(count) + '\t' + str(stab1[1]) + '\t' + str(stab2[1]) + '\t' + str(binding[1]) + '\t' + str(probability) + '\n'
      output.write(to_file)
      output.close()

      shutil.move(new_mutant_name, str(count) + '.pdb')
      shutil.move(old_mutant_name, str(count) + '.wt.pdb')
      prefix = str(count)
      all_kept_mutants.append(new_mutant_name[0:-4])
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
  if exponent > 700:
    return(sys.float_info.max)
  else:
    return(math.exp(exponent))
      
def recode_mutant_pdb(mutation_code, site, prefix):
  recoded_mutant = mutation_code[0] + site + mutation_code[-1]

  new_test = recoded_mutant + '.pdb'
  old_test = recoded_mutant + '.wt.pdb'
  existing = glob.glob(recoded_mutant)

  if len(existing)/2 > 0:
    shutil.move(new_test, new_mutant_name[0:-4] + '_' + str(len(existing)/2) + '.pdb')
    shutil.move(old_test, new_mutant_name[0:-4] + '_' + str(len(existing)/2) + '.wt.pdb')
  
  shutil.copy(prefix + '.pdb', recoded_mutant + '.wt.pdb')  
  shutil.move(foldx.rev_resdict[mutation_code[-1]] + site + '_' + prefix + '.pdb', new_test)
  
  #Remove the unused file that is output from position scan
  old_files = glob.glob('*_' + prefix + '.pdb')
  for a_file in old_files:
    os.remove(a_file)

  return(new_test, old_test)

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

