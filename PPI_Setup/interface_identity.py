#!/usr/bin/python

import foldx, re, shutil, random, os, math, sys, subprocess, glob
import numpy as np
from Bio import *
from Bio.PDB import *
from interface_analysis import *

def main():
  args            =  sys.argv
  in_file         = args[1]
  distance_cutoff = float(sys.argv[2])
  start_structure = args[3]
  ancestral_file  = args[4]

  ancestral_structure1 = capture_pdb(start_structure[0:-4] + '_A.pdb', start_structure, 'A')
  ancestral_structure2 = capture_pdb(start_structure[0:-4] + '_C.pdb', start_structure, 'C')
  
  all_lines = (open(in_file, 'r')).readlines()
  all_data = []

  for a_line in all_lines:
    split = a_line.split('\t')
    split = [thing.strip() for thing in split]
    if split[0] == 'file':
      continue
    else:
      if int(split[1][1:-1]) <= 156:
        all_data.append([split[0], split[1], split[2], split[3], split[4], split[5], split[6], 'A'])
      else:
        all_data.append([split[0], split[1], split[2], split[3], split[4], split[5], split[6], 'C'])
        
  ancestral_lines =  (open(ancestral_file, 'r')).readlines()     
  ancestral_data = []
  
  for a_line in ancestral_lines:
    split = a_line.split('\t')
    split = [thing.strip() for thing in split]
    if split[0] == 'name':
      continue
    else:
      ancestral_data.append([split[0], split[1], split[2], split[3], split[4], split[5]])

  output = open('interface_identity.txt', 'w')
  to_file = 'name\tcount\tidentity\tinterface_identity\tnon_interface_identity\tevolved_interaction\tancestral_interaction\tstability\n'
  output.write(to_file)
  output.close()

  #Make binding comparison using ancestral chain A
  for i in range(0, len(all_data)):
      interval = 5
      score_ob = foldx.Scores()

      a_mutant = all_data[i]

      if float(a_mutant[1][1:-1]) > 156 and i % interval == 0: 
          new_structure2 = capture_pdb(a_mutant[0][0:-4] + '_C.pdb', a_mutant[0], 'C')
          p = PDBParser()
          structure = p.get_structure('temp', a_mutant[0])
          distances = calculate_ca_distance(structure, 'C')
        
          all_sites = (distances[distances[:, 1] <= 999][:, 0]).astype(int)
          identity = calc_identity(ancestral_structure2, new_structure2, all_sites)

          interface_sites = (distances[distances[:, 1] <= distance_cutoff][:, 0]).astype(int)
          interface_identity = calc_identity(ancestral_structure2, new_structure2, interface_sites)
        
          non_interface_sites = (distances[distances[:, 1] > distance_cutoff][:, 0]).astype(int)
          non_interface_identity = calc_identity(ancestral_structure2, new_structure2, non_interface_sites)
        
      elif i % interval == 0:
          new_structure1 = capture_pdb(a_mutant[0][0:-4] + '_A.pdb', a_mutant[0], 'A')
          p = PDBParser()
          structure = p.get_structure('temp', a_mutant[0])
          distances = calculate_ca_distance(structure, 'A')
        
          all_sites = (distances[distances[:, 1] <= 999][:, 0]).astype(int)
          identity = calc_identity(ancestral_structure1, new_structure1, all_sites)
        
          interface_sites = (distances[distances[:, 1] <= distance_cutoff][:, 0]).astype(int)
          interface_identity = calc_identity(ancestral_structure1, new_structure1, interface_sites)
        
          non_interface_sites = (distances[distances[:, 1] > distance_cutoff][:, 0]).astype(int)
          non_interface_identity = calc_identity(ancestral_structure1, new_structure1, non_interface_sites)
        
      else:
          continue

      output = open('interface_identity.txt', 'a')
      print(ancestral_data[i/5][1])
      print(all_data[i][2])
      to_file = ancestral_data[i/5][0] + '\t' + ancestral_data[i/5][1] + '\t' + str(identity) + '\t' + str(interface_identity) + '\t' + str(non_interface_identity) + '\t' + ancestral_data[i/5][3] + '\t' + ancestral_data[i/5][4] + '\t' + ancestral_data[i/5][5] + '\n'
      output.write(to_file)
      output.close()

def capture_pdb(out_name, fn, chain_letter):
  parser = PDBParser()
  structure = parser.get_structure('working_pdb', fn)
  
  writer = PDB.PDBIO()
  writer.set_structure(structure)
  writer.save(out_name, select=SelectChains(chain_letter))
  return(out_name)
    
def calc_identity(fn1, fn2, sites):
  parser = PDBParser()
  structure1 = parser.get_structure('working_pdb1', fn1)
  structure2 = parser.get_structure('working_pdb2', fn2)

  count = 0
  identical = 0

  ppb1 = PPBuilder()
  pp1 = (ppb1.build_peptides(structure1)[0]).get_sequence()

  ppb2 = PPBuilder()
  pp2 = (ppb2.build_peptides(structure2)[0]).get_sequence()
  
  for i in sites:
     if pp1[int(i)] == pp2[int(i)]:
       identical += 1
     count += 1

  return(float(identical)/float(count))

class SelectChains(Select):
  """ Only accept the specified chains when saving. """
  def __init__(self, chain_letters):
    self.chain_letters = chain_letters

  def accept_chain(self, chain):
    return (chain.get_id() in self.chain_letters)

def make_combined_file(files, path):
  with open(path, 'w') as outfile:
    for fname in files:
      with open(fname) as infile:
        for line in infile:
          outfile.write(line)
  subprocess.call("sed -i 's/END//g' " + path, shell=True)
  return(path)

def clean_up(delete_files):
  import glob
  for name in delete_files:
    for files in glob.glob(name):
      if os.path.isfile(files):
        os.remove(files)

def rmWS(string):
    return("".join(string.split()))

#Run main program
if __name__ == '__main__':
   main()
