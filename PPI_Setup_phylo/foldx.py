import subprocess, re, glob, os
import numpy as np

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' } 
            
rev_resdict = { 'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', \
                'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', \
                'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', \
                'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'} 

def makeFoldxBuildModel(mutant, name, output_pdb):
  output = open('individual_list.txt', 'w')
  output.write(mutant + ';')
  output.close()
  
  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<BuildModel>test.out,individual_list.txt;\n'+\
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

def makeFoldxPositionScan(mutant, name, output_pdb):
  output = open('individual_list.txt', 'w')
  output.write(mutant + ';')
  output.close()

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
  
def runFoldxSimpleMutator(mutant, pdbs):
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()
  
  name = 'run_' + str(mutant) + '.foldx'
  makeFoldxPositionScan(mutant, name, 'true')
  command = '/home/agm854/local/bin/foldx64Linux -runfile ' + name
  status = -1
  while status < 0:
    status = subprocess.call(command, shell=True)
  print('\n\nMutation exit status = ' + str(status))

def makeFoldxStability(name):
  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<Stability>Stability.txt;\n'+\
            '<END>#;\n'+\
            '<OPTIONS>FOLDX_optionfile;\n'+\
            '<Temperature>298;\n'+\
            '<R>#;\n'+\
            '<pH>7;\n'+\
            '<IonStrength>0.050;\n'+\
            '<water>-CRYSTAL;\n'+\
            '<metal>-CRYSTAL;\n'+\
            '<VdWDesign>2;\n'+\
            '<OutPDB>false;\n'+\
            '<pdb_hydrogens>;\n'+\
            '<END>#;\n'+\
            '<JOBEND>#;\n'+\
            '<ENDFILE>#;\n'

  output = open(name, 'w')
  output.write(to_file)
  output.close()

def runFoldxStability(name, pdbs):
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()
  
  name = 'run_' + name + '.foldx'
  makeFoldxStability(name)
  command = '/home/agm854/local/bin/foldx64Linux -runfile ' + name
  status = -1
  while status < 0:
    status = subprocess.call(command, shell=True)
  print('\n\nStability exit status = ' + str(status))
  
def makeFoldxAnalyzeComplex(name):
  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<AnalyseComplex>#;\n'+\
            '<END>#;\n'+\
            '<OPTIONS>FOLDX_optionfile;\n'+\
            '<Temperature>298;\n'+\
            '<R>#;\n'+\
            '<pH>7;\n'+\
            '<IonStrength>0.050;\n'+\
            '<water>-CRYSTAL;\n'+\
            '<metal>-CRYSTAL;\n'+\
            '<VdWDesign>2;\n'+\
            '<OutPDB>false;\n'+\
            '<pdb_hydrogens>;\n'+\
            '<END>#;\n'+\
            '<JOBEND>#;\n'+\
            '<ENDFILE>#;\n'

  output = open(name, 'w')
  output.write(to_file)
  output.close()

def runFoldxAnalyzeComplex(name, pdbs):
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()

  name = 'run_' + name + '.foldx'
  makeFoldxAnalyzeComplex(name)
  command = '/home/agm854/local/bin/foldx64Linux -runfile ' + name
  status = -1
  while status < 0:
    status = subprocess.call(command, shell=True)
  print('\n\nAnalyze Complex exit status = ' + str(status))
  
def makeFoldxRepair(name):
  to_file = '<TITLE>FOLDX_runscript;\n'+\
            '<JOBSTART>#;\n'+\
            '<PDBS>#;\n'+\
            '<BATCH>list.txt;\n'+\
            '<COMMANDS>FOLDX_commandfile;\n'+\
            '<RepairPDB>#;\n'+\
            '<END>#;\n'+\
            '<OPTIONS>FOLDX_optionfile;\n'+\
            '<Temperature>298;\n'+\
            '<R>#;\n'+\
            '<pH>7;\n'+\
            '<IonStrength>0.050;\n'+\
            '<water>-CRYSTAL;\n'+\
            '<metal>-CRYSTAL;\n'+\
            '<VdWDesign>2;\n'+\
            '<OutPDB>true;\n'+\
            '<pdb_hydrogens>false;\n'+\
            '<END>#;\n'+\
            '<JOBEND>#;\n'+\
            '<ENDFILE>#;\n'

  output = open(name, 'w')
  output.write(to_file)
  output.close()

def runFoldxRepair(name, pdbs):
  output = open('list.txt', 'w')
  fns = ''
  for pdb in pdbs:
    fns += pdb + '\n'
  output.write(fns)
  output.close()

  name = 'run_' + name + '.foldx'    
  makeFoldxRepair(name)
  subprocess.call('/home/agm854/local/bin/foldx64Linux -runfile ' + name, shell=True)
  command = '/home/agm854/local/bin/foldx64Linux -runfile ' + name
  status = -1
  while status < 0:
    status = subprocess.call(command, shell=True)
  print('\n\nRepair exit status = ' + str(status))

def checkOutputMutator(prefix):
  files = glob.glob('*_' + prefix + '.pdb')
  if len(files) == 2:
    return(True)
  else:
    print('\n\nWrong number of output mutator files.\n')
    return(False)

def checkOutputAnalyzeComplex(prefix):
  indiv_energy1 = 'Indiv_energies_AnalyseComplex_' + prefix + '.fxout'
  indiv_energy2 = 'Indiv_energies_AnalyseComplex_' + prefix + '.wt.fxout'
  inter_energy1 = 'Interaction_AnalyseComplex_' + prefix + '.fxout'
  inter_energy2 = 'Interaction_AnalyseComplex_' + prefix + '.wt.fxout'

  print(indiv_energy1)
  print(indiv_energy2)
  print(inter_energy1)
  print(inter_energy2)

  if os.path.isfile(indiv_energy1) and os.path.isfile(indiv_energy2) and os.path.isfile(inter_energy1) and os.path.isfile(inter_energy2):
    print('\n\nThese files are present: ' + indiv_energy1 + ' ' + indiv_energy2 + ' ' + inter_energy1 + ' ' + inter_energy2 + '\n')
    return(True)
  else:
    print('\n\nSome files are missing from Analyze Complex\n')
    return(False)

class Scores:
  def __init__(self):    
    self.files = []
    self.total_energies = []
    self.ids = []
    self.wts = []
    self.best_mutant = []
    self.group1 = []
    self.group2 = []
    self.clashes1 = []
    self.clashes2 = []
    self.stability = []
    self.interaction_energies = []
    
  def parseAnalyzeComplex(self):
    self.parseFiles('Interaction_AnalyseComplex_')
    if 'wt' not in self.files[0] and len(self.files) > 1:
      new_files = [self.files[1], self.files[0]]
      self.files = new_files
      
    for a_file in self.files:
      lines = open(a_file, 'r').readlines()
      start = False
      for a_line in lines:
        energy_list = a_line.split('\t')
        if not start and energy_list[0] == 'Pdb':
          start = True
        elif start and len(energy_list) > 1:
          (self.ids).append(self.removeWhiteSpace(energy_list[0]))
          (self.group1).append(self.removeWhiteSpace(energy_list[1]))
          (self.group2).append(self.removeWhiteSpace(energy_list[2]))
          (self.clashes1).append(self.removeWhiteSpace(energy_list[3]))
          (self.clashes2).append(self.removeWhiteSpace(energy_list[4]))
          (self.interaction_energies).append(self.removeWhiteSpace(energy_list[5]))

    self.parseFiles('Indiv_energies_')
    if 'wt' not in self.files[0] and len(self.files) > 1:
      new_files = [self.files[1], self.files[0]]
      self.files = new_files
          
    count = 0
    for a_file in self.files:
      lines = open(a_file, 'r').readlines()
      start = False
      stability = []
      for a_line in lines:
        energy_list = a_line.split('\t')
        if not start and energy_list[0] == 'Pdb':
          start = True
        elif start and len(energy_list) > 1 and self.ids[count] == energy_list[0]:
          stability.append(self.removeWhiteSpace(energy_list[2]))
      (self.stability).append(stability)
      count += 1
           
  def parseStability(self):
    self.parseFiles('Stability.txt')
    for a_file in self.files:
      lines = open(a_file, 'r').readlines()
      start = False
      for a_line in lines:
        energy_list = a_line.split('\t')
        if not start and energy_list[0] == '':
          start = True
        elif start and len(energy_list) > 1:
          (self.ids).append(energy_list[0])
          (self.total_energies).append(energy_list[1])
      
  def parseFiles(self, prefix):
    import os
    self.files = []
    raw_files = raw_files = os.listdir('.')
    
    for i in raw_files:
      if i.startswith(prefix):
        (self.files).append(i)
   
  def removeWhiteSpace(self, string):
    return("".join(string.split()))
    
  def findBestMutant(self):
    best_score = np.array(self.energies).argsort()[:1][::-1]
    self.best_mutant = self.ids[best_score]

  def getTotalEnergies(self):
    return(self.total_energies)
  def getWTenergies(self):
    return(self.wts)
  def getIds(self):
    return(self.ids)
  def getBestMutant(self):
    return(self.best_mutant)
  def getGroup1(self):
    return(self.group1)
  def getGroup2(self):
    return(self.group2)
  def getClashes1(self):
    return(self.clashes1)
  def getClashes2(self):
    return(self.clashes2)
  def getInteractionEnergies(self):
    return(self.interaction_energies)
  def getStability1(self):
    return(self.stability[0])
  def getStability2(self):
    return(self.stability[1])

  def cleanUp(self, delete_files):
    import os, glob
    for name in delete_files:
      for files in glob.glob(name):
        if os.path.isfile(files):
          os.remove(files)
    for files in glob.glob('run_*.foldx'):
      if os.path.isfile(files):
        os.remove(files)
    for files in glob.glob('*fxout'):
      if os.path.isfile(files):
        os.remove(files)
    if os.path.isfile('missing.txt'): 
      os.remove('missing.txt')
    if os.path.isfile('Stability.txt'):
      os.remove('Stability.txt')
    if os.path.isfile('list.txt'):
      os.remove('list.txt')
    if os.path.isfile('individual_list.txt'):
      os.remove('individual_list.txt')
    if os.path.isfile('scanning_output.txt'):
      os.remove('scanning_output.txt')
