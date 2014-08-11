import sys, os
from Bio.PDB import *
import numpy as np

def main(): 
    if len( sys.argv ) != 4:
        print '''

        You don't have the right number of arguments.

        Input line should resemble:
        python calculate_distance.py chain distance_cutoff fixed_mutation_file

        '''
    elif len( sys.argv ) ==  4:
        chain           = str(sys.argv[1])
        distance_cutoff = float(sys.argv[2])
        file_name       = str(sys.argv[3])

        pdb = get_final_pdb(file_name)

        p = PDBParser(QUIET=True)
        structure = p.get_structure(pdb[0:-4], pdb)
    
        distances = calculate_ca_distance(structure, chain)
        mutation_sites = find_mutation_sites(file_name)
        interface_sites = distances[distances[:, 1] < distance_cutoff][:, 0]
        
        print(interface_sites)

        print(calculate_fraction_interface(interface_sites, mutation_sites))
    
def get_final_pdb(file_name):
    '''This function has one input.
       1.) A file name in String format for the fixed mutation data from main.py
           Traditionally this file is named final_data.txt, but 
           could be anything the simulator used

       This function will return one output.
       1.) The name of the final pdb structure
    '''

    in_file = open(file_name, 'r')

    return(in_file.readlines()[-1].split('\t')[0])

def calculate_ca_distance(structure, test_chain):
    '''This function has two inputs.  
       1.) A structure object from the biopython parser.
       2.) The chain to measure distances FROM in String format.
     
       This function will return one output.
       1.) A numpy array containing two element lists of the form [Integer, Float] with the site number and 
           the distance to the closest atom in its bound partner.
    '''

    distances = []
    for atom1 in structure.get_atoms():
        full_id = atom1.get_full_id()
        if test_chain == full_id[2] and 'CA' == full_id[4][0]:
            test_position = atom1.get_coord()
            
            tmp_distances = []
            for atom2 in structure.get_atoms():
                ref_id = atom2.get_full_id()
                if test_chain != ref_id[2]:
                    tmp_3d_distance = test_position - atom2.get_coord()
                    tmp_distances.append(np.sqrt(tmp_3d_distance[0]**2 + tmp_3d_distance[1]**2 + tmp_3d_distance[2]**2))
            distances.append([int(full_id[3][1]), np.array(tmp_distances).min()])
            
    return(np.array(distances))

def find_mutation_sites(file_name):
    '''This function has one input.
       1.) A file name in String format for the fixed mutation data from main.py
           Traditionally this file is named final_data.txt, but 
           could be anything the simulator used

       This function will return one output.
       1.) A list of integers of the sites where mutations fixed
    '''
  
    in_file = open(file_name, 'r')
    
    mutation_sites = []

    for a_line in in_file.readlines():
        split = a_line.split('\t')
        if split[0] == 'file':
            continue
        else:
            mutation_sites.append(int(split[1][1:-1]))

    return(mutation_sites)

def calculate_fraction_interface(interface_sites, mutation_sites):
    '''This function has two inputs.
       1.) A list of the interface sites in Integer format
       2.) A list of mutation sites in Integer Format

       This function has one output.
       1.) A float that is the fraction of mutations in the interface
    '''
    
    interface_count = 0
    mutation_count  = 0

    for mutation_site in mutation_sites:
        if mutation_site in interface_sites:
            interface_count += 1
        mutation_count += 1

    return(float(interface_count)/mutation_count)
        

if __name__ == "__main__":
    main()


