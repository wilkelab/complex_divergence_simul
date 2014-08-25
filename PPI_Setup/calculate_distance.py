import sys, os
from Bio.PDB import *
import numpy as np
from interface_analysis import *

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

        p = PDBParser()
        structure = p.get_structure(pdb[0:-4], pdb)
    
        distances = calculate_ca_distance(structure, chain)
        mutation_sites = find_mutation_sites(file_name, chain)
        interface_sites = distances[distances[:, 1] <= distance_cutoff][:, 0]
        non_interface_sites = distances[distances[:, 1] > distance_cutoff][:, 0]

        interface_fraction, non_interface_fraction = calculate_fraction_interface(interface_sites, mutation_sites, len(non_interface_sites))

        print(str(interface_fraction) + '\t' + str(non_interface_fraction))

def find_mutation_sites(file_name, chain):
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
            if chain == 'A' and int(split[1][1:-1]) < 160:
                mutation_sites.append(int(split[1][1:-1]))
            elif chain == 'C' and int(split[1][1:-1]) > 160:
                mutation_sites.appen(int(split[1][1:-1]))

    return(mutation_sites)

if __name__ == "__main__":
    main()


