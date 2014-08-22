import numpy as np

def calculate_fraction_interface(interface_sites, mutation_sites, non_interface_sites):
    '''This function has two inputs.
       1.) A list of the interface sites in Integer format
       2.) A list of mutation sites in Integer Format
       3.) Sequence length of the chain of interest

       This function has two outputs.
       1.) A float that is the fraction of interface mutations per interface site
       2.) A float that is the fraction of non-interface mutations per non-interface site
    '''
    
    interface_count = 0
    non_interface_count  = 0

    for mutation_site in mutation_sites:
        if mutation_site in interface_sites:
            interface_count += 1
        else:
            non_interface_count += 1

    return([float(interface_count)/len(interface_sites), float(non_interface_count)/non_interface_sites])
    
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
        if len(distances) == 0:
            start_num = int(full_id[3][1])
            
        if test_chain == full_id[2] and 'CA' == full_id[4][0]:
            test_position = atom1.get_coord()
            
            tmp_distances = []
            for atom2 in structure.get_atoms():
                ref_id = atom2.get_full_id()
                if test_chain != ref_id[2]:
                    tmp_3d_distance = test_position - atom2.get_coord()
                    tmp_distances.append(np.sqrt(tmp_3d_distance[0]**2 + tmp_3d_distance[1]**2 + tmp_3d_distance[2]**2))
            distances.append([int(full_id[3][1]) - start_num, np.array(tmp_distances).min()])
            
    return(np.array(distances))
    
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
