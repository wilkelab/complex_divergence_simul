import sys, os
from Bio.PDB import *
import numpy as np

def main(): 
    if len( sys.argv ) != 4:
        print '''

        You don't have the right number of arguments.

        '''
    elif len( sys.argv ) ==  4:
        pdb = str(sys.argv[1])
        chain = str(sys.argv[2])
        distance_cutoff = float(sys.argv[3])

    p = PDBParser()
    structure = p.get_structure(pdb[0:-4], pdb)
    
    distances = calculate_ca_distance(structure, chain, distance_cutoff)
    print(distances[distances[:, 1] < distance_cutoff][:, 0])
    
def calculate_ca_distance(structure, test_chain, distance_cutoff):
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
            distances.append([full_id[3][1], np.array(tmp_distances).min()])
            
    return(np.array(distances))
        
if __name__ == "__main__":
    main()