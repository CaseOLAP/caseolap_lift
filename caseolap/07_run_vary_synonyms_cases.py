'''
The purpose of this file is to create case-sensitive variations
of the synonyms. The synonyms will then be queried.
'''
from caseolap._07_vary_synonyms_cases import *


'''
Paths
'''
# Input
entity_dict_path = 'input/id2syns_not_case_varied.json' # entity dict
species = ['human', 'pig', 'mouse', 'rat']  # Species studied. Used for permuting syns.

# Output
case_varied_entites_outpath = 'data/casesensitive_entities.txt' # Case sensitive entity dict
case_varied_entity_dict_path = 'input/id2syns.json' # entity dict



'''
Main code
'''
if __name__ == '__main__':
    # Instantiates the class, loads entity dictionary mapping ID to synonyms
    VSC = VarySynonymsCases(entity_dict_path, case_varied_entites_outpath)
    
    # Adds some more synonyms
    VSC.add_species_syns(species)
    
    # Makes case-sensitive variations of the synonyms
    VSC.gets_final_syns()
    
    # Make a dictionary of the ID to synonym data
    make_id2syns_dict()