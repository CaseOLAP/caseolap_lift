from caseolap._10_make_entity_counts import *
'''
The purpose of this class is to make the entitycount file that maps 
PMIDs to entities to entity counts. This uses a "PMIDs to synonyms to synonym counts"
mapping and an "entity to synonym" mapping to do this. No ElasticSearch querying
is needed in this step.
'''


'''
Parameters
'''
# Input
remove_syns_infile = 'data/remove_these_synonyms.txt'
all_id2syns_path = 'input/id2syns.json'
core_id2syns_path = 'input/core_id2syns.json'
pmid_syn_count_in = 'data/pmid_synonym_counts_2012-2022.json'

# Output
all_entitycount_outfile = 'data/all_entitycount_2012-2022.txt'
core_entitycount_outfile = 'data/core_entitycount_2012-2022.txt'
'''
Main code
'''
if __name__ == '__main__':

    #### Expanded set proteins
    # Instantiate and initialize class
    MEC = MakeEntityCounts(all_id2syns_path, remove_syns_infile, pmid_syn_count_in)
    
    # Makes "pmid->entity->count" from "pmid->syn->count" & "entity->syn" mappings
    MEC.entitycount(all_entitycount_outfile)
    
    
    #### Core Proteins
    # Instantiate and initialize class
    MEC = MakeEntityCounts(core_id2syns_path, remove_syns_infile, pmid_syn_count_in)
    
    # Makes "pmid->entity->count" from "pmid->syn->count" & "entity->syn" mappings
    MEC.entitycount(core_entitycount_outfile)