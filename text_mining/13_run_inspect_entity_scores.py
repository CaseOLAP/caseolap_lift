'''
The purpose of this file is so that the user can see why each entity scored
highly; the user can see the ranked entities.
Get the ranked entities.
Get the ranked synonyms.
'''
from caseolap._13_inspect_entity_scores import *



'''
Parameters
'''
# Input path 
id2syns_in = 'input/id2syns.json'                 # The case-varied entity dict 
caseolap_scores_in = 'result/caseolap.csv'        # The CaseOLAP scores
popular_scores_in = 'result/popularity_score.csv'
distinct_scores_in = 'result/distinctiveness_score.csv'

pmid_syn_count_in = 'data/pmid_synonym_counts_2012-2022.json'  # Counts of each synonym
remove_syns_in = 'data/remove_these_synonyms.txt'    # Syns that were not used
cat2pmids_in = 'data/metadata_category2pmids_2012-2022.json'   # Category->[PMID,...,PMID] 

# Output paths
ranked_syns_out = 'result/ranked_synonyms/ranked_synonyms.txt' # Syns ranked by counts
ranked_ent_caseolap_out = 'result/Ranked Entities/ranked_caseolap_score/ranked_entities.txt'  # Ents ranked by score
ranked_ent_popular_out = 'result/Ranked Entities/ranked_popularity_score/ranked_entities.txt'  # Ents ranked by score
ranked_ent_distinct_out = 'result/Ranked Entities/ranked_distinctiveness_score/ranked_entities.txt'  # Ents ranked by score



'''
Main Code
'''
if __name__ == '__main__':
    # Initialize the class
    IES = InspectEntityScores(caseolap_scores_in, id2syns_in, remove_syns_in,\
                              pmid_syn_count_in, cat2pmids_in)

    '''ranked_synonyms'''
    # Export the ranked synonyms
    IES.get_ranked_synonyms_found(cat2pmids_in, pmid_syn_count_in, ranked_syns_out)
    
    '''Ranked Entities (CaseOLAP Score: Popular + Distinct)'''  
    # Sort entities by CaseOLAP scores
    IES.sort_all_scores(caseolap_scores_in)
    
    # Display the proportion entities found / entities searched for
    IES.prop_entities_found()
    
    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()
    
    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(ranked_ent_caseolap_out, score_component = 'CaseOLAP Score')
    
    
    '''Ranked Entities (Popular)'''
    # Initialize the class
    IES = InspectEntityScores(popular_scores_in, id2syns_in, remove_syns_in,\
                              pmid_syn_count_in, cat2pmids_in)

    IES.get_ranked_synonyms_found(cat2pmids_in, pmid_syn_count_in, ranked_syns_out,
                                  export = False)
    
    # Sort entities by CaseOLAP scores
    IES.sort_all_scores(popular_scores_in)
    
    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()
    
    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(ranked_ent_popular_out, score_component = 'Popularity Score')
    
    
    
    
    '''Ranked Entities (Distinct)'''
    # Initialize the class
    IES = InspectEntityScores(distinct_scores_in, id2syns_in, remove_syns_in,\
                              pmid_syn_count_in, cat2pmids_in)
    
    # Get the ranked synonyms
    IES.get_ranked_synonyms_found(cat2pmids_in, pmid_syn_count_in, ranked_syns_out,\
                                  export = False)
    
    # Sort entities by CaseOLAP scores
    IES.sort_all_scores(distinct_scores_in)
    
    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()
    
    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(ranked_ent_distinct_out, score_component = 'Distinctiveness Score')
    
    
    
    