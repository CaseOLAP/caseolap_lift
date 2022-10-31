'''
The purpose of this file is to create mappings for the PMIDs to their categories.
Categories (e.g., diseases) are identified via MeSH term metadata.
The TextCube is mappings between PMID to category
'''
import json, sys, time, os
# setting path
sys.path.append('..')

from text_mining.caseolap._06_textcube import *



'''
Parameters
'''
# Input data
root_dir = '/caseolap_lift_shared_folder/'
data_dir = os.path.join(root_dir,'data')
meshtree = os.path.join(data_dir,'MeSH/mtrees2022.bin')                 # MeSH Tree ontology
mesh2pmid = os.path.join(data_dir,'mesh2pmid.json') # MeSH to PMID mapping
root_cat = os.path.join(data_dir,'categories.txt')                 # Categories' root MeSH Tree nums
textcube_config = os.path.join(root_dir,'config/textcube_config.json')   # ['CategoryName1',...]

# Output data directories
textcube_pmid2category = os.path.join(data_dir,'textcube_pmid2category.json') # Map PMID to category
textcube_category2pmid = os.path.join(data_dir,'textcube_category2pmid.json') # Map category to PMID
textcube_stat = os.path.join(data_dir,'textcube_stat.txt')                    # Num. documents per category
MeSHterms_percat = os.path.join(data_dir,'meshterms_per_cat.json')            # MeSH *Terms* per category
logfile_path = os.path.join(root_dir,'log/textcube_log.txt')                       # Logs file's messages



'''
Main Code
'''
if __name__ == '__main__':
    logfile = open(logfile_path, 'w') 
    
    # Get category names
    category_names = json.load(open(textcube_config, 'r'))
    print(str(len(category_names)),'categories: ',category_names)

    # Initialize TextCube class
    TC =  TextCube(category_names)
    
    # Get categories' descendant MeSH terms
    TC.descendant_mesh(root_cat, meshtree, MeSHterms_percat, logfile)
    
    # Get categories' PMIDs
    TC.category2pmids_mapping(mesh2pmid, textcube_category2pmid, logfile) 
    
    # Get PMIDs' categories
    TC.pmid2category_mapping(textcube_pmid2category, logfile)                
    
    # Display the number of documents per category
    TC.category_statistics(textcube_stat, logfile)                              

    logfile.close()
