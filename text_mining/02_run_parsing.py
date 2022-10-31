'''
The purpose of this file is to parse the downloaded PubMed documents,
saving their information into a dictionary.
'''
import json, sys, os, time, traceback

# setting path
sys.path.append('..')
from text_mining.caseolap._02_parsing import *
from lxml import etree



'''
Parameters
'''
# Input
root_dir = '/caseolap_lift_shared_folder/'
data_dir = os.path.join(root_dir,'data')
baseline_files = os.path.join(data_dir,'ftp.ncbi.nlm.nih.gov/pubmed/baseline' ) # Documents prior to this year
update_files = os.path.join(data_dir,'ftp.ncbi.nlm.nih.gov/pubmed/updatefiles') # Documents from this year
parsing_config_file = os.path.join(root_dir,'config/parsing_config.json')

# Output 
pubmed_path = os.path.join(data_dir,'pubmed.json')     # The parsed PubMed documents (dictionary)
filestat_path = os.path.join(data_dir,'filestat.json') # Unzipped PubMed download file name : #PMIDs
logfile_path = os.path.join(root_dir,'log/parsing_log.txt') # Logs progress on parsing



'''
Main Code
'''
if __name__ == '__main__':    
    
    # Start time
    t1 = time.time()
    
    # Open the files
    logfile = open(logfile_path, 'w') 
    pubmed = open(pubmed_path, 'w')
    filestat = open(filestat_path, 'w')
    parsing_config = json.load(open(parsing_config_file,'r'))
    print(parsing_config)
    
    # Parse the files (baseline and updatefiles)     
    parse_dir(baseline_files, pubmed, filestat, 'baseline', parsing_config, logfile)
    parse_dir(update_files, pubmed, filestat, 'updatefiles', parsing_config, logfile)
    
    # Close the files
    pubmed.close()
    filestat.close()
    logfile.close()
    t2 = time.time()
    
    # Print info
    print('Parsing finished, results dumped to',pubmed_path)
    print('Total Time: ', str((t2 - t1)/360),'hours')           
    
    