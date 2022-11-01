'''
The purpose of this file is to map MeSH Terms (metadata) to PMIDs
'''
import sys
# setting path
sys.path.append('..')

from text_mining.caseolap._03_mesh2pmid import *
import os


'''
Parameters
'''
### Input
root_dir = '/caseolap_lift_shared_folder/'
data_dir = os.path.join(root_dir,'data')
parsed_data_inputfile = os.path.join(data_dir,"pubmed.json")

### Output 
mesh2pmid_outputfile = os.path.join(data_dir,"mesh2pmid.json" )   # {"MeSH Term":[PMID1,...], ...}
mesh2pmid_statfile = os.path.join(data_dir,"mesh2pmid_stat.json") # {"MeSH Term": #PMIDs, ...}
logFilePath = os.path.join(root_dir,'log/mesh2pmid_log.txt')           # Logs mapping progress



'''
Main Code
'''
if __name__ == '__main__':
    
    # Open log file
    logfile = open(logFilePath, "w")     
    
    # Initialize class
    mesh2PMID = MeSH2PMID() 
    
    # Map MeSH to PMID 
    mesh2PMID.mesh2pmid_mapping(parsed_data_inputfile,mesh2pmid_outputfile,logfile)
    
    # Map MeSH to #PMIDs
    mesh2PMID.mesh2pmid_mapping_stat(mesh2pmid_statfile,logfile)         
    
    # Close log file
    logfile.close()                                                                  
