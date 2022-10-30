'''
The purpose of this file is to map MeSH Terms (metadata) to PMIDs
'''
from caseolap._03_mesh2pmid import *



'''
Parameters
'''
### Input
parsed_data_inputfile = "./data/pubmed.json"

### Output 
mesh2pmid_outputfile = "./data/mesh2pmid.json"    # {"MeSH Term":[PMID1,...], ...}
mesh2pmid_statfile = "./data/mesh2pmid_stat.json" # {"MeSH Term": #PMIDs, ...}
logFilePath = "./log/mesh2pmid_log.txt"           # Logs mapping progress



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