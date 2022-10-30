import time, re, sys, os, json
from collections import defaultdict

# MeSH to PMID object (stored two dictionaries)
class MeSH2PMID(object):
    
    def __init__(self):
        self.mesh2pmid       = dict()  # MeSH to list of PMIDs 
        self.mesh2pmid_count = dict() # MeSH to number of PMIDs
        
    def mesh2pmid_mapping(self, parsed_data_inputfile, mesh2pmid_outputfile, logfile):
        '''
        FUNCTION:
        - Map MeSH to PMIDs and output them as files
        
        PARAMS:
        - parsed_data_inputfile (JSON file): Parsed text files 
        - mesh2pmid_outputfile (JSON file): Stored mesh2pmid mapping
        - logfile (Text file): Stores progress
        '''
        
        # Read in parsed pubmed documents
        with open(parsed_data_inputfile) as fin:
            print('MeSH to PMID mapping is running...')
            start = time.time()
            
            # For each parsed document, map PMID to MeSH
            for num_docs, document in enumerate(fin):
                    
                # Get PMID and MeSH
                pub_info = json.loads(document.strip())
                pmid_mesh = dict()
                pmid_mesh['pmid'] = pub_info.get('PMID', '-1')                   
                pmid_mesh['mesh_heading'] = ' '.join(pub_info['MeshHeadingList'])

                # Map PMID to MeSH
                if pmid_mesh['pmid'] != '-1':
                    for mesh in pub_info['MeshHeadingList']:
                        self.mesh2pmid.setdefault(mesh,[]).append(pmid_mesh['pmid'])

                # Print progress
                if num_docs % 500000 == 0:
                    msg = str(num_docs) +' PMIDs attempted to map to MeSH terms'
                    logfile.write(msg+'\n')
                    print(msg, end = '\r')

            # Export 'MeSH to PMID' mapping table (this loads one line at a time,
            # the textcube file reads it in one line at a time)
            print('Exporting MeSH to PMID mapping')
            mesh2pmid_fout = open(mesh2pmid_outputfile, 'w')
            for mesh, pmids in self.mesh2pmid.items():
                json.dump({mesh:pmids}, mesh2pmid_fout)
                mesh2pmid_fout.write('\n')

            # Print final progress
            print('Finished MeSH to PMID mapping. Time = '+\
                  str(round(time.time()-start)/60)+' minutes')
            
        
    def mesh2pmid_mapping_stat(self, mesh2pmid_countfile, logfile): 
        '''
        FUNCTION:
        - Get count on the numer of MeSH to PMID mappings 
        
        PARAMS:
        - mesh2pmid_countfile (JSON file): Stores counts of MeSH to PMIDs
        - logfile (Text file): Stores progress
        '''
        
        # Printing progress
        msg = 'MeSH to PMID mapping stat is running.'
        logfile.write(msg+'\n')
        print(msg)
        
        # Make mapping
        for mesh, pmids in self.mesh2pmid.items():
            self.mesh2pmid_count[mesh] = len(pmids)
        
        # Export mapping
        json.dump(self.mesh2pmid_count, open(mesh2pmid_countfile,'w'))