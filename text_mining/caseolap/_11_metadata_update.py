import json


class MetadataUpdate():
    
    
    def __init__(self, category_names):
        self.category_names = category_names
        self.pmid2entity2count = dict()
        self.category2pmids = dict()         
        
        # Map Category Name -> Category Number 
        self.category_num2name = dict()
        for category_num, category_name in enumerate(self.category_names):
            self.category_num2name[category_num] = category_name
        print('11. Category Names to Number:', self.category_num2name)

    def update_pmid2entity2count(self, input_file_entitycount, 
                                 output_file_pmid2entity2count, logfile):
        '''
        FUNCTION:
        - This reformats the entitycount file into a nested dictionary so
          it can be used in the CaseOLAP score calculation. 
        
        PARAMS:
        - input_file_entity_count (TXT): Text file with lines for each PMID and 
          its entities and their counts (e.g., 1234567 entity1|4 entity2|3 ...)
        - output_file_pmid2entity2count (JSON): Collected and merged data from the 
          input file will be written here, storing PMID->Entity->Entity Count
        - logfile (TXT): Writes the progress 
        '''
        
        # Print progress
        msg = 'Updating the PMID -to-> Entity Count dictionary'
        logfile.write(msg+'\n'+'='*len(msg)+'\n')
        print(msg)
        
        
        # PMID->Entity->Entity Count Dictionary'''
        with open(input_file_entitycount) as fin:
            for line in fin:
                
                # Read in data
                line = line.strip().split(' ')
                pmid = line[0]
                entities2count_list = line[1:]
                entities2counts = dict()

                for entity2count in entities2count_list:
                    entity2count = entity2count.split('|')
                        
                    # Entity, Count
                    entity = entity2count[0]
                    count = entity2count[1]
                    
                    # Entity->Count
                    entities2counts[entity] = count
                    
                    # PMID->Entity->Count
                    self.pmid2entity2count[pmid] = entities2counts
                            
                # Write to logfile
                logfile.write('PMID: '+str(pmid)+\
                              " Entity Counts: "+str(entities2counts)+'\n')

        # Export
        json.dump(self.pmid2entity2count, open(output_file_pmid2entity2count, 'w') )
                
                


    def map_category2pmid_pmids_with_entities(self, input_file_textcube_pmid2category,\
                         output_file_metadata_category2pmids, logfile):
        '''
        FUNCTION:
        - This maps Category Names to PMIDs for all PMIDs in which an entity of interest
          was found.
         
        PARAMS:
        - input_file_textcube_pmid2category (JSON): This list of lists indicates 
          the PMIDs and the category (number) they belong to (e.g., [['123456',0],...])
          This is for all documents (publications) in the category regardless of if
          they have an entity.
        - output_file_metadata_category2pmids (JSON): This dictionary mapping category
          names to PMIDs has a subset of the above input information; the category-PMID
          mappings are only for PMIDs that had an entity discovered in them.
        - logfile (TXT): Output the progress          
        '''
        
        # Print progress
        msg = 'Updating the Category -to-> Entity-Containing-PMID dictionary'
        logfile.write(msg+'\n'+'='*len(msg)+'\n')
        print(msg)

            
        # PMID to Category Dictionary
        pmids2categories = json.load(open(input_file_textcube_pmid2category))
        for pmid2category in pmids2categories:
            
            # PMID, Category (PMIDs irrespective of entities found, from textcube)
            pmid = pmid2category[0]

            category_num = pmid2category[1]
            if pmid == 'doc_id':
                continue
                
            # PMIDs w/entities found
            if pmid in self.pmid2entity2count:
                
                # PMID->Category (PMIDs w/entities found)
                category_name = self.category_num2name[category_num]
                self.category2pmids.setdefault(category_name,list()).append(pmid)                                
                
                
        # Write results        
        for name in self.category_names:                       
            logfile.write("Category: " +name+" includes"+\
                          str(len(self.category2pmids[name]))+\
                          " documents containing entities.\n")
        
        # Export results
        json.dump(self.category2pmids, open(output_file_metadata_category2pmids, 'w'))       