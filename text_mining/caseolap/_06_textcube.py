import json, sys, time

class TextCube(object):
    def __init__(self,category_names):
        self.category_names = category_names
        self.category2pmids = [set() for _ in range(len(self.category_names))]
        self.concerned_cat = []
        self.pmid2category = []

        
    def descendant_mesh(self, infile_root_cat, infile_meshtree, 
                        outfile_meshterms_percat, logfile):
        '''
        FUNCTION:
        - Starting from the category's root node tree number, 
          find its descendant tree numbers and save their 
          corresponding MeSH terms (names)
          
        PARAMS:
        - infile_root_cat: Input file path where each line is a 
          category's root MeSH tree numbers
        - infile_meshtree: Input file path of the MeSH tree where
          each line is a 'MeSH tree number; MeSH term'
        - outfile_meshterms_percat: Output file path of where the
          list of each category's MeSH terms will be saved
        - logfile: Output file path, writes progress
        '''
        
        # Print update about which step is running
        msg = 'Collecting categories\' subcategory MeSH terms...'
        logfile.write(msg+'\n'+'='*50+'\n')
        print(msg)
        
        
        # Get main category name (list of main MeSH Tree Numbers)
        for line in open(infile_root_cat):
            self.concerned_cat.append(line.strip().split())
           
        # Number of categories
        self.num_cat = len(self.concerned_cat)
        
        
        # MeSH terms in each category
        self.mesh_terms_per_cat = [set() for _ in range(self.num_cat)]
        
        
        # Read in MeSH tree (term; tree number)
        with open(infile_meshtree) as fin:
            for line in fin:
                line = line.strip().split(';')

                # MeSH Term, Tree Number
                mesh_term = line[0]
                tree_number = line[1]

                # For each category
                for cat_num, root_tree_numbers in enumerate(self.concerned_cat):

                    # For each root tree number in the category
                    for root_tree_number in root_tree_numbers:

                        # If this tree number is a descendant of the root
                        if tree_number.startswith(root_tree_number):

                            # Then save its MeSH term
                            self.mesh_terms_per_cat[cat_num].add(mesh_term)

                            
        # Output number of MeSH terms per category
        for cat_num in range(self.num_cat):                    
            logfile.write(self.category_names[cat_num] +\
                          ': includes descendants'+\
                          str(self.mesh_terms_per_cat[cat_num])+'\n')

        # Output MeSH terms per category
        mesh_lists = [list(mesh_terms) for mesh_terms in self.mesh_terms_per_cat]
        json.dump(mesh_lists, open(outfile_meshterms_percat, 'w'))

            
            
            
    def category2pmids_mapping(self, input_file_mesh2pmid, output_file_textcube_category2pmid, logfile):
        '''
        FUNCTION:
        - Map categories to PMIDs (e.g., Disease Category : [PMIDs studying the disease category])

        PARAMS:
        - input_file_mesh2pmid: Input file path, stores the MeSH-PMIDs mappings
        - output_file_textcube_category2pmid: Output file path, stores Category-PMIDs
        - logfile: Output file path, stores progress
        '''

        # Print progress
        msg = 'Textcube category to PMID mapping is being created....'
        logfile.write(msg+'\n'+'='*len(msg)+'\n')       
        print(msg)

        start = time.time()
        for input_index, line in enumerate(open(input_file_mesh2pmid)):

            # Print progress
            if input_index%100 == 0:
                msg = str(input_index) + ' MeSH names analyzed for textcube... ' + \
                      str(round(time.time()-start,4)) + ' seconds'
                logfile.write(msg+'\n')
                print(msg, end='\r')

            # Get one MeSH-PMIDs
            mesh2pmid = json.loads(line.strip())
            mesh = list(mesh2pmid.keys())[0]
            pmids = list(mesh2pmid.values())[0]

            # For each category number
            for cat_num in range(self.num_cat):

                # If MeSH -in-> category
                if mesh in self.mesh_terms_per_cat[cat_num]:

                    # Add: Category -has-> PMIDs
                    self.category2pmids[cat_num] = \
                    self.category2pmids[cat_num].union(set(pmids)) 

        # Switch sets to lists, export                
        for category in range(len(self.category2pmids)):
            self.category2pmids[category] = list(self.category2pmids[category])
        json.dump(self.category2pmids, open(output_file_textcube_category2pmid, 'w'))

        # Print result summary
        for cat_num, name in enumerate(self.category_names):    
            logfile.write(name + " includes " + str(len(self.category2pmids[cat_num])) + " documents. \n")
   


    def pmid2category_mapping(self, output_file_textcube_pmid2category, logfile):
        '''
        FUNCTION:
        - Map PMID to category, list of lists [PMID, category number]
        
        PARAMS:
        - output_file_textcube_category2pmid: Output file path, stores Category-PMIDs
        - logfile: Output file path, stores progress
        '''
        msg = 'Textcube PMID to category mapping is being created....'
        logfile.write(msg+'\n============================================== \n')
        print(msg)
        
        for cat_num in range(self.num_cat):
            for cur_pmid in self.category2pmids[cat_num]:
                self.pmid2category.append([cur_pmid, cat_num])
        
        json.dump(self.pmid2category, open(output_file_textcube_pmid2category, 'w'))
            

            
    def category_statistics(self, outputfile_textcube_stat, logfile): 
        '''
        FUNCTION:
        - Saves the counts of documents (PMIDs) per category.
          This includes duplicate counts of PMIDs in multiple documents

        PARAMS:
        - outputfile_textcube_stat: Output file path, stores textcube summary on categories
        - logfile: Output file path, stores progress
        '''

        # Print update about which step is running
        msg = 'Textcube category statistics is being created....'
        logfile.write(msg+'\n'+'='*50+'\n')
        print(msg)    

        # Open output file to write stats to
        with open(outputfile_textcube_stat, 'w') as fout:
            all_pmids = [] 
            category_count = [0]*self.num_cat

            # Get PMID & Category counts
            for pmid2catnum in self.pmid2category:

                # PMID
                pmid = pmid2catnum[0]
                all_pmids.append(pmid)

                # Category
                cat_num = pmid2catnum[1]
                category_count[cat_num] += 1


            # Write PMID & Category counts
            print('\nCategory: Documents (From All Years)\n'+'='*19)
            fout.write('\nCategory: Documents\n'+'='*19+'\n')
            for cat_num, category_name in enumerate(self.category_names):
                num_documents = str("{:,}".format(category_count[cat_num]))
                msg = str(category_name)+': '+num_documents
                fout.write(msg+'\n')
                print(msg)
                
            num_documents = str("{:,}".format(len(set(all_pmids))))
            msg = 'TOTAL: '+ num_documents
            fout.write(msg+'\n')
            print(msg)

