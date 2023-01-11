import pandas as pd, json, os
from multiprocessing import cpu_count,Process



class VarySynonymsCases(object):
    '''
    Case-vary the good synonym names but not acronyms. 
    (e.g., Peroxidase --> peroxidase, not PreP). 
    I.e., capitalize words when their cases are different due 
    to sentence location (e.g., beginning or not) 
    '''
    
    def __init__(self, entity_dict_path, case_varied_entites_outpath):
        '''
        PARAMS:
        - entity_dict_path: Input file path, where the entity dictionary 
          containing the original synoynsm is
        - case_varied_entites_outpath: Output file path, where the entity 
          dictionary containing case-varied synonyms will go
        '''
        self.id2syns = {}
        self.backup = self.id2syns
        self.backup_id_syn = {}
        self.id2syns = json.load(open(entity_dict_path))  # ID:[syns]
        self.case_varied_entites_path = case_varied_entites_outpath

        
    def add_species_syns(self, species):
        '''
        FUNCTION: 
        Add synonyms (e.g., if one synonym is "P12345 human", add "P12345")
        
        PARAMS: 
        - species = the possible species you are studying.
        '''
        # For each ID
        for ID, synonym_phrases in self.id2syns.items():
            
            # For each synonym
            for synonym_phrase in synonym_phrases:
                synonym_phrase = synonym_phrase.split(' ')
                
                # Check if the synonym ends with a species name 
                if len(synonym_phrase) > 1:
                    if synonym_phrase[-1].lower() in species:  
                        for word in synonym_phrase[:-1]:
                            
                            # Synonyms with 3+ words may be a phrase
                            # e.g., "Protein upregulated in mouse"
                            # instead of "P12345 mouse"
                            if len(synonym_phrase) > 2:
                                print('Maybe remove these synonyms',\
                                      word, 'from', synonym_phrase)
                            else:
                                # Add the ID (e.g., P12345)
                                self.id2syns[ID].append(word)  
    
    
    def vary_case_inds(self, synonym):
        '''
        FUNCTION: 
        Return the locations of the words you want to case-vary 
        e.g., for ['Fake','PrO-10','name'], return [0,2] because 
        PrO-10 shouldn't be case-varied, but 'Fake' and 'name' should 
        be (stored in locations 0 and 2).

        PARAMS: 
        - synonym = synonym you want to case-vary
        '''
        indices_of_words_to_be_case_varied = list()
        
        
        # Check all words in the synonym
        synonym_as_list = synonym.replace('-', ' - ').replace('/','/ ').split(' ')
        
        for word_ind, word in enumerate(synonym_as_list):
            if word == '': 
                continue
                
            # Uppercase word to be case-varied, e.g., Peroxidase, not PrEP
            good_upper_case = (word[0].isupper() and word[1:].replace('/','').islower()) 
            
            # Lowercase word to be case-varied, e.g., perodixase
            good_lower_case = word.replace('/','').islower()  
            
            # If it should be case-varied, store the index
            if (good_upper_case or good_lower_case) and word.replace('/','').isalpha():
                indices_of_words_to_be_case_varied.append(word_ind)
                
        return indices_of_words_to_be_case_varied



    def get_all_case_varied_word_list(self, root1, syn_as_list, 
                                      syn_as_list_opp_case, 
                                      all_final_synoynms, 
                                      current_word_index):
        '''
        FUNCTION: 
        - Intermediate step in case-varying a synonym, often comprised 
         of multiple words. 
         Recursive function that takes two lists and makes all possible 
         permutations them:
         e.g., ['a','b'] + ['A','B'] --> ['a b', 'a B', 'A B', 'A b']
         Each list in a synonym split into its constituent words
        
        PARAMS: 
        - root1 = first word in a (often multi-word) synonym. Like the root 
          in a tree. The first time this recursive function is called, root
          should be l1[0] or l2[0]
        - syn_as_list = list of the original non-case-varied synonym
          e.g., ["fake", "syno9m", "phrase"]
        - syn_as_list_opp_case = list of the synonym with every 
          permitted word case-varied
          e.g., ["Fake", "syno9m", "Phrase"]
          l1 and l2 will be combined in every possible way to populate
          mainset
        - all_final_synoynms = Collects all the final case-varied synonyms
        - current_word_index = the current index in the 
          list holding the synonym phrase (e.g., i = 3 when we will be
          looking at  )
        '''
        # Will be final permuted synonyms, left and right nodes
        left_node = root1.replace(' - ','-') 
        right_node = root1.replace(' - ','-') 
        
        # Look at the next word in the synonym phrase
        current_word_index += 1
        
        # When the recursive process ends, save the new permuted synonym
        if current_word_index >= len(syn_as_list):
            all_final_synoynms.add(left_node)
            all_final_synoynms.add(right_node)
            return
        
        # Recursively add both case-variations of the ith word to 
        # their own left or right branch on the current "tree"
        if(len(left_node)==0):
            print("root1: ",root1)
            print("left_node: ",left_node)
            print("right_node: ",right_node)
            print("syn_as_list: ",str(syn_as_list))
            
        if left_node[-1] == '/':
            left_node += syn_as_list[current_word_index]
        else:
            left_node += ' '+syn_as_list[current_word_index]
        self.get_all_case_varied_word_list(left_node, 
                                           syn_as_list, 
                                           syn_as_list_opp_case, 
                                           all_final_synoynms, 
                                           current_word_index)

        if right_node[-1] == '/':
            right_node += syn_as_list_opp_case[current_word_index]
        else:
            right_node += ' '+syn_as_list_opp_case[current_word_index]
        self.get_all_case_varied_word_list(right_node, 
                                           syn_as_list, 
                                           syn_as_list_opp_case,
                                           all_final_synoynms, 
                                           current_word_index)
            
        
    def vary_case_of_a_phrase(self, syn, relevant_indices, combos):
        '''
        FUNCITON: 
        Takes in a synonym, splits it into its constituent words, 
        and case-varies the appropriate words (i.e. not acronyms).
        "Combos" is the important final variable here that gets modified.
        
        PARAMS: 
        - syn = synonym
        - relevant_indices = indices where the constituent words that \
          should be case-varied are
        - combos = combinations of case-varied words for a synonym 
          (e.g., fake protein, Fake Protein, fake Protein, Fake Protein)
        '''        
        # Split synonym into a list
        syn_as_list = syn.replace('-', ' - ').replace('/','/ ').split(' ')
        if syn_as_list[0] == '':
            syn_as_list = syn_as_list[1:]
        # Split synonym into a list, words are opposite cased 
        syn_as_list_opp_case = list()
        
        for index, word in enumerate(syn_as_list):     
            
            # If this word should be case-varied, do that
            if index in relevant_indices:
                
                # Case-vary
                if word[0].isupper():
                    case_varied_word = word[0].lower() + word[1:]

                else: 
                    case_varied_word = word[0].upper() + word[1:]
                
                # Add case-varied word
                syn_as_list_opp_case.append(case_varied_word)
            
            else:
                # Add original word
                syn_as_list_opp_case.append(word)
        #if syn_as_list[0] == '':
        #    print("This is the original synonym",syn)
        self.get_all_case_varied_word_list(syn_as_list[0], 
                                           syn_as_list, 
                                           syn_as_list_opp_case, 
                                           combos, 0)
        self.get_all_case_varied_word_list(syn_as_list_opp_case[0], 
                                           syn_as_list, 
                                           syn_as_list_opp_case, 
                                           combos, 0)
        
        

    def get_case_varied_syns(self, b_id, ids_and_syns):
        '''
        FUNCTION: Makes case-varied synonyms for a batch of entities. 
        Outputs it to a file. 
        
        PARAMS: 
        - id_and_syns = [{ID: [synonyms,...']},...]
        '''


        # Get and clear the new file location
        temp_outpath = self.case_varied_entites_path[:-4]+'_'+str(b_id)+'.txt'
        open(temp_outpath, 'w')
        
        tot = len(ids_and_syns)
        
        # Case-varies each synonym phrase (e.g., "fake protein")
        for index, id_and_syns in enumerate(ids_and_syns):
            case_varied_syns = set()
            
            # Print progress
            if b_id == 1:
                print('Batch 1 progress:', index+1,'/',tot,end='\r')
            
            # ID (str), synonyms (list)
            ID = list(id_and_syns.keys())[0]
            syns = list(id_and_syns.values())[0]
            
            # Each synonym
            for syn in syns:
                if type(syn) != None:
                    self.vary_case_of_a_phrase(syn, self.vary_case_inds(syn),
                                               case_varied_syns)

            # Writes the entity's new case-varied synonyms to a temp file
            with open(temp_outpath, 'a') as fout:
                case_varied_syns_str = '|'.join(case_varied_syns)
                fout.write(ID+'|'+case_varied_syns_str+'\n')


                    
    def gets_final_syns(self):
        '''
        FUNCTION: 
        - Main case-varying function that calls the other case-varying 
          functions. Makes all case-varied versions of all entities' synonyms
          e.g., takes "Fake prot" and adds "Fake Prot", "fake prot", and 
          "fake Prot" It case-varies the individual words in the synonym 
          unless they are obviously acronyms (e.g., ProT, Pro10, pro-10)
        
        PARAMS: 
        - case_varied_entites_path = where the entities with the 
          finalized synonyms will be stored
        '''
        
        # Use multiprocessing to case-vary the synonyms
        multiprocess_a_dict_of_keys_and_lists_values(
            self.id2syns, 
            self.get_case_varied_syns)
        
        ### Now, the case-varied synonyms are in temp files.
        ### Next, we combine the temp files into one file. 
        
        # Clear the new file location
        open(self.case_varied_entites_path, 'w')
        
        # For each batch
        procs = cpu_count()
        for b_id in range(procs):
            
            # Open the temp file
            temp_path = self.case_varied_entites_path[:-4]\
                        +'_'+str(b_id)+'.txt'
            lines = open(temp_path).readlines()

            # Output and merge each temp file
            with open(self.case_varied_entites_path,'a') as fout:
                for line in lines:
                    fout.write(line)
            
            # Remove the temp file
            os.remove(temp_path)
            
            
            
            
def multiprocess_a_dict_of_keys_and_lists_values(thedict, the_function):
    '''
    FUNCTION:
    - This takes a dictionary whose keys are strings and 
      values are lists of strings and splits it into input
      for separate processes. The processes then output 
      their results to temp files which are then merged. 

    PARARMS:
    - thedict: the dictionary to be split into input for 
      a multiprocessing function
    - the_function: the function that will use the dictionary
      as input for multiprocessing
    '''

    # How many processors can be used
    procs = cpu_count()

    # List of batches for multiprocessing
    batches = [[] for i in range(procs)]   

    # Length of input dictionary 
    tot = len(thedict)

    # Create batches and send to multiprocessing
    for i,(key, values_list) in enumerate(thedict.items()):

        # Add synonym to a batch
        the_values_list = [val for val in values_list if len(val) > 0]
        b_id = i%procs
        batches[b_id].append({key:the_values_list})

    # Create a list of jobs
    print("Running jobs...")
    jobs = []
    for b_id, batch in enumerate(batches):
        jobs.append(Process(target = the_function, \
                            args = [b_id, batch]))

    # Run the jobs
    for j in jobs: j.start()
    for j in jobs: j.join()
    print('Done!')

    
    
    
def make_id2syns_dict(case_sensitive_entities_file='data/casesensitive_entities.txt',
        id2syns_outfile='input/id2syns.json',
        core_proteins_file='../output/core_proteins.txt',
        core_id2syns_outfile='input/core_id2syns.json'):
    '''
    FUNCTION:
    - Make a dictionary of entity IDs (keys) to synonyms (values)
    '''
    id2syns = dict()

    with open(case_sensitive_entities_file) as fin:
        for line in fin:
            line = line.strip('\n').replace('_',' ').split('|')

            # ID, Synonyms
            ID = line[0]
            names = line

            # ID -> Synonyms
            id2syns[ID] = names

    json.dump(id2syns, open(id2syns_outfile,'w'))
    
    
    ### Core proteins
    core_proteins = set()
    with open(core_proteins_file) as fin:
        for line in fin:
            core_proteins.add(line.strip())
    print(len(core_proteins))
    ### Core protein ID -> synonyms
    core_id2syns = dict()
    for core_id in core_proteins:
        if core_id in id2syns:
            synonyms = id2syns[core_id]
            core_id2syns[core_id] = synonyms
        else:
            print("Missing synonyms for %s!!"%(core_id))
    json.dump(core_id2syns, open(core_id2syns_outfile,'w'))


