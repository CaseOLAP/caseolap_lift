import json, os, nltk
#!pip install nltk
#nltk.download('words')
from nltk.corpus import words
from multiprocessing import cpu_count,Process


class ScreenSynonyms(object):
    '''
    ScreenSynonyms screens the entites to see if the entity names/synonyms are 
    poor quality. If they are poor quality (i.e. potentially false positives), 
    ScreenSynonyms picks them out and presents them to the user, storing them 
    in rem_path.
    '''
    
    def __init__(self, id2syns):
        '''
        Initializing object
        '''
        # ID -> Synonyms
        self.id2syns = json.load(open(id2syns)) # Entity dict {ID:[synonyms],...}
        
        # Synonyms
        synonyms = list()
        for syns in list(self.id2syns.values()):
            synonyms += syns
        synonyms = list(set(synonyms))
        self.synonyms = synonyms       # All synonyms
        
        # Suspect synonyms
        self.sus_syns = set()  # Potentially ambiguous synonyms

        
    '''
    Checking synonyms
    '''
    
    def check_if_short(self, syn):
        '''
        FUNCTION: 
        This checks if the word is short
        Note: You may want to change how short is 'too short'
        Too short = 1 or 2 characters or 3 letters

        PARAMS: 
        - syn = the word to be checked, an entity's synonym
        '''
        # If the synonym is multiple words, it's not a short word
        one_word = " " not in syn
        if one_word == False:
            short_word = False
        else:
            short_word = (len(syn) <= 2) or (len(syn) <=3 and syn.isalpha()) 

        return short_word
    
    
    
    def check_if_eng_word(self, syn):
        '''
        FUNCTION: 
        This checks if the word is an English word

        PARAMS: 
        - syn = the word to be checked, an entity's synonym
        '''
        # Default
        eng_word = False

        # Make an ostensibly singular version of plurals that end with s
        if len(syn) > 1 and syn.endswith('s'): singular_syn = syn[:len(syn)-2]
        else: singular_syn = ''

        # If the word is lowercase, ignoring the first letter
        non_1st_low = syn[1:].islower()
        if non_1st_low:

            # If the given word is English
            if syn.lower() in words.words():
                eng_word = True

            # If an ostensibly singular version of the word is English
            elif singular_syn != '' and singular_syn in words.words(): 
                eng_word = True

        return eng_word    
    
    
    
    def check_suspect_synonyms(self, b_id, syns):
        '''
        FUNCTION: 
        This checks a synonym to see if it's too ambiguous. 
        If the synonym is too short or is another English word, it is saved.
        
        PARAMS: 
        - syn = synonym
        - b_id = the batch ID in the threading process. Uses to create temp 
          files which are merged after threading is complete.
        '''
        for syn in syns:
            short = self.check_if_short(syn)
            english = self.check_if_eng_word(syn)
            ambiguous = short or english
            if ambiguous:
                temp_file = 'data/suspect_synonyms_'+str(b_id)+'.txt'
                open(temp_file,'a').write(syn+'\n')

         
        

    def find_suspect_synonyms(self):    
        '''
        FUNCTION: 
        Find synonyms which are potential false positives.
        These synonyms might be acronyms that could stand for other terms. 
        Calls the 'check_suspect_synonyms' function
        '''
        # How many processors can be used
        procs = cpu_count() 
        
        ### Clear temp files if they are left over from last time
        for b_id in range(procs):
            path = 'data/suspect_synonyms_'+str(b_id)+'.txt'
            open(path,'w')
        
        ### FIND THE POTENTIALLY BAD SYNONYMS 
        batch = [[] for i in range(procs)]    
        tot = len(self.synonyms)

        # For each synonym
        for i, syn in enumerate(self.synonyms):
            
            # Add synonym to a batch
            b_id = i%procs
            batch[b_id].append(syn)
            
            # If it's the end of synonym list
            if i+1 == tot:
                print("Running jobs...")
                jobs = []

                # Add each job checking synonyms to a list
                for b_id, syns in enumerate(batch):
                    jobs.append(Process(target = self.check_suspect_synonyms,\
                                        args = [b_id, syns]))
                # Run the jobs
                for j in jobs: j.start()
                for j in jobs: j.join()
        
        
        
    def merge_suspect_syns(self):
        '''
        FUNCTION:
        Merges the temp files of the potentially bad synonyms.
        '''
        procs = cpu_count() # How many processors can be used
        
        # For each temp file
        for b_id in range(procs):
            # Open temp file path
            path = 'data/suspect_synonyms_'+str(b_id)+'.txt'
            with open(path ,'r') as fin:
                
                # Read in each synonym
                for line in fin:
                    synonym = line.strip()

                    # Save potential bad synonym
                    self.sus_syns.add(synonym)   

                # Delete old temp file
                os.remove(path)
    
    
    
    
    
    '''
    Exporting potentially bad synonyms
    '''

    def export_suspect_synonyms(self, eng_path, short_path, rem_path):
        '''
        FUNCTION: 
        This exports the potentially false positive synonyms to a file 
        which the user will manually inspect     

        PARAMS: 
        - eng_path = path to where ambiguous synonyms that are English words are
        - short_path = path to where the ambiguous synonyms that are short words are
        - rem_path = path to where the English and short synonyms are stored. 
          The user manually inspects this file.
        '''
        print("Exporting the suspect synonyms to data/remove_these_synonyms.txt")
        print('Inspect this file when the process is complete\n')
        
        # Remove old files if left over
        open(short_path,'w'); open(eng_path,'w'); open(rem_path,'w')
        
        ### FIND AND SAVE AMBIGUOUS SYNONYMS
        
        # Potentially bad synonym sets (English, short)
        sus_syns_english, sus_syns_short = set(), set()
        # Potentially bad synonyms (English + short)
        sus_syns = list(self.sus_syns) 
        len_ss = len(sus_syns)
        
        # Check each potentially bad synonym
        for i, sus in enumerate(sus_syns):
            print(i,'/',len_ss, end='\r')
            
            # Check if the suspect synonym is ambiguous
            short = self.check_if_short(sus)
            english = self.check_if_eng_word(sus)

            # Save English synonym
            if english:
                sus_syns_english.add(sus)
                
            # Save short synonym
            elif short:
                sus_syns_short.add(sus)
        
        
        ### EXPORT ALL THE POTENTIALLY BAD SYNONYMS
                
        # Export English syns to a file
        sus_syns_english = list(sus_syns_english)
        sus_syns_english.sort()
        for syn in sus_syns_english:
            open(eng_path,'a').write(syn+'\n')
            
        # Export short syns to a file'
        sus_syns_short = list(sus_syns_short)
        sus_syns_short.sort()
        for syn in sus_syns_short:
            open(short_path,'a').write(syn+'\n')

        # Export English and short syns to a file
        bad_syns = [] 
        with open(rem_path,'w') as fout:
            
            # Merge bad syn lists into one
            bad_syns += sus_syns_english + sus_syns_short
            
            # Write to file
            for sus in bad_syns:
                fout.write(sus+'\n')
                
        print('Synonyms in '+rem_path+' will be removed. Check here.')
        
            
