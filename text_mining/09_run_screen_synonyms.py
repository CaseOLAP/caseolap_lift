'''
The purpose of this file is to screen and modify entity synonyms.
- It screens the entity synonyms (checks for potentially ambiguous synonyms).
- It saves the potentially ambiguous synonyms (rem_path). You will look at these
synonyms and determine if they are too ambiguous for you.
- It also may also add some synonyms (add_species_syns).

Suspect = potentially bad = ambiguous
'''
from caseolap._09_screen_synonyms import *


'''
Parameters
'''
# Input
id2syns = 'input/id2syns.json'

# Output
eng_path = 'data/suspect_synonyms_english_words.txt'# Ambiguous syns: English word
short_path = 'data/suspect_synonyms_too_short.txt'  # Ambiguous syns: Short words
rem_path = 'data/remove_these_synonyms.txt'         # Synonyms to remove



'''
Main code
'''
if __name__ == '__main__':
    # Instantiate the class (Load the entity dictionary, synonyms)
    SE = ScreenSynonyms(id2syns)  
    
    # Identifies synonyms to possible remove
    SE.find_suspect_synonyms()   
    
    # Merge the temp files of suspect synonyms
    SE.merge_suspect_syns()

    # Save the synonyms you might want to remove. Inspect and modify this file
    SE.export_suspect_synonyms(eng_path, short_path, rem_path) 

