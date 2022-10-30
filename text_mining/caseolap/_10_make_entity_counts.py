import json

class MakeEntityCounts(object):
    
    def __init__(self, id2syns_path, remove_syns_infile, pmid_syn_count_in):
        '''
        FUNCTION:
        Read in synonyms to remove, pmid2syncount2syn mapping,
        id2syn mapping. Make syn2id mapping, 

        PARAMS:
        - remove_syns_infile: Synonyms to remove and not query.
        - entity_dict_path: entity ID to synonyms mapping
        - pmid_syn_count_out: Output file. PMID -> Synonym -> Synonym counts
        '''
        
        ### Synonyms to not consider in the search
        with open(remove_syns_infile, 'r') as fin:
            self.bad_syns = [syn.strip() for syn in fin.readlines()]
            

        ### Map entity ID -> synonyms
        id2syns = json.load(open(id2syns_path))
        
        ### Map synonyms -> entity IDs
        self.syn2id = dict()
        for ID,syns in id2syns.items():
            for syn in syns:
                if syn not in self.bad_syns:
                    self.syn2id.setdefault(syn,[]).append(ID)

        # Map PMIDs -> synonyms -> synonym counts
        self.pcs = json.load(open(pmid_syn_count_in))

        
        
    def entitycount(self, entitycount_outfile):
        '''
        FUNCTION:
        Writes filtered entitycount.txt:
        PMID EntityID1|Count ... EntityIDn|Count
        
        PARAMS:
        - entitycount_outfile: entitycount output
        '''
        with open(entitycount_outfile,'w') as fout:
            # For each PMID
            for pmid in self.pcs:
                tempd = dict()
                
                # For each good synonym in the publication
                for synonym in self.pcs[pmid]:
                    if synonym not in self.bad_syns:
                        
                        # Get the current synonym count in current publication
                        syn_count = self.pcs[pmid][synonym]

                        # Map the good synonym to its entity ID(s)
                        for entity_id in self.syn2id[synonym]:
                            tempd[entity_id] = tempd.get(entity_id,0) + syn_count 
                
                # Save mappings: PMID entity_ID|count ...
                if len(tempd) > 0:
                    fout.write(pmid+' ')
                    for entity_id in tempd:
                        ent_count = tempd[entity_id]
                        fout.write(entity_id+'|'+str(ent_count)+' ')
                    fout.write('\n')
                
                