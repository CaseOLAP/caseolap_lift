import json, pandas as pd, numpy as np, csv, os

# Note: counts = #hits, the number of times a synonym was found in the text

class InspectEntityScores(object):
    '''
    This class's purpose is to display the ranked discovered synonyms and 
    ranked scoring  entities to the user. These are stored in 'ranked_syns_out' 
    and 'ranked_ent_path'
    '''

    def __init__(self, scores_in, id2syns_in, remove_syns_in, \
                 pmid_syn_count_in, cat2pmids_in):
        '''
        FUNCTION:
        Reads in input.

        PARAMS:
        - scores_in: Input file. CaseOLAP scores csv
        - id2syns_in: Input file. EntityID->Synonyms dict
        - remove_syns_in: Input file. Synonyms not used in the search.
        - pmid_syn_count_in: Input file. PMID->Synonyms->SynonymCounts
        - cat2pmids_in: Input file. Category->PMIDS
        '''

        # {Entity:{Syn1:count1, Syn2:count2, ...},...}
        self.ranked_ents_counts = dict()

        # CaseOLAP score dataframe (entities x categs)
        self.score_dfs = dict()

        # Entity:Score dictionary
        self.ent2score_dicts = dict()

        # Synonyms:count dictionary
        self.syn2count = dict()

        # ID:Synonyms (entity dictionary)
        self.id2syns = json.load(open(id2syns_in))

        # Categories
        scores = pd.read_csv(scores_in)
        cols = list(scores.columns)
        cols.append('TOTAL')
        cols.sort()
        self.cols = cols  # Categories scores (i.e. score df columns)

        # Synonyms not considered in the search
        with open(remove_syns_in, 'r') as fin:
            self.bad_syns = [syn.strip() for syn in fin.readlines()]

        # Get "pmid->synonym->counts" dictionary
        self.pmid2syn2count = json.load(open(pmid_syn_count_in, 'r'))

        # Get Category->PMID
        cat2pmid = json.load(open(cat2pmids_in, 'r'))

        # Get PMID->Category
        pmid2cat = dict()
        for cat, pmids in cat2pmid.items():
            for pmid in pmids:
                pmid2cat[pmid] = pmid2cat.setdefault(pmid, list())
                pmid2cat[pmid].append(cat)
        self.pmid2cat = pmid2cat

    '''
    PART 1: Getting the ranked synonyms
    '''

    def get_ranked_synonyms_found(self, cat2pmids_in, pmid_syn_count_in,
                                  ranked_syns_out, export=True):
        '''
        FUNCTION:
        Find and display the ranked discovered synonyms: the entity names that were
        found the most in the text. This is not the ranked entities.
        This is 1 of 2 end purposes of the file.

        PARAMS:
        - syn_count_path: Input file. The place where the synonym counts are stored.
        - ranked_syns_out: Output file. Where the ranked synonyms are stored.
        - export: Whether to export the ranked synonyms
        '''
        # Initialize synonym->counts
        syn2count = dict()
        for cat in self.cols:
            syn2count[cat] = dict()

        # Get each PMID's category
        for pmid, syn2cnt in self.pmid2syn2count.items():

            # Map category->synonym->counts
            for syn, count in syn2cnt.items():
                if syn in self.bad_syns:
                    continue

                # 'TOTAL' category
                syn2count['TOTAL'][syn] = syn2count['TOTAL'].get(syn, 0) + count

                # PMID's category/categories
                cats = self.pmid2cat[pmid]
                for cat in cats:
                    syn2count[cat][syn] = syn2count[cat].get(syn, 0) + count

        # Sort and export each category's synonym counts
        for cat in self.cols:
            if cat == 'entity':
                continue

            # Sort "synonym->counts" by counts
            syn2count[cat] = dict(sorted(syn2count[cat].items(), \
                                         key=lambda x: x[1], reverse=True))

            # Export category's sorted "synonym->counts"
            if export == True:
                path = ranked_syns_out.split('.txt')[0] + '_' + cat + '.txt'
                # make directory
                path_dir = os.path.dirname(os.path.abspath(path))
                if not (os.path.exists(path_dir)):
                    os.makedirs(path_dir)
                with open(path, 'w') as fout:
                    for syn, count in syn2count[cat].items():
                        fout.write(syn + ': ' + str(count) + '\n')

        self.syn2count = syn2count


    '''
    PART 2: Getting the ranked entities
    '''

    def sort_all_scores(self, scores_in):
        '''
        FUNCTION:
        Sort the entities by CaseOLAP scores (high to low)

        PARAMS:
        - scores_path: Where the CaseOLAP score table is (entities x category).
        '''
        # Read in CaseOLAP scores csv
        scores = pd.read_csv(scores_in)
        cols = scores.columns

        # Create 'total' score category
        scores['TOTAL'] = scores.loc[:, cols[1]:cols[len(cols) - 1]].sum(axis=1)

        # Sort entities by scores (high to low)
        for cat in self.cols:
            self.score_dfs[cat] = scores.sort_values(by=[cat], ascending=False)

            ### Map EntityID->score

        # For each category
        for cat in self.cols:
            if cat == 'entity':
                continue

            # Initialize Entity->Score dictionary
            self.ent2score_dicts[cat] = dict()

            # For each entity (high to low)
            for i in range(0, len(scores)):
                # Entity ID
                ent = self.score_dfs[cat]['entity'].iloc[i]

                # Score
                scre = self.score_dfs[cat][cat].iloc[i]

                # Entity->Score
                self.ent2score_dicts[cat][ent] = scre

    def prop_entities_found(self):
        '''
        FUNCTION:
        Prints the proportion of entities searched for / entities found
        '''
        # Entity IDs/names with scores (i.e. entities found in the text)
        self.ents = list(self.score_dfs['TOTAL']['entity'])

        # Entity IDs/names without scores (i.e. entities not found in the text)
        searched_ents = len(list(self.id2syns.keys()))

        # Print proportion: entities searched for / entities found
        print("{:,}".format(len(self.ents)) + '/' + "{:,}".format(searched_ents) + \
              ' (Found/Searched) Entities \n' + \
              "{:.2f}".format(100 * (len(self.ents)) / searched_ents) + '%')

    def get_entity_syn_counts(self):
        '''
        FUNCTION:
        Saves the ranked entites, their synonyms, and their synonym counts
        in a nested dictionary
        '''
        syn2count = self.syn2count

        for cat in self.cols:
            if cat == 'entity':
                continue

            # Entities ranked by a category
            self.ranked_ents_counts[cat] = dict()
            ents = list(self.score_dfs[cat]['entity'])

            # For each entity
            for ent in ents:
                # Get its synonyms
                syns = [syn for syn in self.id2syns[ent] if syn in syn2count[cat]]

                # Get its synonyms' counts
                ent_counts = {syn: syn2count[cat][syn] for syn in syns}

                # Sort its synonyms' counts
                ent_counts = sorted(ent_counts.items(), key=lambda x: x[1], reverse=True)

                # Save the {Entity:{Syn1:count1, Syn2:count2, ...},...}
                self.ranked_ents_counts[cat][ent] = dict(ent_counts)

    def rank_entites(self, cat, score_df, ranked_ent_path, score_component, print_msg):
        '''
        FUNCTION:
        Prints and writes the ranked entities, their score, their synonyms, and their
        synonym counts for the number of times the synonyms were found in the text.
        This is 1 of 2 end purposes of the file.

        PARAMS:
        - ranked_ent_path: Output file. Where the ranked ranking entities
          will be saved.
        '''
        # make directory
        path_dir = os.path.dirname(os.path.abspath(ranked_ent_path))
        if not (os.path.exists(path_dir)):
            os.makedirs(path_dir)

        with open(ranked_ent_path, 'w') as fout:
            # Print what the format will look like for this document
            output_msg('\nFormat:\n' + \
                       '# [Rank] [ID] | Score [x.xxxx]\n' + \
                       '===============================\n' + \
                       '#counts in all categories: synonym\n\n\n' + \
                       '*********************************\n' + \
                       '******** Ranked Entities ********\n' + \
                       '*********************************', fout, print_msg)

            # Output each entity, its synonyms, and their counts
            for i, (ID, syn_counts) in enumerate(self.ranked_ents_counts[cat].items()):

                # Output header
                score = str("{:,.4f}".format(self.ent2score_dicts[cat][ID]))
                output_msg('\n#' + str(i + 1) + ' ' + ID + ' | ' + score + \
                           ' Total ' + score_component + '\n' + \
                           '=============================================', \
                           fout, print_msg)

                # Output counts: synonym
                for syn, count in syn_counts.items():
                    output_msg(str(count) + ': ' + syn, fout, print_msg)

            ## check if dir exists / make
            #out_dir = os.path.join(self.result_dir,(all_or_core+'_proteins/'))
            #if not os.path.exists(out_dir):
            #    os.makedirs(out_dir)
            path = ranked_ent_path[:-4].split('/')
            path = '/'.join(path[:-1]) + '/dictionary_' + path[-1] + '.json'
            json.dump(self.ranked_ents_counts[cat], \
                      open(path, 'w'))

    def rank_each_category_score(self, ranked_ent_in, score_component):
        '''
        FUNCTION:
        Outputs the top scoring entities for each category.

        PARAMS:
        - ranked_ent_in: Path to the table of ranked entities
        - score_component: Indicate what the score is ranking
          popularity, distinctiveness, or combined (CaseOLAP score)
        '''
        print('\n' + '=' * 25 + score_component + '=' * 25)

        cats = self.cols

        # Call the rank_entities function for each category
        for i, score_df in enumerate(self.score_dfs):
            # Don't sort by 'entity' column. There aren't scores here.
            if cats[i] == 'entity':
                continue

            # Print the entities ranked by total score
            if cats[i] == 'TOTAL':
                print_msg = True
            else:
                print_msg = False

            # Save entities ranked by each category
            path = ranked_ent_in.split('.txt')[0] + '_' + cats[i] + '.txt'
            self.rank_entites(cats[i], score_df, path, score_component, print_msg)


def output_msg(msg, fout, print_msg):
    if print_msg == True:
        print(msg)
    fout.write(msg + '\n')
