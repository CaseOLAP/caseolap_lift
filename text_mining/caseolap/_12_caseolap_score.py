import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, json, os


# Note: "Term Frequency" here doesn't mean relative frequency of a term compared
# to all other terms in the document. It just means the absolute count.

class Caseolap(object):

    def __init__(self, category2pmids, pmid2entity2count, result_dir,
                 categories_path, logfile):

        self.category_names = json.load(open(categories_path))
        self.category2pmids = category2pmids
        self.pmid2entity2count = pmid2entity2count
        self.all_entities = list()
        self.category2pmid2entity2count = dict()
        self.category2entities = dict()
        self.category2entity2count = dict()  # Includes entities with 1+ counts
        self.category2entity2tf = dict()  # Includes entities with 0 counts
        self.category2total_tf = dict()
        self.category2entity2popularity = dict()
        self.category2entity2num_pmids = dict()
        self.category2entity2ntf = dict()
        self.category2entity2ndf = dict()
        self.category2entity2ntf_ndf_ratio = dict()
        self.category2entity2distinctiveness = dict()
        self.category2entity2caseolap = dict()

        self.result_dir = result_dir
        self.result_stat = list()
        self.logfile = logfile

    def df_builder(self, category_quant, file_name):
        '''
        FUNCTION:
        - This takes data and turns it into a dataframe for export.

        PARAMS:
        - file_name: file name
        '''
        flatdata = list()

        # Each entity
        for entity in self.all_entities:
            entity_dict = {'entity': entity}

            # Each category
            for category in self.category_names:
                # Category: Number of entities: entity
                entity_dict[category] = category_quant[category][entity]

            flatdata.append(entity_dict)

            # Process data into a dataframe
        df = pd.DataFrame(flatdata)
        df = df.set_index('entity')
        out_file = os.path.join(self.result_dir,(file_name+'.csv'))
        df.to_csv(out_file)
        return df

    def dump_json(self, data, file_name):
        '''
        FUNCTION:
        - Export a file to a JSON format

        PARAMS:
        - data: Variable to export as a JSON
        - file_name: File name
        '''
        out_file = os.path.join(self.result_dir,(file_name+'.json'))
        with open(out_file, 'w') as fout:
            json.dump(data, fout)

    def print_progress(self, msg):
        '''
        FUNCTION:
        - Print a message and write it to the logfile.

        PARAMS:
        - msg: The message to print
        '''
        self.logfile.write(msg + '\n')
        print(msg)

    def print_categories2pmid(self, all_or_core, dump=False, verbose=False):
        '''
        FUNCTION:
        - Print the categories and number of publications in those categories

        PARAMS:
        - dump: Whether to export the info to a JSON file
        - verbose: Whether to export the info to the screen and a text logfile
        '''
        # Print counts: header
        if verbose:
            self.print_progress('\nCategory: # PMIDs collected\n' + '=' * 27)

        # Iterate through each category and its PMIDs
        total_pmids = set()
        for category, category_pmids in self.category2pmids.items():
            total_pmids = total_pmids.union(category_pmids)

            # Print counts: body
            self.print_progress(category + ': ' + str(len(category_pmids)) + ' PMIDs')
            self.result_stat.append({'Category': category,
                                     'Total pmids collected': len(category_pmids)})

        self.print_progress('Total: ' + str(len(total_pmids)) + ' PMIDs')
        self.result_stat.append({'Category': 'Total',
                                 'Total pmids collected:': len(total_pmids)})

        # Export counts
        if dump:
            # check if dir exists / make
            out_dir = os.path.join(self.result_dir,(all_or_core+'_proteins/'))
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            file_name = all_or_core + '_proteins/' + all_or_core + '_categpmids'
            self.dump_json(self.category2pmids, file_name)

    def map_category2pmid2entity2count(self):
        '''
        FUNCTION:
        - Map category to its PMIDS to its entities to the entity counts
        '''
        # Iterate through each category name and its publication PMIDs
        for category, category_pmids in self.category2pmids.items():

            # For each PMID, map PMID->Entity->Entity_Count
            pmid2entity2count = dict()
            for pmid in category_pmids:
                try:
                    entity2count = self.pmid2entity2count[pmid]
                    pmid2entity2count[pmid] = entity2count
                except:
                    print(pmid)
                    # print(category_pmids)
                    # print(self.pmid2entity2count)

            # Save Category->PMID->Entity->Entity_Count
            self.category2pmid2entity2count[category] = pmid2entity2count

    def get_all_entities(self, prefix, dump=False, verbose=False):
        '''
        FUNCTION:
        - Get all the entities in each category
        - Get all the entities irrespective of category

        PARAMS:
        - dump: Whether to export the info to a JSON file
        - verbose: Whether to export the info to the screen and a text logfile
        '''
        if verbose:
            print('\nCategory: # Entities Found\n' + '=' * 26)

        all_categories_entities = list()

        # Each category
        for category, pmid2ent2cnt in self.category2pmid2entity2count.items():

            # Each PMID in the category
            categorys_entities = list()
            for pmid, ent2cnt in pmid2ent2cnt.items():
                # Each entity in the PMID in the category
                categorys_entities += list(ent2cnt.keys())

                # Save entities in the category
            categorys_entities = list(set(categorys_entities))

            # Save Category->Entities
            self.category2entities[category] = categorys_entities

            # Print counts
            if verbose:
                self.print_progress(category + ' ' + str(len(categorys_entities)))
                self.result_stat.append({'Category Name': category,
                                         ' # Entities': len(categorys_entities)})

            # Save all entities discovered (from this category)
            all_categories_entities += categorys_entities

            # Save all entities discovered (from all categories)
        self.all_entities = list(set(all_categories_entities))

        # Print number of entities
        if verbose:
            self.print_progress('Total entities: ' + str(len(self.all_entities)))

        # Export
        if dump:
            # check if dir exists / make
            out_dir = os.path.join(self.result_dir,(prefix+'_proteins/'))
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            self.dump_json(self.all_entities, file_name=prefix + '_proteins/' + prefix + '_proteins')
            self.dump_json(self.category2entities, file_name=prefix + '_proteins/' + prefix + '_category2proteins')

    def map_entity2pmid_or_tf(self, pmid2entity2count, tf_or_pmid):
        '''
        FUNCTION:
        - Map entity to either term frequency (count) or PMID

        PARAMS:
        - pmid2entity2count (nested dict): Maps PMID to entity to count
        - tf_or_pmid (str): Is either 'tf' or 'pmid'. Indicates whether to
          map entity to TF or PMID
        '''

        map_dict = list()

        # Each PMID
        for pmid, entity2count in pmid2entity2count.items():

            # Each Entity2Count
            for entity, count in entity2count.items():

                # Map entity to term frequency in a PMID or PMID
                if tf_or_pmid == 'tf':
                    map_dict.append({'Entity': entity, 'TF in a PMID': int(count)})
                elif tf_or_pmid == 'pmid':
                    map_dict.append({'Entity': entity, 'PMID': pmid})

        return map_dict

    def map_entity2counts_in_category(self, list_of_dicts, operation):
        '''
        FUNCTION:
        - Start with the entity to the counts in each PMID in the category
        - Map the entity to the counts in the whole category.

        PARAMS:
        - list_of_dicts: List of dicts of {'Entity':'the_entity', 'tf':0}. Turned
          into a dataframe.
        - operation: Whether to take the sum or the count of the entities
        '''
        # Dataframe of entity2count dictionary
        df = pd.DataFrame(list_of_dicts)
        col = df.columns
        df = df.set_index(col[0])

        # Group by entity, take the sum
        if operation == 'sum':
            gdf = df.groupby(col[0]).sum()

        # Group by entity, take the count
        elif operation == 'count':
            gdf = df.groupby(col[0]).count()

        # Pre-mapping
        entity2count = dict()
        entities = list(gdf.index)
        total_counts = list(gdf[col[1]])

        # Map the 'Entity' to 'counts in the category'
        for entity, total_count in zip(entities, total_counts):
            entity2count[entity] = total_count

        return entity2count

    def get_entity_counts_per_category(self):
        '''
        FUNCTION:
        - For each category, get the entity's term frequency (counts)
        '''
        # For each category...
        for category, pmid2entity2count in self.category2pmid2entity2count.items():
            # Map entity to count in each PMID
            entity2counts = self.map_entity2pmid_or_tf(pmid2entity2count, 'tf')

            # Map entity to counts in all PMIDs
            entity2count = self.map_entity2counts_in_category(entity2counts, 'sum')
            self.category2entity2count[category] = entity2count

    def category2entity2tf_finder(self):
        '''
        FUNCTION:
        - Maps category to its entities to their counts. Unlike
          category2pmid2entity2count, category2entity2tf includes entities
          with 0 counts.
        '''
        # Category
        for category, entity2nonzero_count in self.category2entity2count.items():
            entity2all_count = dict()

            # Entity
            for entity in self.all_entities:

                # Count
                try:
                    count = entity2nonzero_count[entity]
                except:
                    count = 0

                # Entity->Count
                entity2all_count[entity] = int(count)

            # Category->Entity->Count
            self.category2entity2tf[category] = entity2all_count

    def calculate_all_popularity_scores(self, all_or_core, dump=False):
        '''
        FUNCTION:
        - Calculate the popularity scores for all entities

        PARAMS:
        - dump: Indicates whether to export the data to a JSON file
        '''

        # Category
        for category, entities2tf in self.category2entity2tf.items():

            # Calculate Entity TF: Total counts of the entities in the category
            tf_all_ents = sum(entities2tf.values())

            # Save entity TF
            self.category2total_tf[category] = tf_all_ents

            # Entity: Calculate entities' popularity score
            entity2popularity = dict()
            for entity, tf_this_ent in entities2tf.items():
                popularity_score = np.log(tf_this_ent + 1) / np.log(tf_all_ents)
                entity2popularity[entity] = popularity_score

                # Entity: Save entities' popularity score
            self.category2entity2popularity[category] = entity2popularity

        # Export category->entity->popularity score
        if dump:
            file_name = all_or_core + '_proteins/' + all_or_core + '_popularity_score'
            self.df_builder(self.category2entity2popularity, file_name)

    def category2entity2num_pmids_finder(self):
        '''
        FUNCTION:
        - In each category, map the entity to the number of PMIDs it is in(?)
        '''
        for category, pmid2entity2count in self.category2pmid2entity2count.items():
            # Map entity to a PMID
            entity2pmid = self.map_entity2pmid_or_tf(pmid2entity2count, 'pmid')

            # Map entity to number of PMIDs
            entity2num_pmids = self.map_entity2counts_in_category(entity2pmid, 'count')

            # Save
            self.category2entity2num_pmids[category] = entity2num_pmids

    def calculate_category2entity2ntf(self):
        '''
        FUNCTION:
        - Calculate normalized term frequency (ntf) for each entity
          in each category.
        '''
        k1 = 1.2
        b = 0.75

        # Category
        for category, entity2tf in self.category2entity2tf.items():

            # Total Number of Hits / Total Number of Entities
            total_num_hits = self.category2total_tf[category]
            total_num_ents = len([tf for tf in entity2tf.values() if tf > 0])
            hit_to_entity_ratio = total_num_hits / total_num_ents

            # Entity
            entity2ntf = dict()
            for entity, tf in entity2tf.items():
                # Calculate NTF
                ntf = (tf * (k1 + 1)) / \
                      (tf + (k1 * (1 - b + (b * (total_num_hits / hit_to_entity_ratio)))))
                entity2ntf[entity] = ntf

            # Save calculation: Category->Entity->NTF
            self.category2entity2ntf[category] = entity2ntf

    def calculate_category2entity2ndf(self):
        '''
        FUNCTION:
        - Calculate normalized document frequency (ndf) for each
          entity in each category. If the entity occurs in the most
          PMIDs (documents), it will have an NDF = 1. If it occurs in
          less documents, it will have an NDF closer to 0.
        '''

        # Category
        for category, entity2num_pmids in self.category2entity2num_pmids.items():
            all_pmid_counts = list()
            entity2ndf = dict()

            # List of #PMIDs for each entity
            all_pmid_counts = [num_pmids for num_pmids in entity2num_pmids.values()]

            # Number of PMIDs for the most popular entity
            max_pmids_any_entity = max(all_pmid_counts)

            # Entity
            for entity in self.all_entities:

                # Calculate entity-#pmids rank (normalized from >0 - 1)
                if entity in self.category2entities[category]:
                    pmids_this_entity = entity2num_pmids[entity]
                    ndf = np.log(1 + pmids_this_entity) / \
                          np.log(1 + max_pmids_any_entity)
                else:
                    ndf = 0

                # Entity->NDF
                entity2ndf[entity] = ndf

                # Save calculation: Category->Entity->NDF
            self.category2entity2ndf[category] = entity2ndf

    def calculate_entity2ntf_ndf_ratio(self):
        '''
        FUNCTION:
        - Map category to entity to ratio of normalized term frequency
          over normalized document frequency
        '''

        # Category
        for category, entity2ntf in self.category2entity2ntf.items():
            entity2ntf_ndf_ratio = dict()

            # Entity
            for entity in self.all_entities:
                # Normalized Term Frequency
                ntf = entity2ntf[entity]

                # Normalized Document Frequency
                ndf = self.category2entity2ndf[category][entity]

                # Normalized Term Frequency / Normalized Document Frequency
                ntf_ndf_ratio = ntf * ndf
                entity2ntf_ndf_ratio[entity] = ntf_ndf_ratio

            # Save
            self.category2entity2ntf_ndf_ratio[category] = entity2ntf_ndf_ratio

    def calculate_all_distinctiveness_scores(self, all_or_core, dump=False):
        '''
        FUNCTION:
        - Calculate the distinctiveness score.

        PARAMS:
        - dump: Indicate whether to export the scores to a JSON
        '''

        ''' Calculate numerator '''
        # Numerator = e^(entity-category ntf/ndf)-1
        # NOTE: Added '-1' so that the entities with zero hits will have a zero score
        category2entity2exp_ntf_ndf_ratio = dict()

        # Category
        for category, entity2ntf_ndf_ratio in self.category2entity2ntf_ndf_ratio.items():
            entity2exp_ntf_ndf_ratio = dict()

            # Entity->Ratio
            for entity, ntf_ndf_ratio in entity2ntf_ndf_ratio.items():
                entity2exp_ntf_ndf_ratio[entity] = np.exp(ntf_ndf_ratio)

            # Category->Entity->Ratio
            category2entity2exp_ntf_ndf_ratio[category] = entity2exp_ntf_ndf_ratio

        ''' Calculate denominator '''
        # Denom = 1 + sum e^(ntf/ndf for this entity in all other categories)
        entity2denominator = dict()

        # Entity
        for entity in self.all_entities:

            # Category
            denominator = 1
            for category in self.category_names:
                denominator += category2entity2exp_ntf_ndf_ratio[category][entity]

            entity2denominator[entity] = denominator

        ''' Calculate score (numerator/denominator) '''
        # Category
        for category, ent2exp_ntf_ndf_ratio in category2entity2exp_ntf_ndf_ratio.items():
            entity2distinctiveness = dict()

            # Entity
            for entity, exp_ntf_ndf_ratio in ent2exp_ntf_ndf_ratio.items():
                # Entity-Category Score
                denominator = entity2denominator[entity]
                entity2distinctiveness[entity] = (exp_ntf_ndf_ratio - 1) / denominator

            # Save
            self.category2entity2distinctiveness[category] = entity2distinctiveness

        # Export
        if dump:
            file_name = all_or_core + '_proteins/' + all_or_core + '_distinctiveness_score'
            self.df_builder(self.category2entity2distinctiveness, file_name)

    def calculate_caseolap_score(self, all_or_core, dump=False):
        '''
        FUNCTION:
        - Calculate the CaseOLAP score of each entity-category pair

        PARAMS:
        - caseolap_name: Name of the CaseOLAP spreadsheet/dataframe
        - dump:
        '''

        # Category
        for category, ent2distinct in self.category2entity2distinctiveness.items():
            entity2caseolap = dict()

            # Entity & Distinctiveness
            for entity, distinctiveness in ent2distinct.items():
                # Popularity
                popularity = self.category2entity2popularity[category][entity]

                # CaseOLAP = Distinctiveness*Popularity
                entity2caseolap[entity] = distinctiveness * popularity

            # Save CaseOLAP Score
            self.category2entity2caseolap[category] = entity2caseolap

        # Export
        if dump:
            file_name = all_or_core + '_proteins/' + all_or_core + '_caseolap'
            self.df_builder(self.category2entity2caseolap, file_name)
            self.dump_json(self.category2entity2caseolap, file_name)
            self.dump_json(self.result_stat, all_or_core + '_proteins/' + all_or_core + '_result_stat')


