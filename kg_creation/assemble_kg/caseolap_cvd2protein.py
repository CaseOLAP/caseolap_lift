import pandas as pd
import os

from tqdm import tqdm


def caseolap2triples(caseolap_csv_file: os.path) -> pd.DataFrame:
    """
    creates triples out of caseolap dataframe
    """
    caseolap_df = pd.read_csv(caseolap_csv_file)
    caseolap_df = caseolap_df
    cvd_type = caseolap_df.columns.to_list()
    cvd_type.remove("entity")
    caseolap_df = caseolap_df.to_dict("records")
    caseolap_kg = pd.DataFrame()

    #iterate throught the rows and create a mini dataframe for each protein
    #after, attach that mini dataframe to the caseolap dataframe
    #so it doesnt crash lol using batches
    for row in caseolap_df:

        #assemble kg
        holder_df = pd.DataFrame()
        holder_df["head"] = [i for i in cvd_type]
        holder_df["relation"] = ["CaseOLAP_score" for i in cvd_type]
        holder_df["tail"] = [row["entity"] for i in cvd_type]
        holder_df["weight"] = [row[i] for i in cvd_type]

        #update master kg
        caseolap_kg = pd.concat([caseolap_kg, holder_df])

    caseolap_kg = caseolap_kg[caseolap_kg["weight"] != 0]
    return caseolap_kg





