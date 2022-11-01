#data handling
import pandas as pd
import numpy as np

#http to reactome
import requests
import json

import os
import re
from tqdm import tqdm

def pathway2triples(protein2pathway: os.path) -> pd.DataFrame:
    """
    given the pathway input file, generate a edge list from protein to pathway
    """
    pathway_edge_list = pd.read_csv(protein2pathway)
    pathway_edge_list = pathway_edge_list.rename(columns = {"h": "head", "r": "relation", "t":"tail"})
    pathway_edge_list["weight"] = 1
    return pathway_edge_list





