#data handling
import pandas as pd
import numpy as np

#http to reactome
import requests
import json

import os
import re
from tqdm import tqdm

def protein2triples(protein2protein: os.path) -> pd.DataFrame:
    """
    given the pathway input file, generate a edge list from protein to protein
    """
    ppi_edge_list = pd.read_csv(protein2protein)
    ppi_edge_list = ppi_edge_list.rename(columns = {"h": "head", "r": "relation", "t":"tail"})
    ppi_edge_list["weight"] = 1
    return ppi_edge_list
