#data handling
import pandas as pd
import numpy as np

#http to reactome
import requests
import json

import os
import re
from tqdm import tqdm

def reactome2reactome(reactome_hierarchy) -> pd.DataFrame:
    reactome_hierarchy_kg = pd.read_csv(reactome_hierarchy, sep = "\t", header = None, names = ["head", "tail"])

    head = [i for i in reactome_hierarchy_kg["head"] if "R-HSA" in i]
    tail = [i for i in reactome_hierarchy_kg["tail"] if "R-HSA" in i]
        
    reactome_hierarchy_kg = pd.DataFrame({"head" : head, "tail" : tail})

    reactome_hierarchy_kg["relation"] = "Reactome_Hierarchy"
    reactome_hierarchy_kg["weight"] = 1
    return reactome_hierarchy_kg