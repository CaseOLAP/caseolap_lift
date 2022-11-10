#stl
import os
import re
from tqdm import tqdm

#data handling
import pandas as pd
import numpy as np
import json



def transcription_protein2triples(PATH_JSON: str) -> pd.DataFrame:
    """
    given the transcription factor json file, create a head relation tail dataframe
    """

    with open(PATH_JSON, 'r') as file:
        tf_dict = json.loads(file.read())

    head = []
    tail = []
    for key in tf_dict:
        for element in tf_dict[key]:
            head.append(key)
            tail.append(element)

    #relationship
    relation = ["transcription_factor" for i in head]

    weight = [1 for i in relation]

    tf_kg = pd.DataFrame({"head" : head, "relation" : relation, "tail" : tail, "weight" : weight})
    return tf_kg