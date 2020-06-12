import sys 
import pandas as pd
import numpy as np
from scipy.stats import hypergeom

if len(sys.argv) < 4:
    print("invalid arguments: \nusage is: {} annotation_table subset_list genome_size".format(sys.argv[0]))
    quit()


proteome_annotation_file = sys.argv[1]
subset_file = sys.argv[2]
proteome_size = int(sys.argv[3])


annotated_proteome = pd.read_csv(proteome_annotation_file, sep = "\t", names= ["pfam_id","seq_id"])

observed_list = set()
with open(subset_file, "r") as subset:
    for line in subset:
        observed_list.add(line.rstrip()[1:])
    

pfam_dict = {}
for pfamid in (annotated_proteome.pfam_id):
    try: 
        pfam_dict[pfamid] += 1 
    except:
        pfam_dict[pfamid] = 1

observed_dataframe = annotated_proteome[annotated_proteome.seq_id.isin(observed_list)]

observed_dict = {}
for pfamid in (observed_dataframe.pfam_id):
    try: 
        observed_dict[pfamid] += 1 
    except:
        observed_dict[pfamid] = 1

universe_dataset = pd.DataFrame.from_dict(pfam_dict, orient = "index", columns=["total_counts"])

obs_dataset = pd.DataFrame.from_dict(observed_dict, orient = "index", columns=["observed_counts"]) 

dataframe = universe_dataset.merge(obs_dataset, how =  "left", left_index=True, right_index=True )

dataframe[dataframe.isna() == True] = 0

N = dataframe.total_counts.sum() + (proteome_size - len(dataframe))
n = dataframe.observed_counts.sum()

dataframe["pvalue"] = hypergeom.sf(dataframe.observed_counts, N, dataframe.total_counts, n)

print(dataframe[dataframe.observed_counts > 1])

