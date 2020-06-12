import sys 
import pandas as pd
import numpy as np
from scipy.stats import hypergeom

####### controllo l'input 
if len(sys.argv) < 4:
    print("invalid arguments: \nusage is: {} annotation_table subset_list genome_size".format(sys.argv[0]))
    quit()

####### caricamento dei parametri 

proteome_annotation_file = sys.argv[1]
subset_file = sys.argv[2]
proteome_size = int(sys.argv[3])

####### lettura della tabella con pfam e id 

annotated_proteome = pd.read_csv(proteome_annotation_file, sep = "\t", names= ["pfam_id","seq_id"])

####### lettura del file con il subset

observed_list = set() # usiamo un set per ottimizzare il tempo di esecuzione

with open(subset_file, "r") as subset: # con questo blocco il file viene aperto in lettura, letto e poi chiuso
    for line in subset:
        observed_list.add(line.rstrip()) # per ogni linea del file leggo l'id corrispondente, tolgo il ritorno a capo e lo aggiungo al set degli osservati

####### creazione delle tabelle di conta dei pfam

pfam_dict = {} # uso un dizionario per raccogliere le conte relative a ciascun id PFAM

for pfamid in (annotated_proteome.pfam_id): # popolazione del dizionario, per ogni elemento
    try:                                    # prova ad incrementare la conta,
        pfam_dict[pfamid] += 1 
    except:                                 # se non esiste allora conta la prima occorrenza
        pfam_dict[pfamid] = 1

####### conversione del dizionario in dataframe

####### creazione del dataframe delle osservazioni

observed_dataframe = annotated_proteome[annotated_proteome.seq_id.isin(observed_list)] # usando il dataframe del proteoma si annotano le sequenze del subset

####### creazione del dataframe delle osservazioni

observed_dict = {}
for pfamid in (observed_dataframe.pfam_id):
    try: 
        observed_dict[pfamid] += 1 
    except:
        observed_dict[pfamid] = 1

universe_dataset = pd.DataFrame.from_dict(pfam_dict, orient = "index", columns=["total_counts"])

obs_dataset = pd.DataFrame.from_dict(observed_dict, orient = "index", columns=["observed_counts"]) 


####### unione dei due dataframe

dataframe = universe_dataset.merge(obs_dataset, how =  "left", left_index=True, right_index=True )

dataframe[dataframe.isna() == True] = 0 # sostituzione dei NaN con 0


####### test ipergeometrico 

N = dataframe.total_counts.sum() + (proteome_size - len(dataframe))
n = dataframe.observed_counts.sum()

dataframe["pvalue"] = hypergeom.sf(dataframe.observed_counts, N, dataframe.total_counts, n)


####### stampa del risultato del test (solo degli elementi che hanno più di una osservazione)

print(dataframe[dataframe.observed_counts > 1] && dataframe[dataframe.pvalue] < 0.05)


