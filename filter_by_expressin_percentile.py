import pandas as pd 
from Bio import SeqIO as SIO
import sys

expr_file = sys.argv[1] 
transcriptome = sys.argv[2] 
percentile = int(sys.argv[3])

expr_table = pd.read_csv(expr_file, sep = "\t")

sorted_expr = expr_table.sort_values(by="NumReads")

sorted_expr["Cumultaive_TPM"] = sorted_expr.TPM.cumsum()

export_set = set(sorted_expr[sorted_expr.Cumultaive_TPM < percentile*1000000/100].Name)

for s in SIO.parse(transcriptome, "fasta"):
    if s.id not in export_set:
        print(">{}\n{}\n".format(s.id, s.seq))
    else:
        pass
