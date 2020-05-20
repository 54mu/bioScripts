from Bio import SeqIO as sio
import sys
filename = sys.argv[1]
with open(".".join(filename.split(".")[:-1])+".tsv", "w") as outfile:
    outfile.write("seq_id\tGC_content\n")
    for seq in sio.parse(filename, "fasta"):
        gc = 100*(str(seq.seq).count("C") + str(seq.seq).count("G"))/len(seq.seq)
        outfile.write("{}\t{}\n".format(seq.id, gc))
