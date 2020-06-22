# pipeline per l'arricchimento di subset

## 1. hmmer del set intero in fasta contro PFAM

```
hmmsearch --domtblout domains.tbl --cpu 16 ~/Shared/Databases/pfam/Pfam-A.hmm Sbr_peptide_set.fa
```

## 2. Estrazione coppie nome-pfam

attenzione che non ci siano coppie ripetute!

```
sed 's/\s\{2,\}/\t/g' domains.tbl | grep -v ^# | cut -f1,4 | sort | uniq
```

## 3. Test vero e proprio

da ottimizzare lo script ma si comporta bene. vuole in ordine: coppie del passaggio precedente, lista del subset, dimensione totale del trascrittoma

```
python ~/Shared/PFAM_enrichment/test.py annotation_table.tsv Gene_Names_Fully_spanned_genes $(grep -c "^>" P.canaliculata_peptideset.fa)
```
crea un file con la parte iniziale del nome uguale a quella della lista subset

## 4. Associazione nomi a pfamid

estrazione delle corrispondenze dalla tabella del primo passaggio

```
sed 's/\s\{2,\}/\t/g' domains.tbl | grep -v ^# | cut -f3,4 | sed 's/ /\t/g' | cut -f2,3 | sort | uniq | awk '{print $2 "\t" $1}' > dict
```
join delle tabelle

```
join -j 1 -o 1.1 1.2 1.3 1.4 2.2 <(sort <(tail -n +2 Gene_Names_Fully_spanned_genes_enrichment_test.tsv)) <(sort dict) | sed 's/ /\t/g'  
```

## 5. eventuale estrazione di sequenze dal fasta iniziale in base a pfam

attenzione che non ci devono essere linee vuote nel fasta (toglierle con sed)

e se ci sono spazi nei nomi del fasta sono da troncare (sempre sed)

```
grep PF18644.1 domains.tbl | sed 's/\s\{2,\}/\t/' | cut -f1 | sort | uniq | xargs -n 1 samtools faidx Achatina_immaculata_peptideset.fa
```

