# Merging Multiple Sequence Alignments Script

This is a simple Python script that merges two multiple sequence alignments if there is a common sequence shared between them. This does not require the SeqIO library to read the fasta files. 


## Running the Script

Requires two alignment inputs in any order and a name of an output file. This is the script with the example files for test


```
python3 alignment_merge_by_index -a1 alignment1.fasta -a2 alignment2.fasta -o merged_alignment.fasta

python3 alignment_merge_by_index --align1 alignment1.fasta --align2 alignment2.fasta --out merged_alignment.fasta
```

I have two example alignments using a few transposable element sequences from [Hemmer et al 2023](https://www.biorxiv.org/content/10.1101/2022.11.25.518008v2) doi: (https://doi.org/10.1101/2022.11.25.518008)

```
python3 alignment_merge_by_index -a1 alignment_Jockey-3_denovo.fasta -a2 alignment_Jockey-3_Y_scaffold6.fasta -o alignment_Jockey-3_Y_scaffold6_denovo.fasta
```

