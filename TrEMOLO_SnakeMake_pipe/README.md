# TrEMOLO Snakemake Pipeline
 
 TrEMOLO is shortened for Transposable Elements MOvement detection using LOng reads.
 
 The original program can be found at the following Github page. 
  
 https://github.com/DrosophilaGenomeEvolution/TrEMOLO 
 
 Citation for TrEMOLO can be found here.
 
 https://www.mdpi.com/2073-4409/9/8/1776
 
 
 I was having a difficult time running the original program so I adopted it into a a Snakemake pipeline.
 
 Please refer to the original paper and program for use of TrEMOLO for copyright purposes.
 
 ## Requirements
 
 Markup : * Snakemake
          * Blast 2.2+
          * Bedtools v2
          *Assemblyitics or RaGOO ouput
          
## Running the script

The tool requires Bedtools and BLAST to be accessible in the path. You need a TE library fasta file to be formated into a BLAST database. This can be done with *makeblastdb* with BLAST tools prior to running. The tool also requires an Assemblytics BED file created through Assemblytics or RaGOO. 

The pipeline requires changes to the *config.yaml* file to specify the sample names to access the FASTA files and Assemblytics BED file, reference FASTA file, and BLAST TE library database (which needs to be indexed).

`
SAMPLES:
    - sample 
    # - sample2 can create a bullet list for multiple samples

refFasta:
    - reference_genome.fasta
     
altdir:
    - /other/directory/ #in case you files are found elsewhere

TEdatabase:
    - specieslib_db

TEfastaIndex:
    - specieslib_db.fai
`

Once these are set up in the *config.yaml* file, run snakemake in the directory and the file should execute.

`snakemake`


## Output

Markup : * Two BED files for the deletions and insertions in the sample relative to the reference genome
         * Two FASTA files with the sequences corresponding to the deletions and insertions
         * Two BLASTn output files for insertions and deltions compared to the TE database
         * A TE sizes file for running the script
         * A tabulated results file with the BLASTn result according to the percentage threshold for identifying TEs

Each line of the tabulated file contains the following

Markup : * The transposon identified by BLAST
         * The location of the insertion relative to the reference genome
         * The Percld, or ID shared of the insertion / deletion shared with the transposon 
         * FragSize, the size of the TE fragment 
         * RefSize, size of the insertion or deletion
         * PercTotal, percentage of the insertion or deletion occupied by the TE size fragment

Note: This will probably need additional filtering. BLAST identifies several TEs within each insertion identified. Most of these can be eliminated with a PercID threshold but the size of the fragment compared to the total fragment size could also be a filtering factor as well. 
          

