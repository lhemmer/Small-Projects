#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Lucas Hemmer
@date: 2/25/2025
@description

alignment_merge_by_index.py will combine two multiple sequence alignments.
Note: There must be AT LEAST ONE common sequence between both alignments that is named the same. There will be an error displayed if this is not true.
The multiple sequence alignments must be in a fasta file format, extensions (.fasta or .fa) should not a concern.

Inputs necessary for running:
    - \-a1 <filename> \--align1 <filename>
        first alignment file
    - \-a2 <filename> \--align2 <filename>
        second alignment file
	- \-o <filename>, --out <filename>
    	Prefix of output files
"""

################################
#### Import Python modules
################################

import argparse




################################
#### Functions
################################

def read_fasta_file(file_path):
    """Reads a FASTA file and returns its contents as a string."""
    with open(file_path, 'r') as file:
        return file.read()
    

def parse_fasta(fasta_str):
    """Parse a FASTA string into a dictionary of sequences."""
    sequences = {}
    lines = fasta_str.strip().split('\n')
    seq_id = None
    seq = []
    for line in lines:
        if line.startswith('>'):
            if seq_id is not None:
                sequences[seq_id] = ''.join(seq)
            seq_id = line[1:]
            seq = []
        else:
            seq.append(line)
    if seq_id is not None:
        sequences[seq_id] = ''.join(seq)
    return sequences


def align_sequences(alignment1, alignment2):
    """Align sequences in alignment1 not found in alignment2 by the position of the shared sequences."""
    # Parse the alignments
    alignment1_seqs = parse_fasta(alignment1)
    alignment2_seqs = parse_fasta(alignment2)

    # Find the shared sequence
    shared_seq_id = None
    for seq_id in alignment1_seqs:
        if seq_id in alignment2_seqs:
            shared_seq_id = seq_id
            break

    if shared_seq_id is None:
        raise ValueError("No shared sequence found between the two alignments.")

    # Align the sequences in alignment1 not found in alignment2
    aligned_sequences1 = {seq_id: list(seq) for seq_id, seq in alignment1_seqs.items()}
    aligned_sequences2 = {seq_id: list(seq) for seq_id, seq in alignment2_seqs.items()}

    shared_seq1 = list(alignment1_seqs[shared_seq_id])
    shared_seq2 = list(alignment2_seqs[shared_seq_id])


    # Adjust other sequences based on updated shared sequences
    i = 0
    while i < len(shared_seq1) or i < len(shared_seq2):
        if i >= len(shared_seq1):
            shared_seq1.append('-')
            for seq_id in aligned_sequences1:
                aligned_sequences1[seq_id].append('-')
        if i >= len(shared_seq2):
            shared_seq2.append('-')
            for seq_id in aligned_sequences2:
                aligned_sequences2[seq_id].append('-')
        if shared_seq1[i] == '-' and shared_seq2[i] != '-':
            shared_seq2.insert(i, '-')
            for seq_id in aligned_sequences2:
                aligned_sequences2[seq_id].insert(i, '-')
        elif shared_seq2[i] == '-' and shared_seq1[i] != '-':
            shared_seq1.insert(i, '-')
            for seq_id in aligned_sequences1:
                aligned_sequences1[seq_id].insert(i, '-')
        i += 1

    # Format the aligned sequences as a FASTA string
    aligned_fasta = []
    for seq_id, seq in aligned_sequences1.items():
        aligned_fasta.append(f">{seq_id}")
        aligned_fasta.append(''.join(seq))
    for seq_id, seq in aligned_sequences2.items():
        if seq_id not in aligned_sequences1:
            aligned_fasta.append(f">{seq_id}")
            aligned_fasta.append(''.join(seq))
    
    return '\n'.join(aligned_fasta)


def write_fasta_file(file_path, fasta_str):
    """Writes a FASTA string to a file."""
    with open(file_path, 'w') as file:
        file.write(fasta_str + '\n')



################################
#### Main
################################

if __name__ == "__main__":
    ## Parameters
    parser = argparse.ArgumentParser(prog='alignment_merge_by_index.py', description='This Program merges two multiple sequence alignments in fasta files that contain a common sequence')
    parser.add_argument('-a1', '--align1', help='The first alignment file in fasta format', required=True)
    parser.add_argument('-a2', '--align2', help='The second alignment file in fasta format', required=True)
    parser.add_argument('-o',  '--out', help='Output file name', required=True)

    ## Parameter check
    args = parser.parse_args()
    ## adding parameters
    align1 = args.align1
    align2 = args.align2
    outFile = args.out
    
    print ("Merging two alignment files")
    print ("Merging two alignment files")
    
    ## Read in Alignment files
    alignment1 = read_fasta_file(align1)
    alignment2 = read_fasta_file(align2)

    ## Merge the alignments
    resulting_alignment = align_sequences(alignment1, alignment2)
    
    ## Write output
    write_fasta_file(outFile, resulting_alignment)

    
    