# Phylogenetic Analysis of CSP haplotypes against NCBI CSP Sequences

Source: NCBI PfCSP Records (n=2564): https://www.ncbi.nlm.nih.gov/nuccore/?term=circumsporozoite%20protein%20plasmodium%20falciparum


## Alignment:

1. Align nucleotides of the 2564 NCBI samples against the 3D7 CSP reference, trimming to CSP length, save as `PfCSP_aligned_n.fasta`
```sh
# mafft --add NCBI_PfCSP_sequences.fasta --keeplength --thread 24 csp_reference.fasta > PfCSP_aligned_n.fasta
mafft --add NCBI_PfCSP_sequences_filtered.fasta --keeplength csp_reference.fasta > PfCSP_aligned_n.fasta

```

2. Translate to amino acids (using AliView) and save as `PfCSP_aligned_aa.fasta`
3. Add CSP Haplotypes to align NCBI AA samples
```sh
mafft --add csp_HaplotypeSeq_NGmergeSNP_filtered_translated.fasta --keeplength PfCSP_aligned_aa.fasta > PfCSP_aligned_aa_wHaplotypes.fasta
```

## Tree Creation

Use 