# GraphHiC

Codes for the manuscript "Improved Hi-C contact matrices via genome graphs: Hardness and Heuristics". Here are the instructions for generating the graph-based Hi-C contact matrices and inferring sample genomes from genome graphs step by step.

## Genome Graph Construction

This step can be omitted if you already have a directed acyclic genome graph. To reproduce the genome graph created in our experiments, 

1. Download the .vcf file "RECOMB2024_invtrans_merged_sort.vcf.gz" containing structural variations of K-562 cancer cell line from https://kilthub.cmu.edu (data is under-review).
2. Index the .vcf.gz file: `tabix -p vcf RECOMB2024_invtrans_merged_sort.vcf.gz`
3. 
