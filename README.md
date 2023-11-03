# GraphHiC

Codes for the manuscript "Improving Hi-C contact matrices using genome graphs". Here are the instructions for generating the graph-based Hi-C contact matrices and inferring sample genomes from genome graphs step by step.

## Genome Graph Construction

This step can be omitted if you already have a directed acyclic genome graph. To reproduce the genome graph created in our experiments, 

1. Download and Install vg-toolkit from https://github.com/vgteam/vg. (we use v1.49)
2. Download the .vcf file "RECOMB2024_invtrans_merged_sort.vcf.gz" containing structural variations of K-562 cancer cell line from https://kilthub.cmu.edu/articles/dataset/RECOMB2024_invtrans_merged_sort_vcf_gz/24481162.
3. Index the .vcf.gz file: `tabix -p vcf RECOMB2024_invtrans_merged_sort.vcf.gz`
4. Download hg19 human reference genome.
5. Construct the genome graph: `vg construct -S -f -r hg19.fa -v RECOMB2024_invtrans_merged_sort.vcf.gz > K562_hg19_RECOMB2024_invtrans.vg`
6. Graph index: `vg index -x K562_hg19_RECOMB2024_invtrans.xg -g K562_hg19_RECOMB2024_invtrans.gcsa -b ./tmp/ K562_hg19_RECOMB2024_invtrans.vg`
7. Separate the graph into different connected components, each component represents one chromosome, for example, for chromosome 1, `vg find -x K562_hg19_RECOMB2024_invtrans.xg -p chr1 -E > K562_hg19_RECOMB2024_invtrans_chr1.vg`

## Hi-C reads mapping

1. Download raw Hi-C reads, e.g. download raw Hi-C reads of sample SRR1658693 from Sequence Read Archive (SRA).
2. Map the reads: `vg map -f SRR1658693.1_1.fastq.gz -x K562_hg19_RECOMB2024_invtrans.xg -g K562_hg19_RECOMB2024_invtrans.gcsa -t 20 > SRR1658693.1_1.gam` and `vg map -f SRR1658693.1_2.fastq.gz -x K562_hg19_RECOMB2024_invtrans.xg -g K562_hg19_RECOMB2024_invtrans.gcsa -t 20 > SRR1658693.1_2.gam`
3. Chunk .gam file according to chromosomes, e.g. for chromosome 1, we can use the commands: `vg chunk -x K562_hg19_RECOMB2024_invtrans.xg -a SRR1658693.1_1.gam -C -p chr1` and `mv chunk_chr1.gam SRR1658693.1_1.chr1.gam`
4. Sort aligned reads in the .gam file according to the IDs, e.g. `vg convert -t 10 K562_hg19_RECOMB2024_invtrans.vg -G SRR1658693.1_1.chr1.gam | sed 's/^SRR1658693.1.//' | sort -nk1 -T ./tmp/ --parallel=10 | sed 's/^/SRR1658693.1./' | gzip > SRR1658693.1_1.chr1.idsort.gaf.gz`

## Graph pruning

1. Compile python module vg_pb2 according to https://github.com/cartoonist/pystream-protobuf/tree/master
2. Install the dependencies: `pip install -r requirements.txt`, we use python 3.7.0.
3. Run graph pruning script for each chromosome, e.g. for chromosome 21, `python graph_pruning.py --graph K562_hg19_RECOMB2024_invtrans_chr21.vg --f1 SRR1658693.1_1.chr21.idsort.gaf.gz --f2 SRR1658693.1_2.chr21.idsort.gaf.gz --prefix SRR1658693 --chr chr21`. All codes are in the folder codes/
4. This step will generate four kinds of files:

   (1) paired_SRR1658693.1_1.chr21.idsort.gaf.gz and paired_SRR1658693.1_2.chr21.idsort.gaf.gz, paired reads that both mapped onto chromosome 21.

   (2) SRR1658693_chr21_pair_out.txt, the information of mapped nodes in the pruned graph for each read pair. It has five columns. Column 1: read pair ID; column 2: the start mapped node of read end 1; column 3: the end mapped node of read end 1, if it is different from the value in column 2, then the read end is aligned across multiple nodes; column 4: the start mapped node of read end 2; column 5: the end mapped node of read end 2.

   (3) SRR1658693_chr21_final_graph.gfa, the graph information after pruning.

   (4) SRR1658693_chr21_M.npy, the graph-based contact matrix.

## Run heuristic dynamic programming algorithm

1. Unzip the file mu0_info.tar.gz, `tar -xf mu0_info.tar.gz`. It contains the necessary information for estimating the $\mu$ function.
2. Run DP algorithm for each chromosome, e.g. for chromosome 21, `python heuristic_DP.py --graph SRR1658693_chr21_final_graph.gfa --matrix SRR1658693_chr21_M.npy --prefix SRR1658693 --chr chr21 --threads 10`.
3. This step will generate three files:

   (1) SRR1658693_chr21_graph_TAD_boundaries.npy, the TAD boundaries called from heuristic DP.

   (2) SRR1658693_chr21_OPT.npy, the OPT values of heuristic DP, each node has one value.

   (3) SRR1658693_chr21_boundaries_raw.npy, the raw boundary information for each node. 

## Reconstruct sample genome

1. Once DP algorithm has been run on chromosomes, infer the genome via: `python genome_inference.py --reference hg19.fa --prefix SRR1658693 --chromosomes chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX`
2. This step will generate two kinds of files:

   (1) SRR1658693_chr*_reconstructed_nodes.npy, for each chromosome, a list of nodes representing the inferred path. 

   (2) SRR1658693_reconstruction.fa, the inferred genome under .fasta format. 
