# Advanced-Programming-Project

The software analyzes mitochondrial genomes. 
Some of the analysis includes:
1. Calculating GC content and sequence statistics;
2. Identifying conserved motifs across sequences;
3. Performing pairwise sequence alignments;
4. Generating and analyzing mitochondrial DNA-derived data, such as motifs and alignment results.

The project includes a friendly web interface.

On the Home page, you can upload the file that gets processes using the class mtDNA_parser

On the next page, you can visualize the names of the genomes present in the file and choose the number of genomes to be analyzed, from 1, 2 or all.

On the 1 genome analysis page, you can select the name of the genome and the type of analysis:
1. GC contennt and length
2. Extract a subsequence(by also imputing the start and stop indexes) and calculate it's GC content and length
3. Search for a motif(input needed) and find out the position index in the genome and occurance count

On the 2 genome analysis page, you can select the 2 genomes and the type of analysis:
1. GC content and length comparison
2. Search for a motif(input needed) and find out the position index in the genome and occurance count
3. Perform pairwise analysis

On the all genome analysis page, you can select the type of analysis:
1. GC content and length comparison
2. Search for a motif(input needed) and find out the position index in the genome and occurance count
