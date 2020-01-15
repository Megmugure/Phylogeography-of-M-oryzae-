# Exploring-the-phylogenetic-evolution-and-geographic-transmission-patterns-of-rice-blast-M.-oryzae-

This is my MSc project repository. I will be documenting all my work here including the data, scripts and tools I use moving forward. My plan for the coming months will also be highlighted in here.

I first run fastqc on my data to assess the quality, I used fastqc/0.11.7. I noticed that there were 8 isolates missing from the 50 initially assumed to be provided.
I then aligned the reads to a reference(MG8) and I used bowtie2/2.3.4.1 and samtools/1.8
Variant calling I used bcftools/1.8 and igvtools/2.3.98
For filtering the vcf files I used vcftools/0.1.15
I created a phylogenetic tree using fasttree/2.1.10













