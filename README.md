# Dynamic Landscape of Human L1 Transposition Revealed with Functional Data Analysis.
This resource provides the perl code used in paper below

-  "Dynamic Landscape of Human L1 Transposition Revealed with Functional Data Analysis". (2020) Molecular Biology and Evolution 37 (12), 3576-3600.  D Chen, MA Cremona, Z Qi et al. 

URL: https://pubmed.ncbi.nlm.nih.gov/32722770/

# Disclaimer
The pipeline has been tailored for the Washington University HTCF computing environment (https://htcf.wustl.edu/docs/) which uses the slurm queueing system (https://slurm.schedmd.com/tutorials.html). No guarantees are made about other systems, setups, configurations, etc.

# Main functions
The main functions are 
- 1) Align the LINE1 (L1) reads to the human genome to idenfiy the integration sites.
- 2) Cluster the L1 integration sites by a custom window size (default: 500 bp).
- 3) de-novo motif discovery of the L1 integration sites to identify potential binding proteins. 
- 4) Motif enrichment anlaysis of specified transcription factor(s).

# Instructions
These scripts are wrapped by a master perl script. Download all the scripts in one folder and run the warpper perl script 'L1_bar_wrapperV3.pl'.  

Usage: `perl L1_bar_wrapperV3.pl <read1.fq> <read2.fq> <barcode> <genome_aligner> <TF_to_scan>`

-  For `<genome_aligner>`: options can be combined without spaces in between (i.e. 12, 123 or 1234); Genome database is hg19.  
 - 1) bowtie2 for read2 single end alingment    
 - 2) bowtie2 for paired end alignment
 - 3) novoalign for read2 single end alignment
 - 4) novoalign for paired end alignment  

-  For `<TF_to_scan>`: underscore delimited transcription factors to scan;

Output name format: prefix is from the first "-"(hyphen) deliminated input name.
