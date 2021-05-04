# Dynamic Landscape of Human L1 Transposition Revealed with Functional Data Analysis.
This resource provides the perl code used in "Dynamic Landscape of Human L1 Transposition Revealed with Functional Data Analysis". Molecular Biology and Evolution 37 (12), 3576-3600.  D Chen, MA Cremona, Z Qi et al. (2020) 

These scripts are used to align LINE1 (L1) insertions back to the human genome and cluster nearby insertions.  
Download all the scripts in one folder and run only the warpper perl script (L1_bar_wrapperV3.pl).  
Usage: <read1.fq> <read2.fq> <barcode> <mapper option> <TF_to_scan> <Mapper option> <TF_to_scan>
For <Mapper option>: options can be combined without spaces in between (i.e. 12, 123 or 1234); Mapping database is hg19.  
 1) bowtie2 for read2 single end alingment    
 2) bowtie2 for paired end alignment
 3) novoalign for read2 single end alignment
 4) novoalign for paired end alignment  
For <TF_to_scan>: underscore delimited transcription factors to scan;  
Output name format: prefix is from the first "-"(hyphen) deliminated input name.
