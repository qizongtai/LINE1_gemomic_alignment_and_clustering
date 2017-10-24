#!/bin/bash
#SBATCH -n 8  
#SBATCH -N 1 
#SBATCH --mem=100000
#SBATCH -o L1CCWT.out # standard out goes here
#SBATCH -e L1CCWT.err # standard error goes here
#SBATCH -J L1CCWT

#set parameters for sbatch.  
#-n = number of cores 
#-N = number of servers
#-mem = max memory per CPU in MB. Unit is M
#-o file to put standard out
#-e  file to put standard error
#-J  job name

#load software into environment
module load perl/5.20.3
module load bowtie2
module load bedtools
module load samtools
#commands at the terminal
perl L1_bar_wrapperV3.pl L1CCWT1to7mix-R1.fastq  L1CCWT1to7mix-R2.fastq WTbarcode-puro.txt 1 SP1
