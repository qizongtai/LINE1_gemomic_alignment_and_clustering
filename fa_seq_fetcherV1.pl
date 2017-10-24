#!/usr/bin/perl -w

use strict;
use warnings;

if (@ARGV != 3) {
    die "input: *.bed file that has the first 3 columns, length_extension_left(bp), length-extension_right(bp),";
}

my $bed = $ARGV[0];
my $length_extend_left = $ARGV[1];
my $length_extend_right = $ARGV[2];
my $fa_path = "/home/qizongtai/chrom_fasta/hg19_singlefile/hg19.fa";
#my $chrom_index = shift (@ARGV);
#my $fa_path = "/home/comp/rmlab/zqi/chrom-fasta/${chrom_index}-singlefile/${chrom_index}.fa"; #{} to delimiting a variable name; can be removed also.
#use above codes if chrom_index is an option in the command line argument.  

my $bed_extend = "$ARGV[0].$ARGV[1]-$ARGV[2]extend";

open (BED_IN, "<", "$bed") or die "bed file input error\n";
open (BED_OUT, ">", "$bed_extend") or die "bed extended file output error\n";
while (<BED_IN>) {
    chomp;    
    my @element = split(/\s+/, $_);
    my $first = $element[0];
    my $second = $element[1] - $length_extend_left;
    my $third = $element[2] + $length_extend_right;
    print BED_OUT "$first\t$second\t$third\n";
    #print OUT "$first\t$second\t$third\t$fourth\t$fifth\t$sixth\n"; #add strand direction (+/-) at the sixth column of the bed-extended file to use "-c" in bedtools getfasta 
}
system ("bedtools getfasta -fi $fa_path -bed $bed_extend -fo $bed_extend.seq -tab"); #this output single-line formated file.
#system ("bedtools getfasta -fi $fa_path -bed $name_input_bed -fo $name_output_seq -tab -c"); #output single-line formated file and consider strand direction(+/-).
#system ("bedtools getfasta -fi $fa_path -bed $name_input_bed -fo $name_output_seq -name"); #this output fasta formated file.

#Two types of error meassages:
#1)Feature (chr1:300-400) beyond the length of chr1 size (203 bp).  Skipping.
#2)WARNING. chromosome (chr3) was not found in the FASTA file. Skipping.
#need to try to merge the unique counts to the file
