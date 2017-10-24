#! /usr/bin/perl -w
#This script processes two paired reads based on index ID
#output: read1 and read 2 fastq files 
#Using "while (!eof(SEQ1) and !eof(SEQ2))" at line 26 to read two files

use strict;
use warnings;

if (@ARGV != 3) {die "input files: read1.fastq; read2.fastq; lineID_for_sorting";}
my $read1 = shift (@ARGV);
my $read2 = shift (@ARGV);
my $lineID_file = shift (@ARGV);

#make sure that the numbers of seq in read1 and 2 are the same 
my $total_R1counter = 0;
my $total_R2counter = 0;
open (SEQ1_check, "<", "$read1") || die 'read1 input error in checking\n';
open (SEQ2_check, "<", "$read2") || die 'read2 input error in checking\n';
while (my $line_id = <SEQ1_check>) { $total_R1counter ++;}
while (my $line_id = <SEQ2_check>) { $total_R2counter ++;}
if ($total_R1counter != $total_R2counter) { die 'the numbers of sequences in read1 and 2 are not the same'}
close(SEQ1_check);
close(SEQ2_check);

#fetch the sequence based on the line_id
open (SEQ1, "<", "$read1") || die 'read1 input error\n';
open (SEQ2, "<", "$read2") || die 'read2 input error\n';
open (IN, "<", "$lineID_file") || die 'lineID input error\n';
open (OUT1, ">", "$lineID_file"."read1.fq") || die 'output1 error\n';
open (OUT2, ">", "$lineID_file"."read2.fq") || die 'output2 error\n';

my $total_read_counter = 0;
my %read1 = ();
my %read2 = ();

#processing the input read1 and read2 as hash
while (!eof(SEQ1) and !eof(SEQ2)) {  #read same line from two files 
    $total_read_counter ++;
    chomp (my $line1 = <SEQ1>);
    chomp (my $line2 = <SEQ1>);
    chomp (my $line3 = <SEQ1>);
    chomp (my $line4 = <SEQ1>);
    chomp (my $lineA = <SEQ2>);
    chomp (my $lineB = <SEQ2>);
    chomp (my $lineC = <SEQ2>);
    chomp (my $lineD = <SEQ2>);
    $read1{$total_read_counter} = [$line1,$line2,$line3,$line4];
    $read2{$total_read_counter} = [$lineA,$lineB,$lineC,$lineD];
}

#read lineID and extract the reads based on the lineID
while (my $line_id = <IN>) {	
    chomp $line_id;
    if (defined $read1{$line_id} and $read2{$line_id}) {
    print OUT1 "${$read1{$line_id}}[0]\n${$read1{$line_id}}[1]\n${$read1{$line_id}}[2]\n${$read1{$line_id}}[3]\n";
    print OUT2 "${$read2{$line_id}}[0]\n${$read2{$line_id}}[1]\n${$read2{$line_id}}[2]\n${$read2{$line_id}}[3]\n"; 
    } else { die 'can not fetch the sequence based on the line_id\n'}
}

close (SEQ1);
close (SEQ2);
close (IN);
close (OUT1);
close (OUT2);

