#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV !=2) {die "input: newly converted bed_file, two-colunm seqID_barcode file";}

my $bed_file = shift @ARGV;
my $seqID_barcode = shift @ARGV;
my $bed_barcode = $bed_file . ".bar";
my $bed_no_barcode = $bed_file . ".nobar";

##read in two-colunm seqID_barcode file into a hash (seqID=>barcode)
print "Adding barcodes to bed_file... \n";
open (ID_BAR, "<", "$seqID_barcode") or die 'input seqID_barcode file error\n';
my %seqID_barcode = ();
while (<ID_BAR>) {
    chomp;
    my ($seqID, $chipID, $barcode, @dump) = split (/\s+/, $_); 
    $seqID =~ /^@(.*)$/; 
    $seqID_barcode{$1} = $barcode;
}
close (ID_BAR);

##read in six-column bed_file and get the barcode for each bed_seqID
open (BED, "<", "$bed_file") or die 'input bed file error\n';
open (BED_BAR, ">", "$bed_barcode") or die 'output bed_barcode file error\n';
open (BED_NO_BAR, ">", "$bed_no_barcode") or die 'output bed_no_barcode file error\n';
my $total_bed_counter = 0;
my $barcode_counter = 0;
my $no_barcode_counter = 0;
while (<BED>) {
    $total_bed_counter ++;
    chomp;
    my ($chr, $start, $end, $id, $score, $strand, @dump) = split (/\s+/, $_); 
    if (defined $seqID_barcode{$id}) {
        $barcode_counter ++;
        print BED_BAR "$chr\t$start\t$end\t$id\t$score\t$strand\t$seqID_barcode{$id}\n";
    } else {
        $no_barcode_counter ++;
        print BED_NO_BAR "$chr\t$start\t$end\t$id\t$score\t$strand\tNA\t$total_bed_counter\n";
    }
}
close (BED);
close (BED_BAR);
close (BED_NO_BAR);

my $percent = sprintf ("%.2f", 100*($barcode_counter/$total_bed_counter));
print "$percent% reads get the corresponding barcodes\n"; 
print "$no_barcode_counter reads do not get the corresponding barcodes\n"; 
