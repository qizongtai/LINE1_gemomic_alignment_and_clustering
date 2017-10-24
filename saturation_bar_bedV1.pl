#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;

if (@ARGV != 4) {
    die 'input: bed.bar or bed.bar.uniqID file,
    points to calculate the saturations (recommend 50),
    length(bp) to merge at pA direction,
    length(bp) to merge at restriction-site direction';
}

#make hash out of bed.bar file: {$line_id} = $line_content 
my $line_id = 0;
my %lineID_seq = ();
print "Calculating staturation for the insertion sites ...\n";
open (IN, "<", "$ARGV[0]") or die 'input bed file error\n';
while (<IN>) {
    chomp;
    $lineID_seq{$line_id} = $_; # from "0" to "line_number -1", has to do this way because the code to generate the random numbers in line27
    $line_id ++; 
}

#1) generate the random numbers for each draw; 50 draws(points) as recommended in argument;
#2) make hash: {number of drawing} -> [seq1, seq2, seq3,....].  
my $points = $ARGV[1];
my %num_seq = ();
my @random_number =();
my $interval = sprintf("%.2f", $line_id/$points); #The rounded number may give a larger lineID than the actual lineID.
for (my $i=1; $i<=$points; $i++) {
    my $range = int($interval * $i);
    #open(NUM, ">", "$ARGV[0]-$i") or die "outfile error\n";
    for (my $j=1; $j<=$range; $j++) {
        my $num = int(rand($line_id));
        push (@random_number, $num); # a random number array between "0" and "$line_id - 1", exactly the range of keys ($line_id) in %lineID_seq
        #print NUM "$num\n";
    }
    foreach my $id (sort {$a <=> $b} @random_number) {
        push (@{$num_seq{$range}}, $lineID_seq{$id}); # number of drawing -> [seq1, seq2, seq3,....]  
    }
}

#calculate the unique insertions for each draw 
open (OUT, ">", "$ARGV[0].saturation") or die 'output saturation file error\n'; 
print OUT "number_of_read\tmerged_unique_insertions\n";
my ($chr, $start, $end, $id, $score, $strand, $barcode, @dump) = ();
my $last_chr; # default for this assignment is "undefined"
my $insertion_counter = 0;
my $length_pA = $ARGV[2];
my $length_restriction = $ARGV[3];
my (@array_start, @array_end, @array_strand, @array_barcode) = ();

foreach my $number (sort { $a <=> $b} keys %num_seq) { #for each draw
    foreach my $line (@{$num_seq{$number}}) { #get the sequence content
	chomp $line;
	($chr, $start, $end, $id, $score, $strand, $barcode, @dump) = split(/\s+/, $line);
	if (!$last_chr) { $last_chr = $chr; } # this argument is used only for the first time 
        elsif ($chr ne $last_chr) { # merge chr by chr; one chr at each time
	    &merge_unique_barcode_bed ($length_pA, $length_restriction);
	    (@array_start, @array_end, @array_strand, @array_barcode) = ();
	    $last_chr = $chr;
        }
        push (@array_start,$start);
	push (@array_end, $end);
	push (@array_strand, $strand);
	push (@array_barcode, $barcode);
    }
    &merge_unique_barcode_bed ($length_pA, $length_restriction); #run "&merge" again for the last chr
    print OUT  "$number\t$insertion_counter\n";
    $insertion_counter = 0;
    $last_chr = undef; #use the return value of undef() function to set a variable to undef; OR use "undef $last_chr";
    (@array_start, @array_end, @array_strand, @array_barcode) = ();
}    

close(OUT);

####################################################################################################################################
sub merge_unique_barcode_bed {
    my ($pA_len, $res_len) = @_;
    if (@array_start == 1 ) {
	$insertion_counter ++; 
    } elsif (@array_start > 1 ) { 
	$insertion_counter ++; # assign the first value, compare from the second value
	for (my $index=1; $index<@array_start; $index++) { 
	    if ($array_strand[$index] eq "+") { # for + strand, fix $start and merge $end
		if ((abs($array_start[$index] - $array_start[$index-1]) <= $res_len) and
		    (abs ($array_end[$index] - $array_end[$index-1]) <= $pA_len) and
		    ($array_barcode[$index] eq $array_barcode[$index-1]) and
		    ($array_strand[$index] eq $array_strand[$index-1]) ) {
		    next;
		} else {$insertion_counter ++;}
	    } elsif ($array_strand[$index] eq "-") { # for - strand, fix $end and merge $start
		if ((abs($array_start[$index] - $array_start[$index-1]) <= $pA_len) and
		    (abs($array_end[$index] - $array_end[$index-1]) <= $res_len) and
		    ($array_barcode[$index] eq $array_barcode[$index-1]) and
		    ($array_strand[$index] eq $array_strand[$index-1]) ) {
		    next;			
		} else { $insertion_counter ++;}
	    } else { print "strand (+ or -)  conditional error.\n";}
	}
    } else { die "\@array_start error.\n";}
}
