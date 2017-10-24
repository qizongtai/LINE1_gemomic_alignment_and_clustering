#!/usr/bin/perl 
use strict;
use warnings;
use List::Util qw(min max);
use Data::Dumper;
$| = 1;

if (@ARGV != 3) {#Both ends are merged: polyA end + restriction-site end
    die "input files: bed.bar.uniqID.uniqPOS.stat,
    length(bp) to merge at pA direction,
    length(bp) to merge at restriction-site direction";
}

my $sort_barcode_bed = $ARGV[0] . "." . "barsort";
my $merge_bed = $ARGV[0] . "." .  $ARGV[1] . "-" . $ARGV[2] . "merge";
my $merge_overview = $ARGV[0] . "." . $ARGV[1] . "-" . $ARGV[2] . "merge.overview";
my $raw_gnashy =  $ARGV[0] . "." . $ARGV[1] . "-" . $ARGV[2] . "merge.rawgnashy";
my $gnashy = $ARGV[0] . "." . $ARGV[1] . "-" . $ARGV[2] . "merge.gnashy";
my $bed_graph = $ARGV[0] . "." . $ARGV[1] . "-" . $ARGV[2] . "merge.bedGraph";

unlink ($merge_bed) if (-e $merge_bed);
unlink ($merge_overview) if (-e $merge_overview);
unlink ($raw_gnashy) if (-e $raw_gnashy);
unlink ($gnashy) if (-e $gnashy);
unlink ($bed_graph) if (-e $bed_graph);

my $length_pA = $ARGV[1];
my $length_restriction = $ARGV[2];

#sort barcode
open (BED, "<", "$ARGV[0]") or die "input *.bed.bar.uniqID.uniqPOS.stat file error\n";
open (SORTBAR_BED, ">", "$sort_barcode_bed") or die "input *.bed.bar.uniqID.uniqPOS.stat file error\n";
my ($chr, $start, $end, $strand, $barcode, $count) = ();
my $last_chr = undef; #undefined variable used in below conditional
my $bed_read_counter = 0;
my @array_all = ();
while (<BED>) {
    chomp;
    $bed_read_counter ++;
    ($chr, $start, $end, $strand, $barcode, $count) = split(/\s+/, $_);
    if (!defined $last_chr) { $last_chr = $chr; } # this argument is used only for the first time 
    elsif ($chr ne $last_chr) { # merge chr by chr; one chr at each time
	#print Dumper \@array_all;
	@array_all = sort {$a->[0] <=> $b->[0] || $a->[3] cmp $b->[3]} @array_all;
	for (my $i=0; $i<@array_all; $i++) {
	    print SORTBAR_BED "$last_chr\t$array_all[$i][0]\t$array_all[$i][1]\t$array_all[$i][2]\t$array_all[$i][3]\t$array_all[$i][4]\n";
	}
	@array_all = ();
	$last_chr = $chr;
    }
    push (@array_all, [$start, $end, $strand, $barcode, $count]);
}
@array_all = sort {$a->[0] <=> $b->[0] || $a->[3] cmp $b->[3]} @array_all;
for (my $i=0; $i<@array_all; $i++) {
    print SORTBAR_BED "$last_chr\t$array_all[$i][0]\t$array_all[$i][1]\t$array_all[$i][2]\t$array_all[$i][3]\t$array_all[$i][4]\n";
}
close (BED);
close (SORTBAR_BED);

open (IN, "<", "$sort_barcode_bed") or die "input *.bed.bar.uniqID.uniqPOS.stat file error\n";
($chr, $start, $end, $strand, $barcode, $count) = ();
$last_chr = undef; #undefined variable used in below conditional
my $barsort_bed_read_counter = 0;
my (@array_start, @array_end, @array_strand, @array_barcode, @array_count)  = ();
while (<IN>) {
    chomp;
    $barsort_bed_read_counter ++;
    ($chr, $start, $end, $strand, $barcode, $count) = split(/\s+/, $_);
    if (!$last_chr) { $last_chr = $chr; } # this argument is used only for the first time 
    elsif ($chr ne $last_chr) { # merge chr by chr; one chr at each time
	&merge_unique_barcode_bed ($length_pA, $length_restriction);
	(@array_start, @array_end, @array_strand, @array_barcode, @array_count)  = ();
	$last_chr = $chr;
    }
    push (@array_start,$start);
    push (@array_end, $end);
    push (@array_strand, $strand);
    push (@array_barcode, $barcode);
    push (@array_count, $count);
}
&merge_unique_barcode_bed ($length_pA, $length_restriction); #run "&merge_by_id" again for the last chr
close IN;

if ($barsort_bed_read_counter != $bed_read_counter) { die "the read number of barcode-sorted bed file is different from the original.\n";}

open (OVERVIEW, ">", "$merge_overview") or die "output merge.overview file error\n";
print OVERVIEW "The number of reads in uniqueID, POS and barcode bedfile is $bed_read_counter.\n";
print OVERVIEW "After merged by $length_pA bp pA direction and $length_restriction at restriction direction, the \"NEW\" unique insertions were shown as below:\n";

system (" wc -l  $merge_bed >> $merge_overview ");
system (" cut -f 1 $merge_bed | sort | uniq -c | sort -nrk1 >> $merge_overview "); # "cut|sort|uniq"; "cut -f1" -> select for the first column; "-c" -> counts for each unique value. 
#system ("wc -l $bedfile_name > chr_list.txt");  #this outputs total counts of unique insertions.
#system ( "sort -u -k1,1 $bedfile > chr_list.txt" ); #"sort" outputs all columns although "-k1,1" sorts on first column. Need "cut -f1" to filter for just first column.

#transfer rawgnashy into GNASHY file 
system (" cut -c 4- $raw_gnashy | grep '^[0-9XY]' | grep -v '[0-9]_g' > $gnashy ");
system (" sed -i 's/X/23/g' $gnashy ");
system (" sed -i 's/Y/24/g' $gnashy ");

#append or prepend a header line to bedGraph file
my $name_mediator = 'name_mediator';
open (my $in, "<",  $bed_graph)  or die "can't read bed graph file: $!";
open (my $out, ">", "$name_mediator") or die "can't write bed graph: $!";
print $out "header\n"; 
while( <$in> ) {print $out $_;}
close ($out);
unlink($bed_graph);
rename("$name_mediator", $bed_graph);

print "Completed!\n";

########################################################################################################################################################
sub merge_unique_barcode_bed {
    my ($pA_len, $res_len) = @_;
    open (OUT, ">>", "$merge_bed") or die "output merge file error\n";
    open (RAW_GNASHY, ">>", "$raw_gnashy") or die "output raw.gnashy file error\n";
    open (BED_GRAPH, ">>", "$bed_graph") or die "output bedGraph file error\n";
    if (@array_start == 1 ) {
	print OUT "$last_chr\t$start\t$end\t$strand\t$barcode\t$count\t1\n";
	if ($strand eq "-" ) {
	    print RAW_GNASHY "$last_chr\t$start\t$count\n";
	    my $start_plus = $start + 1;
	    print BED_GRAPH "$last_chr\t$start\t$start_plus\t1\n";
	} elsif ($strand eq "+"){
	    print RAW_GNASHY "$last_chr\t$end\t$count\n";
	    my $start_plus = $start + 1;
	    print BED_GRAPH "$last_chr\t$end\t$start_plus\t1\n";
	} else {
	    die "strand (+ or -) error \"1\" in generating RAW_GNASHY file.\n";
	} 	
    } elsif (@array_start > 1 ) { # assign the first value, compare from the second value
	my $merge_counter = 1;
	push (my @start_list_temp, $array_start[0]);
	push (my @end_list_temp, $array_end[0]);
	my $unique_counter = $array_count[0];
	my $start_print = 0;
	my $end_print = 0;
	for (my $index=1; $index<@array_start; $index++) { # compare from the second value
	    if ((abs($array_start[$index] - $array_start[$index-1]) <= $res_len) and
		(abs($array_end[$index] - $array_end[$index-1]) <= $pA_len) and
		($array_barcode[$index] eq $array_barcode[$index-1]) and
		($array_strand[$index] eq $array_strand[$index-1]) ) {
	    $merge_counter ++;
	    push (@start_list_temp, $array_start[$index]);
	    push (@end_list_temp, $array_end[$index]);
	    $unique_counter += $array_count[$index];
	    } else {
		#need to print the previous one "$index-1"
		$start_print = min(@start_list_temp);
		$end_print = max(@end_list_temp);
		print OUT "$last_chr\t$start_print\t$end_print\t$array_strand[$index-1]\t$array_barcode[$index-1]\t$unique_counter\t$merge_counter\n";
		if ($array_strand[$index-1] eq "+") { # for + strand, fix $start and merge $end
		    #print RAW_GNASHY "$last_chr\t$start_print\t$unique_counter\n"; # $start_print means closet site to the cutters direction
		    print RAW_GNASHY "$last_chr\t$end_print\t$unique_counter\n"; # $end_print means closest sites to the pA direction
		} elsif ($array_strand[$index-1] eq "-") { # for - strand, fix $end and merge $start
		    #print RAW_GNASHY "$last_chr\t$end_print\t$unique_counter\n"; # $end_print means closet site to the cutters direction
		    print RAW_GNASHY "$last_chr\t$start_print\t$unique_counter\n"; # $start_print means closest sites to the pA direction
		} else { die "strand (+ or -) error in generating RAW_GNASHY file.\n";}
		my $start_plus = $start + 1;
		print BED_GRAPH "$last_chr\t$start\t$start_plus\t$merge_counter\n";
		#initalize the values for the next comparison as line92-97
		$unique_counter = 0;
		$start_print = 0;
		$end_print = 0;
		$merge_counter = 1;
		@start_list_temp = ();
		@end_list_temp = ();
		#process the "current index" value for the next comparison  
		push (@start_list_temp, $array_start[$index]);
		push (@end_list_temp, $array_end[$index]);
		$unique_counter = $array_count[$index];
	    }
	    if ($index == $#array_start) { # deal with the last element in the loop
		$start_print = min(@start_list_temp);
		$end_print = max(@end_list_temp);	    
		print OUT "$last_chr\t$start_print\t$end_print\t$array_strand[$index]\t$array_barcode[$index-1]\t$unique_counter\t$merge_counter\n";
		if ($array_strand[$index] eq "-" ) {
		    print RAW_GNASHY "$last_chr\t$start_print\t$unique_counter\n";
		    my $start_plus = $start + 1;
		    print BED_GRAPH "$last_chr\t$start\t$start_plus\t$merge_counter\n";
		} elsif ($array_strand[$index] eq "+"){
		    print RAW_GNASHY "$last_chr\t$end_print\t$unique_counter\n";
		    my $start_plus = $start + 1;
		    print BED_GRAPH "$last_chr\t$start\t$start_plus\t$merge_counter\n";
		} else {
		    die "strand (+ or -) error in generating RAW_GNASHY file.\n";
		} 
	    }
	}
    } else {
	die "\@array_start error.\n";
    }
    close OUT;
    close RAW_GNASHY;    
}
    

