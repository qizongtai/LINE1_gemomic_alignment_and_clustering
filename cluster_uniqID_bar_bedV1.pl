#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);
use Data::Dumper;
$| = 1;

if (@ARGV != 3) {die "input files: 'bed.bar' or 'bed.bar.uniqID', length(bp) to cluster, minimal read number to call an unique barcode (insertion)"; }

my $cluster_bed = $ARGV[0] . "." . $ARGV[1] . "cluster";
my $cluster_overview = $ARGV[0] . "." . $ARGV[1] . "cluster.overview";	
my $raw_gnashy =  $ARGV[0] . "." . $ARGV[1] . "cluster.rawgnashy";
my $gnashy = $ARGV[0] . "." . $ARGV[1] . "cluster.gnashy";
my $more_1barcode = $ARGV[0] . "." . $ARGV[1] . "cluster.more1bar";
my $more_1barcode_more_reads = $ARGV[0] . "." . $ARGV[1] . "cluster.more1bar" . "." . "$ARGV[2]reads";
unlink ($cluster_bed) if (-e $cluster_bed);
unlink ($cluster_overview) if (-e $cluster_overview);
unlink ($raw_gnashy) if (-e $raw_gnashy);
unlink ($gnashy) if (-e $gnashy);
unlink ($more_1barcode) if (-e $more_1barcode);
unlink ($more_1barcode_more_reads) if (-e $more_1barcode_more_reads);

my $length = $ARGV[1];
my $read_cutoff = $ARGV[2];
my $total_read_counter = 0;
my $cluster_counter = 0;
my $more_1barcode_counter = 0;
my $more_1barcode_more_reads_counter = 0;

open (IN, "<", "$ARGV[0]") or die "input bed.uniquePOS.stat.merge file error\n";
open (OVERVIEW, ">", "$cluster_overview") or die "output cluster.overview file error\n";
my ($chr, $start, $end, $id, $score, $strand, $barcode, @dump) = ();
my (@array_start, @array_end, @array_strand, @array_barcode)  = ();
my $last_chr; # undef this variable
while (<IN>) {
	chomp;
	$total_read_counter ++;
	($chr, $start, $end, $id, $score, $strand, $barcode, @dump) = split(/\s+/,$_);
	if (!$last_chr) { $last_chr = $chr; } # this argument is used only for the first time 
        elsif ($chr ne $last_chr) { # cluster chr by chr; one chr at each time
	    &cluster_bed_unique ($length);
	    (@array_start, @array_end, @array_strand, @array_barcode)  = ();
	    $last_chr = $chr;
        }
        push(@array_start,$start); 
	push(@array_end, $end);
	push(@array_strand, $strand);
	push(@array_barcode, $barcode);
}
&cluster_bed_unique ($length); #run "&cluster_by_id" again for the last chr
close IN;

my $percent_morebar = sprintf ("%.2f", 100*($more_1barcode_counter/$cluster_counter));
my $percent_onebar = 100 - $percent_morebar;
my $onebar_counter = $cluster_counter - $more_1barcode_counter;
my $percent_morebar_cutoff = sprintf ("%.2f", 100*($more_1barcode_more_reads_counter/$cluster_counter));
my $percent_onebar_cutoff = 100 - $percent_morebar_cutoff;
my $onebar_cutoff_counter = $cluster_counter - $more_1barcode_more_reads_counter;
print OVERVIEW "The total number of reads in the initial merged unique bedfile is $total_read_counter.\n";
print OVERVIEW "Clustered by $length bp at both ends, $cluster_counter clusters were identified.\n";
print OVERVIEW "Out of these clusters, $percent_onebar % ($onebar_counter) has one barcode and $percent_morebar % ($more_1barcode_counter) has more than one barcode.\n";
print OVERVIEW "If requiring at least $read_cutoff reads to call an unique barcode, ";
print OVERVIEW "$percent_onebar_cutoff % ($onebar_cutoff_counter) has one barcode and $percent_morebar_cutoff % ($more_1barcode_more_reads_counter) has more than one barcode.\n";
print OVERVIEW "These clusters were shown as below:\n";
close(OVERVIEW);

system (" wc -l  $cluster_bed >> $cluster_overview ");
system (" cut -f 1 $cluster_bed | sort | uniq -c | sort -nrk1 >> $cluster_overview "); # "cut|sort|uniq"; "cut -f1" -> select for the first column; "-c" -> counts for each unique value. 
#system ("wc -l $bedfile_name > chr_list.txt");  #this outputs total counts of unique insertions.
#system ( "sort -u -k1,1 $bedfile > chr_list.txt" ); #"sort" outputs all columns although "-k1,1" sorts on first column. Need "cut -f1" to filter for just first column.

system (" cut -c 4- $raw_gnashy | grep '^[0-9XY]' | grep -v '[0-9]_g' > $gnashy ");
system (" sed -i 's/X/23/g' $gnashy ");
system (" sed -i 's/Y/24/g' $gnashy ");

print "$total_read_counter\t$cluster_counter\t$more_1barcode_counter\t$more_1barcode_more_reads_counter\n";

###############################################################################################################################################
sub cluster_bed_unique {
    my ($len) = @_;
    open (OUT, ">>", "$cluster_bed") or die "output cluster file error\n";
    open (GNASHY, ">>", "$raw_gnashy") or die "output raw.gnashy file error\n";
    open (MORE1BAR, ">>", "$more_1barcode") or die "output more_1barcode file error\n";
    open (MORE1BAR_MOREREADS, ">>", "$more_1barcode_more_reads") or die "output more_1barcode_more_reads file error\n";  
    if (@array_start == 1 ) { #if there is only one read in chr, no comparison is needed.
	$cluster_counter ++;
	print OUT "$last_chr\t$start\t$end\t$strand\t1\n";
	print GNASHY "$last_chr\t$start\t1\n";
    } elsif (@array_start > 1 ) { 
	my $start_print = 0;
	my $end_print = 0;
	my $uniq_barcode_counter = 0;
	my %barcode_read = ();
	my $read_counter = 0; #count merging times (number of reads) for each cluster, in the output;
	#the first value was assigned
	push (my @start_list_temp, $array_start[0]);
	push (my @end_list_temp, $array_end[0]);
	$barcode_read{$array_barcode[0]} = 1;
	$read_counter = 1;
	for (my $index=1; $index<@array_start; $index++) { # start from the second value; 
	    if ((abs($array_start[$index] - $array_start[$index-1]) <= $len) and
		(abs($array_end[$index] - $array_end[$index-1]) <= $len) ) {
		push (@start_list_temp, $array_start[$index]);
		push (@end_list_temp, $array_end[$index]);
		$barcode_read{$array_barcode[$index]} ++; 
		$read_counter ++; 
	    } else {
		$cluster_counter ++;
		$start_print = min(@start_list_temp);
		$end_print = max(@end_list_temp);
		$uniq_barcode_counter = keys %barcode_read;
		#output all clusters
		print OUT "$last_chr\t$start_print\t$end_print\t$array_strand[$index-1]\t$read_counter\t$uniq_barcode_counter\n"; #need to print the previous one "$index-1"
		print GNASHY "$last_chr\t$start_print\t$uniq_barcode_counter\n";
		#output clusters with more than one barcode
		if ($uniq_barcode_counter > 1) {
			$more_1barcode_counter ++;
			print MORE1BAR "$last_chr\t$start_print\t$end_print\t$array_strand[$index-1]\t$read_counter\t$uniq_barcode_counter";
			foreach my $uniq_barcode (sort {$a cmp $b} keys %barcode_read) {
				print MORE1BAR "\t$uniq_barcode\t$barcode_read{$uniq_barcode}"
			}
			print MORE1BAR "\n";
		}
		#output clusters with more than one barcode; the barcodes was counted when at least read_cutoff reads exist.
		my %barcode_morereads = (); #set up a new hash: {barcode} -> numbe of reads (barcode has more than read_cutoff reads)
		foreach my $uniq_barcode (sort {$a cmp $b} keys %barcode_read) { 
			if ($barcode_read{$uniq_barcode} > $read_cutoff) {
				$barcode_morereads{$uniq_barcode} = $barcode_read{$uniq_barcode};
			}
		}
		my $uniq_barcode_morereads_counter = keys %barcode_morereads; #the same as line105-113
		if ($uniq_barcode_morereads_counter > 1) {
			$more_1barcode_more_reads_counter ++;
			print MORE1BAR_MOREREADS "$last_chr\t$start_print\t$end_print\t$array_strand[$index-1]\t$read_counter\t$uniq_barcode_morereads_counter";
			foreach my $uniq_barcode (sort {$a cmp $b} keys %barcode_morereads) {
				print MORE1BAR_MOREREADS "\t$uniq_barcode\t$barcode_morereads{$uniq_barcode}"
			}
			print MORE1BAR_MOREREADS "\n";			
		}		
		#initalize the the values for the next comparison as line 75-82
		$start_print = 0;
		$end_print = 0;
		$read_counter = 0;
		$uniq_barcode_counter = 0;
		$uniq_barcode_morereads_counter = 0;
		@start_list_temp = ();
		@end_list_temp = ();
		%barcode_read = ();
		%barcode_morereads = ();
		#assign the "current index" for the next comparison
		push (@start_list_temp, $array_start[$index]);
		push (@end_list_temp, $array_end[$index]);
		$barcode_read{$array_barcode[$index]} = 1;
		$read_counter = 1
	    }
	    if ($index == $#array_start) { # deal with the last element in the loop
		$cluster_counter ++; #same as line92-129
		$start_print = min(@start_list_temp); 
		$end_print = max(@end_list_temp);
		$uniq_barcode_counter = keys %barcode_read;
		#output all clusters
		print OUT "$last_chr\t$start_print\t$end_print\t$array_strand[$index]\t$read_counter\t$uniq_barcode_counter\n";
		print GNASHY "$last_chr\t$start_print\t$uniq_barcode_counter\n";
		#output clusters with more than one barcode
		if ($uniq_barcode_counter > 1) {
			$more_1barcode_counter ++;
			print MORE1BAR "$last_chr\t$start_print\t$end_print\t$array_strand[$index]\t$read_counter\t$uniq_barcode_counter";
			foreach my $uniq_barcode (sort {$a cmp $b} keys %barcode_read) {
				print MORE1BAR "\t$uniq_barcode\t$barcode_read{$uniq_barcode}"
			}
			print MORE1BAR "\n";
		}	
		#output clusters with more than one barcode; the barcodes was counted when at least read_cutoff reads exist.
		my %barcode_morereads = (); #set up a new hash: {barcode} -> numbe of reads (barcode has more than read_cutoff reads)
		foreach my $uniq_barcode (sort {$a cmp $b} keys %barcode_read) {
			if ($barcode_read{$uniq_barcode} > $read_cutoff) {
				$barcode_morereads{$uniq_barcode} = $barcode_read{$uniq_barcode};
			}
		}
		my $uniq_barcode_morereads_counter = keys %barcode_morereads; #the same as line105-113
		if ($uniq_barcode_morereads_counter > 1) {
			$more_1barcode_more_reads_counter ++;
			print MORE1BAR_MOREREADS "$last_chr\t$start_print\t$end_print\t$array_strand[$index]\t$read_counter\t$uniq_barcode_morereads_counter";
			foreach my $uniq_barcode (sort {$a cmp $b} keys %barcode_morereads) {
				print MORE1BAR_MOREREADS "\t$uniq_barcode\t$barcode_morereads{$uniq_barcode}"
			}
			print MORE1BAR_MOREREADS "\n";			
		}
	    }
	}
    } else {
	die "\@array_start error.\n";
    }
    close OUT;
    close GNASHY;
}
    

