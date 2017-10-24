#!/usr/bin/perl #need to change $pwm_filepath to screen more TFs
use strict;
use warnings;
use Data::Dumper;
my $starttime = (time);

my $usage = "script.pl -fa <fasta.formatted.sequence.to.search> OR -txt <sigleline.formatted.sequence.to.search> -tf <tf.name>\n";
my $pwm_filepath = "all.pwms.frequencies"; #my $pwm_filepath = "sp1.pwm";
my $seq_input;
my $tf_query;
my $file_flag;
my $score_cutoff;

if (@ARGV){
    ($seq_input, $tf_query, $file_flag, $score_cutoff) = &processarguments(@ARGV);
} else { &noarguments; }

#reads in a previously formatted position weight matrices file as a hash
open (IN, "<", "$pwm_filepath") || die "cannot find pwm file: $pwm_filepath.\n";
my %tfs = (); #tf_name -> 1 (a random number to take up the value space in a hash) 
my %pwm = (); #tf_name_base(A, C, G or T)-> pwm_values (score for each pos as an array)
while (<IN>) {
    chomp;
    my @line = split(/\|/, $_);
    my $tf_name_base = $line[0]; #print $key."\n"; #exit;
    my $length = length ($tf_name_base);
    my $tf_name = substr($tf_name_base, 0, ($length -2)); #print $substr."\n";
    if (!$tfs{$tf_name}) { $tfs{$tf_name} = 1; }
    my @values = split(/\s+/, $line[1]);
    $pwm{$tf_name_base} = [@values]; 
}
close IN;

#scan the sequences by pwm
open (IN, "<", "$seq_input") || die "couldn't open input file\n$usage\n";
open (RAW, ">", "$seq_input.$tf_query.logscores") || die "cannot open raw output file\n";
print RAW "Chromosome\tStart\tEnd\tTF\tTF_Start\tTF_End\tLogScore\tPlotscore\tStrand\n";
open (COUNT, ">", "$seq_input.$tf_query.tfbs") || die "cannot open tf-PWM-counter file\n";
print COUNT "Chromosome\tStart\tEnd\tTF\tCount\n";
open (WEIGHT, ">", "$seq_input.$tf_query.logscores.weighted") || die "cannot open weighted output file\n";
print WEIGHT "Chromosome\tStart\tEnd\tTF\tTF_Start\tStrand\tWeighted_Score\n";
open (OVERVIEW, ">", "$seq_input.$tf_query.overview") || die "cannot open overview output file\n";

my @tf_name_database = sort{ $a cmp $b } keys (%tfs); #print Dumper \@tf_name_database; print "$tf_query\n";
my $chr = '';
my $seq = '';
my $seq_start = 0; # the start coordiante of a sequence in the genome
my $seq_end = 0; # the end coordiante of a sequence in the genome
my %backhash = (A => 0.29, G => 0.21, T => 0.29, C => 0.21);
my $denom = log(2); 
my $strand;
my $found_seq_counter = 0;
my $unfound_seq_counter = 0;
my $seq_counter = 0;
my $sum_tf_bs_counter = 0; #summed counts of potential tf binding sites.
my $sum_pos_seq_counter = 0; #summed counts of sequence having at least one potential tf binding site.

#determine the file type and call the subroutine to scan the sequences
if ($file_flag == 1) { #for fasta file, the headers don't have $seq_start and $seq_end cornidates)
    while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^>/) {
	    $seq_counter ++;
	    if (!$chr) { $chr = substr($line, 1);} #only executed for the first time
	    &scan; #do sub for the previous sequence
	    $chr = substr($line, 1);
	    $seq = ''; 	   		   
	} else {$seq = $seq.$line;}
	if (eof(IN)) { &scan; }		
    }
} elsif ($file_flag == 2) {#for single-line formatted file, the headers have start and end cornidates) 
    while (my $line = <IN>) {
	chomp $line;
	$seq_counter ++;
	if ($line =~ /^chr/) {
	    $found_seq_counter ++;
	    ($chr, $seq_start, $seq_end, $seq) = split (/\t+|-{1}|:{1}/, $line);
	    &scan; 
	    $chr = '';
	    $seq = '';
	    $seq_start = 0;
	    $seq_end = 0;
	} else { $unfound_seq_counter ++; }
	
    }	
} else { die "Error: format of the sequence file to search is not available in this perl script"; }

if ($file_flag == 1) {
    my $percent_pos_seq = sprintf("%.3f", 100* $sum_pos_seq_counter/$seq_counter);
    print OVERVIEW "There are $seq_counter sequences.\n";
    print OVERVIEW "There are $sum_tf_bs_counter specified $tf_query binding sites and $sum_pos_seq_counter ($percent_pos_seq %) sequences that has at least one $tf_query binding site.\n\n";
    print "$seq_counter\t$sum_pos_seq_counter\t$sum_tf_bs_counter\n"; 
} elsif ($file_flag == 2){
    my $percent_pos_seq = sprintf ("%.3f", 100* ($sum_pos_seq_counter/$seq_counter));
    my $percent_foundseq = sprintf ("%.3f", 100* ($found_seq_counter/$seq_counter)); #$seq_counter=$unfound_seq_counter+$found_seq_counter
    print OVERVIEW "There are $seq_counter sequences, out of these:\n";
    print OVERVIEW "$found_seq_counter ($percent_foundseq %) positions fetch their corresponding sequences and have them scanned.\n";
    print OVERVIEW "$unfound_seq_counter positions don't fetch their corresponding sequences.\n\n";
    print OVERVIEW "There are $sum_tf_bs_counter specified $tf_query binding sites and $sum_pos_seq_counter ($percent_pos_seq %) sequences that has at lease one $tf_query binding site.\n\n";
    print "$seq_counter\t$sum_pos_seq_counter\t$sum_tf_bs_counter\n";
}

my $stoptime = (time);
my $elapsedtime = $stoptime - $starttime;
print OVERVIEW "This program took $elapsedtime second to finish\n";

close IN; close RAW; close WEIGHT; close COUNT; close OVERVIEW;

#########################################################################################################################################################
#This subroutine processes the arguments to find the few transcription factors to test#
sub processarguments {
    my @args = @_;
    my $args_counter = 0;
    foreach my $argument (@args){
	if ($argument =~ /-fa/) { ##fasta formated file
	    $seq_input = $args[$args_counter + 1]; #get the variable after '-fa'
	    $file_flag = 1;
	} elsif ($argument =~ /-txt/) { ##single-line formated file
	    $seq_input = $args[$args_counter + 1]; #get the variable after '-txt'
	    $file_flag = 2;
	} elsif ($argument =~ /-tf/) {
	    $tf_query = $args[$args_counter + 1];  #get the variable after '-tf'
	} elsif ($argument =~ /-cf/) {
	    $score_cutoff = $args[$args_counter + 1];#get the variable after '-cf'
	}
	$args_counter ++;
    }
    return ($seq_input, $tf_query, $file_flag, $score_cutoff);
}

sub noarguments {
    print "Error: no arguments\n";
    print "Argument options include:\n";
    print "-fa <fasta formatted sequence file to search>\n";
    print "OR\n";
    print "-txt <sigleline formatted sequence file to search>\n";
    print "-tf <underscore separated list of transcription factors to look for>\n";
    print $usage."\n";
    die;
}

sub scan { #scan one sequence at a time using a window corresponding to the length of each pwm
    my $length = length ($seq); 
    my $tf_bs_counter = 0; #summed counts of a single pwm for each sequence 
    my $total_tf_bs_counter = 0; #summed counts of all possible pwm for each sequence
    my $match = 0;
    my $undefined_base = 0;
    for (my $e=0; $e<@tf_name_database; $e++) {
	if ($tf_name_database[$e] =~ /$tf_query/i) {
	    $match ++;
	    my @sites;
	    my @strands;
	    my $outname = $tf_name_database[$e];
	    my $key = $tf_name_database[$e]."_A"; #print "$key\n";
	    #unless (defined @{$pwm{$key}}) { print ERRORS "$key is undefined\n"; }
	    my $motif_length = @{$pwm{$key}} ; #print "$motif_length\n"; exit;
	    my $end = $length - $motif_length;
	    for (my $pos_start=0; $pos_start <= $end; $pos_start++) { #loop through the given sequence, stopping after only the length of motif is left
		my $pos_end = $pos_start + $motif_length;
		my $logscore = -1000;
		my $sub_seq = substr( $seq, $pos_start, $motif_length); #substr(EXPR,OFFSET,LENGTH)
		my $comp_sub_seq = $sub_seq;
		$comp_sub_seq =~ tr/AGTC/TCAG/;
		my $rev_comp_sub = reverse($comp_sub_seq);			
		my @sub_seqs = ( $sub_seq , $rev_comp_sub ); #in order to check both strands...
		foreach my $string (@sub_seqs) {
		    my @bases = split(//,$string);
		    my $score = 1; #sets the probability score to 1, will rapidly decrease, usually
		    my $background_score = 1;
		    $strand = "+";  
		    for (my $i=0; $i<@bases; $i++ ) {#now find the probability of the base at each position for the selected TF
			my $base = uc $bases[$i];
			my $key = $tf_name_database[$e]."_".$base; #print "$key";
			my $base_probability = ${$pwm{$key}}[$i]; #print "$base_probability"; #pulls the probability of that base at that position
			if (!defined $base_probability ){
			    #print "undefined base_probability\n$key\t$string\n"; #print into STDOUT, cause problem if use `` to extract outputs in STDOUT (line100)
			    $undefined_base ++;
			    last;
			} #print ERRORS "key:".$key." is undefined for $i th base of string $string\n";								
			if ($base_probability == 0)  { $base_probability = 0.001; }
			$score = ($score * $base_probability);
			$background_score = ($background_score * ($backhash{$base}));#also have to calculate a background probability value for log odds calculation
			if ($i == (@bases - 1)) {
			    my $newlogscore = sprintf("%.3f", ((log($score/$background_score))/$denom)); #print $newlogscore."\n";
			    if ($newlogscore > $logscore) {
				$logscore = $newlogscore;  # assign the highest score of the two complementary sequences
				if ($string eq $rev_comp_sub ) { $strand = "-";} else { $strand = "+"; }  #sets the strand as forward or reverse
			    }
			}
		    }
		}
		if ( $logscore > $score_cutoff) {
		    $tf_bs_counter ++;
		    my $plotscore = ($logscore**2)*5; 
		    print RAW "$chr\t$seq_start\t$seq_end\t$tf_name_database[$e]\t$pos_start\t$pos_end\t$logscore\t$plotscore\t$strand\t$tf_bs_counter\n";
		    push(@sites, $pos_start);
		    push(@strands, $strand);
		}
	    }
	    print COUNT "$chr\t$seq_start\t$seq_end\t$tf_name_database[$e]\t$tf_bs_counter\n";
	    if ($tf_bs_counter >= 1) {$sum_pos_seq_counter ++;} #count the sequence that has at least 1 TF binding site.
	    $sum_tf_bs_counter +=  $tf_bs_counter; #count the total TF binding sites for all input sequences.
	    $tf_bs_counter = 0; 
	    #this is a weighting system for multiple sites within a specific window of sequence
	    for (my $i = 0; $i < (scalar(@sites)-1); $i++) {
		my $gene_score = 1;
		my $pos_start = $sites[$i];
		if (scalar @sites > 1) {
		    my $s = $strands[$i];
		    my $k = $i + 1;		
		    my $count = 1;
		    while (($k < scalar @sites) and (($sites[$k] - $sites[$i]) < 50)) {
			$k ++;
			$count ++; #print $count."\n";
		    }
		if ($count > 1){ $gene_score = ($gene_score + ($count ** $count)); }
		print WEIGHT "$chr\t$seq_start\t$seq_end\t$tf_name_database[$e]\t$pos_start\t$s\t$gene_score\n";	
		}
	    } 
	    last; #use "last" to find the first matched pwm and then stop looping
	} else {
	    if ($e == (@tf_name_database-1) and $match = 0) { die "no matched $tf_query was found in $pwm_filepath"; }
	    next; 
	}
    }
}
