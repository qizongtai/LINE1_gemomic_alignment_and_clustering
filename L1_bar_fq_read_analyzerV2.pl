#! /usr/bin/perl -w
#------------------------------------------------------------------------------------------------------------
# Two parameters can be adjusted in the script:
# 1) the minimal polyA length (25bp)
# 2) the sequencing errors allowed within the polyA sequence (5bp)  
#------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
use Data::Dumper;

use constant CGT_MAX => 5; ## the max seq error within the ployA length
use constant POLYA_MIN => 25; ## the min polyA length

if (@ARGV != 3) { die "input: read1.fastq, read2.fastq, barcode"; }

my $read1 = shift (@ARGV);
my $read2 = shift (@ARGV);
my $bar_file = shift (@ARGV);
my $phix_path = "/scratch/rmlab/ref/phiX.fa"; #"~/indices/ref-seq/phiX.fa" doesn't work;
my @name_input = split (/-/,$read1); #get prefix for output files from the input fastq names.
my $name_output = $name_input[0];

unlink ("$name_output-phiX-read1ID") if (-e "$name_output-phiX-read1ID");
unlink ("$name_output-unknown-read1ID") if (-e "$name_output-unknown-read1ID");
unlink ("$name_output-bar-taq-polyTonly-read2ID") if (-e "$name_output-bar-taq-polyTonly-read2ID");
unlink ("$name_output-bar-msp-polyTonly-read2ID") if (-e "$name_output-bar-msp-polyTonly-read2ID");
unlink ("$name_output-phiX-read2ID") if (-e "$name_output-phiX-read2ID");
unlink ("$name_output-unknown-read2ID") if (-e "$name_output-unknown-read2ID");

open (SEQ1, "<", "$read1") || die 'read1 file error';
open (SEQ2, "<", "$read2") || die 'read2 file error';
open (PHIX, "<", "$phix_path") || die 'phix file error';
open (BAR, "<", "$bar_file") || die 'barcode file error';

##processing barcode file
my @barcode = ();
my %barcode_length = ();
while (my $line = <BAR>) {
    if ($line =~ /^\w/) {	
    chomp $line;
    push (@barcode, $line);
    $barcode_length {uc ($line)} = length ($line);
    }
}

##processing the input phix sequence
my $phix = "";
while (my $line = <PHIX>) {
    if ($line =~ /^\w/) {	
    chomp $line;
    $phix .= uc ($line);
    }
}

##output files for processed reads and statistics
open (OUT1, ">", "$name_output-polyA-length.txt") || die 'output1 error\n';
open (OUT2, ">", "$name_output-polyA-del-all-read1.fq") || die 'output2 error\n';
open (OUT3, ">", "$name_output-polyT-length.txt") || die 'output3 error\n';
open (OUT4, ">", "$name_output-polyT-del-all-read2.fq") || die 'output4 error\n';
open (OUT5, ">", "$name_output-read.overview") || die 'output5 error\n';
open (OUT6, ">", "$name_output-R1length-polyA-del.txt") || die 'output6 error\n';
open (OUT7, ">", "$name_output-R2length-polyT-del.txt") || die 'output7 error\n';
open (OUT8, ">", "$name_output-bar-distribution.txt") || die 'output8 error\n';

##output files for read1 lined ID (line number)
open (R1IDA, ">", "$name_output-plasmid-read1ID") || die 'outputA1 error\n';
open (R1IDB, ">", "$name_output-polyA-read1ID") || die 'outputB1 error\n'; # this include the IDs from polyA only reads
open (R1IDC, ">", "$name_output-polyA-only-read1ID") || die 'outputC1 error\n';
open (R1IDD, ">", "$name_output-polyA-no3bp-read1ID") || die 'outputD1 error\n';
open (R1IDE, ">", "$name_output-polyA-only-no3bp-read1ID") || die 'outputE1 error\n';
open (my $R1IDF, ">>", "$name_output-phiX-read1ID") || die 'outputF1 error\n';
open (my $R1IDG, ">>", "$name_output-unknown-read1ID") || die 'outputG1 error\n';
    
##processing the input read1
my ($total_R1counter, $phix_R1counter, $plasmid_R1counter, $polyA_3bp_R1counter, $polyAonly_3bp_R1counter,
    $polyA_no3bp_R1counter, $polyAonly_no3bp_R1counter, $unknown_3bp_R1counter, $unknown_no3bp_R1counter) = (0) x 9;

##variables to find the overlapping lineIDs (=line number) of qualified R1 and R2
my @polyAID = ();
my @cutterID = ();

while (my $line1 = <SEQ1>) {
    chomp $line1;
    chomp (my $line2 = <SEQ1>);
    chomp (my $line3 = <SEQ1>);
    chomp (my $line4 = <SEQ1>);
    $total_R1counter ++;
    if ($line2 =~ /^CAAC/) { ##with 4bp confirmation seq
	if ($line2 =~ /^CAACAACAATTGCATTCATT/) {
	    $plasmid_R1counter ++;
	    print OUT2 "$line1\n$line2\n$line3\n$line4\n";
	    print R1IDA "$total_R1counter\n";
	} else {
	    my ($target, $polyA_len, $cgt_len) = &removePolyA($line2, 0, 0);		
	    if ($target ne "") {
		$polyA_3bp_R1counter ++;
		print OUT1 "$polyA_len\n";
		print OUT2 "$line1\n$target\n$line3\n"; ##print only polyA-trimmed read1
		my $length_target = length($target);
		if ($line4 =~ /(.{$length_target})$/) {print OUT2 "$1\n";}
		print OUT6 "$length_target\n";
		print R1IDB "$total_R1counter\n";
		push (@polyAID, $total_R1counter);
	    } elsif ($polyA_len > (length($line2)-40) ) {
		$polyAonly_3bp_R1counter ++;
		print OUT1 "$polyA_len\n";
		print OUT2 "$line1\n\n$line3\n\n";#remove $line2 & line4
		print OUT6 "0\n";
		print R1IDB "$total_R1counter\n";
		print R1IDC "$total_R1counter\n";
		push (@polyAID, $total_R1counter);
	    } else {
		print OUT2 "$line1\n$line2\n$line3\n$line4\n";
		($phix_R1counter, $unknown_3bp_R1counter) =
		&phiX_seq($total_R1counter, $phix_R1counter, $unknown_3bp_R1counter, $line2, $R1IDF, $R1IDG);
	    }
	}
    } else { ##no 4bp confirmation seq
	my ($target, $polyA_len, $cgt_len) = &removePolyA($line2, 0, 0);		
	if ($target ne "") {
	    $polyA_no3bp_R1counter ++;
	    #print OUT1 "$polyA_len\n";
	    print OUT2 "$line1\n$target\n$line3\n"; ##print only polyA-trimmed read1
	    my $length_target = length($target);
	    if ($line4 =~ /(.{$length_target})$/) {print OUT2 "$1\n";}
	    #print OUT6 "$length_target\n";
	    print R1IDD "$total_R1counter\n";
	} elsif ($polyA_len > (length($line2)-40) ) {
	    $polyAonly_no3bp_R1counter ++;
	    #print OUT1 "$polyA_len\n";
	    print OUT2 "$line1\n$line2\n$line3\n$line4\n";
	    #print OUT6 "0\n";
	    print R1IDE "$total_R1counter\n";
	} else {
	    print OUT2 "$line1\n$line2\n$line3\n$line4\n";  
	    ($phix_R1counter, $unknown_no3bp_R1counter) =
	    &phiX_seq($total_R1counter, $phix_R1counter, $unknown_no3bp_R1counter, $line2, $R1IDF, $R1IDG);
	}
    }
}

##output files for read2 lined ID (line number)
open (R2IDA, ">", "$name_output-plasmid-read2ID") || die 'output2A error\n';
open (R2IDB, ">", "$name_output-bar-taq1-read2ID") || die 'output2B error\n';
open (R2IDC, ">", "$name_output-bar-msp1-read2ID") || die 'output2C error\n';
open (R2IDX, ">", "$name_output-bar-bfa1-read2ID") || die 'output2X error\n'; #newly added 
open (R2IDY, ">", "$name_output-bar-nocutter-read2ID") || die 'output2Y error\n'; #newly added
open (R2IDD, ">", "$name_output-bar-cutter-read2ID") || die 'output2D error\n';
open (my $R2IDE, ">>", "$name_output-bar-taq-polyTonly-read2ID") || die 'output2E error\n';
open (my $R2IDF, ">>", "$name_output-bar-msp-polyTonly-read2ID") || die 'output2F error\n';
open (my $R2IDP, ">>", "$name_output-bar-bfa-polyTonly-read2ID") || die 'output2P error\n'; #newly added
open (my $R2IDQ, ">>", "$name_output-bar-nocutter-polyTonly-read2ID") || die 'output2Q error\n'; #newly added
open (R2IDG, ">", "$name_output-nobar-read2ID") || die 'output2G error\n';
open (my $R2IDH, ">>", "$name_output-phiX-read2ID") || die 'output2H error\n';
open (my $R2IDI, ">>", "$name_output-unknown-read2ID") || die 'output2I error\n';
open (ID_BAR, ">", "$name_output-ID-barcode") || die 'output2J error\n';

##processing the input read2
my ($total_R2counter, $phix_R2counter, $plasmid_R2counter,
    $taq_counter, $msp_counter, $bfa_counter, $barcode_confseq_nocutter_counter,
    $taq_polyT_R2counter, $msp_polyT_R2counter, $bfa_polyT_R2counter, $nocutter_polyT_R2counter,
    $taq_no_polyT_R2counter, $msp_no_polyT_R2counter, $bfa_no_polyT_R2counter, $nocutter_no_polyT_R2counter,
    $taq_polyTonly_R2counter, $msp_polyTonly_R2counter, $bfa_polyTonly_R2counter, $nocutter_polyTonly_R2counter,
    $barcode_noconfseq_nocutter_R2counter,
    $nobarcode_R2counter, $unknown_R2counter) = (0) x 22;
my %seqID_barcode =();
my %barcode_seqID =();
my $barcode_size = keys %barcode_length; #get the size or number of keys in hash

while (my $line1 = <SEQ2>) {
    chomp $line1;
    chomp (my $line2 = <SEQ2>);
    chomp (my $line3 = <SEQ2>);
    chomp (my $line4 = <SEQ2>);
    $total_R2counter ++;
    my $nobarcode_counter = 0;
    foreach my $barcode (keys %barcode_length) {
	if ($line2 =~ /^($barcode)(.*)$/) {
	    $seqID_barcode {$line1} = $barcode;
	    push (@{$barcode_seqID {$barcode}}, $line1);
	    print ID_BAR "$line1\t$barcode\n";
	    my $line2_del_bar = $2; #$line2_bar = substr($line2, $barcode_length{$barcode}-1);
	    if ($line2_del_bar =~ /^AGCAGTGTTTAAACTAGTCGACCGGTAATAGTAATCAATTACGGG/
		or $line2_del_bar =~ /^AGCAGTGTTTAAACTAGTCGACCGGTAATAGTAATCAATTACG/
		or $line2_del_bar =~ /^AGCAGTGTTTAAACTAGTCGACCGGTAATAGTAATCAATTA/) {
		$plasmid_R2counter ++; # uncut or Taq1 or Msp1 or Bfa1
		print OUT4 "$line1\n$line2\n$line3\n$line4\n";
		print R2IDA "$total_R2counter\n";
		last;
	    } elsif ($line2_del_bar =~ /^AGCAGTGTTTAAACTAGTCGA(CCGG)(.*)$/) {
		$msp_counter ++;
		push (@cutterID, $total_R2counter);
		print R2IDB "$total_R2counter\n";
		print R2IDD "$total_R2counter\n";
		my $read2_msp_cutter = "CCGG";
		my $read2_msp_genomic = $1 . $2;
		my $rc_read2_msp_genomic = &rc_seq($read2_msp_genomic);
		my ($target, $polyA_len, $cgt_len) = &removePolyA ($rc_read2_msp_genomic, 0, 0);
		my %args = ('line1' => $line1, 'line2' => $line2, 'line3' => $line3, 'line4' => $line4,
			    'barcode' => $barcode, 'read2_cutter' => $read2_msp_cutter, 'read2_genomic' => $read2_msp_genomic,
			    'target' => $target, 'polyA_len' => $polyA_len, 'R2ID_polyTonly' => $R2IDF,
			    'total_R2counter' => $total_R2counter,
			    'polyT_R2counter' => $msp_polyT_R2counter,
			    'polyTonly_R2counter' => $msp_polyTonly_R2counter,
			    'no_polyT_R2counter' => $msp_no_polyT_R2counter);
		($msp_polyT_R2counter, $msp_polyTonly_R2counter, $msp_no_polyT_R2counter) = &read2_seq (\%args);
		last;
	    } elsif ($line2_del_bar =~ /^AGCAGTGTTTAAACTAG(TCGA)(.*)$/) {
		$taq_counter ++;
		push (@cutterID, $total_R2counter);
		print R2IDC "$total_R2counter\n";
		print R2IDD "$total_R2counter\n";
		my $read2_taq_cutter = "TCGA";
		my $read2_taq_genomic = $1 . $2;
		my $rc_read2_taq_genomic = &rc_seq ($read2_taq_genomic);
		my ($target, $polyA_len, $cgt_len) = &removePolyA ($rc_read2_taq_genomic, 0, 0);
		my %args = ('line1' => $line1, 'line2' => $line2, 'line3' => $line3, 'line4' => $line4, 
			    'barcode' => $barcode, 'read2_cutter' => $read2_taq_cutter, 'read2_genomic' => $read2_taq_genomic,
			    'target' => $target, 'polyA_len' => $polyA_len, 'R2ID_polyTonly' => $R2IDE,
			    'total_R2counter' => $total_R2counter,
			    'polyT_R2counter' => $taq_polyT_R2counter,
			    'polyTonly_R2counter' => $taq_polyTonly_R2counter,
			    'no_polyT_R2counter' => $taq_no_polyT_R2counter);
		($taq_polyT_R2counter, $taq_polyTonly_R2counter, $taq_no_polyT_R2counter) = &read2_seq (\%args);
		last;
	    } elsif ($line2_del_bar =~ /^AGCAGTGTTTAAA(CTAG)(.*)$/) {
		$bfa_counter ++;
		push (@cutterID, $total_R2counter);
		print R2IDX "$total_R2counter\n";
		print R2IDD "$total_R2counter\n";
		my $read2_bfa_cutter = "CTAG";
		my $read2_bfa_genomic = $1 . $2;
		my $rc_read2_bfa_genomic = &rc_seq ($read2_bfa_genomic);
		my ($target, $polyA_len, $cgt_len) = &removePolyA ($rc_read2_bfa_genomic, 0, 0);
		my %args = ('line1' => $line1, 'line2' => $line2, 'line3' => $line3, 'line4' => $line4, 
			    'barcode' => $barcode, 'read2_cutter' => $read2_bfa_cutter, 'read2_genomic' => $read2_bfa_genomic,
			    'target' => $target, 'polyA_len' => $polyA_len, 'R2ID_polyTonly' => $R2IDP,
			    'total_R2counter' => $total_R2counter,
			    'polyT_R2counter' => $bfa_polyT_R2counter,
			    'polyTonly_R2counter' => $bfa_polyTonly_R2counter,
			    'no_polyT_R2counter' => $bfa_no_polyT_R2counter);
		($bfa_polyT_R2counter, $bfa_polyTonly_R2counter, $bfa_no_polyT_R2counter) = &read2_seq (\%args);
		last;
	    } elsif ($line2_del_bar =~ /^AGCAGTGTTTAAA(.*)$/) {
		$barcode_confseq_nocutter_counter ++;
		push (@cutterID, $total_R2counter); #this is optional, use # if do not want to add this seq into the qualified R2
		print R2IDY "$total_R2counter\n";
		print R2IDD "$total_R2counter\n";
		my $read2_nocutter = "";
		my $read2_nocutter_genomic = $1;
		my $rc_read2_nocutter_genomic = &rc_seq ($read2_nocutter_genomic);
		my ($target, $polyA_len, $cgt_len) = &removePolyA ($rc_read2_nocutter_genomic, 0, 0);
		my %args = ('line1' => $line1, 'line2' => $line2, 'line3' => $line3, 'line4' => $line4, 
			    'barcode' => $barcode, 'read2_cutter' => $read2_nocutter, 'read2_genomic' => $read2_nocutter_genomic,
			    'target' => $target, 'polyA_len' => $polyA_len, 'R2ID_polyTonly' => $R2IDQ,
			    'total_R2counter' => $total_R2counter,
			    'polyT_R2counter' => $nocutter_polyT_R2counter,
			    'polyTonly_R2counter' => $nocutter_polyTonly_R2counter,
			    'no_polyT_R2counter' => $nocutter_no_polyT_R2counter);
		($nocutter_polyT_R2counter, $nocutter_polyTonly_R2counter, $nocutter_no_polyT_R2counter) = &read2_seq (\%args);
		last;
	    } else {
		$barcode_noconfseq_nocutter_R2counter ++;
		print OUT4 "$line1\n$line2\n$line3\n$line4\n";
		($phix_R2counter, $unknown_R2counter) =
		&phiX_seq($total_R2counter, $phix_R2counter, $unknown_R2counter, $line2, $R2IDH, $R2IDI);
		last;
	    } 
	} else {
	    $nobarcode_counter ++;
	    if ($nobarcode_counter == $barcode_size) {
		$nobarcode_R2counter ++;
		print OUT4 "$line1\n$line2\n$line3\n$line4\n"; #content
		print R2IDG "$total_R2counter\n"; #ID
		 ($phix_R2counter, $unknown_R2counter) =
		 &phiX_seq($total_R2counter, $phix_R2counter, $unknown_R2counter, $line2, $R2IDH, $R2IDI); 
	    }
	}
    }	
}

##output the number of reads assigned to each barcode
foreach my $barcode_id (@barcode) {
    if (!defined $barcode_seqID{$barcode_id}) {
	my $num = 0;
	print OUT8 "$barcode_id\t$num\n";
    } else {
	my $num = @{$barcode_seqID{$barcode_id}};
	print OUT8 "$barcode_id\t$num\n";
    }
}

##output the overlap ID between polyA_read1ID and cutter_read2ID 
open (OVERLAP, ">", "$name_output-polyA-cutter-ID") || die 'outputA3 error\n';
my @polyA_cutterID = ();
my %union = ();
my %isect = ();
foreach my $e (@polyAID,@cutterID) {
    $isect{$e}++ if ($union{$e});    
    $union{$e}++;
}

#alternative method 1
#foreach my $e (@polyAID) { $union{$e}++ }
#foreach my $e (@cutterID) {
#   if ( $union{$e} ) { $isect{$e} = ++ }
#   $union{$e} = ++;
#}

@polyA_cutterID = sort {$a <=> $b} keys %isect;
foreach my $element (@polyA_cutterID) {
   print OVERLAP "$element\n";
}

##output the statistical summary for each category of the reads
print OUT5 "the number of read1 is $total_R1counter\n";
print OUT5 "the number of phiX sequence is $phix_R1counter\n";
print OUT5 "the number of plasmid sequence is $plasmid_R1counter\n";
print OUT5 "the number of reads with 4bp, polyA and genomic sequences is $polyA_3bp_R1counter\n";
print OUT5 "the number of reads with 4bp and polyA only is $polyAonly_3bp_R1counter\n";
print OUT5 "the number of reads with 4bp but whitout polyA is $unknown_3bp_R1counter\n";
print OUT5 "the number of reads without 4bp but with polyA and genomic sequences is $polyA_no3bp_R1counter\n";
print OUT5 "the number of reads without 4bp but with polyA only is $polyAonly_no3bp_R1counter\n";
print OUT5 "the number of reads without 4bp or polyA is $unknown_no3bp_R1counter\n\n";

print OUT5 "the number of read2 is $total_R2counter\n";
print OUT5 "the number of phiX sequence is $phix_R2counter\n";
print OUT5 "the number of plasmid sequence is $plasmid_R2counter\n";
print OUT5 "the number of reads with barcode+confseq+Msp1 is $msp_counter\n";
print OUT5 "\t with genomic sequence and polyT is $msp_polyT_R2counter\n";
print OUT5 "\t with genomic sequence but not polyT is $msp_no_polyT_R2counter\n";
print OUT5 "\t with polyT ONLY but no genomic sequence is $msp_polyTonly_R2counter\n";
print OUT5 "the number of reads with barcode+confseq+TaqI site is $taq_counter\n";
print OUT5 "\t with genomic sequence and polyT is $taq_polyT_R2counter\n";
print OUT5 "\t with genomic sequence but not polyT is $taq_no_polyT_R2counter\n";
print OUT5 "\t with polyT ONLY but no genomic sequence is $taq_polyTonly_R2counter\n";
print OUT5 "the number of reads with barcode+confseq+Bfa1 site is $bfa_counter\n";
print OUT5 "\t with genomic sequence and polyT is $bfa_polyT_R2counter\n";
print OUT5 "\t with genomic sequence but not polyT is $bfa_no_polyT_R2counter\n";
print OUT5 "\t with polyT ONLY but no genomic sequence is $bfa_polyTonly_R2counter\n";
print OUT5 "the number of reads with barcode+confseq but no cutter is $barcode_confseq_nocutter_counter\n";
print OUT5 "\t with genomic sequence and polyT is $nocutter_polyT_R2counter\n";
print OUT5 "\t with genomic sequence but not polyT is $nocutter_polyTonly_R2counter\n";
print OUT5 "\t with polyT ONLY but no genomic sequence is $nocutter_no_polyT_R2counter\n";
print OUT5 "the number of reads with barcode but no confseq cutter site is $barcode_noconfseq_nocutter_R2counter\n";
print OUT5 "the number of reads with no barcode is $nobarcode_R2counter\n";
print OUT5 "the number of unknown sequence is $unknown_R2counter\n";

close (SEQ1); close (SEQ2); close (PHIX); close (BAR);
close (OUT1); close (OUT2); close (OUT3); close (OUT4); close (OUT5); close (OUT6); close (OUT7); close (OUT8); 
close (R1IDA); close (R1IDB); close (R1IDC); close (R1IDD); close (R1IDE); close ($R1IDF); close ($R1IDG); 
close (R2IDA); close (R2IDB); close (R2IDC); close (R2IDD); close ($R2IDE); close ($R2IDF); close (R2IDG); close ($R2IDH); close ($R2IDI);
close (R2IDX); close (R2IDY); close ($R2IDP); close ($R2IDQ); #newly added
close (ID_BAR); close (OVERLAP);

####################################################################################################################
sub rc_seq {
    my $seq = shift @_;
    my $rc_read = reverse $seq;   ##Get the reverse sequence###
    $rc_read =~ tr/ATGC/TACG/;    ##Get the complement of the sequence###
    return $rc_read;
}

sub removePolyA() {
    my ($seq, $polyA, $cgt) = @_;  ##$polyA=polyA counter ; $cgt=mismatch counter  
    if ($seq =~ /^(.*?)(A{6,})(.*)/) { ##six or more consecutive As are considered as a posible ployA
	if ($polyA) {
	    my $aNum = $1 =~ tr/A/A/; ## for situation like "TAAAAAC"
	    $cgt += (length($1) - $aNum);
	    $polyA += (length($2) + $aNum);
	    if ($cgt > CGT_MAX) { ##if the sequence error > 5, stop and check the previous ployA
		$polyA -= (length($2) + $aNum); 
		if ($polyA < POLYA_MIN){
		    return ("", 0, 0); ##in previous version:($seq, $polyA, $cgt) = &removePolyA($seq, 0, 0)
		} else {
		    return ($seq, $polyA, $cgt);
		}
	    }					
	} else {
	    $polyA = length($2);
	}
	($seq, $polyA, $cgt) = &removePolyA($3, $polyA, $cgt)
    } elsif ($polyA >= POLYA_MIN) { 
	return ($seq, $polyA, $cgt);
    } else {
	return ("", 0, 0);
    }	
}

sub read2_seq {
    my %args = %{$_[0]};
    if ($args{target} ne "") { 
	$args{polyT_R2counter} ++;
	print OUT3 "$args{polyA_len}\n";
	my $cutter_polyT_del_read2 = &rc_seq($args{target});
	print OUT4 "$args{line1}\n$cutter_polyT_del_read2\n$args{line3}\n";##print cutter and polyT trimmed read2
	my $length_barcode = length($args{barcode});
	my $length_cutter = length($args{read2_cutter});
	my $length_genomic = length($cutter_polyT_del_read2);
	if ($args{line4} =~ /^.{$length_barcode}.{$length_cutter}(.{$length_genomic})/) { print OUT4 "$1\n";}
	print OUT7 "$length_genomic\n";
    } elsif ($args{polyA_len} > 100) {
	$args{polyTonly_R2counter} ++;
	print OUT3 "$args{polyA_len}\n";
	print OUT4 "$args{line1}\n$args{read2_genomic}\n$args{line3}\n";
	my $length_barcode = length($args{barcode});
	my $length_cutter = length($args{read2_cutter});
	my $length_genomic = length($args{read2_genomic});
	if ($args{line4} =~ /^.{$length_barcode}.{$length_cutter}(.{$length_genomic})/) { print OUT4 "$1\n";}
	print OUT7 "0\n";
	print {$args{R2ID_polyTonly}} "$args{total_R2counter}\n";
    } else {
	$args{no_polyT_R2counter} ++;
	print OUT3 "0\n";
	print OUT4 "$args{line1}\n$args{read2_genomic}\n$args{line3}\n";
	my $length_barcode = length($args{barcode});
	my $length_cutter = length($args{read2_cutter});
	my $length_genomic = length($args{read2_genomic});
	if ($args{line4} =~ /^.{$length_barcode}.{$length_cutter}(.{$length_genomic})/) { print OUT4 "$1\n";}
	print OUT7 "$length_genomic\n";
    }
    return ($args{polyT_R2counter}, $args{polyTonly_R2counter}, $args{no_polyT_R2counter});
}

sub phiX_seq {
    my ($total_counter, $phix_counter, $unknown_counter, $tr_line2, $fh_phixID, $fh_unknowID) = @_;
    $tr_line2 =~ s/N/./g;
    my $rc_tr_line2 = &rc_seq($tr_line2);
    if ($tr_line2 =~/^(.{25})(.{25})(.{25})(.{25})(.{25})(.*)$/) {
	my $i = 0;
	my @match = ($1, $2, $3, $4, $5, $6);
	for (my $j = 0; $j <= 5; $j ++) { if ($phix =~ /$match[$j]/){$i++;} }				
	if ($i >= 3 ) {
	    $phix_counter ++;
	    print $fh_phixID "$total_counter\n";
	} else {
	    $rc_tr_line2 =~/^(.{25})(.{25})(.{25})(.{25})(.{25})(.*)$/ ;
	    my $i = 0;
	    my @match = ($1, $2, $3, $4, $5, $6);
	    for (my $j = 0; $j <= 5; $j ++) { if ($phix =~ /$match[$j]/){$i++;} }				
	    if ($i >= 3 ) {
		$phix_counter ++;
		print $fh_phixID "$total_counter\n";
	    } else {
		$unknown_counter ++;
		print $fh_unknowID "$total_counter\n";
	    }
	}
    }
    return ($phix_counter, $unknown_counter);
}
