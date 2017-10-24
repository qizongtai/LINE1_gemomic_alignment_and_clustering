#! /usr/bin/perl -w
#----------------------------------------------------------#
# Pepline V2 for analysis of L1 unique insertion sites     #
#                                                          #
# Author: Zongtai Qi                                       #
# Mitra Lab, Washington Univ in St Louis.                  #
# Send all comments to qizongtai@wustl.edu                 #
# All Rights Reserved.                                     #
# Copyright (C) 2015                                       #
#----------------------------------------------------------#
use strict;
use warnings;
use IO::CaptureOutput qw(capture_exec);

my $usage = &usage_info;
if (@ARGV != 5) { die "$usage";}

my ($read1, $read2, $barcode, $map, $tf) = @ARGV;

my @map = split(//, $map);
my @tf = split(/\_/, $tf);
my %mapper_option =(
    "1" => "bowtie_read2",
    "2" => "bowtie_paired",
    "3" => "novo_read2",
    "4" => "novo_paired",
);

my @name_input = split (/-/,$read1);
my $name_output = $name_input[0];
my $map_option = '';
my $hg19_bowtie = "/scratch/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome";
my $hg19_novo = "/scratch/rmlab/ref/novoalign_indexes/hg19/hg19.ndx";
my $hg19_cyto = "~/indices/cytoplot/hg/cytoband_hg19_2.txt";

#step1-2:read filtering and fetching
&step1to2 ($name_output);
#step3-12:read mapping and analyzing
foreach my $option (@map) {
    $map_option = $mapper_option{$option};
    #Sept3: only one of the step3 commands is excuted based on mapper option by a "if" loop.
    my $step31 = "bowtie2 -x $hg19_bowtie $name_output-polyA-cutter-IDread2.fq -S $name_output-$map_option.sam";
    my $step32 = "bowtie2 -I 0 -X 1000 --dovetail -x $hg19_bowtie -1 $name_output-polyA-cutter-IDread1.fq -2 $name_output-polyA-cutter-IDread2.fq -S $name_output-$map_option.sam";
    my $step33 = "novoalign -r None -e 4 -d $hg19_novo -f $name_output-polyA-cutter-IDread2.fq -o SAM > $name_output-$map_option.sam";
    my $step34 = "novoalign -r None -e 4 -d $hg19_novo -f $name_output-polyA-cutter-IDread1.fq $name_output-polyA-cutter-IDread2.fq -o SAM > $name_output-$map_option.sam";
    if ($option == 1) {
	print "$step31\n";
	my ($stdout, $stderr, $success, $exit_code) = capture_exec ($step31);
	open(OVERVIEW, ">>", "$name_output-read.overview") or die 'overview file output error';
	print OVERVIEW "\n$stdout\n$stderr\n";
	print "$stdout\n$stderr\n";
	print "------read mapping(bwotie-single-end) done------\n\n";
	&step4to7 ($name_output, $map_option);
	&step8 ($name_output, $map_option);
	&step9 ($name_output, $map_option);
	&step10to12 ($name_output, $map_option, \@tf);
	&cat_file;
	next;
    }
    if ($option == 2) {
	print "$step32\n";
	my ($stdout, $stderr, $success, $exit_code) = capture_exec ($step32);
	open(OVERVIEW, ">>", "$name_output-read.overview") or die 'overview file output error';
	print OVERVIEW "\n$stdout\n$stderr\n";
	print "$stdout\n$stderr\n";
	print "------read mapping(bwotie-paired-end) done------\n\n";
	&step4to7 ($name_output, $map_option);
	&step8 ($name_output, $map_option);
	&step9 ($name_output, $map_option);
	&step10to12 ($name_output, $map_option, \@tf);
	&cat_file;
	next;
    }
    if ($option == 3) {
	print "$step33\n";
	my ($stdout, $stderr, $success, $exit_code) = capture_exec ($step33);
	open(OVERVIEW, ">>", "$name_output-read.overview") or die 'overview file output error';
	print OVERVIEW "\n$stdout\n$stderr\n";
	print "$stdout\n$stderr\n";
	print "------read mapping(novoalign-single-end) done------\n\n";
	&step4to7 ($name_output, $map_option);
	&step8 ($name_output, $map_option);
	&step9 ($name_output, $map_option);
	&step10to12 ($name_output, $map_option, \@tf);
	&cat_file;
	next;
    }
    if ($option == 4) {
	print "$step34\n";
	my ($stdout, $stderr, $success, $exit_code) = capture_exec ($step34);
	open(OVERVIEW, ">>", "$name_output-read.overview") or die 'overview file output error';
	print OVERVIEW "\n$stdout\n$stderr\n";
	print "$stdout\n$stderr\n";
	print "------read mapping(novoalign-paired-end) done------\n\n";
	&step4to7 ($name_output, $map_option);
	&step8 ($name_output, $map_option);
	&step9 ($name_output, $map_option);
	&step10to12 ($name_output, $map_option, \@tf);
	&cat_file;
	next;
    }
}
print "-----------all completed!-----------\n\n\n\n";

##############################################################################################################################################################
sub step1to2 {
    my $name_output = $_[0] ;
    my $step11 = "perl L1_bar_fq_read_analyzerV2.pl $read1 $read2 $barcode";
    #my $step11 = "perl L1_pT_bar_fq_read_analyzerV1.pl $read1 $read2 $barcode";
    my $step12 = "Rscript length_plotV1.R -a $name_output-polyA-length.txt -b $name_output-R1length-polyA-del.txt -c 250 -o $name_output-R1.pdf";
    my $step13 = "Rscript length_plotV1.R -a $name_output-polyT-length.txt -b $name_output-R2length-polyT-del.txt -c 250 -o $name_output-R2.pdf";
    my $step2 = "perl fq_read_fetcher_pairV3.pl $name_output-polyA-del-all-read1.fq $name_output-polyT-del-all-read2.fq $name_output-polyA-cutter-ID";
    print "\n$step11\n";
    system ( "$step11" );
    print "\n$step12\n";
    system ( "$step12" );
    print "\n$step13\n";
    system ( "$step13" );
    print "--------read analysis done--------\n\n";
    print "$step2\n";
    system ( "$step2" );
    print "--------read fetching done--------\n\n";
}

sub step4to7 {
    my ($name_output, $map_option) = @_ ; 
    my $step41 = "samtools view $name_output-$map_option.sam -S -b -q 10 -o $name_output-$map_option.bam"; ##-q 10 (mapping quality cut off=10)
    my $step42 = "samtools sort $name_output-$map_option.bam -o $name_output-$map_option-sorted.bam";
    my $step43 = "samtools index $name_output-$map_option-sorted.bam";
    my $step44 = "bedtools bamtobed -i $name_output-$map_option-sorted.bam > $name_output-$map_option-sorted.bed"; ## changed to "bedtools bamtobed" due to new version of bedtools
    my $step5 = "perl bar_bedV1.pl $name_output-$map_option-sorted.bed $name_output-ID-barcode";
    my $step61 = "perl uniqID_uniqPOS_bar_bedV1.pl $name_output-$map_option-sorted.bed.bar";
    my $step62 = "perl merge_uniq_bar_bedV2.pl $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat 180 8"; ## (180 -> pA direction) (8 -> restriction-site direction)
    my $step71 = "perl saturation_bar_bedV1.pl $name_output-$map_option-sorted.bed.bar.uniqID 50 100 50";
    my $step72 = "Rscript saturation_plotV1.R -i $name_output-$map_option-sorted.bed.bar.uniqID.saturation -o $name_output-$map_option-sorted.bed.bar.uniqID.saturation.pdf";

    print "$step41\n";
    system ( "$step41" );
    print "$step42\n";
    system ( "$step42" );
    print "$step43\n";
    system ( "$step43" );
    print "$step44\n";
    system ( "$step44" );    
    print "--------SamtoBed conversion done--------\n\n";
    
    print "$step5\n";
    system ( "$step5" );
    print "--------barcode attachment done--------\n\n";
    
    print "$step61\n";
    system ( "$step61" );
    print "$step62\n";
    system ( "$step62" );
    print "-------unique&merge insertion done--------\n\n";
    
    print "$step71\n";
    system ( "$step71" );
    print "$step72\n";
    system ( "$step72" );
    print "--------saturation analysis done--------\n\n";    
}

sub step8 {
    my ($name_output, $map_option) = @_ ;
    open(CLUSTER_MERGE, ">", "$name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge.clusterstat") || die "output clusterstat1 error\n";
    open(CLUSTER_UNIQID, ">", "$name_output-$map_option-sorted.bed.bar.uniqID.clusterstat") || die "output clusterstat2 error\n";
    print CLUSTER_MERGE "cluster_len\treads\tclusters\tmore_1barcode\tmore_1barcode_more_reads\n";
    print CLUSTER_UNIQID "cluster_len\treads\tclusters\tmore_1barcode\tmore_1barcode_more_reads\n";
    for (my $i=100; $i<=2000; $i+=100) {
	my $step81 = "perl cluster_merge_uniq_bar_bedV2.pl $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge $i 5"; # $i->distance to cluster; 5->minimal reads to call an uniq barcode
	my $step82 = "perl cluster_uniqID_bar_bedV1.pl $name_output-$map_option-sorted.bed.bar.uniqID $i 5"; # $i->distance to cluster; 5->minimal reads to call an uniq barcode
	print "$step81\n";
	my $cluster_stat1 = `$step81`;
	print "$step82\n";
	my $cluster_stat2 = `$step82`;
	print CLUSTER_MERGE "$i*2\t$cluster_stat1";
	print CLUSTER_UNIQID "$i*2\t$cluster_stat2";
    }
    print "---------cluster analysis done---------\n\n";
}

sub step9 {
    my ($name_output, $map_option) = @_ ;
    my $step91 = "Rscript cyto_plotV2.R -i $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge -c $hg19_cyto -o $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.100-50merge.pdf";
    my $step92 = "Rscript cyto_plotV2.R -i $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge.500cluster -c $hg19_cyto -o $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.100-50merge.500cluster.pdf";
    print "$step91\n";    
    system ( "$step91" );
    print "$step92\n";    
    system ( "$step92" );
    print "--------chr_distribution plot done--------\n\n";
}

sub step10to12 {
    my ($name_output, $map_option) = ($_[0], $_[1]); 
    my @tf = @{$_[2]};
    open (SCAN_ALL, ">", "$name_output-$map_option-sorted.bed.bar.uniqID.500cluster.scanstat") || die "output scanstat1 error\n";
    open (SCAN_MORE1BAR, ">", "$name_output-$map_option-sorted.bed.bar.uniqID.500cluster.more1bar.5reads.scanstat") || die "output scanstat2 error\n";
    print SCAN_ALL "extension\ttf_scan\tcutoff_score\tseq_counter\tsum_pos_seq_counter\tsum_tf_bs_counter\tcutoff_score\tseq_counter\tsum_pos_seq_counter\tsum_tf_bs_counter\n";
    print SCAN_MORE1BAR "extension\ttf_scan\tcutoff_score\tseq_counter\tsum_pos_seq_counter\tsum_tf_bs_counter\tcutoff_score\tseq_counter\tsum_pos_seq_counter\tsum_tf_bs_counter\n";
    for (my $i=0; $i<=500; $i+=100) {
	if ($i == 0) {
	    print SCAN_ALL "$i";
	    print SCAN_MORE1BAR "$i";
	} else {
	    print SCAN_ALL "\n$i";
	    print SCAN_MORE1BAR "\n$i";
	}
	my $step101 = "perl fa_seq_fetcherV1.pl $name_output-$map_option-sorted.bed.bar.uniqID.500cluster $i $i";
	my $step102 = "perl fa_seq_fetcherV1.pl $name_output-$map_option-sorted.bed.bar.uniqID.500cluster.more1bar.5reads $i $i";
	#my $step101 = "perl fa_seq_fetcherV1.pl $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge.500cluster. $i $i";
	#my $step102 = "perl fa_seq_fetcherV1.pl $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge.500cluster.more1bar.5reads $i $i";
	print "$step101\n";    
	system ( "$step101" );
	print "$step102\n";    
	system ( "$step102" );
	print "---------sequence retrieval done---------\n";
	foreach my $tf_scan (@tf) {
	    print SCAN_ALL "\t$tf_scan";
	    print SCAN_MORE1BAR "\t$tf_scan";
	    my @cutoff = (10, 12);
	    foreach my $cutoff (@cutoff) {
		my $step111 = "perl pwm_scannerV1.pl -txt $name_output-$map_option-sorted.bed.bar.uniqID.500cluster.$i-${i}extend.seq -tf $tf_scan -cf $cutoff";
		my $step112 = "perl pwm_scannerV1.pl -txt $name_output-$map_option-sorted.bed.bar.uniqID.500cluster.more1bar.5reads.$i-${i}extend.seq -tf $tf_scan -cf $cutoff";
		#my $step111 = "perl pwm_scannerV1.pl -txt $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge.500cluster.$i-${i}extend.seq -tf $tf_scan -cf $cutoff";
		#my $step112 = "perl pwm_scannerV1.pl -txt $name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge.500cluster.more1bar.5reads.$i-${i}extend.seq -tf $tf_scan -cf $cutoff";
		print "$step111\n";   
		my $scan_stat1 = `$step111`; chomp $scan_stat1;
		print "$step112\n";    
		my $scan_stat2 = `$step112`; chomp $scan_stat2;
		print SCAN_ALL "\t$cutoff\t$scan_stat1";
		print SCAN_MORE1BAR "\t$cutoff\t$scan_stat2";
	    }
	}
	print "-----------TF PMWs scan done-----------\n\n";
    }
    close (SCAN_ALL); close (SCAN_MORE1BAR);
}

sub cat_file {
    my $cat_step = "head -n -0 $name_output-read.overview " .
		   "$name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.overview " .
		   "$name_output-$map_option-sorted.bed.bar.uniqID.clusterstat " .
                   "$name_output-$map_option-sorted.bed.bar.uniqID.500cluster.scanstat " .
		   "$name_output-$map_option-sorted.bed.bar.uniqID.500cluster.more1bar.5reads.scanstat " .
		   "$name_output-$map_option-sorted.bed.bar.uniqID.uniqPOS.stat.100-50merge.clusterstat " .
		   "> $name_output.overview";
    print "$cat_step\n";
    system ( "$cat_step" );
    print "-------merge files into $name_output.overview-------\n\n";
}

sub usage_info {
    my $usage_info = '
Usage: <read1.fq> <read2.fq> <barcode> <mapper option> <TF_to_scan>; 
<Mapper option>: options can be combined without spaces in between(i.e. 12, 123 or 1234). Mapping database is hg19. 
 1) bowtie2 for read2 single end alingment    
 2) bowtie2 for paired end alignment
 3) novoalign for read2 single end alignment
 4) novoalign for paired end alignment
<TF_to_scan>: underscore delimited transcription factors to scan; 
Output name: prefix is from the first "-"(hyphen) deliminated input name;
';
    return $usage_info;
}

