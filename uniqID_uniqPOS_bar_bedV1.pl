#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(min max);
$| = 1;

if (@ARGV !=1) {die "input: bed.bar or bed.sortbar file";}

my $bed = $ARGV[0];
my $bed_uniqID = $ARGV[0] . ".uniqID";
my $bed_uniqPOS = $ARGV[0] . ".uniqID.uniqPOS";
my $bed_uniqPOS_stat = $ARGV[0] . ".uniqID.uniqPOS.stat";
my $rawgnashy = $ARGV[0] . ".uniqID.uniqPOS.rawgnashy";
my $gnashy = $ARGV[0] . ".uniqID.uniqPOS.gnashy";
my $same_id_error = $ARGV[0] . ".iderror";
my $overview = $ARGV[0] . ".uniqID.uniqPOS.overview";
unlink ($bed_uniqID) if (-e $bed_uniqID);
unlink ($bed_uniqPOS) if (-e $bed_uniqPOS);
unlink ($bed_uniqPOS_stat) if (-e $bed_uniqPOS_stat);
unlink ($rawgnashy) if (-e $rawgnashy);
unlink ($gnashy) if (-e $gnashy);
unlink ($same_id_error) if (-e $same_id_error);
unlink ($overview) if (-e $overview);

############################################################################################################################
print "Merging pair-aligned reads by illumina ID ... ";
open (FH, "<", "$bed") or die 'input bed file error\n';  # (before was) open(FH,$bed);
my $total_read_counter = 0;
my $merge_id_couter = 0; # same as below "$total_paired_unpaired_counter" and generated in merging ID process
my $unpaired_counter = 0;
my $paired_counter = 0;
my $same_id_error_counter = 0;
my $more_id_error_counter = 0;
my ($chr, $start, $end, $id, $score, $strand, $barcode, @dump1) = ();
my %hash = ();
my @ids = ();
my $last_chr;

while (<FH>) {
    $total_read_counter ++;
    chomp;
    ($chr, $start, $end, $id, $score, $strand, $barcode, @dump1) = split(/\t+/,$_);
    if (!$last_chr) { $last_chr = $chr; } # this argument is used only for the first time 
    elsif ($chr ne $last_chr) { # chr by chr; one chr at each time
        &bed_uniqID;
        %hash = ();
        @ids = ();
        $last_chr = $chr;
    }
    #Important: the following codes need to be after "&bed_uniqID" to deal with the last chr. 
    my ($id_sn, $id_no) = (0, 0);
    if ($id =~ /\//) { 
        ($id_sn, $id_no) = split('/', $id); # the paired read has id_no (1 or 2) wiht "/" separator.
    } else {
        $id_sn = $id;# the single read doesn't have id_no and set the id_no = 0.
    }
    push(@ids, $id_sn) if (!$hash{$id_sn});# use array to preserve the order of the hash below
    push(@{$hash{$id_sn}}, "$start\t$end\t$id_no\t$score\t$strand\t$barcode");
}
&bed_uniqID; #run "&bed_uniqID" again for the last chr
close FH;
print "Completed!\n";

#############################################################################################################################
print "Merging reads by start pos, end pos, strand and barcode... ";
open(FH, "<", "$bed_uniqID");
my $total_paired_unpaired_counter = 0; #same as above "$merge_id_couter" and generated in merging POS process
my $unique_pos_counter = 0;
my @pos_strand_bar = ();
($chr, $start, $end, $id, $score, $strand, $barcode, @dump1) = ();
%hash = ();
$last_chr = '';

while (<FH>) {
    $total_paired_unpaired_counter++;
    chomp;
    ($chr, $start, $end, $id, $score, $strand, $barcode, @dump1) = split(/\t+/, $_);
    if (!$last_chr) { $last_chr = $chr; }
    elsif ($chr ne $last_chr) {
        &bed_uniqPOS;
        %hash = ();
        @pos_strand_bar = ();
        $last_chr = $chr;
    }
    #Important: the following codes need to be after "&bed_uniqPOS" to deal with the last chr.
    if (!$hash{"$start:$end:$strand:$barcode"}) { #use non-defined of a hash key to get the unique positions. Note: only the first unique element get stored!
        push(@pos_strand_bar, "$start:$end:$strand:$barcode"); #store the first unique positions into a list.
        $hash{"$start:$end:$strand:$barcode"}{'content'} = $_; #store the content from the very first element.
        # $hash{"$start:$end:$strand"} = {} hash reference
        # $hash{"$start:$end:$strand"} -> {'content'} = $_;
        # ${$hash{"$start:$end:$strand"}}{'content'} = $_;
    }
    $hash{"$start:$end:$strand:$barcode"}{'count'} ++;
    # use defined hash key to count the repeated values;
    # $hash{"$start:$end:$strand"}++ is wrong as $hash{"$start:$end:$strand"} is hashref already.
}
&bed_uniqPOS; #run "&bed_uniqPOS" again for the last chr
close FH;
print "Completed!\n";
print "The total number of reads is $total_read_counter:\n";
print "$unpaired_counter are single reads (from single-end aligment) and $paired_counter are paired reads (from pair-end aligment).\n";
print "Merged by unique readID assigned by illumina seq (only paired reads are merged), $total_paired_unpaired_counter reads are left.\n";
print "Merged by unique chr posistion, $unique_pos_counter reads are left \n\n ";

open (OVERVIEW, ">", "$overview") or die 'output overview file error\n'; 
print OVERVIEW "The total number of reads is $total_read_counter:\n";
print OVERVIEW "$unpaired_counter are single reads (from single-end aligment) and $paired_counter are paired reads (from pair-end aligment).\n";
print OVERVIEW "Merged by unique readID assigned by illumina seq (only paired reads are merged), $total_paired_unpaired_counter reads are left.\n";
print OVERVIEW "Merged by unique chr posistion, $unique_pos_counter reads are left and are distributed as follows:\n\n";
close (OVERVIEW);

system (" cut -f 1 $bed_uniqPOS | sort | uniq -c >> $overview "); # "cut|sort|uniq"; "cut -f1" -> select for the first column; "-c" -> counts for each unique value. 
#system ("wc -l $bedfile_name > chr_list.txt");  #this outputs total counts of unique insertions.
#system ( "sort -u -k1,1 $bedfile > chr_list.txt" ); #"sort" outputs all columns although "-k1,1" sorts on first column. Need "cut -f1" to filter for just first column.

#this sub can append the overview info onto name-read.overview
#&append_info ($overview, $ARGV[0]);

system (" cut -c 4- $rawgnashy | grep '^[0-9XY]' | grep -v '[0-9]_g' > $gnashy ");  # from the first 4 characters to the end (omit the end position).
system (" sed -i 's/X/23/g' $gnashy ");
system (" sed -i 's/Y/24/g' $gnashy ");

################################################################################################################################
sub bed_uniqID {
    open(OUT, ">>", "$bed_uniqID");
    open(SAME_ID_ERROR, ">>", "$same_id_error");
    my ($start, $end, $id_no, $score, $strand, $barcode) = ();
    foreach my $id_sn (@ids) { # use array to preserve the order of a hash
        $merge_id_couter ++;
        ID_CONDITION: {
        if (@{$hash{$id_sn}} == 1) {
            # seqID ($id_sn) only has one corresponding seq, meaning seq was from single end mapping
            # the single read doesn't have $id_no and $id_no was set to 0 in the above codes $id_no = 0
            $unpaired_counter ++; 
            ($start, $end, $id_no, $score, $strand, $barcode) = split("\t", ${$hash{$id_sn}}[0]);
            print OUT "$last_chr\t$start\t$end\t$id_sn/$id_no\t$score\t$strand\t$barcode\n";
        } elsif (@{$hash{$id_sn}} == 2) {
            # seqID ($id_sn) has two corresponding seq, meaning seq was from paired end mapping
            # the paired read has either id_no = 1 or 2
            $paired_counter ++; 
            my $id_to_bar = '';
            my @pos_all = ();
            foreach my $data (@{$hash{$id_sn}}) {
                my ($start, $end, $id_no, $score, $strand, $barcode) = split("\t", $data);
                push (@pos_all, ($start, $end));
                if ($id_no == 0) { # id_no=0 means signle end mapping which shouldn't have two seqs with the same seqID. An error!
                    $same_id_error_counter ++;
                    $paired_counter --;
                    print SAME_ID_ERROR "$last_chr\t$start\t$end\t$id_sn/$id_no\t$score\t$strand\t$barcode\n";
                    last ID_CONDITION;
                } elsif ($id_no == 1) { # for a megred paired read, keep the info from "$id_no = 1" and discard that = 2. 
                    $id_to_bar = "$id_sn/$id_no\t$score\t$strand\t$barcode";
                }
            }
            $start = min(@pos_all);
            $end = max(@pos_all);
            print OUT "$last_chr\t$start\t$end\t$id_to_bar\n";
        } else {
            $more_id_error_counter ++;
            foreach my $data (@{$hash{$id_sn}}) {
                my ($start, $end, $id_no, $score, $strand, $barcode) = split("\t", $data);
                print SAME_ID_ERROR "$last_chr\t$start\t$end\t$id_sn/$id_no\t$score\t$strand\t$barcode\n";
            }
        }
        }
    }
    close OUT;
}

sub bed_uniqPOS {
    open(OUT, ">>$bed_uniqPOS");
    open(STAT, ">>$bed_uniqPOS_stat");
    open(GNASHY, ">>$rawgnashy");
    foreach my $pos_strand_bar (@pos_strand_bar) { # use array to preserve the order of a hash
        $unique_pos_counter ++;
        print OUT "$hash{$pos_strand_bar}{'content'}\n";
        my ($start, $end, $strand, $barcode) = split(':',$pos_strand_bar);
        print STAT "$last_chr\t$start\t$end\t$strand\t$barcode\t$hash{$pos_strand_bar}{'count'}\n";
        print GNASHY "$last_chr\t$start\t$hash{$pos_strand_bar}{'count'}\n";
    }
    close STAT;
    close OUT;
    close GNASHY;
}

sub append_info {
    my ($info, $target_name) = @_;
    my @name_input = split (/-/, $target_name);
    my $name_output = $name_input[0];
    open (IN, "<", "$info") or die $!;
    open (ADD, ">>", "$name_output-read.overview") or die $!;
    while (<IN>) { print ADD $_; }
    close (IN);
    close (ADD);
}

