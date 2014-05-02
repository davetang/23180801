#!/usr/bin/perl

#This program is free software: you can redistribute it and/or modify it under the terms of the GNU
#General Public License as published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This package is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with this program. If
#not, see <http://www.gnu.org/licenses/>

#Written by Dave Tang <davetingpongtang@gmail.com>
#Thu Feb 23 2012
#Updated Thu, Mar 07 2013 for nanoCAGE 2.0 protocol.

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('f:s:g:e:h', \%opts);

if ($opts{'h'} || !exists $opts{'f'} || !exists $opts{'s'} || !exists $opts{'g'} || !exists $opts{'e'}){
   usage();
}

my $infile = $opts{'f'};
my $spacer = $opts{'s'};
my $edit_threshold = $opts{'e'};
my $genome = $opts{'g'};

#check spacer sequence, translate into uppercase and append GGG
die "Spacer $spacer not recognised; please enter nucleotide spacer sequence\n" unless $spacer =~ /^[AGTC]+$/i;
$spacer =~ tr/acgt/ACGT/;
$spacer =~ s/$/GGG/;
my $spacer_length = length($spacer);

if ($infile !~ /[sb]am$/){
   die "Is $infile a sam or bam file\n";
}

#read hg19 genome into memory
my %genome = read_genome($genome);
warn "Genome stored into memory\n";

#open file handle for removed tags
my $random_string = random_string(6);
my $basename = $infile;
$basename =~ s/\.[sb]am$//;
my $removed = $basename . "_nw_${edit_threshold}_${random_string}_removed.sam";
open(REMOVED,'>',$removed) || die "Could not open $removed for writing: $!\n";
my $outfile = $basename . "_nw_${edit_threshold}_${random_string}_filtered.sam";
open(OUT,'>',$outfile) || die "Could not open $outfile for writing: $!\n";

#store number of removed tags
my $removed_tally = '0';

#read sam file
if ($infile =~ /\.sam$/){
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
} elsif ($infile =~ /\.bam$/){
   open(IN,'-|',"samtools view -h $infile") || die "Could not open $infile: $!\n";
}
while(<IN>){
   chomp;
   if (/^@/){
      #hack for files produced by internal pipeline
      if (/^\@RG/){
         s/BC:[ACGT]{6}\t//;
      }
      #output header information into file handles
      print OUT "$_\n";
      print REMOVED "$_\n";
      next;
   }
   #HWI-EAS420_0001:1:51:12051:15569#0/1    0       chr1    3004760 2       36M     *       0       0       ACAGTGGGGCACCTTGGGCACGGAGTTGGCAGACAC    *       NM:i:4  MD:Z:0G2C1A2T27 XP:Z:$$$$$$((((((~~~~~~~~~~~~''''''EECA?H
   my ($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,$tag,$vtype,$value) = split;
   next if $flag & 0x0004; #skip unmapped
   my $strand = '+';
   $strand = '-' if $flag & 0x0010;

   my $length = length($seq);
   my $end = $pos + $length;

   my $upstream = '';
   if ($strand eq '+'){
      $upstream = substr($genome{$rname}, $pos -1 - $spacer_length, $spacer_length);
   } else {
      $upstream = substr($genome{$rname}, $end -1, $spacer_length);
      $upstream =~ tr/ACGT/TGCA/;
      $upstream =~ tr/acgt/tgca/;
      $upstream = reverse($upstream);
   }
   $upstream = uc($upstream);

   my $edit = align($spacer,$upstream);

   #if the number of edits is greater than the threshold, remove this tag
   if ($edit <= $edit_threshold){
      my $last_three = substr($upstream, $spacer_length - 3, 3);
      #only remove if at least 2 of the last 3 nucleotides are guanosines
      if ($last_three =~ /GG$/ || $last_three =~ /G[ACGT]G$/ || $last_three =~ /GG[ACGT]$/){
         ++$removed_tally;
         print REMOVED "$_\n";
      } else {
         print OUT "$_\n";
      }
   } else {
      print OUT "$_\n";
   }

}
close(IN);

#convert to bam

foreach my $file_to_bam ($outfile, $removed){
   my $outfile_bam = $file_to_bam;
   $outfile_bam =~ s/sam$/bam/;
   my $convert_bam = "samtools view -bS -T $genome $file_to_bam > $outfile_bam";
   system($convert_bam);

   #sort bam
   my $outfile_bam_sorted = $outfile_bam;
   $outfile_bam_sorted =~ s/\.bam$/_sorted/;
   my $sort_bam = "samtools sort $outfile_bam $outfile_bam_sorted";
   system($sort_bam);

   #remove sam and unsorted bam file
   unlink($file_to_bam);
   unlink($outfile_bam);
}

warn "$removed_tally entries removed at threshold $edit_threshold\n";

#align spacer to upstream sequence
sub align {
   my ($seq1,$seq2) = @_;
   my $edit = '0';

   my $match = '1';
   my $mismatch = '-1';
   my $gap = '-1';

   my @matrix;
   $matrix[0][0]{score}   = 0;
   $matrix[0][0]{pointer} = "none";
   for(my $j = 1; $j <= length($seq1); $j++) {
      $matrix[0][$j]{score}   = $gap * $j;
      $matrix[0][$j]{pointer} = "left";
   }
   for (my $i = 1; $i <= length($seq2); $i++) {
      $matrix[$i][0]{score}   = $gap * $i;
      $matrix[$i][0]{pointer} = "up";
   }

   # fill
   for(my $i = 1; $i <= length($seq2); $i++) {
      for(my $j = 1; $j <= length($seq1); $j++) {
         my ($diagonal_score, $left_score, $up_score);

         # calculate match score
         my $letter1 = substr($seq1, $j-1, 1);
         my $letter2 = substr($seq2, $i-1, 1);
         if ($letter1 eq $letter2) {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $match;
         } else {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $mismatch;
         }

         # calculate gap scores
         $up_score   = $matrix[$i-1][$j]{score} + $gap;
         $left_score = $matrix[$i][$j-1]{score} + $gap;

         # choose best score
         if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
               $matrix[$i][$j]{score}   = $diagonal_score;
               $matrix[$i][$j]{pointer} = "diagonal";
            } else {
               $matrix[$i][$j]{score}   = $left_score;
               $matrix[$i][$j]{pointer} = "left";
            }
         } else {
            if ($up_score >= $left_score) {
               $matrix[$i][$j]{score}   = $up_score;
               $matrix[$i][$j]{pointer} = "up";
            } else {
               $matrix[$i][$j]{score}   = $left_score;
               $matrix[$i][$j]{pointer} = "left";
            }
         }
      }
   }
   # trace-back
   my $align1 = "";
   my $align2 = "";

   # start at last cell of matrix
   my $j = length($seq1);
   my $i = length($seq2);

   while (1) {
      last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix

      if ($matrix[$i][$j]{pointer} eq "diagonal") {
         $align1 .= substr($seq1, $j-1, 1);
         $align2 .= substr($seq2, $i-1, 1);
         $i--;
         $j--;
      } elsif ($matrix[$i][$j]{pointer} eq "left") {
         $align1 .= substr($seq1, $j-1, 1);
         $align2 .= "-";
         $j--;
      } elsif ($matrix[$i][$j]{pointer} eq "up") {
         $align1 .= "-";
         $align2 .= substr($seq2, $i-1, 1);
         $i--;
      }
   }
   $align1 = reverse $align1;
   $align2 = reverse $align2;

   for (my $i=0; $i<length($align1); ++$i){
      if (substr($align1,$i,1) ne substr($align2,$i,1)){
         ++$edit;
      }
   }
   return($edit);
}

#store genome into %genome
sub read_genome {
   my ($infile) = @_;
   warn "Reading genome into memory...\n";
   my %genome = ();
   my $current = '';
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      if (/^>(\S*)/){
         $current = $1;
         warn("Processing $current\n");
      } else {
         $genome{$current} .= $_;
      }
   }
   close(IN);
   return(%genome);
}

sub random_string {
   my ($length) = @_;

   my @char = (0 .. 9);
   push(@char, qw/A B C D E F G H I J K L M N O P Q R S T U V W X Y Z/);
   push(@char, qw/a b c d e f g h i j k l m n o p q r s t u v w x y z/);

   my $l = scalar(@char);
   my $random_string = '';
   for (1 .. $length){
      my $rand = int(rand($l));
      $random_string .= $char[$rand];
   }
   return($random_string);
}

sub usage {
print STDERR <<EOF;
Usage: $0 -f file -s spacer -g genome.fa -e edit distance

Where:   -f infile.bam          input sam or bam file
         -s TATA                nucleotide spacer sequence
         -g hg19.fa             fasta file containing genome sequence
         -e 1                   the number of edits at which to remove reads
         -h                     this helpful usage message

EOF
exit();
}

__END__
