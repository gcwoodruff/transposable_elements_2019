#!/usr/bin/perl
#from https://gist.github.com/avrilcoghlan/4037d6b8cca32eaf48b0 
use strict;
use warnings;

my $fasta = $ARGV[0];
my $gff = $ARGV[1];

# read in the gff file to find which sequences to take:
my %TAKE = ();
my $num_to_take = 0;
open(GFF,"$gff");
while(<GFF>)
{  
   my $line = $_;
   chomp $line;
   if (substr($line,0,1) ne '#')
   {  
      # seq0    LTRharvest      LTR_retrotransposon     227145  239120  .       -       .       ID=LTR_retrotransposon1;Parent=repeat_region 1;ltr_similarity=98.06;seq_number=0
      # seq0    LTRharvest      LTR_retrotransposon     1847908 1862843 .       -       .       ID=LTR_retrotransposon2;Parent=repeat_region 2;ltr_similarity=87.74;seq_number=0
      # seq0    LTRharvest      LTR_retrotransposon     3847365 3852453 .       +       .       ID=LTR_retrotransposon3;Parent=repeat_region 3;ltr_similarity=96.35;seq_number=0
      my @temp = split(/\t+/,$line);
     my  $seq = $temp[0]; # eg. seq0
     my  $feature = $temp[2];
      if ($feature eq 'LTR_retrotransposon')
      {  
        my  $start = $temp[3];
        my  $end = $temp[4];
         $seq = $seq."_".$start."_".$end;
         $TAKE{$seq} = 1;
         $num_to_take++;
      }
   }
}
close(GFF);

# read in the fasta file of sequences, and print out those we want to take:
my $take = 0;
my $num_found = 0;
open(FASTA,"$fasta");
while(<FASTA>)
{  
  my  $line = $_;
   chomp $line;
   if (substr($line,0,1) eq '>') # >chromosome:WBcel235:I:1:15072434:1 chromosome I (dbseq-nr 0) [65157,67865]
   {
     my  @temp = split(/\s+/,$line);
     my  $seqno = $temp[$#temp-1]; # eg. 0)
      chop($seqno); # eg. 0
     my  $pos = $temp[$#temp]; # eg. [65157,67865]
      $pos = substr($pos,1,length($pos)-2);
      @temp = split(/\,/,$pos);
      my $start = $temp[0];
      my $end = $temp[1];
      my $seq = "seq".$seqno."_".$start."_".$end;
      if ($TAKE{$seq}) { $take = 1; $num_found++;} else { $take = 0;}
   }    
   if ($take == 1) { print "$line\n";}
}
close(FASTA);
if ($num_found != $num_to_take) { print STDERR "ERROR: num_found $num_found num_to_take $num_to_take\n"; exit;}
print STDERR "FINISHED\n";
