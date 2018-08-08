#!/bin/env perl
#usage:perl ExtractSeqFromList.pl gene_name.ls all_gene.fasta



open (List,@ARGV[0]) || die ("could not open file");
open (FASTA , @ARGV[1]);
my @list= <List>;
$/ = '>';
my @fasta = <FASTA> ;
close FASTA;
 $/ = '';
foreach my $pattern (@list){
    chomp $pattern;
      foreach (@fasta){
          if(/^>$/){
              next;
          }
          s/>//;
          s/(\w+)/>\1/;
          if ($_=~/$pattern/){
            print "$_"; 
            last;
         }
      }
		}
	

