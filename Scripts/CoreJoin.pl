#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my $sKey;
open (DATA,"total");
while(<DATA>){
   chomp;
   if(/(>.*)/){
      $sKey = $1; # keep the key
    }else{
    $hash{$sKey}.= $_; #combine all strings for the same key
    

    }
}

foreach(keys %hash){
 
 print "$_\n$hash{$_}\n";
}


