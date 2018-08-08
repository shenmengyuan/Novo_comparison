#! /usr/bin/perl -w
use strict;
my($lines,@lines,$line_nums,@b,$j);

open(ID,"count_all");
@lines = <ID>;
$line_nums = @lines;
foreach $_(@lines){
	chomp;
	@b=split/\s/;
	$b[0]=~s/://;
	#print "$b[0]\n";
	open(FW,">$b[0].txt");
	for ($j=1;$j<@b ;$j++) {
        	print FW "$b[$j]\n";
	system("~/software/seqkit grep -f $b[0].txt ../total.faa >$b[0].fa")
	}
	close(FW);
}
