#!/bin/perl -w
use strict;
use warnings;



sub help{
	if(grep {$_ eq "-h"} @ARGV){
		print "##################################################################################\n\n";
		print "	This script can be used to prepare the family file used in i-ADhore tool\n";
		print "	Usage: perl write_families.pl MAAs.table > families.csv\n";
		print "	Input: MAAs.table\n";
		print "	Output: families.csv\n";
		print "	output with two columns as: (1)	MAAs_ID;(2)	HOM_ID\n;"
		#print "###################################################################################\n\n";
		}
	}

&help();

my ($MAAs_table)=@ARGV;

open IN, $MAAs_table; # MAAs.table file
while (<IN>){
my @all=split /\t+/, $_;
print "$all[1].$all[0]\tHOM$all[0]\n";
}
close IN;

