#!/bin/perl -w
use strict;
use warnings;



sub help{
	if(grep {$_ eq "-h"} @ARGV){
		print "############################################################################################\n\n";
		print "	This script can be used to prepare the species_scaffold_lst/* files used in i-ADhore tool \n";
		print "	Usage:	perl write_scaffold_lst.pl MAAs_gff(or genes_gff) species \n";
		print "	Input:	MAAs_gff, species \n";
		print "	Output:	species_scaffold_lst/*scaffold.lst \n";
		print "###################################################################################\n\n";
		}
	}

&help();

####shift: spec_gene.gff OR spec_MAAs.gff
(my $MAAs_gff, my $spec)=@ARGV;
open IN, $MAAs_gff;
my $last_scaffold="";
while (<IN>){
	chomp($_);
    my @all=split /\s+/, $_;
    my $gene_id="$spec.$all[0]";
	if ($last_scaffold ne $all[2]){
		###scaffold_lst formate: Cucumber1+
		open OUT , ">$spec\_scaffold_lst/$all[2]\.lst";
		print OUT  "$gene_id$all[5]\n";
		close OUT;
	}else {
		open OUT, ">>$spec\_scaffold_lst/$last_scaffold\.lst";
		print OUT  "$gene_id$all[5]\n";
		close OUT;
		}
	$last_scaffold=$all[2];
	}
close IN;



