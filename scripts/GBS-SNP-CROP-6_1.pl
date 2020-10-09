#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw( min max );
use File::Basename;

my ($file,$type);

GetOptions(
'b=s' => \$file,		# file
'p=s' => \$type,			    	# string
) or die "Error in options\n";

my $sttime = time;
$file = basename($file,  ".mpileup");

###################################
# Identifying both SNPs and indels 
###################################
if ($type eq "indel") {
	 
	  my $mpileup_input = join (".", "$file", "mpileup");
		my $count_out = join (".", "$file","count","txt");
		my $ref_file = join (".", "$file","ref","txt");

		open my $PILEUP, "<", "$mpileup_input" or die "Can't load $mpileup_input file";
		open my $OUT1, ">", "$count_out" or die "Can't initialized $count_out ouput file";
		open my $REF, ">", "$ref_file" or die "Unable to identify $ref_file\n";

		while (<$PILEUP>) {
	
			my @input1 = split("\t", $_);
			my $ref = $input1[2];
			my $algn = $input1[4];		
			$algn =~ s/\^.\././g;
			
			my @positions;
			my @sizes;
			my @Indels;
			while ( $algn =~ /[\.|,]{1}[\+|-]{1}(\d+)/g  ) { 
				push @positions, pos($algn);
				push @sizes, $1;
			}

			my $indices = scalar @positions - 1;

			for (my $k = $indices; $k >= 0; $k--) {
				my $indel = substr ( $algn, $positions[$k] - length($sizes[$k]) - 1, 1 + length($sizes[$k]) + $sizes[$k]);
				$indel = uc($indel);
				push @Indels, $indel;
				my $start = substr ( $algn, 0, $positions[$k] - length($sizes[$k]) - 2 );
				my $end = substr ( $algn, $positions[$k] + $sizes[$k] );
				$algn = join ("", "$start", "$end" );
			}

			# Count Indels frequency and store a hash with Indel type and count	
			my %Indel_cnt;
			$Indel_cnt{$_}++ foreach @Indels;
		
			my @IndelType = sort {$Indel_cnt{$b} <=> $Indel_cnt{$a}} keys %Indel_cnt;
			my @IndelCount = @Indel_cnt{@IndelType};
	
			# Work on a nucleotides specific string
			$algn =~ s/\$//g;   
			$algn =~ s/\*//g;
		
			my $uc_ref = uc $ref;
			my $lc_ref = lc $ref;
			$algn =~ s/\./$uc_ref/g;
			$algn =~ s/\,/$lc_ref/g;
		
			my @bases = split(//, $algn);
			@bases = grep /\S/, @bases;
		
			my $A = 0; my $a = 0; my $xA = 0;
			my $C = 0; my $c = 0; my $xC = 0;	
			my $G = 0; my $g = 0; my $xG = 0;
			my $T = 0; my $t = 0; my $xT = 0;
		
			for (my $x=0; $x<scalar(@bases);$x++) {
				if ($bases[$x] =~ /A/){
					$A++;
				}
				if ($bases[$x] =~ /a/){
					$a++;
				}
				if ($bases[$x] =~ /C/){
					$C++;
				}
				if ($bases[$x] =~ /c/){
					$c++;
				}
				if ($bases[$x] =~ /G/){
					$G++;
				}
				if ($bases[$x] =~ /g/){
					$g++;
				}
				if ($bases[$x] =~ /T/){
					$T++;
				}
				if ($bases[$x] =~ /t/){
					$t++;
				}
			}
		
			if ($A >= $a) {
				$xA = $A;
			} else {
				$xA = $a
			}
		
			if ($C >= $c) {
				$xC = $C;
			} else {
				$xC = $c
			}
		
			if ($G >= $g) {
				$xG = $G;
			} else {
				$xG = $g
			}
		
			if ($T >= $t) {
				$xT = $T;
			} else {
				$xT = $t
			}
		
			if (scalar @IndelCount == 0) { 	
				print $OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$xA","$xC","$xG","$xT","_","_","_","_"),"\n";
			
				if ( (($xA + $xC) * ($xG + $xT)) > 0 or (($xA + $xG) * ($xC + $xT)) > 0 ) {
					print $REF join ("\t", $input1[0], $input1[1], $input1[2]),"\n";
				}
				next;
		
			} elsif (scalar @IndelCount == 1) {
				print $OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$xA","$xC","$xG","$xT","$IndelCount[0]","$IndelType[0]","_","_"),"\n";
				print $REF join ("\t", $input1[0], $input1[1], $input1[2]),"\n";
					
			} elsif (scalar @IndelCount > 1) {
				print $OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$xA","$xC","$xG","$xT","$IndelCount[0]","$IndelType[0]","$IndelCount[1]","$IndelType[1]"),"\n";
				print $REF join ("\t", $input1[0], $input1[1], $input1[2]),"\n";
			
			} else {
				next;
			}
		}
		close $PILEUP;
		close $OUT1;
		close $REF;

	print "DONE.\n";

#########################
# Identifying SNPs only
#########################

} elsif ($type eq "snp") {
  
    my $mpileup_input = join (".", "$file", "mpileup");
		my $count_out = join (".", "$file","count","txt");
		my $ref_file = join (".", "$file","ref","txt");

		open my $PILEUP, "<", "$mpileup_input" or die "Can't load $mpileup_input file";
		open my $OUT1, ">", "$count_out" or die "Can't initialized $count_out ouput file";
		open my $REF, ">", "$ref_file" or die "Unable to identify $ref_file\n";
	
		while (<$PILEUP>) {	
			my @input1 = split("\t", $_);
			my $ref = $input1[2];
			my $algn = $input1[4];		
			$algn =~ s/\^.\././g;
			
			my @positions;
			my @sizes;
			my @Indels;
			while ( $algn =~ /[\.|,]{1}[\+|-]{1}(\d+)/g  ) { 
				push @positions, pos($algn);
				push @sizes, $1;
			}

			my $indices = scalar @positions - 1;

			for (my $k = $indices; $k >= 0; $k--) {
				my $indel = substr ( $algn, $positions[$k] - length($sizes[$k]) - 1, 1 + length($sizes[$k]) + $sizes[$k]);
				$indel = uc($indel);
				push @Indels, $indel;
				my $start = substr ( $algn, 0, $positions[$k] - length($sizes[$k]) - 2 );
				my $end = substr ( $algn, $positions[$k] + $sizes[$k] );
				$algn = join ("", "$start", "$end" );
			}

			$algn =~ s/\$//g;   
			$algn =~ s/\*//g;
		
			my $uc_ref = uc $ref;
			my $lc_ref = lc $ref;
			$algn =~ s/\./$uc_ref/g;
			$algn =~ s/\,/$lc_ref/g;
		
			my @bases = split(//, $algn);
			@bases = grep /\S/, @bases;
		
			my $A = 0; my $a = 0; my $xA = 0;
			my $C = 0; my $c = 0; my $xC = 0;	
			my $G = 0; my $g = 0; my $xG = 0;
			my $T = 0; my $t = 0; my $xT = 0;
		
			for(my $x=0; $x<scalar(@bases);$x++) {
				if ($bases[$x] =~ /A/){
					$A++;
				}
				if ($bases[$x] =~ /a/){
					$a++;
				}
				if ($bases[$x] =~ /C/){
					$C++;
				}
				if ($bases[$x] =~ /c/){
					$c++;
				}
				if ($bases[$x] =~ /G/){
					$G++;
				}
				if ($bases[$x] =~ /g/){
					$g++;
				}
				if ($bases[$x] =~ /T/){
					$T++;
				}
				if ($bases[$x] =~ /t/){
					$t++;
				}
			}
		
			if ($A >= $a) {
				$xA = $A;
			} else {
				$xA = $a
			}
		
			if ($C >= $c) {
				$xC = $C;
			} else {
				$xC = $c
			}
		
			if ($G >= $g) {
				$xG = $G;
			} else {
				$xG = $g
			}
		
		if ($T >= $t) {
				$xT = $T;
			} else {
				$xT = $t
			}
			
			print $OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$xA","$xC","$xG","$xT"),"\n";
		
			if ( (($xA + $xC) * ($xG + $xT)) > 0 or (($xA + $xG) * ($xC + $xT)) > 0 ) {
				print $REF join ("\t", $input1[0], $input1[1], $input1[2]),"\n";
			}
		}
		close $PILEUP;
		close $OUT1;
		close $REF;

	print "DONE.\n";
}


FINAL:
exit;
