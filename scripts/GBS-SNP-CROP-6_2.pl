#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw( min max );
use File::Basename;

my ($posFile, $countList,$output_file,$type);

$output_file = 'GSC.DiscoveryMatrix.txt';

GetOptions(
'in=s' => \$posFile,	     	# file
'count=s' => \$countList,   # file
'out=s' => \$output_file,		# file
'p=s' => \$type 			    	# string
) or die "Input options wrong\n";

my $sttime = time;
#$posFile = basename($posFile);
#$countList = basename($countList);

###################################
# Identifying both SNPs and indels 
###################################
if ($type eq "indel") {
	print "\nCreating a comprehensive master matrix, with genotype-specific alignment summaries,for all putative variants positions (INDEL) ...\n";

	print "\nOpen posFile ... ";
	open my $POS, "<", "$posFile" || die "Can't load file $!";
	print "\nDONE! \nOpen countList ...";
	open my $LIST, "<","$countList" || die "Can't load file $!";
	print "\nDONE! \nOpen output_file ...";
	open my $DEST, ">", "$output_file" || die "Can't load file $!";
	print "\nDONE! \n";

	my %posHash;

	while (<$POS>){
		chomp;
		my @input4 = split("\t", $_);
		my $chr_pos_ref = join("\t", "$input4[0]", "$input4[1]", "$input4[2]");
		if ( $posHash{$chr_pos_ref} ) {
			next;
		} else {
			$posHash{$chr_pos_ref} = $chr_pos_ref;
			next;
		}
	}
	close $POS;

	while ( my $fileName = <$LIST> ) {
		chomp $fileName;
		open my $GENO_BASE_COUNT_FILE, "<", "$fileName" or die "Can't load file $!";

		my %genoHash;

		while ( my $line = <$GENO_BASE_COUNT_FILE>){
			my @input5 = split ("\t", $line);
			chomp @input5;
			my $chr_pos_ref = join ("\t","$input5[0]", "$input5[1]","$input5[2]");
			$genoHash{$chr_pos_ref} = $input5[3];
			@input5 = ();
		}

		foreach my $chr_pos_ref ( keys %posHash ){
			if ( $genoHash{$chr_pos_ref} ) {
				$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "$genoHash{$chr_pos_ref}");
			} else {
				$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "_,_,_,_,_,_,_,_,");
			}
		}

		%genoHash = ();
		close $GENO_BASE_COUNT_FILE;
	}
	close $LIST;

	foreach my $key ( sort {(split /\t/, $a)[0] cmp (split /\t/, $b)[0] || (split /\t/, $a)[1] <=> (split /\t/, $b)[1]} keys %posHash ) {
		print $DEST "$posHash{$key}\n";
	}
	close $DEST;
	print "DONE.\n";
	
#########################
# Identifying SNPs only
#########################

} elsif ($type eq "snp") {
	
  print "Preparing to create master matrix (SNP).\n";

	open my $POS, "<", "$posFile" or die "Can't load file $!";
	open my $LIST, "<","$countList" or die "Can't load file $!";
	open my $DEST, ">", "$output_file" or die "Can't load file $!";

	my %posHash;
	print "\nCreating a comprehensive master matrix, with genotype-specific alignment summaries,for all putative variants positions ... \n";

	while (<$POS>){
		chomp;
		my @input4 = split("\t", $_);
		my $chr_pos_ref = join("\t", "$input4[0]", "$input4[1]", "$input4[2]");
		if ( $posHash{$chr_pos_ref} ) {
			next;
		} else {
			$posHash{$chr_pos_ref} = $chr_pos_ref;
			next;
		}
	}
	close $POS;

	while ( my $fileName = <$LIST> ) {
		chomp $fileName;
		open my $GENO_BASE_COUNT_FILE, "<", "$fileName" or die "Can't load file $!";

    print "Open $GENO_BASE_COUNT_FILE";

		my %genoHash;

		while ( my $line = <$GENO_BASE_COUNT_FILE>){
			my @input5 = split ("\t", $line);
			chomp @input5;
			my $chr_pos_ref = join ("\t","$input5[0]", "$input5[1]","$input5[2]");
			$genoHash{$chr_pos_ref} = $input5[3];
			@input5 = ();
		}

		foreach my $chr_pos_ref ( keys %posHash ){
			if ( $genoHash{$chr_pos_ref} ) {
				$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "$genoHash{$chr_pos_ref}");
			} else {
				$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "_,_,_,_");
			}
		}

		%genoHash = ();
		close $GENO_BASE_COUNT_FILE;
	}
	close $LIST;

	foreach my $key ( sort {(split /\t/, $a)[0] cmp (split /\t/, $b)[0] || (split /\t/, $a)[1] <=> (split /\t/, $b)[1]} keys %posHash ) {
		print $DEST "$posHash{$key}\n";
	}
	close $DEST;
	print "DONE.\n";
}

FINAL:
exit;
