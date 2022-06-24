#!/usr/bin/perl -w

##SNV caller which generates file.pgsnp based on mpileup data on locs of interest with following:
##		Chr
##		Pos
##		Pos+1 (why?)
##		Mogelijke NTs op deze pos
##		Hoeveelheid NTs
##		Reads per NT called
##		Allele_score (just a ‘0’ for every NT variant)

use strict;

my $pileup_file = shift @ARGV;
my $vp = shift @ARGV;
$vp //= 0;

my $minimal_count = 0.15; #SNV frequency threshold
my $minimal_coverage = 25; #Minimal coverage at pos to determine SNV

open FILE, $pileup_file or die "Cannot open file: $!";


my @records;
my ($min,$max,$chrom_store);

while(<FILE>){
	chomp;
	my @data = split /\t/;
	my ($chrom,$pos,$nuc) = ($data[0], $data[1],$data[2]);
	my $i = 3;
	my ($coverage, $string);
	while(defined $data[$i]){
		$coverage += $data[$i];
		if($coverage>0){
			$string .= $data[$i+1];
		}else{
			$string = ""
		}
		$i+=3;
	}
	
	$min = $pos if not defined $min or $pos < $min;
	$max = $pos if not defined $max or $pos > $max;
	$chrom_store = $chrom;
	$nuc = uc($nuc);
	my @indels;
	@indels = identifyIndels($string);
	@indels = (); #empty indels to not do them; see if this makes the data better.
	$string = removeIndels($string);
	
	
	my $indelslength = @indels;
	my $indelcount = 0;
	

	#some strings contain stars which are a place holder for something and count in the coverage
	#but we cannot do anything with them so we remove them
	my $star_count=0;
	$star_count++ while /\*/g;
	$coverage -= $star_count;
	
	#If SNV, determine some values for thresholds 
	if($string =~ /[ACGTacgt]/ or $indelslength>0){

		my $cnt = 0;
		my ($fw,$rev) = (0,0);
		$fw++ while $string =~ /[ACGT]/g;
		$rev++ while $string =~ /[acgt]/g;
		foreach (@indels){
			if ($_ =~ /[-+]\d*([ACGT])*/){
				$fw +=1;
			}elsif($_ =~ /[-+]\d*([acgt])*/){
				$rev +=1;
			}
		}
		my $strandedness = 0;
		if($fw > $rev){
			$strandedness = $rev/($fw+$rev);
		}else{
			$strandedness = $fw/($rev+$fw);
		}	
		$cnt = $fw + $rev;
	
		next if $coverage == 0 or $strandedness < 0;
		$cnt /= $coverage;

	
		#if SNV is over threshold; count the amount of times each NT in the string occurs and store in dictionary 
		if(($cnt >= $minimal_count and $cnt <= 1-$minimal_count and $coverage >= $minimal_coverage) or $pos == $vp){ #change to VP
			my %count_nt;

			$count_nt{uc($1)}++ while $string =~ /([ACGTacgt])/g;
			foreach(@indels){
				$count_nt{uc($_)}++
			}
			
			my $ref_cnt = 0;
			$ref_cnt++ while $string =~ /[.,]/g;
			$count_nt{$nuc} = $ref_cnt;

			#when terminal reads are too abundant, remove them 
			my ($terminal_ref,%terminal_nt);
			$terminal_nt{uc($1)}++ while $string =~ /([ACGTacgt])\$/g;
			$terminal_ref++ while $string =~ /[.,]\$/g;
			$terminal_nt{$nuc} = $terminal_ref;
			for(keys %terminal_nt){
				if(defined $terminal_nt{$_} and $terminal_nt{$_}/$count_nt{$_} > 0.2){
					$count_nt{$_} -= $terminal_nt{$_};
					$coverage -= $terminal_nt{$_};
				}
			}	

			#create the records for in the snp file
			my @nuc = keys %count_nt;
			my @high_nuc;
			my @count_array;
			for(@nuc){
				#minimal coverage > 5% to be included in the data (include VP at all times, just not if the coverage is below 1%)
				if(($count_nt{$_}/$coverage > $minimal_count and $count_nt{$_}/$coverage < 1-$minimal_count) or ($pos == $vp and $count_nt{$_}/$coverage > 0.01)){
					push @high_nuc, $_;
					push @count_array, $count_nt{$_};
				}
			}

			if((scalar @high_nuc == 1 and $high_nuc[0] eq $nuc) and $pos != $vp){
				next;
			}	

			my $nuc_entry = join("/",@high_nuc);
			
			my $count_entry = join(",", @count_array);
			
			my @zeroes;
			my $allele_count = scalar @high_nuc;
			push @zeroes, 0 for(1..$allele_count);
			my $allele_score = join(",", @zeroes);
			push @records, join("\t", ($chrom,$pos,$pos+1, $nuc_entry, $allele_count, $count_entry, $allele_score));

		}	
	}	
}

print "browser position $chrom_store:$min-$max\n"; 

print join("\n", @records), "\n";


sub identifyIndels{
	my $string = shift;
	my @indels;
	while( $string =~ /[+-](\d+)/g ){
		my $len = $1;
		$string =~ s/([+-]\d+.{$len})//;
		push(@indels, $1);
	}
	return @indels;
}	

sub removeIndels{
	my $string = shift;
	my @indels;
	while( $string =~ /[+-](\d+)/g ){
		my $len = $1;
		$string =~ s/([+-]\d+.{$len})//;
		push(@indels, $1);
	}
	return $string;
	
}	

