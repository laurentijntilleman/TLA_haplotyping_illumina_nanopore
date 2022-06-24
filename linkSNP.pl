#!/usr/bin/perl -w
## Takes the file with all the separate mentions of reads with SNVs
## When a read is mentioned multiple times, the SNVs in this read are reported as a link, and it's counted how often this link occurs
## Output file is a file with the linked SNVs (NT) with its POS and the count of how often.
## Didn't go through the code in detail how it's done yet..
use strict;

my $snp_read_file = shift @ARGV or die "supply list of sorted reads with SNVlinks\n";


my $prev_id = "";
my $prev_record;

my @snp_list;
my %snp_count;
open SNP, $snp_read_file or die "Cannot open file: $!";
while(<SNP>){
	chomp;
	my @data = split /\t/;
	my $read_id = $data[0];
	$read_id =~ s/#.*//;
	
	#if the read ids are mentioned multiple times - it contains multiple heterozygous SNVs -  push it on the same list
	if($read_id eq $prev_id){
		push @snp_list, $_;
	}else{
		#do something to extract the snps and print them
		if(scalar @snp_list > 1){
			my @linked_snps = parseSNPs(@snp_list);
			for my $snp_link( @linked_snps ){
				$snp_count{join("\t", @$snp_link)}++;
				@$snp_link[(2,3,0,1)] = @$snp_link;
				$snp_count{join("\t", @$snp_link)}++;
			}	

		}	
		@snp_list = ();
		push @snp_list, $_;
	}	
	$prev_id = $read_id;
}	

print join("\t", ($_, $snp_count{$_})), "\n" for keys %snp_count; #Don't delete,t he output line!

sub parseSNPs{
	my @list = @_;
	my @links;
	for my $snp1_record( @list ){
		for my $snp2_record ( @list ){
			last if $snp1_record eq $snp2_record;
			my @snp1 = split /\t/, $snp1_record;
			my @snp2 = split /\t/, $snp2_record;
			my ($read_pos1,$seq1,$snp_pos1,$snp1,$cigar1) = ($snp1[2],$snp1[4],$snp1[5],$snp1[6],$snp1[3]);
			my ($read_pos2,$seq2,$snp_pos2,$snp2,$cigar2) = ($snp2[2],$snp2[4],$snp2[5],$snp2[6],$snp2[3]);
			next if $snp_pos1 == $snp_pos2;

			my $found_snp1;
			my $found_snp2;

			if ($snp1 !~ /[-+]/){
				$found_snp1 = getSNPpos($read_pos1,$snp_pos1,$seq1,$cigar1); 
			}else{
				$found_snp1 = getINDELpos($read_pos1,$snp_pos1,$seq1,$cigar1, $snp1); 
				#Pass to similar function as SNP, for indel. Pass $snp1/$snp2 to it as well. 
				#Remember, the pos where it happens is +1!
				#Check at the end whether the SNP is identical to the indel or alternative (reference most likely). 
				#If its an 'in' (+), return the 'indel SNP' if it matches the indel and the alternative if it matches the alternative
				#If its an 'del' (-), return the alternative if it matches the indel, and the indel if it does not match 
			}
			if ($snp2 !~ /[-+]/){
				$found_snp2 = getSNPpos($read_pos2,$snp_pos2,$seq2,$cigar2);
			}else{
				$found_snp2 = getINDELpos($read_pos2,$snp_pos2,$seq2,$cigar2, $snp2);
				#Pass to similar function as SNP, for indel. Pass $snp1/$snp2 to it as well. Check whether it's the indel or the ref and return that.
			}
			if(defined $found_snp1 and defined $found_snp2 and $snp1 =~ /[$found_snp1]/ and $snp2 =~ /[$found_snp2]/){
				push @links, [$found_snp1, $snp_pos1, $found_snp2, $snp_pos2];
			}								
		}
	}
	return @links;
}	



sub getSNPpos{
	my ($read_pos,$snp_pos,$seq,$cigar) = @_; 
	my $del_pos = $read_pos;
	my $shift_seq = 0;
	#deal with insertions in the reference
	while( $cigar =~ /(\d+)M(\d+)D/g){
		$del_pos += $1;
		$shift_seq += $2 if $del_pos < $snp_pos;
	}
	#print($snp_pos);
	#print("\n");
	#print($read_pos);
	#print("\n");
	#print($shift_seq);
	#print("\n");

	
	my $pos_in_read = $snp_pos-$read_pos-$shift_seq;
	#print($pos_in_read);
	#print("\n");
	#print("\n");
	my $match_len = 0;
	while($cigar =~ /(\d+)M/g){
		$match_len += $1;
	}	
	#remove unmatched part from seq
	if($cigar=~/^(\d+)S/){
		substr($seq, 0, $1, "");
	}
	#remove the end of the sequence
	if($cigar =~ /(\d+)S^/){
		substr($seq, -$1, $1, "");
	}	
	#remove the insertions in the sequence
	while($cigar =~ /(\d+)M(\d+)I/g){
		substr($seq, $1, $2, "");
	}
	if($pos_in_read < $match_len and $match_len > 0 and $pos_in_read < length($seq) ){
		return substr($seq, $pos_in_read, 1); #This returns the NT in the read at the SNV position, fails to report +/- for indels. Also does not take larger indels into account. So right now, with indels making it to this point, this is actually faulty until I fix this properly
	}else{
		return 'N';
	}	
}	


sub getINDELpos{
	my ($read_pos,$snp_pos,$seq,$cigar,$snp_nts) = @_; 

	#split the snp in the indel and ref and draw info on in (+) or del (-), amount of NTs and which NTs
	my @snp_nt_split = split("/", $snp_nts);
	my $indel;
	my $ref_nt;
	if ($snp_nt_split[0] =~ /[-+]/){
		$indel = $snp_nt_split[0];
		$ref_nt = $snp_nt_split[1];
	}else{
		$indel = $snp_nt_split[1];
		$ref_nt = $snp_nt_split[0];
	}
	my $in_or_del;
	my $indel_length;
	my $indel_nt;
		
	if ($indel =~ /([-+])(\d*)([ACGTacgt])/){
		$in_or_del = $1;
		$indel_length = $2;
		$indel_nt = $3;
	}
	
	#correct the position for where to find the change in the sequence
	$snp_pos += 1;

	my $snp_in_read;

	if ($cigar =~ /(\d*)M(\d*)[DI]/){
		if ($read_pos+$1 == $snp_pos and $indel_length == $2){
			$snp_in_read = $indel;	
		}else{
			$snp_in_read = $ref_nt;
		}
	}else{
		$snp_in_read = $ref_nt;
	}
	return $snp_in_read;
}
