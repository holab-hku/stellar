#!/usr/bin/perl
use strict; use warnings;

my $threads = 16;

my ($BAM, $fastq, $outputFastq) = ($ARGV[0], $ARGV[1], $ARGV[2]);

print("Inside Bam to fastq converter\n");

open my $IN, "samtools view -@ $threads $BAM |" or die "can't open $BAM\n";

my %read_ids;

my $read = "";
while (<$IN>) {
	my @line = split("\t", $_);
	$read = "@" . $line[0];
	$read_ids{$read}++;
	
}
close $IN;

open my $OUT, ">$outputFastq" or die "can't open $outputFastq\n";


open my $READ, "gunzip -c $fastq |" or die "can't open $fastq\n";
my $count = 0;
my $getline = 0;
while (my $line = <$READ>) {
	if($getline == -1){
		if($count == 2){
			$count = 0;
			$getline = 0;
		}
		else{
			$count+=1;
		}
	}
	elsif($getline == 1){
		print $OUT $line;
		if($count == 2){
			$count = 0;
			$getline = 0;
		}
		else{
			$count+=1;
		}
	}
	else{
		my @reads = split(' ',$line);
		if(exists($read_ids{$reads[0]})){
			$getline = 1;
			print $OUT $line;
		}
		else{
			$getline = -1;
		}
	}
	
}
close $READ;
close $OUT;
system("gzip $outputFastq");
