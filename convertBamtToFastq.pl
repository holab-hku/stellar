#!/usr/bin/perl
use strict; use warnings;

#The number of threads to be used. The more the faster
my $threads = 16;

### Takes the following:
# BAM: The bam file produced by STAR
# fastq: The original fastq file that was used to generate the STAR BAM file
# outputFastq: The path+name of the fastq file that should be produced
my ($BAM, $fastq, $outputFastq) = ($ARGV[0], $ARGV[1], $ARGV[2]);

print("Inside Bam to fastq converter\n");

#Opening the given BAM file
open my $IN, "samtools view -@ $threads $BAM |" or die "can't open $BAM\n";

my %read_ids;

#Getting the read ids of the reads in the BAM file and putting these into the %read_ids dictionary 
my $read = "";
while (<$IN>) {
	my @line = split("\t", $_);
	$read = "@" . $line[0];
	$read_ids{$read}++;
	
}
close $IN;

#Opening the output_fastq which shall be the output of this script
open my $OUT, ">$outputFastq" or die "can't open $outputFastq\n";

#Reading the given original fastq file
open my $READ, "gunzip -c $fastq |" or die "can't open $fastq\n";

my $count = 0;
my $getline = 0;

#Copying the reads from original fastq to the output_fastq depending on their presence inside the BAM file
#Also, implements the logic to skip/place 4 lines in fast as each read is represented as 4 lines in fastq
while (my $line = <$READ>) {
	if($getline == -1){ #Skip the next 3 lines
		if($count == 2){
			$count = 0;
			$getline = 0;
		}
		else{
			$count+=1;
		}
	}
	elsif($getline == 1){ #Print the next 3 lines
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
		my @reads = split(' ',$line); #Splitting the line in fastq by space as the delimiter
		if(exists($read_ids{$reads[0]})){ #See if this read id is present in the provided BAM file. If yes, print this and the next 3 lines to output_fastq
			$getline = 1;
			print $OUT $line;
		}
		else{ #Skip the next three lines and check the next read in the fastq
			$getline = -1;
		}
	}
	
}
close $READ;
close $OUT;

#Zipping the fastq file that is just produced
system("gzip $outputFastq");
