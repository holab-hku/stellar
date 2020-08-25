#!/usr/bin/perl
use strict; use warnings;
use POSIX;
#use List::Util 'shuffle';

 

my $genome_dir = "/home/msnaveed/sra_local_repo/chimera_index/v3";
my $white_list = "/home/msnaveed/sra_local_repo/10x_genomics/3M-february-2018.txt";
#737K-august-2016.txt

my $species_output_dir = "";

my $human_gtf = "/home/msnaveed/sra_local_repo/chimera_genome/human_genome/v3/*.gtf";
my $mouse_gtf = "/home/msnaveed/sra_local_repo/chimera_genome/altered_mouse_genome/v3/*.gtf";

my $system_cmd;

my ($Fastq1, $Fastq2, $base_output_dir) = ($ARGV[0], $ARGV[1], $ARGV[2]);

my $datetime = localtime();  
print "Current date and time according to the system - $base_output_dir : $datetime\n"; 
### Comment from here


$system_cmd = system("mkdir $base_output_dir");


$species_output_dir = "$base_output_dir/processFiles/";
$system_cmd = system("mkdir $species_output_dir");

$system_cmd = system("./runSTARAnalysis.sh $genome_dir $white_list $Fastq1 $Fastq2 $species_output_dir > Log_STARAnalysis_$base_output_dir.txt");
	

$system_cmd = system("perl classify_barcodes_pipeline.pl $mouse_gtf $human_gtf $base_output_dir stellar");


for my $specie (qw(processFiles)) {
	
	$species_output_dir = "$base_output_dir/$specie/analysis";
	my $currentLabelledFile = "$species_output_dir/$specie" . "_barcodes_classification.csv";

	open my $BARCODE, "<$currentLabelledFile" or die "can't open $currentLabelledFile\n";

	my %barcodes;
	my $header = <$BARCODE>;
	my @header_ar = split(",", $header);
	chomp(@header_ar);

	while (<$BARCODE>) {
		chomp;
		my @line = split(",", $_);
		my $barcode = $line[0];
		my $species = $line[1];
		$barcodes{$species}{$barcode}+=1;

	}

	if(exists $barcodes{"mouse"}){


		my $currentSamFile = "$species_output_dir/Aligned.sortedByCoord.out.bam";

		my $prem_currentSamFile = "$species_output_dir/MouseFilteredAligned.out.sam";

		#Writing the header of the sam file to the filteredSamFile
		$system_cmd = system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		
		#Opening the filteredSamFile for writing
		open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#Opening the currentSamFile for reading
		open my $IN, "samtools view -@ 32 $currentSamFile |" or die "can't open $currentSamFile\n";

		while (<$IN>) {

			my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
			
			if ($barcode && (exists $barcodes{"mouse"}{$barcode})) {
				print $OUT $_;
			}
		
		}
		close $IN;
		close $OUT;

		$system_cmd = system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam");
		$system_cmd = system("rm $prem_currentSamFile");


		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Fastq1 $species_output_dir/mouse_R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Fastq2 $species_output_dir/mouse_R1.fastq";
		system("$sys1");
		system("$sys2");


	}

	if(exists $barcodes{"human"}){

		my $currentSamFile = "$species_output_dir/Aligned.sortedByCoord.out.bam";

		my $prem_currentSamFile = "$species_output_dir/HumanFilteredAligned.out.sam";

		#Writing the header of the sam file to the filteredSamFile
		$system_cmd = system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		
		#Opening the filteredSamFile for writing
		open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#Opening the currentSamFile for reading
		open my $IN, "samtools view -@ 32 $currentSamFile |" or die "can't open $currentSamFile\n";

		while (<$IN>) {

			my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
			
			if ($barcode && (exists $barcodes{"human"}{$barcode})) {
				print $OUT $_;
			}
		
		}
		close $IN;
		close $OUT;

		$system_cmd = system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam");
		$system_cmd = system("rm $prem_currentSamFile");

		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Fastq1 $species_output_dir/human_R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Fastq2 $species_output_dir/human_R1.fastq";
		#$system_cmd = system("$sys1" . "&" . "$sys2");

		system("$sys1");
		system("$sys2");


	}
	
	
}


$datetime = localtime();  
print "Running final Star - $base_output_dir : $datetime\n"; 

my $results_dir = "$base_output_dir/results";
my $exp_dir = "";
my $mouse_genome_dir = "/home/msnaveed/sra_local_repo/chimera_genome/mouse_index/v3";
my $human_genome_dir = "/home/msnaveed/sra_local_repo/chimera_genome/human_index/v3";

system("mkdir $results_dir");



for my $exp (qw(h m)) {
	system("mkdir $results_dir/$exp");
	$exp_dir = "$results_dir/$exp/";
	if($exp eq "h"){
		if (-e "$base_output_dir/processFiles/analysis/human_R1.fastq.gz") {
			
			$system_cmd = system("./runFinalSTAR.sh $human_genome_dir $white_list '$base_output_dir/processFiles/analysis/human_R2.fastq.gz' '$base_output_dir/processFiles/analysis/human_R1.fastq.gz' $exp_dir > Human_Log_STARAnalysis.txt");
		}
		else{
			print ("h fastq doesn't exist\n");
		}
	}
	if($exp eq "m"){
		if (-e "$base_output_dir/processFiles/analysis/mouse_R1.fastq.gz") {
			$system_cmd = system("./runFinalSTAR.sh $mouse_genome_dir $white_list '$base_output_dir/processFiles/analysis/mouse_R2.fastq.gz' '$base_output_dir/processFiles/analysis/mouse_R1.fastq.gz' $exp_dir > Mouse_Log_STARAnalysis.txt");
		}
		else{
			print ("m fastq doesn't exist\n");
		}
	}

}


$datetime = localtime();  
print "Current date and time according to the system - $base_output_dir : $datetime\n"; 