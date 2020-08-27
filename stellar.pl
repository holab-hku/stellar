#!/usr/bin/perl
use strict; use warnings;
use POSIX;

#Getting arguments from command line
my ($Fastq1, $Fastq2, $chimera_genome, $barcode_whitelist, $gtf_human, $gtf_mouse, $mouse_index, $human_index, $threshold, $base_output_dir) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4], $ARGV[5], $ARGV[6], $ARGV[7], $ARGV[8], $ARGV[9]);

#The chimera genome index directory: Mouse and Human index concatenated
my $genome_dir = $chimera_genome;#"/home/msnaveed/sra_local_repo/chimera_index/v3";

#The barcode whitelist for 10x data
#737K-august-2016.txt to be used for v2 chemistry
#3M-february-2018.txt to be used for v3 chemistry
my $white_list = $barcode_whitelist;#"/home/msnaveed/sra_local_repo/10x_genomics/3M-february-2018.txt"; 


my $species_output_dir = "";

#This is the latest human annotation.gtf file
my $human_gtf = $gtf_human;#"/home/msnaveed/sra_local_repo/chimera_genome/human_genome/v3/*.gtf";

#This is the latest mouse annotation.gtf file
my $mouse_gtf = $gtf_mouse; #"/home/msnaveed/sra_local_repo/chimera_genome/altered_mouse_genome/v3/*.gtf";

my $system_cmd;


my $datetime = localtime();  
print "Checkpoint 1 - Current date and time according to the system - $base_output_dir : $datetime\n"; 


#Creating a base output directory where all the output results would be directed
$system_cmd = system("mkdir $base_output_dir");

#species_output_dir is the path of the directory that holds the peliminary STAR results
$species_output_dir = "$base_output_dir/processFiles/";
$system_cmd = system("mkdir $species_output_dir");

#Run preliminary STAR on the raw fastq files input
$system_cmd = system("./runSTARAnalysis.sh $genome_dir $white_list $Fastq1 $Fastq2 $species_output_dir > $base_output_dir/Log_STARAnalysis_$base_output_dir.txt");
	
#It calls a perl scrip classify_barcodes_pipeline.pl that uses the preliminary STAR results to produce a labelled csv files classifying each barcode to either human, mouse or unspecified.
$system_cmd = system("perl classify_barcodes_pipeline.pl $mouse_gtf $human_gtf $base_output_dir stellar $threshold");

#It used the previously generated labelled csv file to break down the preliminary STAR's sam file and produces separate human and mouse bam files 
for my $specie (qw(processFiles)) {

	$species_output_dir = "$base_output_dir/$specie/analysis";
	
	#$currentLabelledFile is the path to labelled file
	my $currentLabelledFile = "$species_output_dir/$specie" . "_barcodes_classification.csv";

	open my $BARCODE, "<$currentLabelledFile" or die "can't open $currentLabelledFile\n";

	#Reading the labelled file and putting the barcodes in the dictionary as per their label(human/mouse/unspecified)
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

	#For the mouse barcodes in the labelled file
	if(exists $barcodes{"mouse"}){

		#BAM file produced by the preliminary STAR
		my $currentSamFile = "$species_output_dir/Aligned.sortedByCoord.out.bam";

		#New mouse BAM file
		my $prem_currentSamFile = "$species_output_dir/MouseFilteredAligned.out.sam";

		#Writing the header of the sam file to the $prem_currentSamFile
		$system_cmd = system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		
		#Opening the prem_currentSamFile for writing
		open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#Opening the currentSamFile for reading
		open my $IN, "samtools view -@ 32 $currentSamFile |" or die "can't open $currentSamFile\n";

		#Copying the barcodes present in labelled file and in the $currentSamFile to $prem_currentSamFile
		while (<$IN>) {

			my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
			
			if ($barcode && (exists $barcodes{"mouse"}{$barcode})) {
				print $OUT $_;
			}
		
		}
		close $IN;
		close $OUT;

		#Converting filtered SAM file to BAM file
		$system_cmd = system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam");
		#Removing the SAM File as we have not the binary BAM version
		$system_cmd = system("rm $prem_currentSamFile");

		#Converting the newly created filtered BAM file to fastq using the perl script convertBamtToFastq.pl
		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Fastq1 $species_output_dir/mouse_R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Fastq2 $species_output_dir/mouse_R1.fastq";

		#Running the conversion of BAM to fastq in parallel because we need to individualy convert for R1 and R2. Converting in parallel would save us a lot of time
		#Logic: Creating two child processes and letting these complete
		foreach my $i (1, 2) {
			my $pid = fork();
			if ($pid==0) { # child
				if($i==1){
					exec("$sys1"); #Exec shall exit the child process
				}
				if($i==2){
					exec("$sys2");	#Exec shall exit the child process
				}
				die "Exec $i failed: $!\n";
			} elsif (!defined $pid) {
				warn "Fork $i failed: $!\n";
			}
		}

		1 while wait() >= 0;

		#system("$sys1");
		#system("$sys2");


	}

	#For the human barcodes in the labelled file
	if(exists $barcodes{"human"}){

		#BAM file produced by the preliminary STAR
		my $currentSamFile = "$species_output_dir/Aligned.sortedByCoord.out.bam";

		#New human BAM file
		my $prem_currentSamFile = "$species_output_dir/HumanFilteredAligned.out.sam";

		#Writing the header of the sam file to the $prem_currentSamFile
		$system_cmd = system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		
		#Opening the $prem_currentSamFile for writing the human barcodes
		open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#Opening the currentSamFile for reading
		open my $IN, "samtools view -@ 32 $currentSamFile |" or die "can't open $currentSamFile\n";

		#Copying the barcodes present in labelled file and in the $currentSamFile to $prem_currentSamFile
		while (<$IN>) {

			my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
			
			if ($barcode && (exists $barcodes{"human"}{$barcode})) {
				print $OUT $_;
			}
		
		}
		close $IN;
		close $OUT;

		#Converting filtered SAM file to BAM file
		$system_cmd = system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam");
		#Removing the SAM File as we have not the binary BAM version
		$system_cmd = system("rm $prem_currentSamFile");

		#Converting the newly created filtered BAM file to fastq using the perl script convertBamtToFastq.pl
		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Fastq1 $species_output_dir/human_R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Fastq2 $species_output_dir/human_R1.fastq";
		

		#Running the conversion of BAM to fastq in parallel because we need to individualy convert for R1 and R2. Converting in parallel would save us a lot of time
		#Logic: Creating two child processes and letting these complete
		foreach my $i (1, 2) {
			my $pid = fork();
			if ($pid==0) { # child
				if($i==1){
					exec("$sys1"); #Exec shall exit the child process
				}
				if($i==2){
					exec("$sys2"); #Exec shall exit the child process
				}
				die "Exec $i failed: $!\n";
			} elsif (!defined $pid) {
				warn "Fork $i failed: $!\n";
			}
		}

		1 while wait() >= 0;

		#system("$sys1");
		#system("$sys2");


	}
	
	
}


$datetime = localtime();  
print "Checkpoint 2 - $base_output_dir : $datetime\n"; 


### THe below part is the running of final STAR after the processing done above


my $results_dir = "$base_output_dir/results";
my $exp_dir = "";

#Paths to human and mouse genome index directories
my $mouse_genome_dir = $mouse_index; #"/home/msnaveed/sra_local_repo/chimera_genome/mouse_index/v3";
my $human_genome_dir = $human_index; #"/home/msnaveed/sra_local_repo/chimera_genome/human_index/v3";

#Creating the results directory
system("mkdir $results_dir");

#Running STAR command on the human and mouse fastq individually that are produced by the processing done earlier
for my $exp (qw(h m)) {

	#Creating the h or m directory for the final STAR results
	system("mkdir $results_dir/$exp");

	$exp_dir = "$results_dir/$exp/";
	if($exp eq "h"){
		if (-e "$base_output_dir/processFiles/analysis/human_R1.fastq.gz") { #Check if the human fastq files even exist
			
			$system_cmd = system("./runFinalSTAR.sh $human_genome_dir $white_list '$base_output_dir/processFiles/analysis/human_R2.fastq.gz' '$base_output_dir/processFiles/analysis/human_R1.fastq.gz' $exp_dir > $results_dir/Human_Log_STARAnalysis.txt");
		}
		else{ #Printed when the initial input raw fastq files doesn't contain any human reads
			print ("h fastq doesn't exist\n");
		}
	}
	if($exp eq "m"){ #Check if the mouse fastq files even exist
		if (-e "$base_output_dir/processFiles/analysis/mouse_R1.fastq.gz") {
			$system_cmd = system("./runFinalSTAR.sh $mouse_genome_dir $white_list '$base_output_dir/processFiles/analysis/mouse_R2.fastq.gz' '$base_output_dir/processFiles/analysis/mouse_R1.fastq.gz' $exp_dir > $results_dir/Mouse_Log_STARAnalysis.txt");
		}
		else{ #Printed when the initial input raw fastq files doesn't contain any mouse reads
			print ("m fastq doesn't exist\n");
		}
	}

}

#Final check point that can help us determine the total run time for the script
$datetime = localtime();  
print "Finished processing - Current date and time according to the system - $base_output_dir : $datetime\n"; 