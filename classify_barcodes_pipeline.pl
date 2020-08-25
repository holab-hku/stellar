#!/usr/bin/perl
use strict; use warnings;

my $MIN_COUNT = 1;

die"
Usage: output_percent_mouse_genes_per_barcode.pl <Path/to/gencode.vM24.unique_gene_names.gtf> <Path/to/gencode.v33.unique_gene_names.gtf> <base_output_dir> <type> <threshold>
" unless @ARGV==5;

### Arguments
# mouse_gtf: Mouse gtf file
# human_gtf: Human gtf file
# path: Path to the directory containing the filtered gene files
# type: In general, this determines the file structure. It could be chimera, stellar or any other value
# threshold: The value on the basis of which a particular barcode is classified as either mouse/human/unspecified. For eg if it is 80% then a barcode must have gene pct >=80% for either human/mouse for it to be classified otherwise it'll be regarded as 'unspecified'
my ($mouse_gtf, $human_gtf, $path, $type, $threshold) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]);

#Checking if correct files have been input
if ($mouse_gtf !~ /.*gencode.vM24.unique_gene_names.gtf/) {
	die "mouse gtf isn't gencode.vM24.unique_gene_names.gtf\n";
} elsif ($human_gtf !~ /.*gencode.v33.unique_gene_names.gtf/) {
	die "human gtf isn't gencode.v33.unique_gene_names.gtf\n";
}

#Reading the gene names from gtf files
my %all_genes;
#get whole set of mouse genes and human genes
for my $species (qw(mouse human)) {
	my $IN;
	if ($species eq "mouse") {
		open $IN, "<$mouse_gtf" or die "can't open $mouse_gtf\n";
	} elsif ($species eq "human") {
		open $IN, "<$human_gtf" or die "can't open $human_gtf\n";
	}
	while (<$IN>) {
		if ($_ =~ /^#/) {
			next;
		}
	
		my $gene_name;
		if ($_ =~ /gene_name \"(.*?)\"/) {
			$gene_name = $1;
		}
		$all_genes{$species}{$gene_name}++;
	}
	
	close $IN;
	
}

### The following three variants are supported:
# chimera: Goes through all variants of experiment 5s(mix with different proportions of human & mouse) and outputs a label file for each exp
# stellar: Goes through one directory processFiles
# <else>: Goes through human and mouse directory
#
# In general, it takes the matrix.csv file from the STAR results and output a labelled csv file. The three types above are designed with respect to the directory structure


if($type eq "chimera"){
	for my $exp (qw(5.1 5.2 5.3)) {
		my $exp_matrix = "$path/$exp/analysis/matrix.csv";
		# my $mouse_output_file = "exp$exp/exp$exp" . "_mouse_barcodes.csv";
		# my $human_output_file = "exp$exp/exp$exp" . "_human_barcodes.csv";
		my $output_file = "$path/$exp/analysis/" . $exp . "_barcodes_classification.csv";
		
		open my $IN, "<$exp_matrix" or die "can't open $exp_matrix\n";
		my $header = <$IN>;
		my @header_ar = split(",", $header);
		chomp(@header_ar);
		
		my %count_expressed_barcodes_gene;
		
		while (<$IN>) {
			chomp;
			my @line = split(",", $_);
			my $gene = $line[0];
			for (my $i=1; $i < @line; $i++) {
				if ($line[$i] >= $MIN_COUNT) {
					for my $species (qw(mouse human)) {
						if (exists $all_genes{$species}{$gene}) {
							$count_expressed_barcodes_gene{$header_ar[$i]}{$species} += $line[$i];
						}
					}
				}
			}
		}
		close $IN;
		
		# open my $HUMAN, ">$human_output_file" or die "can't open $human_output_file\n";
		# print $HUMAN "Barcode\n";

		# open my $MOUSE, ">$mouse_output_file" or die "can't open $mouse_output_file\n";
		# print $MOUSE "Barcode\n";

		open my $OUT, ">$output_file" or die "can't open $output_file\n";
		print $OUT "Barcode,Species\n";

		for my $barcode (keys %count_expressed_barcodes_gene) {
			for my $species (qw(mouse human)) {
				if (!exists $count_expressed_barcodes_gene{$barcode}{$species}) {
					$count_expressed_barcodes_gene{$barcode}{$species} = 0;
				}
			}
			
			my $tot_genes;
			for my $species (qw(mouse human)) {
				$tot_genes += $count_expressed_barcodes_gene{$barcode}{$species};
			}
			
			my $mouse_genes_pct = ($count_expressed_barcodes_gene{$barcode}{"mouse"} / $tot_genes) * 100;
			if($mouse_genes_pct>=$threshold){
				print $OUT $barcode . ",mouse\n";
			}
			elsif($mouse_genes_pct<=100-$threshold){
				print $OUT $barcode . ",human\n";
			}
			else{
				print $OUT $barcode . ",unspecified\n";
			}
			
		}
		
		close $OUT;
		# close $MOUSE;
		# close $HUMAN;

	}


}
elsif($type eq "stellar"){
	for my $exp (qw(processFiles)) {
		my $exp_matrix = "$path/$exp/analysis/matrix.csv";
		# my $mouse_output_file = "exp$exp/exp$exp" . "_mouse_barcodes.csv";
		# my $human_output_file = "exp$exp/exp$exp" . "_human_barcodes.csv";
		my $output_file = "$path/$exp/analysis/" . $exp . "_barcodes_classification.csv";
		
		open my $IN, "<$exp_matrix" or die "can't open $exp_matrix\n";
		my $header = <$IN>;
		my @header_ar = split(",", $header);
		chomp(@header_ar);
		
		my %count_expressed_barcodes_gene;
		
		while (<$IN>) {
			chomp;
			my @line = split(",", $_);
			my $gene = $line[0];
			for (my $i=1; $i < @line; $i++) {
				if ($line[$i] >= $MIN_COUNT) {
					for my $species (qw(mouse human)) {
						if (exists $all_genes{$species}{$gene}) {
							$count_expressed_barcodes_gene{$header_ar[$i]}{$species} += $line[$i];
						}
					}
				}
			}
		}
		close $IN;
		
		# open my $HUMAN, ">$human_output_file" or die "can't open $human_output_file\n";
		# print $HUMAN "Barcode\n";

		# open my $MOUSE, ">$mouse_output_file" or die "can't open $mouse_output_file\n";
		# print $MOUSE "Barcode\n";

		open my $OUT, ">$output_file" or die "can't open $output_file\n";
		print $OUT "Barcode,Species\n";

		for my $barcode (keys %count_expressed_barcodes_gene) {
			for my $species (qw(mouse human)) {
				if (!exists $count_expressed_barcodes_gene{$barcode}{$species}) {
					$count_expressed_barcodes_gene{$barcode}{$species} = 0;
				}
			}
			
			my $tot_genes;
			for my $species (qw(mouse human)) {
				$tot_genes += $count_expressed_barcodes_gene{$barcode}{$species};
			}
			
			my $mouse_genes_pct = ($count_expressed_barcodes_gene{$barcode}{"mouse"} / $tot_genes) * 100;
			if($mouse_genes_pct>=$threshold){
				print $OUT $barcode . ",mouse\n";
			}
			elsif($mouse_genes_pct<=100-$threshold){
				print $OUT $barcode . ",human\n";
			}
			else{
				print $OUT $barcode . ",unspecified\n";
			}
			
		}
		
		close $OUT;
		# close $MOUSE;
		# close $HUMAN;

	}
}
else{
	for my $exp (qw(human mouse)) {
		my $exp_matrix = "$path/$exp/analysis/matrix.csv";
		# my $mouse_output_file = "exp$exp/exp$exp" . "_mouse_barcodes.csv";
		# my $human_output_file = "exp$exp/exp$exp" . "_human_barcodes.csv";
		my $output_file = "$path/$exp/analysis/" . $exp . "_barcodes_classification.csv";
		
		open my $IN, "<$exp_matrix" or die "can't open $exp_matrix\n";
		my $header = <$IN>;
		my @header_ar = split(",", $header);
		chomp(@header_ar);
		
		my %count_expressed_barcodes_gene;
		
		while (<$IN>) {
			chomp;
			my @line = split(",", $_);
			my $gene = $line[0];
			for (my $i=1; $i < @line; $i++) {
				if ($line[$i] >= $MIN_COUNT) {
					for my $species (qw(mouse human)) {
						if (exists $all_genes{$species}{$gene}) {
							$count_expressed_barcodes_gene{$header_ar[$i]}{$species} += $line[$i];
						}
					}
				}
			}
		}
		close $IN;
		
		# open my $HUMAN, ">$human_output_file" or die "can't open $human_output_file\n";
		# print $HUMAN "Barcode\n";

		# open my $MOUSE, ">$mouse_output_file" or die "can't open $mouse_output_file\n";
		# print $MOUSE "Barcode\n";

		open my $OUT, ">$output_file" or die "can't open $output_file\n";
		print $OUT "Barcode,Species\n";

		for my $barcode (keys %count_expressed_barcodes_gene) {
			for my $species (qw(mouse human)) {
				if (!exists $count_expressed_barcodes_gene{$barcode}{$species}) {
					$count_expressed_barcodes_gene{$barcode}{$species} = 0;
				}
			}
			
			my $tot_genes;
			for my $species (qw(mouse human)) {
				$tot_genes += $count_expressed_barcodes_gene{$barcode}{$species};
			}
			
			my $mouse_genes_pct = ($count_expressed_barcodes_gene{$barcode}{"mouse"} / $tot_genes) * 100;
			if($mouse_genes_pct>=$threshold){
				print $OUT $barcode . ",mouse\n";
			}
			elsif($mouse_genes_pct<=100-$threshold){
				print $OUT $barcode . ",human\n";
			}
			else{
				print $OUT $barcode . ",unspecified\n";
			}
			
		}
		
		close $OUT;
		# close $MOUSE;
		# close $HUMAN;

	}

}
