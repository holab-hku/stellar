# Stellar

Stellar is a software written in perl that takes in raw fastq files and uses STAR aligner to divide the mixed human+mouse reads in the raw fastq files into two individual human and mouse sets while finally aligning these inidividual sets using STAR. Currently, it is only supported on 10X genomics data.

## Usage?
<ol>
  <li> <strong>Bring prelimiary files from the storage directory to your personal account:</strong> <code>cp -r /storage/holab/datamart {your personalised path}</code></li>
  <li> <strong>Put <i>STARsolo</i> and <i>samtools</i> to path:</strong> It is recommended to fetch the latest version of these two softwares to get the enhanced experience but a copy of executables is present in the storage directory:
    <code>/storage/holab/datamart/STAR</code> & <code>/storage/holab/datamart/samtools-1.10</code><br/><br/>
    <i>To put these into path:</i><br/>
    - Go to the executable directory and then do <code>export PATH="$PATH:`pwd`"</code><br/>
    - Test if it was succesfully added to the PATH variable: <code>STAR --version</code>
   
  </li>
  
  <li> <strong>Create a directory and clone this repository:</strong> <code>git clone </code> </li>
  <li> <strong>Running the main script:</strong></li><br>
</ol>

  ```ruby
  perl stellar.pl {R2 fastq file} {R1 fastq file} {path_to_chimera_index_files} {path_to_whitelist_barcode} {path_to_human_gtf} {path_to_mouse_gtf} {path_to_mouse_index_files} {path_to_human_index_files} {threshold} {output_directory} 

  Required:
    -{R2 fastq file}: Path to the R2 fastq file for 10x data
    -{R1 fastq file}: Path to the R1 fastq file for 10x data
    -{path_to_chimera_index_files}: Path to the chimera genome index - human + mouse mixed
    -{path_to_whitelist_barcode}: Path to the barcode whitelist file which is needed for 10X data. Depends on the chemistry(v2/v3) used to produce the 10X data
    -{path_to_human_gtf}: Path to human annotation gtf file
    -{path_to_mouse_gtf}: Path to mouse annotation gtf file
    -{path_to_mouse_index_files}: Path to mouse genome index file - only mouse
    -{path_to_human_index_files}: Path to human genome index file - only human
    -{threshold}: The threshold on the basis of which a particular barcode is classified as either mouse/human/unspecified. For eg if it is 80% then a barcode must have gene pct >=80% for either human/mouse for it to be classified otherwise it'll be regarded as 'unspecified' 
    -{output_directory}: The path of the directory where all the results of the current instance are put.
      
  ```
  
  ### Result Structure:
  This should be the tree view of the resulting directories/files after succesffuly running this script.
  ```bash
    {output_directory}
      |
      |__ processFiles
      |
      |__ results
      |    |
      |    |__ h: STAR results for the human reads present in the raw data input
      |    |
      |    |__ m: STAR results for the mouse reads present in the raw data input
      |
      |__ Log_STARAnalysis_{output_directory}.txt: The log file of the initial STAR run
  ```
  
  ### Example Run
  ```ruby
  nohup perl stellar.pl ../fastq/L44TX3-1_S37_L004_R2_001.fastq.gz ../fastq/L44TX3-1_S37_L004_R1_001.fastq.gz /home/msnaveed/sra_local_repo/chimera_index/v3 /home/msnaveed/sra_local_repo/10x_genomics/3M-february-2018.txt /home/msnaveed/sra_local_repo/chimera_genome/human_genome/v3/*.gtf /home/msnaveed/sra_local_repo/chimera_genome/altered_mouse_genome/v3/*.gtf /home/msnaveed/sra_local_repo/chimera_genome/mouse_index/v3 /home/msnaveed/sra_local_repo/chimera_genome/human_index/v3 80 L44TX3-1_S37 &
  ```
  #### Note
  - Make sure Rscript and Perl is also installed and is already on the PATH variable.
  - Use <code>nohup ... &</code> so that the command doesn't break upon logging off and continue to run in the backgroun even after signing off.
  - Within each script, code is commented for convenience
  - To learn on how to create genome index files visit the holab page[https://holab-hku.github.io/10X-workshop/processing-10x-rna-seq-data.html]

## File Hierarchy

### Stellar Files
  - **stellar.pl**: The main script that contains the logic of running stellar
  - **convert_to_genes_cells_matrices.r**: Converts dense matrix to sparse using R package Seurat. The output file is matrix.csv
  - **runSTARAnalysis.sh**: Runs STAR command and then runs convert_to_genes_cells_matrices.r on it to produce a matrix.csv
  - **classify_barcodes_pipeline.pl**: This script uses matrix.cvs to create a labelled file classifying each barcode in the raw input data as either human/mouse or unspecified
  - **convertBamtToFastq.pl**: Converts a bam file to fastq file. Need to run this individually for each R1 and R2 for 10X data
  - **runFinalSTAR.sh**: Runs a STAR command individually on each human and mouse dataset produced from original dataset
  
    <i>Note: To know more about the arguments of these files take, reads the header comments within each file.</i>
  
### DataMart Files
Present in the holab storage directory: <code>/storage/holab/datamart</code>
  - **gencode.v33.unique_gene_names.gtf**: Human gtf fie with unique gene names
  - **gencode.vM24.unique_gene_names.gtf**: Mouse gtf with unique gene names
  - **gencode.v33.vM24.concat.unique_gene_names.gtf**: Concatenated human + mouse gtf file with unique gene names
  - **human_index**: Indexed files generated using human genome
  - **mouse_index**: Indexed files generated using mouse genome
  - **chimera_index**: Indexed files generated using human + mouse genome. Concatenated gtf file and two fasta files were used
  - **STAR**: Executable directory of STAR version 2.7.3a
  - **samtools-1.10**: Executable directory of samtools version samtools 1.10

## Time Complexity
The time it takes to run the script is directly related to the size of the fastq files fed. On running the script on R1(7.2GB) and R2(7.7GB) fastq files, it took aronund 4 hours to complete the whole process.


## Caution
This stellar version was tested on v3 chemistry 10X dataset If there's any error thrown by STAR command, it is likely that some parameters might need to be changed within the command. These changes should be made in the STAR commands present in runSTARAnalysis.sh & runFinalSTAR.sh.
