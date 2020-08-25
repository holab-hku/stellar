# example run: ./runFinalSTAR.sh $genome_dir $white_list 'R2.fastq.gz' 'R1.fastq.gz' $exp_dir

STAR --genomeDir $1 --soloType Droplet --soloCBwhitelist $2 --runThreadN 16 --readFilesCommand zcat --readFilesIn $3 $4 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMattrRGline ID:$3 --soloBarcodeReadLength 0 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $5
