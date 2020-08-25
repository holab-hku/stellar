STAR --genomeDir $1 --soloType Droplet --soloCBwhitelist $2 --runThreadN 16 --readFilesCommand zcat --readFilesIn $3 $4 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMattrRGline ID:$3 --soloBarcodeReadLength 0 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $5
mkdir ${5}analysis
mv ${5}Solo.out/Gene/filtered/* ${5}analysis
gzip -r ${5}analysis/*
Rscript convert_to_genes_cells_matrices.r ${5}analysis
mv ${5}Aligned.sortedByCoord.out.bam ${5}analysis
