conda activate chippy

genome="ara_refs/TAIR10_chr_count.txt"
normalizeby=10000000 # scaling factor

while read line; do

    samp=$(echo $line | cut -d ";" -f1)
    readSize=$(echo $line | cut -d ";" -f2)
    fragSize=$(echo $line | cut -d ";" -f3)
    extend=`awk -v f=${fragSize} -v r=${readSize} 'BEGIN{print f-r}'`
    orig_bam="./mapped/${samp}.sorted.bam"
    # convert BAM to BED format
    bamToBed -i ${orig_bam} > ./mapped/${samp}.bed

    # sort the BED file
    sort -k 1,1 ./mapped/${samp}.bed > ./mapped/${samp}.sorted.bed

    # extend the reads
    slopBed -i ./mapped/${samp}.sorted.bed -l 0 -r ${extend} -s -g ${genome} > ./mapped/${samp}.extend.bed

    # normalize signal and output BEDGRAPH
    totreads=`samtools view -c ${orig_bam}`
    scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

    genomeCoverageBed -i ./mapped/${samp}.extend.bed -g ${genome} -bg -scale $scaling | awk 'BEGIN{OFS=FS="\t"} {$4=sprintf("%.2f",$4)}{print}' > ./mapped/${samp}.bg

    # compress BEDGRAPH to BIGWIG FORMAT
    bedGraphToBigWig ./mapped/${samp}.bg ${genome} ./mapped/${samp}.bw

done < ./qc_out/chipqc_all/chip_readsize_fragsize.csv

conda deactivate