conda activate chippy

# 1. get down-sampled reads into bedgraph format
mkdir mapped/downsampled

genome="ara_refs/TAIR10_chr_count.txt"
working_dir="mapped/pooled"
normalizeby=10000000 # scaling factor

while read line; do
    
    replicate=$(echo $line | cut -d ";" -f1)
    samp=$(echo $line | cut -d ";" -f4)
    echo $replicate
    echo $samp

    readSize=$(echo $line | cut -d ";" -f2)
    fragSize=$(echo $line | cut -d ";" -f3)
    extend=`awk -v f=${fragSize} -v r=${readSize} 'BEGIN{print f-r}'`
    orig_bam="mapped/${replicate}.sorted.bam"
    
    if [[ ! -f ${working_dir}/${samp}.bigwig ]]; then
        # convert BAM to BED format
        bamToBed \
            -i ${orig_bam} > ${working_dir}/${replicate}.bed

        # sort the BED file
        sort -k 1,1 ${working_dir}/${replicate}.bed > ${working_dir}/${replicate}.sorted.bed && rm ${working_dir}/${replicate}.bed

        # extend the reads
        slopBed -i ${working_dir}/${replicate}.sorted.bed -l 0 -r ${extend} -s -g ${genome} > ${working_dir}/${replicate}.extend.bed && rm ${working_dir}/${replicate}.sorted.bed

        # normalize signal and output BEDGRAPH
        totreads=`samtools view -c ${orig_bam}`
        scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

        genomeCoverageBed -i ${working_dir}/${replicate}.extend.bed -g ${genome} -bg -scale $scaling | awk 'BEGIN{OFS=FS="\t"} {$4=sprintf("%.2f",$4)}{print}' > ${working_dir}/${replicate}.bg && rm ${working_dir}/${replicate}.extend.bed
    
    fi

done < ./qc_out/chipqc_all/chip_readsize_fragsize.csv


# 2. return average signal of two reps
dos2unix ./qc_out/chipqc_all/chip_readsize_fragsize.csv
my_array=( $(cut -d ";" -f4 ./qc_out/chipqc_all/chip_readsize_fragsize.csv | sort -u) )

for samp in "${my_array[@]}"; do
    echo "starting averaging for genotype ${samp}"
    if [[ ! -f ${working_dir}/${samp}.bigwig ]]; then

        rep1=$(awk -F ";" -v pat="$samp""$" '$4~pat {print $1}' ./qc_out/chipqc_all/chip_readsize_fragsize.csv | awk 'NR==1{print $0}')
        rep2=$(awk -F ";" -v pat="$samp""$" '$4~pat {print $1}' ./qc_out/chipqc_all/chip_readsize_fragsize.csv | awk 'NR==2{print $0}')

        echo "replicate1 is ${rep1}"
        echo "replicate2 is ${rep2}"

        bedtools unionbedg -i ${working_dir}/${rep1}.bg ${working_dir}/${rep2}.bg > ${working_dir}/${samp}.unionbg && rm ${working_dir}/${rep1}.bg ${working_dir}/${rep2}.bg

        awk 'BEGIN{OFS="\t"} $1!="ChrC" && $1!="ChrM"{ avg=($4+$5)/2; print $1,$2,$3,avg}' ${working_dir}/${samp}.unionbg > ${working_dir}/${samp}.bedgraph && rm ${working_dir}/${samp}.unionbg

        bedGraphToBigWig ${working_dir}/${samp}.bedgraph ${genome} ${working_dir}/${samp}.bigwig && rm ${working_dir}/${samp}.bedgraph

    fi

done

conda deactivate chippy