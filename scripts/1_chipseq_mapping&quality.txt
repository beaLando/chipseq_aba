# CREATE NECESSARY DIRECTORIES
mkdir ./fastq/raw
mkdir ./fastq/trimmed
mkdir ./mapped
mkdir ./qc_out

conda activate chippy

# QUALITY CHECK & CLEAN READS
fastqc ./fastq/raw/*.fq.gz -o ./qc_out -t 6 -noextract

for f1 in ./fastq/raw/*1.fq.gz
do
echo "working with file $f1"

dir="./fastq/trimmed"
f2=${f1%%1.fq.gz}"2.fq.gz"
logf=${dir}/$(basename -s .fq.gz $f1)
f1p=${dir}/$(basename -s .fq.gz $f1)_P.filtered.fq.gz
f1u=${dir}/$(basename -s .fq.gz $f1)_U.filtered.fq.gz 
f2p=${dir}/$(basename -s .fq.gz $f2)_P.filtered.fq.gz 
f2u=${dir}/$(basename -s .fq.gz $f2)_U.filtered.fq.gz 

fastp -i $f1 -I $f2 -q 20 -u 40 -l 40 --detect_adapter_for_pe --thread 6 -o $f1p -O $f2p --unpaired1 $f1u --unpaired2 $f2u -j ${logf}.json -h ${logf}.html

done

fastqc ./fastq/trimmed/*P.filtered.fq.gz -o ./qc_out -t 6 -noextract


# MAP TO REFERENCE, MARK DUPLICATES, and GENERATE BAM FILES
genome=./ara_refs/GCF_000001735.4_TAIR10.1_genomic.fna.gz
bowtie2-build -f $genome Arabidopsis_Genome

for fqgz in `ls ./fastq/trimmed/*1_P.filtered.fq.gz`  
    do

    dir="./fastq/trimmed"
    base=$(basename $fqgz "_1_P.filtered.fq.gz")
    echo "working with sample $base"

    fq1=${dir}/${base}_1_P.filtered.fq
    fq2=${dir}/${base}_2_P.filtered.fq
    sam=./mapped/${base}.sam
    bam=./mapped/${base}.bam
    bam_rg=./mapped/${base}.readgroup.bam
    bam_filt=./mapped/${base}.filtered.bam
    bam_dup=./mapped/${base}.dupmrkd.bam
    bam_dupQ=./qc_out/${base}.dupmrkdQC.txt
    bam_sorted=./mapped/${base}.sorted.bam
    bam_sorted_nodups=./mapped/${base}.noDups.sorted.bam

    # Alignment to reference genome
    gunzip -k ${dir}/${base}_1_P.filtered.fq.gz ${dir}/${base}_2_P.filtered.fq.gz
    bowtie2 -x Arabidopsis_Genome -1 $fq1 -2 $fq2 -S $sam --threads 6
    
    # Filter out reads that did not map, mapped in multiple places, had mapping quality lower than 30
    samtools sort $sam -o $bam
    picard AddOrReplaceReadGroups -I $bam -O $bam_rg --RGID rg1 --RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM $base
    samtools index $bam_rg
    samtools view --threads 3 -F 772 -q 30 -b $bam_rg NC_003070.9 NC_003071.7 NC_003074.8 NC_003075.7 NC_003076.8 | samtools sort - -o $bam_filt
    samtools index $bam_filt
    
    # Mark duplicates with picard
    picard MarkDuplicates -I $bam_filt -O $bam_dup -M $bam_dupQ --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES false --ASSUME_SORTED true
    
    # Create bam with duplicates filtered out and produce final files
    samtools sort $bam_dup -o $bam_sorted
    samtools index $bam_sorted
    samtools view --threads 3 -F 1796 -q 30 -b $bam_sorted | samtools sort - -o $bam_sorted_nodups
    samtools index $bam_sorted_nodups

    rm $fq1
    rm $fq2
    rm $bam_rg
    rm $bam_dup

    done


# PRODUCE QC REPORT for TRIMMING & MAPPING
conda deactivate chippy
conda activate mapping #conda could not install qualimap in chippy env., not sure why, but there is some sort of incompatibility, so I run it separately in another environment
DISPLAY=:0 #for no problem with qualimap

for fqgz in `ls ./fastq/trimmed/*1_P.filtered.fq.gz` #around 1min per sample 
    do

    dir="./fastq/trimmed"
    base=$(basename $fqgz "_1_P.filtered.fq.gz")
    echo "working with sample $base"

    bam=./mapped/${base}.bam
    bam_sorted=./mapped/${base}.sorted.bam
    bam_stats=./qc_out/${base}_qualimap
    bam_stats_clean=./qc_out/${base}_clean_qualimap

    qualimap bamqc -bam $bam -c -nt 6 -outdir $bam_stats
    qualimap bamqc -bam $bam_sorted -c -nt 6 -outdir $bam_stats_clean

    done

conda deactivate mapping

conda activate chippy
multiqc -ip -f ./qc_out/*_fastqc.zip ./qc_out/*.dupmrkdQC.txt #./mapping/*_qualimap ##donwload qualimap results separately
conda deactivate