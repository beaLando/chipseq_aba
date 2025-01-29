#!/bin/bash
greenscreen_regions="ara_refs/arabidopsis_greenscreen_20inputs.bed"

# 1. mask reads that overlap Greenscreen
# create output directory for reads
# that do and do not overlap
# greenscreen regions
mkdir -p mapped/greenscreen_mask
mkdir -p mapped/greenscreen_regions

while read line; do
    samp=`echo $line | cut -d "," -f1`
    if [[ "$samp" != "SampleID" ]]; then
       index_bam="mapped/greenscreen_mask/${samp}.sorted.bam.bai"
	if [[ ! -f "$index_bam" ]]; then
            bedtools intersect -ubam -v -a mapped/${samp}.sorted.bam -b ara_refs/arabidopsis_greenscreen_20inputs.bed > mapped/greenscreen_regions/${samp}.sorted.bam
            bedtools intersect -ubam -a mapped/${samp}.sorted.bam -b ara_refs/arabidopsis_greenscreen_20inputs.bed > mapped/greenscreen_mask/${samp}.sorted.bam 
    	    samtools index mapped/greenscreen_mask/${samp}.sorted.bam
        fi
    fi
done < sampleIDs.csv

# 2. Run ChIPQC
conda activate Rrr

## create an output directory for ChIPQC results
mkdir -p qc_out/chip_gsMask

## To know chromosome lengths
cut -f1-2 ./ara_refs/GCF_000001735.4_TAIR10.1_genomic.fna.fai #can save to file

## Fix gff file
### Remove header
./ara_refs/GCF_000001735.4_TAIR10.1_genomic.gff #manually remove header and save as _noheader.gff
### Remove hash symbol and one entry where start is bigger than end
cat ./ara_refs/GCF_000001735.4_TAIR10.1_genomic_noheader.gff | sed -e 's/#//g' | awk -F $'\t' '($5 > $4){print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS="\t" > ./ara_refs/GCF_000001735.4_TAIR10.1_genomic_noheaderFilt.gff
### Remove two genes that give problems
cat ./ara_refs/GCF_000001735.4_TAIR10.1_genomic_noheaderFilt.gff | awk '!/GeneID:36335702/' OFS="\t" | awk '!/GeneID:36335684/' OFS="\t" | cat > ./ara_refs/GCF_000001735.4_TAIR10.1_genomic_noheaderFilt1.gff
### Remove chloroplast and mithocondrial genomes
cat ./ara_refs/GCF_000001735.4_TAIR10.1_genomic_noheaderFilt1.gff | awk '!/NC_037304.1/' OFS="\t" | awk '!/NC_000932.1/' OFS="\t" | cat > ./ara_refs/GCF_000001735.4_TAIR10.1_genomic_noheaderFilt2.gff

## Run
Rscript ./scripts/ChIPQC_samples.R

conda deactivate Rrr
