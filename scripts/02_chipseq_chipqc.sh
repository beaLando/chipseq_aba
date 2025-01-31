#!/bin/bash

## NB: First check how you created bigwig files on bam alignment based on greenscreen (extendReads and which file?). In case re-create them.
## See here for more infos: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html

# GFF formatting
## First need to create a matrix with scores based on the read density values in the bigWig files for each gene/genomic region (in BED format!!):
### Going from GFF to BED is surprisingly difficult, because GFF files are nasty
#### First, I fix original GFF file, then I convert it to GTF, using agat:
conda activate gff_gtf
agat_convert_sp_gxf2gxf.pl --gff ara_refs/GCF_000001735.4_TAIR10.1_genomic.gff -o ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gff
agat_convert_sp_gff2gtf.pl --gff ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gff -o ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf

#### I also save bed file with only +/-1000bp around TSS genes
awk 'BEGIN{OFS="\t"}{if($3=="gene"){$3="tss";if($7 == "+"){$5=$4}else{$4=$5};print $0}}' ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gff | awk {'print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7'} | bedops --range -1000:1000 --everything - > ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat_TSS.bed
conda deactivate

# MASK READS that OVERLAP GREENSCREEN
# create output directory for reads that do and do not overlap greenscreen regions
greenscreen_regions="ara_refs/arabidopsis_greenscreen_20inputs.bed"

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


# QUALITY CHECKS IN DEEPTOOLS
conda activate chippy

## FRAGMENT SIZES
bamPEFragmentSize --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --samplesLabel cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --blackListFileName ara_refs/arabidopsis_greenscreen_20inputs.bed \
  --plotFileFormat pdf \
  --table mapped/fragment_sizes.txt \
  -o mapped/fragment_sizes.pdf -p 6

## GENOME COVERAGE
plotCoverage --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --blackListFileName ara_refs/arabidopsis_greenscreen_20inputs.bed \
  --outRawCounts mapped/coverage.txt \
  --plotFileFormat pdf \
  --plotFile mapped/coverage.pdf -p 6

## ENRICHMENT
### With tair gtf
plotEnrichment --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --blackListFileName ara_refs/arabidopsis_greenscreen_20inputs.bed \
  --BED ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf \
  --extendReads 150 \
  --plotFileFormat pdf \
  -o mapped/enrichment.pdf -p 6

### Only with promoter regions
plotEnrichment --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --blackListFileName ara_refs/arabidopsis_greenscreen_20inputs.bed \
  --BED ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat_TSS.bed \
  --extendReads 150 \
  --plotFileFormat pdf \
  -o mapped/enrichment_tss.pdf -p 6

## FINGERPRINTs
plotFingerprint --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --extendReads 150  --binSize=1000 --plotFile ./mapped/fingerprints.pdf \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./mapped/fingerprint.log

plotFingerprint --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  --extendReads 150  --binSize=1000 --plotFile ./mapped/fingerprints_cocit.pdf \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 \
  -p 6 &> ./mapped/fingerprint_cocit.log

plotFingerprint --bamfiles mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --extendReads 150  --binSize=1000 --plotFile ./mapped/fingerprints_cocitaba.pdf \
  --labels cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./mapped/fingerprint_cocitaba.log

## CORRELOGRAM SAMPLES
### Mapped BAMs
multiBamSummary bins --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --outFileName ./mapped/multiBamArray.npz --binSize=5000 \
  --extendReads=150 --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./mapped/multiBamSummary.log

plotCorrelation --corData ./mapped/multiBamArray.npz \
  --plotFile ./mapped/multiBamArray_correl.pdf --outFileCorMatrix ./mapped/multiBamArray_correl_mx.txt \
  --whatToPlot heatmap --corMethod spearman

## BIGWIGs & MATRIX plots
### If you want to recreate them:
#### For each sample bigwig:
bamCoverage --bam mapped/Unknown_CE349-003R0001.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0001.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 150 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0001_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0003.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0003.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 150 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0003_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0005.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0005.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 150 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0005_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0007.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0007.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 150 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0007_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0002.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0002.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 150 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0002_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0004.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0004.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 150 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0004_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0006.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0006.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 150 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0006_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0008.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 150 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0008_bamCoverage.sorted.log

#### For average of bigwig by genotype and input/gfp
bigwigAverage -b mapped/Unknown_CE349-003R0001.sorted.bw mapped/Unknown_CE349-003R0003.sorted.bw -o mapped/Unknown_CE349-003R000_cocitInput_AVG.sorted.bw
bigwigAverage -b mapped/Unknown_CE349-003R0002.sorted.bw mapped/Unknown_CE349-003R0004.sorted.bw -o mapped/Unknown_CE349-003R000_cocit_AVG.sorted.bw
bigwigAverage -b mapped/Unknown_CE349-003R0005.sorted.bw mapped/Unknown_CE349-003R0007.sorted.bw -o mapped/Unknown_CE349-003R000_cocitabaInput_AVG.sorted.bw
bigwigAverage -b mapped/Unknown_CE349-003R0006.sorted.bw mapped/Unknown_CE349-003R0008.sorted.bw -o mapped/Unknown_CE349-003R000_cocitaba_AVG.sorted.bw

#### For bigwig of chip normalized against input:
bamCompare -b1 mapped/Unknown_CE349-003R0002.sorted.bam \
-b2 mapped/Unknown_CE349-003R0001.sorted.bam \
-o mapped/bamcompare_cocitRep1.bw \
--scaleFactorsMethod SES \
--centerReads \
--binSize 10 \
-p 6 2> mapped/bamcompare_cocitRep1.log

bamCompare -b1 mapped/Unknown_CE349-003R0004.sorted.bam \
-b2 mapped/Unknown_CE349-003R0003.sorted.bam \
-o mapped/bamcompare_cocitRep2.bw \
--scaleFactorsMethod SES \
--centerReads \
--binSize 10 \
-p 6 2> mapped/bamcompare_cocitRep2.log

bamCompare -b1 mapped/Unknown_CE349-003R0006.sorted.bam \
-b2 mapped/Unknown_CE349-003R0005.sorted.bam \
-o mapped/bamcompare_cocitabaRep1.bw \
--scaleFactorsMethod SES \
--centerReads \
--binSize 10 \
-p 6 2> mapped/bamcompare_cocitabaRep1.log

bamCompare -b1 mapped/Unknown_CE349-003R0008.sorted.bam \
-b2 mapped/Unknown_CE349-003R0007.sorted.bam \
-o mapped/bamcompare_cocitabaRep2.bw \
--scaleFactorsMethod SES \
--centerReads \
--binSize 10 \
-p 6 2> mapped/bamcompare_cocitabaRep2.log


#### Then I actually compute matrix by gene (usually it is done by transcript)
computeMatrix reference-point -S mapped/Unknown_CE349-003R0002.sorted.bw mapped/Unknown_CE349-003R0004.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o macs2_out/matrix_cocit_TSS.gz -p 6 --outFileSortedRegions macs2_out/matrix_cocit_TSS.bed --transcriptID gene --transcript_id_designator gene_id
computeMatrix reference-point -S mapped/Unknown_CE349-003R0006.sorted.bw mapped/Unknown_CE349-003R0008.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o macs2_out/matrix_cocitaba_TSS.gz -p 6 --outFileSortedRegions macs2_out/matrix_cocitaba_TSS.bed --transcriptID gene --transcript_id_designator gene_id 

computeMatrix reference-point -S mapped/Unknown_CE349-003R000_cocit_AVG.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o macs2_out/matrix_cocit_pooled_TSS.gz -p 6 --outFileSortedRegions macs2_out/matrix_cocit_pooled_TSS.bed --transcriptID gene --transcript_id_designator gene_id
computeMatrix reference-point -S mapped/Unknown_CE349-003R000_cocitaba_AVG.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o macs2_out/matrix_cocitaba_pooled_TSS.gz -p 6 --outFileSortedRegions macs2_out/matrix_cocitaba_pooled_TSS.bed --transcriptID gene --transcript_id_designator gene_id 

#### Then you can do plotting
plotHeatmap -m macs2_out/matrix_cocit_TSS.gz \
-out macs2_out/matrix_cocit_TSS_profile-heatmap.pdf \
--colorList white,black \
--zMin 0 --zMax 2.5 \
--missingDataColor 0.5 \
--plotFileFormat pdf 

plotHeatmap -m macs2_out/matrix_cocitaba_TSS.gz \
-out macs2_out/matrix_cocitaba_TSS_profile-heatmap.pdf \
--colorList white,black  \
--zMin 0 --zMax 2.5 \
--missingDataColor 0.5 \
--plotFileFormat pdf 

plotHeatmap -m macs2_out/matrix_cocit_pooled_TSS.gz \
-out macs2_out/matrix_cocit_pooled_TSS_profile-heatmap.pdf \
--colorList white,black \
--zMin 0 --zMax 2.5 \
--missingDataColor 0.5 \
--plotFileFormat pdf 

plotHeatmap -m macs2_out/matrix_cocitaba_pooled_TSS.gz \
-out macs2_out/matrix_cocitaba_pooled_TSS_profile-heatmap.pdf \
--colorList white,black \
--zMin 0 --zMax 2.5 \
--missingDataColor 0.5 \
--plotFileFormat pdf 


conda deactivate
