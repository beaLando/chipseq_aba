#!/bin/bash

## NB: First check how you created bigwig files on bam alignment based on greenscreen (extendReads and which file?). In case re-create them.
## See here for more infos: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html

# CHECK FRAGMENT SIZE
## MACS2:
### In our case it should be between 100-200bp since fragments were of that length before sequencing (indeed, in our case, read length and fragment length overlap)
conda activate chippy

mkdir -p ./qc_out/mapped

while read line; do

    samp=`echo $line | cut -d "," -f1`
    if [[ "$samp" != "SampleID" ]]; then

       sorted_nodups="mapped/${samp}.noDups.sorted.bam" #need deduped bams
       out_R="${samp}_predictdPlot.R" #need deduped bams
       
       macs2 predictd -i $sorted_nodups --outdir ./qc_out/mapped --rfile $out_R -g 101274395 -m 2 50

    fi

done < sampleIDs.csv # I should actually run it on immunoprec. samples only, but here i'm running it on inputs

## DEEPTOOLS:
bamPEFragmentSize --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --samplesLabel cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --blackListFileName ara_refs/arabidopsis_greenscreen_20inputs.bed \
  --plotFileFormat pdf \
  --table ./qc_out/mapped/fragment_sizes.txt \
  -o ./qc_out/mapped/fragment_sizes.pdf -p 6


# QUALITY CHECKS IN DEEPTOOLS
## GENOME COVERAGE
plotCoverage --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --blackListFileName ara_refs/arabidopsis_greenscreen_20inputs.bed \
  --outRawCounts ./qc_out/mapped/coverage.txt \
  --plotFileFormat pdf \
  --plotFile ./qc_out/mapped/coverage.pdf -p 6

## ENRICHMENT
### With tair gtf
plotEnrichment --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --blackListFileName ara_refs/arabidopsis_greenscreen_20inputs.bed \
  --BED ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf \
  --extendReads 170 \
  --plotFileFormat pdf \
  -o ./qc_out/mapped/enrichment.pdf -p 6

### Only with promoter regions
plotEnrichment --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --blackListFileName ara_refs/arabidopsis_greenscreen_20inputs.bed \
  --BED ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat_TSS.bed \
  --extendReads 170 \
  --plotFileFormat pdf \
  -o ./qc_out/mapped/enrichment_tss.pdf -p 6

## FINGERPRINTs
plotFingerprint --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --extendReads 170  --binSize=1000 --plotFile ./qc_out/mapped/fingerprints.pdf \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./qc_out/mapped/fingerprint.log

plotFingerprint --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  --extendReads 170  --binSize=1000 --plotFile ./qc_out/mapped/fingerprints_cocit.pdf \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 \
  -p 6 &> ./qc_out/mapped/fingerprint_cocit.log

plotFingerprint --bamfiles mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --extendReads 170  --binSize=1000 --plotFile ./qc_out/mapped/fingerprints_cocitaba.pdf \
  --labels cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./qc_out/mapped/fingerprint_cocitaba.log

## CORRELOGRAM SAMPLES
### Mapped BAMs
multiBamSummary bins --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --outFileName ./qc_out/mapped/multiBamArray.npz --binSize=5000 \
  --extendReads=170 --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./qc_out/mapped/multiBamSummary.log

plotCorrelation --corData ./qc_out/mapped/multiBamArray.npz \
  --plotFile ./qc_out/mapped/multiBamArray_correl.pdf --outFileCorMatrix ./qc_out/mapped/multiBamArray_correl_mx.txt \
  --whatToPlot heatmap --corMethod spearman

## BIGWIGs & MATRIX plots
### If you want to recreate them:
#### For each sample bigwig:
bamCoverage --bam mapped/Unknown_CE349-003R0001.sorted.bam \
  --outFileName ./qc_out/mapped/Unknown_CE349-003R0001.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 101274395 --extendReads 170 \
  --binSize 50 -p 6 2> ./qc_out/mapped/Unknown_CE349-003R0001_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0003.sorted.bam \
  --outFileName ./qc_out/mapped/Unknown_CE349-003R0003.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 101274395 --extendReads 170 \
  --binSize 50 -p 6 2> ./qc_out/mapped/Unknown_CE349-003R0003_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0005.sorted.bam \
  --outFileName ./qc_out/mapped/Unknown_CE349-003R0005.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 101274395 --extendReads 170 \
  --binSize 50 -p 6 2> ./qc_out/mapped/Unknown_CE349-003R0005_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0007.sorted.bam \
  --outFileName ./qc_out/mapped/Unknown_CE349-003R0007.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 101274395 --extendReads 170 \
  --binSize 50 -p 6 2> ./qc_out/mapped/Unknown_CE349-003R0007_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0002.sorted.bam \
  --outFileName ./qc_out/mapped/Unknown_CE349-003R0002.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 101274395 --extendReads 170 \
  --binSize 50 -p 6 2> ./qc_out/mapped/Unknown_CE349-003R0002_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0004.sorted.bam \
  --outFileName ./qc_out/mapped/Unknown_CE349-003R0004.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 101274395 --extendReads 170 \
  --binSize 50 -p 6 2> ./qc_out/mapped/Unknown_CE349-003R0004_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0006.sorted.bam \
  --outFileName ./qc_out/mapped/Unknown_CE349-003R0006.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 101274395 --extendReads 170 \
  --binSize 50 -p 6 2> ./qc_out/mapped/Unknown_CE349-003R0006_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --outFileName ./qc_out/mapped/Unknown_CE349-003R0008.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 101274395 --extendReads 170 \
  --binSize 50 -p 6 2> ./qc_out/mapped/Unknown_CE349-003R0008_bamCoverage.sorted.log

#### For average of bigwig by genotype and input/gfp
bigwigAverage -b ./qc_out/mapped/Unknown_CE349-003R0001.sorted.bw ./qc_out/mapped/Unknown_CE349-003R0003.sorted.bw -o ./qc_out/mapped/Unknown_CE349-003R000_cocitInput_AVG.sorted.bw
bigwigAverage -b ./qc_out/mapped/Unknown_CE349-003R0002.sorted.bw ./qc_out/mapped/Unknown_CE349-003R0004.sorted.bw -o ./qc_out/mapped/Unknown_CE349-003R000_cocit_AVG.sorted.bw
bigwigAverage -b ./qc_out/mapped/Unknown_CE349-003R0005.sorted.bw ./qc_out/mapped/Unknown_CE349-003R0007.sorted.bw -o ./qc_out/mapped/Unknown_CE349-003R000_cocitabaInput_AVG.sorted.bw
bigwigAverage -b ./qc_out/mapped/Unknown_CE349-003R0006.sorted.bw ./qc_out/mapped/Unknown_CE349-003R0008.sorted.bw -o ./qc_out/mapped/Unknown_CE349-003R000_cocitaba_AVG.sorted.bw

#### For bigwig of chip normalized against input:
bamCompare -b1 mapped/Unknown_CE349-003R0002.sorted.bam \
-b2 mapped/Unknown_CE349-003R0001.sorted.bam \
-o ./qc_out/mapped/bamcompare_cocitRep1.bw \
--scaleFactorsMethod SES \
--centerReads \
--binSize 10 \
-p 6 2> ./qc_out/mapped/bamcompare_cocitRep1.log

bamCompare -b1 mapped/Unknown_CE349-003R0004.sorted.bam \
-b2 mapped/Unknown_CE349-003R0003.sorted.bam \
-o ./qc_out/mapped/bamcompare_cocitRep2.bw \
--scaleFactorsMethod SES \
--centerReads \
--binSize 10 \
-p 6 2> ./qc_out/mapped/bamcompare_cocitRep2.log

bamCompare -b1 mapped/Unknown_CE349-003R0006.sorted.bam \
-b2 mapped/Unknown_CE349-003R0005.sorted.bam \
-o ./qc_out/mapped/bamcompare_cocitabaRep1.bw \
--scaleFactorsMethod SES \
--centerReads \
--binSize 10 \
-p 6 2> ./qc_out/mapped/bamcompare_cocitabaRep1.log

bamCompare -b1 mapped/Unknown_CE349-003R0008.sorted.bam \
-b2 mapped/Unknown_CE349-003R0007.sorted.bam \
-o ./qc_out/mapped/bamcompare_cocitabaRep2.bw \
--scaleFactorsMethod SES \
--centerReads \
--binSize 10 \
-p 6 2> ./qc_out/mapped/bamcompare_cocitabaRep2.log


#### Then I actually compute matrix by gene
computeMatrix reference-point -S ./qc_out/mapped/Unknown_CE349-003R0002.sorted.bw ./qc_out/mapped/Unknown_CE349-003R0004.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o ./qc_out/mapped/matrix_cocit_TSS.gz -p 6 --outFileSortedRegions ./qc_out/mapped/matrix_cocit_TSS.bed --transcriptID gene --transcript_id_designator gene_id
computeMatrix reference-point -S ./qc_out/mapped/Unknown_CE349-003R0006.sorted.bw ./qc_out/mapped/Unknown_CE349-003R0008.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o ./qc_out/mapped/matrix_cocitaba_TSS.gz -p 6 --outFileSortedRegions ./qc_out/mapped/matrix_cocitaba_TSS.bed --transcriptID gene --transcript_id_designator gene_id 

computeMatrix reference-point -S ./qc_out/mapped/Unknown_CE349-003R000_cocit_AVG.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o ./qc_out/mapped/matrix_cocit_pooled_TSS.gz -p 6 --outFileSortedRegions ./qc_out/mapped/matrix_cocit_pooled_TSS.bed --transcriptID gene --transcript_id_designator gene_id
computeMatrix reference-point -S ./qc_out/mapped/Unknown_CE349-003R000_cocitaba_AVG.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o ./qc_out/mapped/matrix_cocitaba_pooled_TSS.gz -p 6 --outFileSortedRegions ./qc_out/mapped/matrix_cocitaba_pooled_TSS.bed --transcriptID gene --transcript_id_designator gene_id 

#### Then you can do plotting
plotHeatmap -m ./qc_out/mapped/matrix_cocit_TSS.gz \
-out ./qc_out/mapped/matrix_cocit_TSS_profile-heatmap.pdf \
--colorList white,black \
--zMin 0 --zMax 3.8 \
--missingDataColor 0.5 \
--plotFileFormat pdf 

plotHeatmap -m ./qc_out/mapped/matrix_cocitaba_TSS.gz \
-out ./qc_out/mapped/matrix_cocitaba_TSS_profile-heatmap.pdf \
--colorList white,black  \
--zMin 0 --zMax 3.8 \
--missingDataColor 0.5 \
--plotFileFormat pdf 

plotHeatmap -m ./qc_out/mapped/matrix_cocit_pooled_TSS.gz \
-out ./qc_out/mapped/matrix_cocit_pooled_TSS_profile-heatmap.pdf \
--colorList white,black \
--zMin 0 --zMax 3.8 \
--missingDataColor 0.5 \
--plotFileFormat pdf 

plotHeatmap -m ./qc_out/mapped/matrix_cocitaba_pooled_TSS.gz \
-out ./qc_out/mapped/matrix_cocitaba_pooled_TSS_profile-heatmap.pdf \
--colorList white,black \
--zMin 0 --zMax 3.8 \
--missingDataColor 0.5 \
--plotFileFormat pdf 

conda deactivate
