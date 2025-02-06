# OVERLAP in PEAKS BETWEEN REPLICATES OR IDR-peaks
mkdir -p qc_out/macs2

## NB: files must be sorted by coordinates (chromosome and start, but narrowpeak files are already sorted)
## Also, I need bed format, not narrowpeak:
cat macs2_out/gsMask/Unknown_CE349-003R0002_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/gsMask/Unknown_CE349-003R0002_peaks.bed
cat macs2_out/gsMask/Unknown_CE349-003R0004_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/gsMask/Unknown_CE349-003R0004_peaks.bed
cat macs2_out/gsMask/Unknown_CE349-003R0006_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/gsMask/Unknown_CE349-003R0006_peaks.bed
cat macs2_out/gsMask/Unknown_CE349-003R0008_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/gsMask/Unknown_CE349-003R0008_peaks.bed

cat macs2_out/gsMask/cocit_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/gsMask/cocit_peaks.bed
cat macs2_out/gsMask/cocitaba_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/gsMask/cocitaba_peaks.bed

cat macs2_out/IDR/cocit_idr01_gsMask.peaks | tr -s '\t' | cut -f-6 > macs2_out/IDR/cocit_idr01_gsMask.bed
cat macs2_out/IDR/cocitaba_idr01_gsMask.peaks | tr -s '\t' | cut -f-6 > macs2_out/IDR/cocitaba_idr01_gsMask.bed

## Intersections between replicates
intersectBed -a ./macs2_out/gsMask/Unknown_CE349-003R0002_peaks.bed -b ./macs2_out/gsMask/Unknown_CE349-003R0004_peaks.bed > qc_out/macs2/003R0002_003R0004_intersect.bed
intersectBed -a ./macs2_out/gsMask/Unknown_CE349-003R0006_peaks.bed -b ./macs2_out/gsMask/Unknown_CE349-003R0008_peaks.bed > qc_out/macs2/003R0006_003R0008_intersect.bed

## Intersections between peaks called with pooling and respective replicates
intersectBed -a macs2_out/gsMask/cocit_peaks.bed -b ./macs2_out/gsMask/Unknown_CE349-003R0002_peaks.bed > qc_out/macs2/cocit_pooled_003R0002_intersect.bed
intersectBed -a macs2_out/gsMask/cocit_peaks.bed -b ./macs2_out/gsMask/Unknown_CE349-003R0004_peaks.bed > qc_out/macs2/cocit_pooled_003R0004_intersect.bed
intersectBed -a macs2_out/gsMask/cocitaba_peaks.bed -b ./macs2_out/gsMask/Unknown_CE349-003R0006_peaks.bed > qc_out/macs2/cocitaba_pooled_003R0006_intersect.bed
intersectBed -a macs2_out/gsMask/cocitaba_peaks.bed -b ./macs2_out/gsMask/Unknown_CE349-003R0008_peaks.bed > qc_out/macs2/cocitaba_pooled_003R0008_intersect.bed

## Intersections between peaks called with pooling and peaks called with IDR
intersectBed -a macs2_out/gsMask/cocit_peaks.bed -b macs2_out/IDR/cocit_idr01_gsMask.bed > qc_out/macs2/cocit_pooled_idr01_gsMask.bed
intersectBed -a macs2_out/gsMask/cocitaba_peaks.bed -b macs2_out/IDR/cocitaba_idr01_gsMask.bed > qc_out/macs2/cocitaba_pooled_idr01_gsMask.bed

# QUALITY CHECKs on CALLED PEAK REGIONS
## CORRELOGRAMs for mapped BAMs considering only peak regions
conda activate chippy

multiBamSummary BED-file --BED ./macs2_out/IDR/cocit_idr01_gsMask.bed \
      --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
      mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
      mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
      mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
      --outFileName ./qc_out/macs2/multiBamArray_IDRpeaks.npz \
      --extendReads=170 -p 6 \
      --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2

plotCorrelation --corData ./qc_out/macs2/multiBamArray_IDRpeaks.npz \
      --plotFile ./qc_out/macs2/multiBamArray_IDRpeaks_correl.pdf --outFileCorMatrix ./qc_out/macs2/multiBamArray_IDRpeaks_correl_mx.txt \
      --whatToPlot heatmap --corMethod pearson --plotNumbers --removeOutliers

## ENRICHMENT in PEAKS
plotEnrichment --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  --BED ./macs2_out/IDR/cocit_idr01_gsMask.bed \
  --extendReads 170 \
  --plotFileFormat pdf \
  -o ./qc_out/macs2/enrichment_at_peaks.pdf -p 6

conda deactivate

## See here for interpretation of figures and further quality checks with creation of pseudoreplicates:
## https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

## Another way to check quality is to check overlap in peaks between replicates using bedtools
## You can also visually inspect peaks: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/11_qualitative_assessment_IGV.html