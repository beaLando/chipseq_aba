# IDR & INTERSECTION

conda activate idr_chip

## Make IDR output directory
mkdir macs2_out/IDR

## Sort peak by -log10(p-value)
sort -k8,8nr macs2_out/Unknown_CE349-003R0002_peaks.narrowPeak > macs2_out/Unknown_CE349-003R0002_peaks_sortedbylP.narrowPeak 
sort -k8,8nr macs2_out/Unknown_CE349-003R0004_peaks.narrowPeak > macs2_out/Unknown_CE349-003R0004_peaks_sortedbylP.narrowPeak
sort -k8,8nr macs2_out/Unknown_CE349-003R0006_peaks.narrowPeak > macs2_out/Unknown_CE349-003R0006_peaks_sortedbylP.narrowPeak
sort -k8,8nr macs2_out/Unknown_CE349-003R0008_peaks.narrowPeak > macs2_out/Unknown_CE349-003R0008_peaks_sortedbylP.narrowPeak

## IDR
### NB: for now I run it on peak files as they are and see how it goes but, ideally, I should re-run macs2 with less stringent options and then run idr

idr --samples macs2_out/Unknown_CE349-003R0002_peaks_sortedbylP.narrowPeak macs2_out/Unknown_CE349-003R0004_peaks_sortedbylP.narrowPeak \
--input-file-type narrowPeak \
--rank q.value \
--output-file macs2_out/IDR/cocit_idr \
--plot \
--log-output-file macs2_out/IDR/cocit.idr.log

idr --samples macs2_out/Unknown_CE349-003R0006_peaks_sortedbylP.narrowPeak macs2_out/Unknown_CE349-003R0008_peaks_sortedbylP.narrowPeak \
--input-file-type narrowPeak \
--rank q.value \
--output-file macs2_out/IDR/cocitaba_idr \
--plot \
--log-output-file macs2_out/IDR/cocitaba.idr.log

intersectBed -v -wa -a macs2_out/IDR/cocit_idr -b ara_refs/arabidopsis_greenscreen_20inputs.bed > macs2_out/IDR/cocit_idr_gsMask.bed
intersectBed -v -wa -a macs2_out/IDR/cocitaba_idr -b ara_refs/arabidopsis_greenscreen_20inputs.bed > macs2_out/IDR/cocitaba_idr_gsMask.bed
awk -F"\t" -v q=10 'BEGIN{OFS="\t"} $9>=q && $1!="ChrC" && $1!="ChrM"{print}' macs2_out/IDR/cocitaba_idr_gsMask.bed > macs2_out/IDR/cocitaba_idr_gsMask10.bed
awk -F"\t" -v q=10 'BEGIN{OFS="\t"} $9>=q && $1!="ChrC" && $1!="ChrM"{print}' macs2_out/IDR/cocit_idr_gsMask.bed > macs2_out/IDR/cocit_idr_gsMask10.bed

## OVERLAP in PEAKS
### NB: files must be sorted by coordinates (chromosome and start, but narrowpeak files are already sorted)
### Also, I need bed format, not narrowpeak:
cat macs2_out/Unknown_CE349-003R0002_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/Unknown_CE349-003R0002_peaks.bed
cat macs2_out/Unknown_CE349-003R0004_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/Unknown_CE349-003R0004_peaks.bed
cat macs2_out/Unknown_CE349-003R0006_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/Unknown_CE349-003R0006_peaks.bed
cat macs2_out/Unknown_CE349-003R0008_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/Unknown_CE349-003R0008_peaks.bed

### Intersections between replicates
intersectBed -a ./macs2_out/Unknown_CE349-003R0002_peaks.bed -b ./macs2_out/Unknown_CE349-003R0004_peaks.bed > macs2_out/IDR/003R0002_003R0004_intersect.bed
intersectBed -a ./macs2_out/Unknown_CE349-003R0006_peaks.bed -b ./macs2_out/Unknown_CE349-003R0008_peaks.bed > macs2_out/IDR/003R0006_003R0008_intersect.bed

### No-intersection between replicates
intersectBed -a ./macs2_out/Unknown_CE349-003R0002_peaks.bed -b ./macs2_out/Unknown_CE349-003R0004_peaks.bed -v > macs2_out/IDR/003R0002_003R0004_NOintersect.bed
intersectBed -a ./macs2_out/Unknown_CE349-003R0006_peaks.bed -b ./macs2_out/Unknown_CE349-003R0008_peaks.bed -v > macs2_out/IDR/003R0006_003R0008_NOintersect.bed
intersectBed -a ./macs2_out/Unknown_CE349-003R0004_peaks.bed -b ./macs2_out/Unknown_CE349-003R0002_peaks.bed -v > macs2_out/IDR/003R0004_003R0002_NOintersect.bed
intersectBed -a ./macs2_out/Unknown_CE349-003R0008_peaks.bed -b ./macs2_out/Unknown_CE349-003R0006_peaks.bed -v > macs2_out/IDR/003R0008_003R0006_NOintersect.bed

## Intersections and no-intersections between peaks called with pooling and respective replicates
cat macs2_out/cocit_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/cocit_peaks.bed
cat macs2_out/cocitaba_peaks.narrowPeak | tr -s '\t' | cut -f-6 > macs2_out/cocitaba_peaks.bed

intersectBed -a macs2_out/cocit_peaks.bed -b ./macs2_out/Unknown_CE349-003R0002_peaks.bed > macs2_out/IDR/cocit_003R0002_intersect.bed
intersectBed -a macs2_out/cocit_peaks.bed -b ./macs2_out/Unknown_CE349-003R0004_peaks.bed > macs2_out/IDR/cocit_003R0004_intersect.bed
intersectBed -a macs2_out/cocitaba_peaks.bed -b ./macs2_out/Unknown_CE349-003R0006_peaks.bed > macs2_out/IDR/cocitaba_003R0006_intersect.bed
intersectBed -a macs2_out/cocitaba_peaks.bed -b ./macs2_out/Unknown_CE349-003R0008_peaks.bed > macs2_out/IDR/cocitaba_003R0008_intersect.bed

intersectBed -a macs2_out/cocit_peaks.bed -b ./macs2_out/Unknown_CE349-003R0002_peaks.bed -v > macs2_out/IDR/cocit_003R0002_NOintersect.bed
intersectBed -a macs2_out/cocit_peaks.bed -b ./macs2_out/Unknown_CE349-003R0004_peaks.bed -v > macs2_out/IDR/cocit_003R0004_NOintersect.bed
intersectBed -a macs2_out/cocitaba_peaks.bed -b ./macs2_out/Unknown_CE349-003R0006_peaks.bed -v > macs2_out/IDR/cocitaba_003R0006_NOintersect.bed
intersectBed -a macs2_out/cocitaba_peaks.bed -b ./macs2_out/Unknown_CE349-003R0008_peaks.bed -v > macs2_out/IDR/cocitaba_003R0008_NOintersect.bed

conda deactivate

## See here for interpretation of figures and further quality checks with creation of pseudoreplicates:
## https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

## Another way to check quality is to check overlap in peaks between replicates using bedtools
## You can also visually inspect peaks: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/11_qualitative_assessment_IGV.html