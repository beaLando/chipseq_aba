# HEATMAPS & CO.
## NB: First check how you created bigwig files on bam alignment based on greenscreen (extendReads and which file?). In case re-create them.
## See here for more infos: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html

conda activate chippy

## FINGERPRINTs
plotFingerprint --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --extendReads 250  --binSize=1000 --plotFile ./mapped/fingerprints.pdf \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./mapped/fingerprint.log

plotFingerprint --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  --extendReads 250  --binSize=1000 --plotFile ./mapped/fingerprints_cocit.pdf \
  --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 \
  -p 6 &> ./mapped/fingerprint_cocit.log

plotFingerprint --bamfiles mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --extendReads 250  --binSize=1000 --plotFile ./mapped/fingerprints_cocitaba.pdf \
  --labels cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./mapped/fingerprint_cocitaba.log

## CORRELOGRAMs
### Mapped BAMs
multiBamSummary bins --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
  mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
  mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
  mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --outFileName ./mapped/multiBamArray.npz --binSize=5000 \
  --extendReads=250 --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2 \
  -p 6 &> ./mapped/multiBamSummary.log

plotCorrelation --corData ./mapped/multiBamArray.npz \
  --plotFile ./mapped/multiBamArray_correl.pdf --outFileCorMatrix ./mapped/multiBamArray_correl_mx.txt \
  --whatToPlot heatmap --corMethod spearman

### Mapped BAMs considering only peak regions
multiBamSummary BED-file --BED ./macs2_out/IDR/cocit_idr_gsMask10.bed \
      --bamfiles mapped/Unknown_CE349-003R0001.sorted.bam mapped/Unknown_CE349-003R0003.sorted.bam  \
      mapped/Unknown_CE349-003R0002.sorted.bam mapped/Unknown_CE349-003R0004.sorted.bam \
      mapped/Unknown_CE349-003R0003.sorted.bam mapped/Unknown_CE349-003R0005.sorted.bam \
      mapped/Unknown_CE349-003R0006.sorted.bam mapped/Unknown_CE349-003R0008.sorted.bam \
      --outFileName ./mapped/multiBamArray_peaks.npz \
      --extendReads=250 -p 6 \
      --labels cocit_input_rep1 cocit_input_rep2 cocit_rep1 cocit_rep2 cocitaba_input_rep1 cocitaba_input_rep2 cocitaba_rep1 cocitaba_rep2

plotCorrelation --corData ./mapped/multiBamArray_peaks.npz \
      --plotFile ./mapped/multiBamArray_peaks_correl.pdf --outFileCorMatrix ./mapped/multiBamArray_peaks_correl_mx.txt \
      --whatToPlot heatmap --corMethod pearson --plotNumbers --removeOutliers

## BIGWIGs
### If you want to recreate them:
#### For normal bigwig:
bamCoverage --bam mapped/Unknown_CE349-003R0001.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0001.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 250 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0001_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0003.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0003.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 250 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0003_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0005.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0005.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 250 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0005_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0007.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0007.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 250 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0007_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0002.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0002.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 250 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0002_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0004.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0004.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 250 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0004_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0006.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0006.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 250 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0006_bamCoverage.sorted.log

bamCoverage --bam mapped/Unknown_CE349-003R0008.sorted.bam \
  --outFileName mapped/Unknown_CE349-003R0008.sorted.bw \
  --normalizeUsing RPGC --effectiveGenomeSize 119482012 --extendReads 250 \
  --binSize 50 -p 6 2> mapped/Unknown_CE349-003R0008_bamCoverage.sorted.log

#### Calculate average of bigwig by genotype and input/gfp
##### Method above
bigwigAverage -b mapped/Unknown_CE349-003R0001.sorted.bw mapped/Unknown_CE349-003R0003.sorted.bw -o mapped/Unknown_CE349-003R000_cocitInput_AVG.sorted.bw
bigwigAverage -b mapped/Unknown_CE349-003R0002.sorted.bw mapped/Unknown_CE349-003R0004.sorted.bw -o mapped/Unknown_CE349-003R000_cocit_AVG.sorted.bw
bigwigAverage -b mapped/Unknown_CE349-003R0005.sorted.bw mapped/Unknown_CE349-003R0007.sorted.bw -o mapped/Unknown_CE349-003R000_cocitabaInput_AVG.sorted.bw
bigwigAverage -b mapped/Unknown_CE349-003R0006.sorted.bw mapped/Unknown_CE349-003R0008.sorted.bw -o mapped/Unknown_CE349-003R000_cocitaba_AVG.sorted.bw

##### Method greenscreen
bigwigAverage -b mapped/Unknown_CE349-003R0001.bw mapped/Unknown_CE349-003R0003.bw -o mapped/Unknown_CE349-003R000_cocitInput_AVG.bw
bigwigAverage -b mapped/Unknown_CE349-003R0002.bw mapped/Unknown_CE349-003R0004.bw -o mapped/Unknown_CE349-003R000_cocit_AVG.bw
bigwigAverage -b mapped/Unknown_CE349-003R0005.bw mapped/Unknown_CE349-003R0007.bw -o mapped/Unknown_CE349-003R000_cocitabaInput_AVG.bw
bigwigAverage -b mapped/Unknown_CE349-003R0006.bw mapped/Unknown_CE349-003R0008.bw -o mapped/Unknown_CE349-003R000_cocitaba_AVG.bw

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

conda deactivate


## PLOTTING
### First need to create a matrix with scores based on the read density values in the bigWig files for each gene/genomic region (in BED format!!):
#### Going from GFF to BED is surprisingly difficult, because GFF files are nasty
#### First, I fix original GFF file, then I convert it to GTF, using agat:
conda activate gff_gtf
agat_convert_sp_gxf2gxf.pl --gff ara_refs/GCF_000001735.4_TAIR10.1_genomic.gff -o ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gff
agat_convert_sp_gff2gtf.pl --gff ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gff -o ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf
conda deactivate

#### Then I actually compute matrix by gene (usually it is done by transcript)
conda activate chippy
computeMatrix reference-point -S mapped/Unknown_CE349-003R0002.sorted.bw mapped/Unknown_CE349-003R0004.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o macs2_out/matrix_cocit_TSS.gz -p 6 --outFileSortedRegions macs2_out/matrix_cocit_TSS.bed --transcriptID gene --transcript_id_designator gene_id
computeMatrix reference-point -S mapped/Unknown_CE349-003R0006.sorted.bw mapped/Unknown_CE349-003R0008.sorted.bw -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o macs2_out/matrix_cocitaba_TSS.gz -p 6 --outFileSortedRegions macs2_out/matrix_cocitaba_TSS.bed --transcriptID gene --transcript_id_designator gene_id 

computeMatrix reference-point -S mapped/pooled/Unknown_CE349-003R000_gfpcocit.bigwig -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o macs2_out/matrix_cocit_pooled_TSS.gz -p 6 --outFileSortedRegions macs2_out/matrix_cocit_pooled_TSS.bed --transcriptID gene --transcript_id_designator gene_id
computeMatrix reference-point -S mapped/pooled/Unknown_CE349-003R000_gfpcocitaba.bigwig -R ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf -b 2000 -a 2000 --referencePoint TSS --skipZeros -o macs2_out/matrix_cocitaba_pooled_TSS.gz -p 6 --outFileSortedRegions macs2_out/matrix_cocitaba_pooled_TSS.bed --transcriptID gene --transcript_id_designator gene_id 

### Then you can do plotting (I guess I can also do it in R, once I have the matrix above, it's just a heatmap and density plot after all):
plotHeatmap -m macs2_out/matrix_cocit_TSS.gz \
-out macs2_out/matrix_cocit_TSS_profile-heatmap.png \
--colorList white,black \
--zMin 0 --zMax 2.1  

plotHeatmap -m macs2_out/matrix_cocitaba_TSS.gz \
-out macs2_out/matrix_cocitaba_TSS_profile-heatmap.png \
--colorList white,black \
--zMin 0 --zMax 2.1  

plotHeatmap -m macs2_out/matrix_cocit_pooled_TSS.gz \
-out macs2_out/matrix_cocit_pooled_TSS_profile-heatmap.png \
--colorList white,black

plotHeatmap -m macs2_out/matrix_cocitaba_pooled_TSS.gz \
-out macs2_out/matrix_cocitaba_pooled_TSS_profile-heatmap.png \
--colorList white,black

### I can also look at enrichment in specific region (e.g. for FT > I need to create BED file for this // all genes close to significan peaks // ...):
#### I try for FT
computeMatrix scale-regions \
-R ara_refs/FT.bed \
-S mapped/Unknown_CE349-003R0002.sorted.bw mapped/Unknown_CE349-003R0004.sorted.bw mapped/Unknown_CE349-003R0006.sorted.bw mapped/Unknown_CE349-003R0008.sorted.bw \
--skipZeros \
-p 6 \
-a 2000 -b 2000 \
-o macs2_out/matrix_FT_cocit_cocitaba_TSS_binding_sites.gz

#--regionBodyLength 2000 \


plotProfile -m macs2_out/matrix_FT_cocit_cocitaba_TSS_binding_sites.gz \
-out macs2_out/matrix_FT_cocit_cocitaba_TSS_binding_sites.png \
--perGroup  --plotTitle "" \
--samplesLabel "cocit-Rep1" "cocit-Rep2" "cocitaba-Rep1" "cocitaba-Rep2" \
-T "Binding sites at FT"  -z "" \
--startLabel "" \
--endLabel "" \
--colors red red darkblue darkblue

conda deactivate