# LIBRARIES
library(readxl)
library(ChIPseeker)
library(GeneOverlap)
library(VennDiagram)
library(clusterProfiler)
library(gprofiler2) #!!! 4 TOMORROW: you can try to use gprofiler2 instead of clusterProfiler !!! > use code you wrote for rnaseq paolo
library(org.At.tair.db)
library(AnnotationDbi)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(tidyverse)


# DATA
## Annotation TAIR10
txdb <- TxDb.Athaliana.BioMart.plantsmart51

## De Los Reyes
### Summary methods
# The potential input of photoperiod in this network was determined by carrying out a genome-wide analysis 
# of ChIP-seq experiments in 35S:CO and co-10 7-DAG (days after germination) seedlings around Zeitgeber time (ZT) 14 
# under long-day conditions (LD) using specific CO antibodies (Serrano-Bueno et al., 2020). 
# The analysis of CO ChIP-seq data identified 3214 peaks and 2418 putative target genes (Supplemental Table 1).

# Finally, target gene selection was carried out using PeakAnnotator (v1.3), part of the PeakAnalyzer program (Salmon-Divon et al., 2010), 
# and an ad hoc R script developed in this study according to the criterion of the nearest downstream gene.
# For the CO ChIP-seq analysis, peaks detected in the 35S:CO experiment were filtered based on the peaks found in co-10 immunoprecipitation.

### 35SCO peaks
co35s <- read.table("./data/delosReyes2024/GSE222657_35SCO_peaks.narrowpeak")
str(co35s)
colnames(co35s) <- c("chr", "start", "end", "peak_name", "score", "strand", "signal", "pv", "qv", "summit")
co35s <- makeGRangesFromDataFrame(co35s, keep.extra.columns = TRUE)

peakAnno.co35s <- annotatePeak(co35s, TxDb=txdb, tssRegion=c(-1000, 1000), verbose = TRUE)
plotAnnoBar(peakAnno.co35s)
plotDistToTSS(peakAnno.co35s, title = "Distribution of transcription factor-binding loci \n relative to TSS")

### co-10 peaks
co10 <- read.table("./data/delosReyes2024/GSE222657_co10_peaks.narrowpeak")
str(co10)
colnames(co10) <- c("chr", "start", "end", "peak_name", "score", "strand", "signal", "pv", "qv", "summit")
co10 <- makeGRangesFromDataFrame(co10, keep.extra.columns = TRUE)

peakAnno.co10 <- annotatePeak(co10, TxDb=txdb, tssRegion=c(-1000, 1000), verbose = TRUE)
plotAnnoBar(peakAnno.co10)
plotDistToTSS(peakAnno.co10, title = "Distribution of transcription factor-binding loci \n relative to TSS")

### filtered peaks
filt <- read.table("./data/delosReyes2024/GSE222657_filtered_peaks.narrowpeak")
str(filt)
colnames(filt) <- c("chr", "start", "end", "peak_name", "score", "strand", "signal", "pv", "qv", "summit")
filt <- makeGRangesFromDataFrame(filt, keep.extra.columns = TRUE)

peakAnno.filt <- annotatePeak(filt, TxDb=txdb, tssRegion=c(-1000, 1000), verbose = TRUE)
plotAnnoBar(peakAnno.filt)
plotDistToTSS(peakAnno.filt, title = "Distribution of transcription factor-binding loci \n relative to TSS")

#### is the non-intersecting portion of co35s peaks with co10, the same as "filtered peaks"?
filt.maybe <- subsetByOverlaps(co35s, co10, invert = TRUE)
length(intersect(filt$peak_name, filt.maybe$peak_name)) #yes, there is only one location less, but overlap is full

### targets CO
chip.co.trgts <- read.table("./data/delosReyes2024/TableS1_COtargets_chipseq.txt")
str(chip.co.trgts)
colnames(chip.co.trgts) <- "geneID"

#### are these CO targets the same as the genes falling next to filtered peaks?
##### not clear...they used a custom script to annotate peaks close to genes, so it is possible that they get a slightly different annotations/number of annotations
##### we still have an overlap of 1902 genes between the two lists, which I would say suggests the target list is derived from the filtered genes
filt.trgts <- unique(peakAnno.filt@anno$geneId)
length(filt.trgts) #2781
length(unique(chip.co.trgts$geneID)) #2417
length(intersect(unique(chip.co.trgts$geneID), filt.trgts)) #1902

#### Rna-seq data
rna.35sco.up <- read_xlsx("./data/delosReyes2024/1-s2.0-S1674205224001862-mmc3.xlsx", skip = 1, sheet = 1)
str(rna.35sco.up)

rna.suc2co.up <- read_xlsx("./data/delosReyes2024/1-s2.0-S1674205224001862-mmc3.xlsx", skip = 1, sheet = 2)
str(rna.suc2co.up)

rna.co.dwn <- read_xlsx("./data/delosReyes2024/1-s2.0-S1674205224001862-mmc3.xlsx", skip = 1, sheet = 3)
rna.co.dwn

## Ours
### SUC2CO Peaks (IDR-checked, q>=10 filtered, greenscreen masked)
samplefiles <- list.files("output/macs2_out/IDR", pattern = "_gsMask10_chrFix.bed", full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("cocit", "cocitaba")

suc2co <- read.table("./output/macs2_out/IDR/cocit_idr_gsMask10_chrFix.bed")
str(suc2co)
colnames(suc2co) <- c("chr", "start", "end", "peak_name", "score", "strand", "V7", "V8", "V9", "summit", paste0("V", 11:20))
suc2co$peak_name <- as.character(1:nrow(suc2co))
suc2co <- makeGRangesFromDataFrame(suc2co, keep.extra.columns = TRUE)

peakAnno.suc2co <- annotatePeak(suc2co, TxDb=txdb, tssRegion=c(-1000, 1000), verbose = TRUE)
plotAnnoBar(peakAnno.suc2co)
plotDistToTSS(peakAnno.suc2co, title = "Distribution of transcription factor-binding loci \n relative to TSS")

peakAnno.suc2co

# Annotated peaks generated by ChIPseeker
# 4733/4733  peaks were annotated
# Genomic Annotation Summary:
#   Feature   Frequency
# 9           Promoter 64.37777308
# 4             5' UTR  0.04225650
# 3             3' UTR  4.26790619
# 1           1st Exon  0.08451299
# 7         Other Exon  0.71836045
# 2         1st Intron  0.21128248
# 8       Other Intron  0.88738644
# 6 Downstream (<=300)  1.11979717
# 5  Distal Intergenic 28.29072470

prop.anno.suc2co <- data.frame(Feature = c("Promoter", "UTRs", "UTRs", 
                                           "Other (introns, exons, downstream)", "Other (introns, exons, downstream)", "Other (introns, exons, downstream)",
                                           "Other (introns, exons, downstream)", "Other (introns, exons, downstream)", "Distal Intergenic"),
                               Percentage = c(64.37777308, 0.04225650, 4.26790619,
                                              0.08451299, 0.71836045, 0.21128248,
                                              0.88738644, 1.11979717, 28.29072470)) %>%
  group_by(Feature) %>%
  summarise(Percentage = sum(Percentage)) %>%
  ungroup() %>%
  as.data.frame()

ggplot(prop.anno.suc2co, 
       aes(x = "", y = Percentage, fill = Feature)) +
  geom_col(col = "black") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("black", "#666666", "#CCCCCC", "white")) +
  theme_void(base_size = 18)

## Retrieve peak annotation data frames
keytypes(txdb) #GENEID
columns(org.At.tair.db)

suc2co_annot <- data.frame(peakAnno.suc2co@anno)
co35s_annot <- data.frame(peakAnno.co35s@anno)
co10_annot <- data.frame(peakAnno.co10@anno)

### Get IDs
ids.suc2co <- unique(suc2co_annot$geneId)
ids.co35s <- unique(co35s_annot$geneId)
ids.co10 <- unique(co10_annot$geneId)

### Return the gene symbol for the set of IDs
annotations_db.suc2co <- AnnotationDbi::select(org.At.tair.db,
                                               keys = ids.suc2co,
                                               columns = c("ENTREZID", "GENENAME"),
                                               keytype = "TAIR")

annotations_db.co35s <- AnnotationDbi::select(org.At.tair.db,
                                               keys = ids.co35s,
                                               columns = c("ENTREZID", "GENENAME"),
                                               keytype = "TAIR")

annotations_db.co10 <- AnnotationDbi::select(org.At.tair.db,
                                               keys = ids.co10,
                                               columns = c("ENTREZID", "GENENAME"),
                                               keytype = "TAIR")

### Change IDs to character type to merge
annotations_db.suc2co$ENTREZID <- as.character(annotations_db.suc2co$ENTREZID)
annotations_db.co35s$ENTREZID <- as.character(annotations_db.co35s$ENTREZID)
annotations_db.co10$ENTREZID <- as.character(annotations_db.co10$ENTREZID)


# OVERLAPS 
## in peaks between our suc2co and 35sco
### with co-10 peaks
sum(countOverlaps(suc2co, co35s, ignore.strand = TRUE)) #219
over.suc2.co35 <- subsetByOverlaps(suc2co, co35s, ignore.strand = TRUE)

### without co-10 peaks
sum(countOverlaps(suc2co, co10, ignore.strand = TRUE)) #1 only overlap
sum(countOverlaps(co35s, co10, ignore.strand = TRUE)) #102 overlaps

filt.suc2 <- subsetByOverlaps(suc2co, co10, invert = TRUE, ignore.strand = TRUE)
filt.co35 <- subsetByOverlaps(co35s, co10, invert = TRUE, ignore.strand = TRUE) #this is same as filt.maybe I saved above
sum(countOverlaps(filt.suc2, filt.co35, ignore.strand = TRUE)) #219
over.suc2.co35.filt <- subsetByOverlaps(filt.suc2, filt.co35, ignore.strand = TRUE)

sum(countOverlaps(over.suc2.co35, over.suc2.co35.filt, ignore.strand = TRUE)) #219 - even after removing co-10 peaks, the overlap in peaks stays the same

filt.co35.FIXintersect <- subsetByOverlaps(filt.co35, over.suc2.co35.filt, ignore.strand = TRUE, invert = TRUE)
filt.co35.FIXintersect <- c(filt.co35.FIXintersect, over.suc2.co35.filt)

### overlap significance based on ath genes
over.suc2.co35.test <- newGeneOverlap(unique(filt.suc2$peak_name), unique(filt.co35.FIXintersect$peak_name), genome.size = 27000)
over.suc2.co35.test <- testGeneOverlap(over.suc2.co35.test)
print(over.suc2.co35.test) #not sig.

### Venn Diagram
grid.newpage()
venn.plot <- draw.pairwise.venn(area1 = length(unique(suc2co$peak_name)), #with co-10 peaks
                                area2 = length(unique(co35s$peak_name)), 
                                cross.area = 219,
                                fill = c("black", "white"),
                                c("SUC2CO chip", "CO35S chip"))
grid.newpage()
venn.plot <- draw.pairwise.venn(area1 = length(unique(filt.suc2$peak_name)), #without co-10 peaks
                                area2 = length(unique(filt.co35$peak_name)), 
                                cross.area = 219, 
                                fill = c("black", "white"),
                                c("SUC2CO chip", "CO35S chip"))
grid.newpage()


## in genes falling close to peaks between our suc2co and 35sco
### with co-10 peaks
length(intersect(ids.suc2co, ids.co35s)) #485
over.suc2.co35.gene.test <- newGeneOverlap(ids.suc2co, ids.co35s, genome.size = 27000)
over.suc2.co35.gene.test <- testGeneOverlap(over.suc2.co35.gene.test)
print(over.suc2.co35.gene.test) #sig.

### without co-10 peaks
ids.suc2co.filt <- ids.suc2co[!(ids.suc2co %in% ids.co10)]
ids.co35s.filt <- ids.co35s[!(ids.co35s %in% ids.co10)]

over.suc2.co35.gene.filt.test <- newGeneOverlap(ids.suc2co.filt, ids.co35s.filt, genome.size = 27000)
over.suc2.co35.gene.filt.test <- testGeneOverlap(over.suc2.co35.gene.filt.test)
print(over.suc2.co35.gene.filt.test) #sig.

### Venn Diagram
#### Chip genes only
length(intersect(ids.suc2co, ids.co35s))

grid.newpage()
draw.pairwise.venn(area1 = length(ids.suc2co),
                   area2 = length(ids.co35s),
                   fill = c("black", "white"),
                   cross.area = 485,
                   c("SUC2CO chip", "CO35S chip"))
grid.newpage()

length(intersect(ids.suc2co.filt, ids.co35s.filt))

grid.newpage()
draw.pairwise.venn(area1 = length(ids.suc2co),
                   area2 = length(ids.co35s),
                   fill = c("black", "white"),
                   cross.area = 449,
                   c("SUC2CO chip", "CO35S chip"))
grid.newpage()

#### ALL up DEGs and peaks
venn.plot <- draw.quad.venn(area1 = length(ids.suc2co), #4069 #with co-10 peaks
                                area2 = length(ids.co35s), #2826
                                area3 = length(unique(rna.suc2co.up$TAIR)), #1765
                                area4 = length(unique(rna.35sco.up$TAIR)), #637
                                n12 = length(intersect(ids.suc2co, ids.co35s)), #485
                                n13 = length(intersect(ids.suc2co, unique(rna.suc2co.up$TAIR))), #244
                                n14 = length(intersect(ids.suc2co, unique(rna.35sco.up$TAIR))), #106
                                n23 = length(intersect(ids.co35s, unique(rna.suc2co.up$TAIR))), #370
                                n24 = length(intersect(ids.co35s, unique(rna.35sco.up$TAIR))), #363
                                n34 = length(intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR))), #275
                                n123 = length(intersect(ids.suc2co, intersect(ids.co35s, unique(rna.suc2co.up$TAIR)))), #61
                                n124 = length(intersect(ids.suc2co, intersect(ids.co35s, unique(rna.35sco.up$TAIR)))), #62
                                n134 = length(intersect(ids.suc2co, intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR)))), #41
                                n234 = length(intersect(ids.co35s, intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR)))), #152
                                n1234 = length(intersect(intersect(ids.suc2co, ids.co35s), intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR)))), #22
                                fill = c("blue", "red", "gray", "goldenrod"),
                                c("SUC2CO chip", "CO35S chip", "SUC2CO Up-DEGs", "CO35S Up-DEGs"))
grid.newpage()
venn.plot <- draw.quad.venn(area1 = length(ids.suc2co.filt), #with co-10 peaks
                            area2 = length(ids.co35s.filt), 
                            area3 = length(unique(rna.suc2co.up$TAIR)), 
                            area4 = length(unique(rna.35sco.up$TAIR)),
                            n12 = length(intersect(ids.suc2co.filt, ids.co35s.filt)), 
                            n13 = length(intersect(ids.suc2co.filt, unique(rna.suc2co.up$TAIR))), 
                            n14 = length(intersect(ids.suc2co.filt, unique(rna.35sco.up$TAIR))), 
                            n23 = length(intersect(ids.co35s.filt, unique(rna.suc2co.up$TAIR))),
                            n24 = length(intersect(ids.co35s.filt, unique(rna.35sco.up$TAIR))), 
                            n34 = length(intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR))),
                            n123 = length(intersect(ids.suc2co.filt, intersect(ids.co35s.filt, unique(rna.suc2co.up$TAIR)))), 
                            n124 = length(intersect(ids.suc2co.filt, intersect(ids.co35s.filt, unique(rna.35sco.up$TAIR)))), 
                            n134 = length(intersect(ids.suc2co.filt, intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR)))), 
                            n234 = length(intersect(ids.co35s.filt, intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR)))), 
                            n1234 = length(intersect(intersect(ids.suc2co.filt, ids.co35s.filt), intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR)))), 
                            fill = c("blue", "red", "gray", "goldenrod"),
                            c("SUC2CO chip", "CO35S chip", "SUC2CO Up-DEGs", "CO35S Up-DEGs"))
grid.newpage()

intersect(intersect(intersect(ids.suc2co, ids.co35s), intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR))), unique(rna.co.dwn$TAIR))
intersect(intersect(intersect(ids.suc2co.filt, ids.co35s.filt), intersect(unique(rna.suc2co.up$TAIR), unique(rna.35sco.up$TAIR))), unique(rna.co.dwn$TAIR))


# ENRICHMENT
## Chip X Chip
over.suc2co.ci35s.filt.gene <- intersect(ids.suc2co.filt, ids.co35s.filt)
all.suc2co.ci35s.filt.gene <- unique(c(ids.suc2co.filt, ids.co35s.filt))

### With all chip genes as background
# chip.bp <- enrichGO(gene = over.suc2co.ci35s.filt.gene, #no enrichment
#                    keyType = "TAIR", 
#                    OrgDb = org.At.tair.db, 
#                    ont = "BP", 
#                    pAdjustMethod = "BH", 
#                    qvalueCutoff = 0.1, 
#                    readable = TRUE,
#                    universe = all.suc2co.ci35s.filt.gene)
# 
# dotplot(chip.bp, showCategory = 50)
# goplot(chip.bp)

# chip.cc <- enrichGO(gene = over.suc2co.ci35s.filt.gene, #no enrichment
#                    keyType = "TAIR",
#                    OrgDb = org.At.tair.db,
#                    ont = "CC",
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.1,
#                    readable = TRUE,
#                    universe = all.suc2co.ci35s.filt.gene)
# 
# dotplot(chip.cc, showCategory = 50)
# goplot(chip.cc)

chip.mf <- enrichGO(gene = over.suc2co.ci35s.filt.gene, #enriched for transcription-related stuff
                   keyType = "TAIR",
                   OrgDb = org.At.tair.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.1,
                   readable = TRUE,
                   universe = all.suc2co.ci35s.filt.gene)

dotplot(chip.mf, showCategory = 50)
goplot(chip.mf)

## Chip X Up DEGs SUC2CO
over.suc2co.degs.filt.gene <- intersect(ids.suc2co.filt, unique(rna.suc2co.up$TAIR))
all.suc2co.degs.filt.gene <- unique(c(ids.suc2co.filt, unique(rna.suc2co.up$TAIR)))

### Check overlap with circadian genes
circadian.genes.dlr <- c("AT1G65480", #FT
                     "AT5G24470", #PRR5
                     "AT2G46790", #PRR9
                     "AT3G09600", #RVE8
                     "AT1G01060", #LHY
                     "AT2G46830") #CCA1

intersect(ids.suc2co.filt, circadian.genes.dlr) #FT is among peaks
intersect(over.suc2co.degs.filt.gene, circadian.genes.dlr) #FT is among DEGs + Peaks

### Statistical test for overlap
length(over.suc2co.degs.filt.gene) #242
over.suc2.degs.test <- newGeneOverlap(ids.suc2co.filt, 
                                      unique(rna.suc2co.up$TAIR), 
                                      genome.size = 27000)

over.suc2.degs.test <- testGeneOverlap(over.suc2.degs.test)
print(over.suc2.degs.test) #not sig.

grid.newpage()
draw.pairwise.venn(area1 = length(ids.suc2co.filt),
                   area2 = length(unique(rna.suc2co.up$TAIR)),
                   fill = c("black", "white"),
                   cross.area = 242,
                   c("SUC2CO chip", "SUC2CO DEGs"))
grid.newpage()

### Enrichment
rna.bp <- enrichGO(gene = over.suc2co.degs.filt.gene, #enrichment for stress response
                   keyType = "TAIR",
                   OrgDb = org.At.tair.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.1,
                   readable = TRUE,
                   universe = all.suc2co.degs.filt.gene)

dotplot(rna.bp, showCategory = 50)
goplot(rna.bp)

# rna.cc <- enrichGO(gene = over.suc2co.degs.filt.gene, #no enrichment
#                    keyType = "TAIR",
#                    OrgDb = org.At.tair.db,
#                    ont = "CC",
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.1,
#                    readable = TRUE,
#                    universe = all.suc2co.degs.filt.gene)
# 
# dotplot(rna.cc, showCategory = 50)
# goplot(rna.cc)

rna.mf <- enrichGO(gene = over.suc2co.degs.filt.gene, #only 1 cat enriched (oxidoreductase activity) and only 9 genes
                    keyType = "TAIR",
                    OrgDb = org.At.tair.db,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.1,
                    readable = TRUE,
                    universe = all.suc2co.degs.filt.gene)

dotplot(rna.mf, showCategory = 50)
#goplot(rna.mf)

## Chip X Up DEGs 35SCO
### Statistical test for overlap
over.co35s.degs.filt.gene <- intersect(ids.co35s.filt, unique(rna.35sco.up$TAIR))
all.co35s.degs.filt.gene <- unique(c(ids.co35s.filt, unique(rna.35sco.up$TAIR)))

### Statistical test for overlap
length(over.co35s.degs.filt.gene) #242
over.co35s.degs.test <- newGeneOverlap(ids.co35s.filt, 
                                      unique(rna.35sco.up$TAIR), 
                                      genome.size = 27000)

over.co35s.degs.test <- testGeneOverlap(over.co35s.degs.test)
print(over.co35s.degs.test) #not sig.

grid.newpage()
draw.pairwise.venn(area1 = length(ids.co35s.filt),
                   area2 = length(unique(rna.35sco.up$TAIR)),
                   fill = c("black", "white"),
                   cross.area = 242,
                   c("35SCO chip", "35SCO DEGs"))
grid.newpage()

### Enrichment
rna.s.bp <- enrichGO(gene = over.co35s.degs.filt.gene, #enrichment for stress response and light
                   keyType = "TAIR",
                   OrgDb = org.At.tair.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.1,
                   readable = TRUE,
                   universe = all.co35s.degs.filt.gene)

dotplot(rna.s.bp, showCategory = 20)
goplot(rna.s.bp)

# rna.s.cc <- enrichGO(gene = over.co35s.degs.filt.gene, #no enrichment
#                    keyType = "TAIR",
#                    OrgDb = org.At.tair.db,
#                    ont = "CC",
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.1,
#                    readable = TRUE,
#                    universe = all.co35s.degs.filt.gene)
# 
# dotplot(rna.s.cc, showCategory = 20)
# goplot(rna.s.cc)

# rna.s.mf <- enrichGO(gene = over.co35s.degs.filt.gene, #no enrichment
#                    keyType = "TAIR",
#                    OrgDb = org.At.tair.db,
#                    ont = "MF",
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.1,
#                    readable = TRUE,
#                    universe = all.co35s.degs.filt.gene)
# 
# dotplot(rna.s.mf, showCategory = 20)
# goplot(rna.s.mf)

## Chip X Up DEGs all
over.chip.degs.filt.gene <- intersect(unique(c(ids.suc2co.filt, ids.co35s.filt)), unique(c(rna.suc2co.up$TAIR, rna.35sco.up$TAIR)))
all.chip.degs.filt.gene <- unique(c(ids.suc2co.filt, ids.co35s.filt, rna.suc2co.up$TAIR, rna.35sco.up$TAIR))

rna.a.bp <- enrichGO(gene = over.chip.degs.filt.gene, #enrichment for aba response, stress response, and light response
                     keyType = "TAIR",
                     OrgDb = org.At.tair.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1,
                     readable = TRUE,
                     universe = all.chip.degs.filt.gene)

dotplot(rna.a.bp, showCategory = 20)
goplot(rna.a.bp)

# rna.a.cc <- enrichGO(gene = over.chip.degs.filt.gene, #no enrichment
#                    keyType = "TAIR",
#                    OrgDb = org.At.tair.db,
#                    ont = "CC",
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.1,
#                    readable = TRUE,
#                    universe = all.chip.degs.filt.gene)
# 
# dotplot(rna.a.cc, showCategory = 20)
# goplot(rna.a.cc)

# rna.a.mf <- enrichGO(gene = over.chip.degs.filt.gene, #no enrichment
#                    keyType = "TAIR",
#                    OrgDb = org.At.tair.db,
#                    ont = "MF",
#                    pAdjustMethod = "BH",
#                    qvalueCutoff = 0.1,
#                    readable = TRUE,
#                    universe = all.chip.degs.filt.gene)
# 
# dotplot(rna.a.mf, showCategory = 20)
# goplot(rna.a.mf)
