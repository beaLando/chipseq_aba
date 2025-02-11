# LIBRARIES
library(Gviz)
library(org.At.tair.db)
library(AnnotationDbi)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(tidyverse)
library(xlsx)

options(ucscChromosomeNames=FALSE)

# IMPORT ANNOTATION FILE FOR TAIR 10
txdb <- TxDb.Athaliana.BioMart.plantsmart51
seqlevels(txdb) <- c("NC_003070.9", "NC_003071.7", "NC_003074.8", "NC_003075.7", "NC_003076.8", "NC_037304.1", "NC_000932.1")
annot.track <- GeneRegionTrack(txdb, 
                               fill = "white", 
                               col = "black", 
                               name = "TAIR10 genes")

axisTrack <- GenomeAxisTrack()

plotTracks(c(axisTrack, annot.track), 
           chromosome = "NC_003070.9", from = 24327000, to = 24337000)


# IMPORT target regions BED FILE
targets.ft <- read.xlsx("./data/FTprimers.xlsx", 
                         sheetIndex=2)

str(targets.ft)
targets.ft <- targets.ft[-c(nrow(targets.ft)),]

targets.ft <- targets.ft %>%
  mutate(Chr = case_when(Chr=="1" ~ "NC_003070.9",
                         Chr=="2" ~ "NC_003071.7",
                         Chr=="3" ~ "NC_003074.8",
                         Chr=="4" ~ "NC_003075.7",
                         Chr=="5" ~ "NC_003076.8",
                         Chr=="Mt" ~ "NC_037304.1",
                         Chr=="Pt" ~ "NC_000932.1")) %>%
  as.data.frame()

targets.track <- GeneRegionTrack(GRanges(seqnames = targets.ft$Chr,
                                         ranges = IRanges(targets.ft$Start, targets.ft$End), 
                                        strand = "+"),
                                 fill = "white",
                                 col = "black",
                                 name = "Target Regions on FT")

plotTracks(c(axisTrack, annot.track, targets.track), 
           #sizes = c(0.1,0.10,0.05),
           #ylim=c(0, 5), #activate for not scaled
           chromosome = "NC_003070.9", from = 24328000, to = 24337000)
