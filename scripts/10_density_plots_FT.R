library(Gviz)
library(org.At.tair.db)
library(AnnotationDbi)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(tidyverse)

options(ucscChromosomeNames=FALSE)

# IMPORT BMP-scaled PEAK PROFILES for FT REGION
### SUC2CO
#### Peaks
cocit <- DataTrack(range = "./output/mapped/Unknown_CE349-003R000_cocit_AVG.sorted.bw",
                   genome = "tair10",
                   type = "hist",
                   #window = -1, 
                   windowSize = 5,
                   col ='red',
                   fill = "red",
                   #alpha = 0.7,
                   name = "SUC2CO",
                   transformation = function(x) { log(x + 0.00001) },
                   ylim=c(0,2)
                   )

plotTracks(cocit, chromosome = "NC_003070.9", from = 24330000, to = 24335000)

#### Inputs
cocit.input <- DataTrack(range = "./output/mapped/Unknown_CE349-003R000_cocitInput_AVG.sorted.bw", 
                         genome = "tair10", 
                         type = "hist", 
                         #window = -1, 
                         windowSize = 5,
                         col ='black',
                         fill = "black",
                         #alpha = 0.7,
                         name = "SUC2CO",
                         transformation = function(x) { log(x + 0.00001) },
                         ylim=c(0,2)
                         )

cocit2 <- OverlayTrack(trackList = list(cocit.input, cocit))
plotTracks(cocit2, 
           chromosome = "NC_003070.9", from = 24330000, to = 24335000)

### SUC2CO-aba
#### Peaks
cocitaba <- DataTrack(range = "./output/mapped/Unknown_CE349-003R000_cocitaba_AVG.sorted.bw",
                   genome = "tair10",
                   type = "hist",
                   #window = -1, 
                   windowSize = 5,
                   col ='blue',
                   fill = "blue",
                   #alpha = 0.7,
                   transformation = function(x) { log(x + 0.00001) },
                   name = "SUC2CO-aba",
                   ylim=c(0,2)
)

plotTracks(cocitaba, chromosome = "NC_003070.9", from = 24330000, to = 24335000)

#### Inputs
cocitaba.input <- DataTrack(range = "./output/mapped/Unknown_CE349-003R000_cocitabaInput_AVG.sorted.bw", 
                         genome = "tair10", 
                         type = "hist", 
                         #window = -1, 
                         windowSize = 5,
                         col ='black',
                         fill = "black",
                         #alpha = 0.7,
                         transformation = function(x) { log(x + 0.00001) },
                         name = "SUC2CO-aba",
                         ylim=c(0,2)
)

cocitaba2 <- OverlayTrack(trackList = list(cocitaba.input, cocitaba))
plotTracks(cocitaba2, 
           chromosome = "NC_003070.9", from = 24330000, to = 24335000)

plotTracks(c(cocit2, cocitaba2), 
           chromosome = "NC_003070.9", from = 24330000, to = 24335000,
           ylim=c(0.09,1.5))



# IMPORT ANNOTATION FILE FOR TAIR 10
txdb <- TxDb.Athaliana.BioMart.plantsmart51
seqlevels(txdb) <- c("NC_003070.9", "NC_003071.7", "NC_003074.8", "NC_003075.7", "NC_003076.8", "NC_037304.1", "NC_000932.1")
annot.track <- GeneRegionTrack(txdb, 
                               fill = "white", 
                               col = "black", 
                               name = "TAIR10 genes")

plotTracks(c(cocit2, cocitaba2, annot.track), 
           chromosome = "NC_003070.9", from = 24327000, to = 24337000)

axisTrack <- GenomeAxisTrack()

plotTracks(c(axisTrack, cocit2, cocitaba2, annot.track), 
           chromosome = "NC_003070.9", from = 24327000, to = 24337000)


# IMPORT IDR BED FILE (ours) + IMPORT de los Reyes BED FILEù
suc2co.pks <- read.delim("./output/macs2_out/IDR/cocit_idr_gsMask10.bed", 
                           header=FALSE)

suc2co_aba.pks <- read.delim("./output/macs2_out/IDR/cocitaba_idr_gsMask10.bed", 
                           header=FALSE)

s35co.pks <- read.delim("./data/delosReyes2024/GSE222657_35SCO_peaks.narrowpeak", 
                           header=FALSE)

co10.pks <- read.delim("./data/delosReyes2024/GSE222657_co10_peaks.narrowpeak", 
                           header=FALSE)

str(suc2co.pks)
str(suc2co_aba.pks)
str(s35co.pks)
str(co10.pks)

s35co.pks <- s35co.pks %>%
  mutate(V1 = case_when(V1=="1" ~ "NC_003070.9",
                        V1=="2" ~ "NC_003071.7",
                        V1=="3" ~ "NC_003074.8",
                        V1=="4" ~ "NC_003075.7",
                        V1=="5" ~ "NC_003076.8",
                        V1=="Mt" ~ "NC_037304.1",
                        V1=="Pt" ~ "NC_000932.1")) %>%
  as.data.frame()

co10.pks <- co10.pks %>%
  mutate(V1 = case_when(V1=="1" ~ "NC_003070.9",
                        V1=="2" ~ "NC_003071.7",
                        V1=="3" ~ "NC_003074.8",
                        V1=="4" ~ "NC_003075.7",
                        V1=="5" ~ "NC_003076.8",
                        V1=="Mt" ~ "NC_037304.1",
                        V1=="Pt" ~ "NC_000932.1")) %>%
  as.data.frame()


suc2co.track <- GeneRegionTrack(GRanges(seqnames = suc2co.pks$V1,
                                        ranges = IRanges(suc2co.pks$V2, suc2co.pks$V3), 
                                        strand = "*"),
                                fill = "white",
                                col = "black",
                                name = "SUC2CO peaks")

suc2co_aba.track <- GeneRegionTrack(GRanges(seqnames = suc2co_aba.pks$V1,
                                        ranges = IRanges(suc2co_aba.pks$V2, suc2co_aba.pks$V3), 
                                        strand = "*"),
                                fill = "white",
                                col = "black",
                                name = "SUC2CO-aba peaks")

s35co.track <- GeneRegionTrack(GRanges(seqnames = s35co.pks$V1,
                                        ranges = IRanges(s35co.pks$V2, s35co.pks$V3), 
                                        strand = "*"),
                                fill = "white",
                                col = "black",
                                name = "35SCO peaks")

co10.track <- GeneRegionTrack(GRanges(seqnames = co10.pks$V1,
                                        ranges = IRanges(co10.pks$V2, co10.pks$V3), 
                                        strand = "*"),
                                fill = "white",
                                col = "black",
                                name = "co-10 peaks")

plotTracks(c(axisTrack, annot.track, cocit2, cocitaba2, suc2co.track, suc2co_aba.track, s35co.track, co10.track), 
           sizes = c(0.1,0.10,0.15,0.15,0.05,0.05,0.05,0.05),
           chromosome = "NC_003070.9", from = 24328000, to = 24337000)

plotTracks(c(axisTrack, annot.track, cocit2, cocitaba2, suc2co.track, suc2co_aba.track, s35co.track, co10.track), 
           sizes = c(0.1,0.10,0.15,0.15,0.05,0.05,0.05,0.05),
           type = "l",
           ylim=c(0, 1.5),
           lwd = 2,
           chromosome = "NC_003070.9", from = 24328000, to = 24337000)
