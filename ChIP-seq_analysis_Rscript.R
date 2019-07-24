#============================================#
       ## ChIP-seq analysis in R ##
#============================================#

# Load ChIPpeakAnno
library(ChIPpeakAnno)

# Read the input H3K27ac ChIP-seq data (the output peak table files from MACS2)
#setwd("C:/ChIP-seq_data/MACS2/ChIP-seq")
Control.df <- read.table("Control.txt", header=T)
head(Control.df)

KD.df <- read.table("KD.txt", header = T)
head(KD.df)


## Convert peaks into a GRanges object
# Load the GenomicRanges
library(GenomicRanges)

# Create a GRanges object for Control.df
Control <- GRanges(Control.df$chr, IRanges(Control.df$start, Control.df$end),
                    strand = "*")
head(Control)

# Add more data into the Control GRanges object
names(Control) <- paste("CTRL_peak", 1:nrow(Control.df), sep = "-")
score(Control) <- Control.df$X.log10.qvalue.
head(Control)
tail(Control)
length(Control)


# Create a GRanges object for KD.df
KD <- GRanges(KD.df$chr, IRanges(KD.df$start, KD.df$end), strand = "*")
head(KD)

# Add more data into the Control GRanges object
names(KD) <- paste("KD_peak", 1:nrow(KD.df), sep = "-")
score(KD) <- KD.df$X.log10.qvalue.
head(KD)
tail(KD)
length(KD)



#--------------------------------------------------#
     ## Basic analysis of ChIP-seq peaks ##
#--------------------------------------------------#

## Analysis of mean, median, and max size of th peaks
# The function width() returns the sizes of all ranges in a GRanges object as vector. 
Control_peak_size <- width(Control)
summary(Control_peak_size)

KD_peak_size <- width(KD)
summary(KD_peak_size)


## Remove peaks with size larger than 10K bp
# A GRanges object can be indexed using the [] operator similar than vectors. 
Control <- Control[width(Control) <= 5000]
KD <- KD[width(KD) <= 5000]


## Analysis of the distribution of peak sizes
hist(width(Control), xlab="Control peak size", col="gray", ylim=c(0,12000), 
     breaks = 50)
hist(width(KD), xlab="Knockdown peak size", col="gray", ylim=c(0,14000),
     breaks = 50)


## Analysis of the distribution of ChIP-seq peak p-values
# extract the -log_10 transformed p-values from the score column
Control_Pvals <- score(Control)
KD_Pvals <- score(KD)

hist(Control_Pvals, xlab="-log10(p-value)", col="gray", xlim=c(0,800), 
     ylim=c(0,60000), breaks = 200)
hist(KD_Pvals, xlab="-log10(p-value)", col="gray", xlim=c(0,800), 
     ylim=c(0,60000), breaks = 200)


## Compare the peaks of Control and Knockdown

# Barplot analysis of the number of peaks
bplot <- barplot(c(length(Control), length(KD)), names=c("Control", "Knockdown"),
                 width=0.2, xlim=c(0,0.8), ylim=c(0,80000))
bplot
text(bplot, c(length(Control), length(KD)), labels=c(length(Control), 
            length(KD)), pos=1)



## Analysis of Control peaks overlapping Knockdown peaks
# compute overlaps
overlap_hits <- findOverlaps(Control, KD)
overlap_hits


# get the subsets of overlapping peaks
ovlp1 <- subsetByOverlaps(Control, KD)
Control
ovlp1
length(ovlp1)

ovlp2 <- subsetByOverlaps(KD, Control)
KD
ovlp2
length(ovlp2)

# compute the percentage
length(ovlp1) / length(Control) * 100
length(ovlp2) / length(KD) * 100


## Venn-diagram analysis of the overlap of Control and Knockdown peaks
# Use the function "setdiff" to get the subset of peaks that are unique to Control or Knockdown
Control_uniq <- Control[!names(Control) %in% names(ovlp1),]
KD_uniq <- KD[!names(KD) %in% names(ovlp2),]

Control_uniq
length(Control_uniq)
KD_uniq
length(KD_uniq)



# plot unique and overlap peaks of Control and Knockdown
#install Rtools (https://cran.r-project.org/bin/windows/Rtools/)
#install.packages("remotes")
#library(remotes)
#remotes::install_github("js229/Vennerable")
library(Vennerable)

venn <- Venn(SetNames=c("Control", "Knockdown"),
             Weight=c(
               "10"=length(Control_uniq),
               "11"=length(ovlp1),
               "01"=length(KD_uniq)
             )
          )
plot(venn)



#-------------------------------------------------------#
     ## Functional annotation of ChIP-seq peaks ##
#-------------------------------------------------------#

# Create a TxDb object from TxDb.Hsapiens.UCSC.hg38.knownGene
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("RSQLite")
# BiocManager::install("curl")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


## Analyzing overlap of ChIP-seq peaks with genomic features like introns, exons, 5'UTR
# Use the ChIPpeakAnno package. It has the function "assignChromosomeRegion" for 
# calculating these overlaps by using a annotation database.

library(ChIPpeakAnno)

# Analyzing the overlap of peaks with genomic features 
Control.features <- assignChromosomeRegion(Control, TxDb=txdb, nucleotideLevel=FALSE) 
KD.features <- assignChromosomeRegion(KD, TxDb=txdb, nucleotideLevel=FALSE) 


# show the results 
Control.features
KD.features


# Plot the percentages of peaks overlapping each feature as pie-chart
pie(Control.features$percentage, col = c("purple", "violetred1", "green3", "yellow", 
          "cyan", "white", "red"), main="Control_H3K27ac")
pie(KD.features$percentage, col = c("purple", "violetred1", "green3", "yellow", 
                                         "cyan", "white", "red"), main="Knockdown_H3K27ac")



# Plot the percentages of peaks overlapping each feature as barplots
bplot2 <- barplot(Control.features$percentage, ylab="percent", ylim=c(0,80), 
                  main="Control_H3K27ac") 
text(bplot2, Control.features$percentage, signif(Control.features$percentage, 4), pos=1)

bplot3 <- barplot(KD.features$percentage, ylab="percent", ylim=c(0,80),
                  main="Knockdown_H3K27ac") 
text(bplot2, KD.features$percentage, signif(KD.features$percentage, 4), pos=1)


# Identifying the subset of ChIP-seq peaks that overlap gene promoter regions
# get all genes from the data base as GRanges object 
genes <- genes(txdb)
# take the the region arround the gene start as promoter (Use "promoters" function 
# from the IRanges package)
prom <- promoters(genes, upstream=2000, downstream=200) 
prom

# get only those ChIP-seq peaks that overlap gene promoter regions
CTRLatProm <- subsetByOverlaps(Control, prom)
KDatProm <- subsetByOverlaps(KD, prom)

CTRLatProm
KDatProm

# subset size 
length(CTRLatProm)
length(KDatProm)

# percent of all ChIP-seq peaks 
length(CTRLatProm) / length(Control) * 100
length(KDatProm) / length(KD) * 100


# Identifying the subset of genes that have a peak in their promoter regions using 
# the "findOverlaps" function from the IRanges package and write their gene IDs 
# into an output file

# Search for overlap between ChIP-seq peaks and promoters
CTRLatProm.Hits = findOverlaps(Control, prom)
CTRLprom = genes[subjectHits(CTRLatProm.Hits)]

KDatProm.Hits = findOverlaps(KD, prom)
KDprom = genes[subjectHits(KDatProm.Hits)]


# take only unique ids 
Control.gene.ids <- unique(names(CTRLprom))
KD.gene.ids <- unique(names(KDprom))

# write names to an output file 
write.table(Control.gene.ids, file="Control_H3K27ac_genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(KD.gene.ids, file="Knockdown_H3K27ac_genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


## Calculate the average distances between ChIP-seq peaks and the nearest transciption start sites (TSS)
# Make TSS as the only start coordinate of genes. 
tss <- resize(genes, width=1, fix="start")



# calculate the distances from peaks to tss 
CTRL.dist <- distanceToNearest(Control, tss)
KD.dist <- distanceToNearest(KD, tss)

CTRL.dist
KD.dist

# mcols() can returns the metadata column
mcols(CTRL.dist)
mcols(KD.dist)


# calculate the average distance in kb
CTRL_dist <- mcols(CTRL.dist)[,1]
mean(CTRL_dist)*10^-3

KD_dist <- mcols(KD.dist)[,1]
mean(KD_dist)*10^-3


## Identify the peaks nearest to genes within a range of 5 kb
CTRL.nes <- subset(CTRL.dist, mcols(CTRL.dist)[,1] <=5000)
KD.nes <- subset(KD.dist, mcols(KD.dist)[,1] <=5000)

CTRL.nes
KD.nes

# extract the indece of genes
CTRL.genes <- genes[subjectHits(CTRL.nes)]
KD.genes <- genes[subjectHits(KD.nes)]

CTRL.genes
KD.genes

# extract the gene id names from the GRanges object
CTRL.gene.ids <- names(CTRL.genes)
KD.gene.ids <- names(KD.genes)


# get only unique gene ids
CTRL.gene.ids <- unique(CTRL.gene.ids)
KD.gene.ids <- unique(KD.gene.ids)


# Save gene ids names
write.table(CTRL.gene.ids, file="Control_H3K27ac_genes_2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(KD.gene.ids, file="Knockdown_H3K27ac_genes_2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Save Granges objects as BED files
library(rtracklayer)

export.bed(Control, "Control_peaks.bed") 
export.bed(KD, "Knockdown_peaks.bed")






