#=============================================================================#
  ## Differential gene expression analysis of RNA-seq datasets using R ##
#=============================================================================#


#Load required library packages
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(cluster)
library(ComplexHeatmap)
library(dendextend)
library(dplyr)
library(tidyverse)



# Load raw counts of each sample and create an object to store these raw counts of each sample
smoc2_fibrosis1 <- read.table("GSM2260469_SMOC2_UUO_1.rawcounts.txt", header = F, sep = "\t")
head(smoc2_fibrosis1)
nrow(smoc2_fibrosis1)   # check the number of rows in a dataframe.
smoc2_fibrosis2 <- read.table("GSM2260470_SMOC2_UUO_2.rawcounts.txt", header = F, sep = "\t")
smoc2_fibrosis3 <- read.table("GSM2260471_SMOC2_UUO_3.rawcounts.txt", header = F, sep = "\t")
smoc2_fibrosis4 <- read.table("GSM2260472_SMOC2_UUO_4.rawcounts.txt", header = F, sep = "\t")
smoc2_normal1 <- read.table("GSM2260466_SMOC2_normal_1.rawcounts.txt", header = F, sep = "\t")
smoc2_normal3 <- read.table("GSM2260467_SMOC2_normal_3.rawcounts.txt", header = F, sep = "\t")
smoc2_normal4 <- read.table("GSM2260468_SMOC2_normal_4.rawcounts.txt", header = F, sep = "\t")

table(smoc2_fibrosis1[,1]==smoc2_normal1[,1])   # check whether the Esembl gene ids are same between two different sample datasets.
table(smoc2_fibrosis1[,2]==0)


# Combine individual sample raw counts together to make a data frame
smoc2_rawcounts <- data.frame()
smoc2_rawcounts <- cbind(smoc2_normal1[,2], smoc2_normal3[,2], smoc2_normal4[,2], smoc2_fibrosis1[,2], smoc2_fibrosis2[,2], smoc2_fibrosis3[,2], smoc2_fibrosis4[,2])
head(smoc2_rawcounts)
nrow(smoc2_rawcounts)


# Add Esembl gene ids into rawcount dataframe (as rownames)
rownames(smoc2_rawcounts) <- smoc2_fibrosis1[,1]
head(smoc2_rawcounts)
colnames(smoc2_rawcounts) <- c("smoc2_normal1", "smoc2_normal3", "smoc2_normal4", "smoc2_fibrosis1", "smoc2_fibrosis2", "smoc2_fibrosis3", "smoc2_fibrosis4")

# Load Esembl gene ID and their corresponding genenames and create an object for them
Mouse_ids.genenames <- read.table("Esembl_ID_genename.txt", header = T, sep = "\t")
head(Mouse_ids.genenames)
nrow(Mouse_ids.genenames)


# Create genotype vector
genotype <- c("smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe")

# Create condition vector
condition <- c("normal", "normal", "normal", "fibrosis", "fibrosis", "fibrosis", "fibrosis")

# Create data frame
smoc2_metadata <- data.frame(genotype, condition)

# Assign the row names of the data frame
rownames(smoc2_metadata) <- c("smoc2_normal1", "smoc2_normal3", "smoc2_normal4", "smoc2_fibrosis1", "smoc2_fibrosis2", "smoc2_fibrosis3", "smoc2_fibrosis4")
head(smoc2_metadata)
smoc2_metadata


## DESeq2 analysis ##

dds <- DESeqDataSetFromMatrix(countData = smoc2_rawcounts[rowSums(smoc2_rawcounts)>0,], colData = smoc2_metadata, design = ~condition)
head(dds)
DESeq2.results <- results(DESeq(dds))
head(DESeq2.results)
nrow(DESeq2.results)

DESeq2.filter <- DESeq2.results[(abs(DESeq2.results$log2FoldChange)>1) & !is.na(DESeq2.results$padj),]
head(DESeq2.filter)
nrow(DESeq2.filter)
DESeq2.filter2 <- DESeq2.filter[(DESeq2.filter$padj<0.05),]
nrow(DESeq2.filter2)

DESeq2.filter2.sorted <- DESeq2.filter2[order(DESeq2.filter2$log2FoldChange, decreasing = T),]
head(DESeq2.filter2.sorted)
tail(DESeq2.filter2.sorted)
nrow(DESeq2.filter2.sorted)
table(DESeq2.filter2.sorted$padj<0.05)



## Adding gene annotation using merge() function ##

# merge function works fast.

?merge
merge(df1, df2, by.x='colname in df1', by.y='colname in df2')

df1 <- data.frame(Gene_ID=rownames(DESeq2.filter2.sorted), DESeq2.filter2.sorted[,1:6])
head(df1)

df2 <- Mouse_ids.genenames
head(df2)

df <- merge(df1, df2, by.x="Gene_ID", by.y="Gene.stable.ID")[,c(1,8,2:7)]
df.sorted <- df[order(df$log2FoldChange, decreasing = T),]
head(df)
nrow(df)

# Sort df
df.sorted <- df[order(df$log2FoldChange, decreasing = T),]
rownames(df.sorted) <- NULL
head(df.sorted)
tail(df.sorted)
nrow(df.sorted)


# Save the analyzed results

write.csv(df.sorted, "Normal_vs_Fibrosis_DESeq2.csv", row.names = F)
write.table(df.sorted, "Normal_vs_Fibrosis_DESeq2.txt", quote = F, sep = "\t", row.names = F)



#Load required library packages
library(edgeR)
library(limma)
library(GenomicFeatures)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(GO.db)
library(BiasedUrn)
library(cluster)
library(ComplexHeatmap)


## Filtering to remove lowly expressed genes ##

# use the cpm function from the edgeR library (M D Robinson, McCarthy, and Smyth 2010) to generate the
# CPM values and then filter. Note that by converting to CPMs we are normalising for the different sequencing
# depths for each sample.

# Obtain CPMs
myCPM <- cpm(smoc2_rawcounts)
# Have a look at the output
head(myCPM)


# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)


# Summary of how many TRUEs there are in each row
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
smoc2_rawcounts.keep <- smoc2_rawcounts[keep,]
summary(keep)
dim(smoc2_rawcounts.keep)
head(smoc2_rawcounts.keep)
nrow(smoc2_rawcounts.keep)

smoc2_rawcounts.keep <- data.frame(Gene_ID=rownames(smoc2_rawcounts.keep), smoc2_rawcounts.keep[,1:7])
rownames(smoc2_rawcounts.keep) <- NULL


smoc2.df <- merge(smoc2_rawcounts.keep, Mouse_ids.genenames, by.x="Gene_ID",
                  by.y="Gene.stable.ID")
head(smoc2.df)
nrow(smoc2.df)

# Use the following script to remove duplication under the "Gene.name" column
smoc2.df <- smoc2.df %>% distinct(Gene.name, .keep_all = TRUE)

rownames(smoc2.df) <- smoc2.df[,9]
smoc2.df <- data.frame(smoc2.df[,2:8])
smoc2_rawcounts.keep <- smoc2.df
head(smoc2_rawcounts.keep)



# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],smoc2_rawcounts[,1])

plot(myCPM[,1],smoc2_rawcounts[,1],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5)

plot(myCPM[,2],countdata[,2],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, h=10, col = "blue")


## Convert counts to DGEList object ##

# Next create a DGEList object. This is an object used by edgeR to store count data. It has a number of slots
# for storing various parameters about the data.
y <- DGEList(smoc2_rawcounts.keep)

# have a look at y
y


# See what slots are stored in y
names(y)


# Library size information is stored in the samples slot
y$samples



## Quality control ##

# First, we can check how many reads we have for each sample in the y .
y$samples$lib.size


# We can also plot the library sizes as a barplot to see whether there are any major discrepancies between the samples more easily.

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size, names=colnames(y), las=2)
title("Barplot of library sizes")


# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
head(logcounts)
nrow(logcounts)

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")


# Multidimensional scaling plots
plotMDS(y)

levels(smoc2_metadata$condition)
str(smoc2_metadata)

## Let's choose red for fibrosis and green for normal
col.condition <- c("red","green")[smoc2_metadata$condition]
col.condition
data.frame(smoc2_metadata$condition,col.condition)

# Redo the MDS with cell type colouring
plotMDS(y,col=col.condition)

# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("red","green"),legend=levels(smoc2_metadata$condition))
# Add a title
title("condition")


# MDS plots using different PCs
char.condition <- c(1, 4)[smoc2_metadata$condition]   # assign a number to each level. The number represents the point character (1:circle; 4:cross)
char.condition
plotMDS(y, dim=c(1,2), col=col.condition, pch=char.condition, cex=2)
legend("topleft",fill=c("red","green"),legend=levels(smoc2_metadata$condition))
title("condition")

plotMDS(y, dim=c(2,3), col=col.condition, pch=char.condition, cex=2)
legend("topleft",fill=c("red","green"),legend=levels(smoc2_metadata$condition))
title("condition")

plotMDS(y, dim=c(3,4), col=col.condition, pch=char.condition, cex=2)
legend("topleft",fill=c("red","green"),legend=levels(smoc2_metadata$condition))
title("condition")



## Hierarchical clustering with heatmaps ##

# We estimate the variance for each row in the logcounts matrix
head(logcounts)
var_genes <- apply(logcounts, 1, var)
head(var_genes)



# Get the gene names for the top 100 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var, 20)


# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
nrow(highly_variable_lcpm)


# Plot heatmap using Heatmap (ComplexHeatmap package)
?Heatmap
fontsize <- 0.5
Heatmap(highly_variable_lcpm,
        row_names_side = "left",
        column_names_side = "top",
        row_names_gp = gpar(cex=fontsize),
        column_names_gp = gpar(cex=0.7),
        row_dend_width = unit(3, "cm"),
        column_dend_height = unit(1.5, "cm"))




## Normalisation for composition bias ##

# Apply normalisation to DGEList object
y <- calcNormFactors(y)
y$samples
str(y$samples)


# Plot logcounts
par(mfrow=c(1,2))
plotMD(logcounts,column = 1)
abline(h=0,col="grey")
plotMD(logcounts,column = 5)
abline(h=0,col="grey")

# Plot normalized logcounts
par(mfrow=c(1,2))
plotMD(y,column = 1)
abline(h=0,col="grey")
plotMD(y,column = 5)
abline(h=0,col="grey")


par(mfrow=c(1,2))
plotMD(logcounts,column = 1)
abline(h=0,col="grey")
plotMD(y,column = 1)
abline(h=0,col="grey")



## Differential expression with limma-voom ##

group <- smoc2_metadata$condition
group
group <- factor(group)

# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design
colnames(design) <- c("fibrosis", "normal")

# Voom transform the data #
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
v

# What is contained in this object?
names(v)
v$targets
str(v$design)
y$samples
v$targets==y$samples
dim(v$weights)
class(v)


# We can repeat the box plots for the normalised data to compare to before normalisation. The expression values in
# v$E are already log2 values so we don’t need to log-transform.
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")

## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")

## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")



# Fit the linear model
fit <- lmFit(v)
names(fit)
colnames(fit$design) <- c("fibrosis", "normal")
fit$design
str(fit$design)


# The comparison of interest can be specified using the makeContrasts function.
cont.matrix <- makeContrasts(FIB.vs.N=fibrosis - normal,levels=design)
cont.matrix

# Now we can apply the contrasts matrix to the fit object to get the statistics and estimated parameters of our
# comparison that we are interested in. Here we call the contrasts.fit function in limma.
colnames(fit$coefficients) <- c("fibrosis", "normal")
head(fit$coefficients)

fit.cont <- contrasts.fit(fit, cont.matrix)


# The final step is to call the eBayes function, which performs empirical Bayes shrinkage on the variances, and
# estimates moderated t-statistics and the associated p-values.
fit.cont <- eBayes(fit.cont)


# Check the dimensions of the fit object
dim(fit.cont)
fit.cont
names(fit.cont)

# We can use the limma decideTests function to generate a quick summary of DE genes for the contrasts.
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
head(summa.fit)
head(row.names(summa.fit))
a <- row.names(summa.fit)
head(a)



## This will give all genes with abs(logFC)>=1 and adj p value < 0.05 and be ordered in a decreasing manner
FIB.vs.N <- topTable(fit.cont, coef=1, sort.by="p", number = nrow(fit.cont))
FIB.vs.N <- FIB.vs.N[abs(FIB.vs.N$logFC)>=1&FIB.vs.N$adj.P.Val<0.05, ]
FIB.vs.N <- FIB.vs.N[order(FIB.vs.N$logFC, decreasing = TRUE),]
head(FIB.vs.N)
FIB.vs.N <- data.frame(Gene_symbol=rownames(FIB.vs.N), FIB.vs.N[,1:6])
rownames(FIB.vs.N) <- NULL

write.csv(FIB.vs.N, file = "Fibrosis_vs_Normal_DE_analysis.csv", row.names = F)



## Plots after testing for DE ##

# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"FIB.vs.N"], values = c(-1, 1))

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=30,names=rownames(fit.cont$coefficients))




## Testing relative to a threshold (TREAT) ##

fit.treat <- treat(fit.cont,lfc=1)
res.treat <- decideTests(fit.treat)
summary(res.treat)
head(res.treat)

topTable(fit.treat,coef=1,sort.by="p")



# Notice that much fewer genes are highlighted in the MAplot
par(mfrow=c(1,2))
plotMD(fit.treat,coef=1,status=res.treat[,"FIB.vs.N"], values=c(-1,1))
abline(h=0,col="grey")

volcanoplot(fit.treat,coef=1,highlight=30,names=rownames(fit.treat$coefficients))