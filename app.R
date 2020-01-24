#----------------------------------------------------------#
       ## Shiny App for DESeq2 Analysis of RNA-seq ##
#----------------------------------------------------------#

# Download app.R and fibrosis.RData for the demonstration of this shiny App. 
# Put both files in a same folder.

# Load required R packages for this App
library(shiny)
library(shinythemes)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(annotables)
library(DT)


# Prepare metadata and rawcounts from RNA-seq and save them together as a RData file.
# Put RData together with app.R in a same folder before run this App.
# Load saved RData containing metadata and rawcounts.
load("fibrosis.RData")


# Assign the condition
colnames(metadata) <- "condition"
cdn <- colnames(metadata)
conditions <- sapply(unique(metadata$condition), toString)
cdn_1 <- conditions[1]
cdn_2 <- conditions[2]

# Create color schemes
brewer_pal <- rownames_to_column(brewer.pal.info, var= "color_pal")


# The following code is for generating DESeq2 analyzed data
dge <- DESeqDataSetFromMatrix(countData = rawcounts,
                                 colData =  metadata,
                                 design = ~ condition)
dge <- estimateSizeFactors(dge)
dge_normalized <- counts(dge, normalized = TRUE)
vsd <- vst(dge, blind = TRUE)
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)

dge2 <- DESeqDataSetFromMatrix(countData = rawcounts,
                                 colData =  metadata,
                                 design = ~ condition)
de <- DESeq(dge2)
res <- results(de, 
               contrast = c(cdn, cdn_2, cdn_1), 
               alpha = 0.05)
res <- lfcShrink(de, 
                 contrast =  c(cdn, cdn_2, cdn_1),
                 res = res)
res <- results(de, 
               contrast = c(cdn, cdn_2, cdn_1), 
               alpha = 0.05, 
               lfcThreshold = 0.32)
res <- lfcShrink(de, 
                 contrast = c(cdn, cdn_2, cdn_1), 
                 res = res)



## The following code is the Shiny script for displaying interactive results

# Define UI for DESeq2 application
ui <- fluidPage(

    # Application title
    titlePanel("Fibrosis RNA-Seq Data Analysis using DESeq2"),
    
    theme = shinythemes::shinytheme("sandstone"),

    # Sidebar with select inputs for threshold, gene annotation, color scheme and the number of top-rank genes
    sidebarLayout(
        sidebarPanel(
            selectInput(
                inputId = "padj_th",
                label = "Threshold for Adjusted P Value",
                choices = c(0.05, 0.01, 0.005, 0.001)
            ),
            selectInput(
                inputId = "annotation",
                label = "Select Species for Gene Annotation",
                choices = c("Human build 38", "Human build 37", "Mouse", "Rat", "Chicken", "Worm", "Fly")
            ),
            selectInput(
                inputId = "brewer_pal",
                label = "Select a Color Scheme",
                choices = unique(brewer_pal$color_pal)
            ),
            sliderInput(
                inputId = "gene_n",
                label = "The Number of Top Rank Genes",
                min = 1, max = 100, value = 10
            )
        ),

        # The following code is for the output of plots and tables in UI
        mainPanel(
           tabsetPanel(
               tabPanel("Sample Heatmap", plotOutput("s_heatmap")),
               tabPanel("PCA Plot", plotOutput("pca")),
               tabPanel("Dispersion Plot", plotOutput("disp")),
               tabPanel("MA Plot", plotOutput("ma")),
               tabPanel("Volcano Plot", plotOutput("volcano")),
               tabPanel("Hierarchical Heatmap", plotOutput("h_heatmap")),
               tabPanel("DE Table", DT::DTOutput("table"))
           )
        )
    )
)

# Define server logic required to generate plots and tables
server <- function(input, output) {
    res_all <- reactive({
        annot <- function (ann) {
            if (ann == "Human build 38") {
                grch38
            } else if (ann == "Human build 37") {
                grch37
            } else if (ann == "Mouse") {
                grcm38
            } else if (ann == "Rat") {
                rnor6
            } else if (ann == "Chicken") {
                galgal4
            } else if (ann == "Worm") {
                wbcel235
            } else {
                bdgp6
            }
        }
        gene_annot <- annot(input$annotation)               
        df <- data.frame(res) %>% rownames_to_column(var = "ensgene")
        left_join(x = df,
                  y = gene_annot[, c("ensgene", "entrez", "symbol", "description")],
                  by = "ensgene")
    })
    
    output$s_heatmap <- renderPlot({
        pheatmap(vsd_cor, annotation = select(metadata, condition))  ### Sample Heatmap
    })
    output$pca <- renderPlot({
        plotPCA(vsd, intgroup=cdn)  ### PCA Plot
    })
    output$disp <- renderPlot({
        plotDispEsts(de)   ### Dispersion Plot
    })
    output$ma <- renderPlot({
        plotMA(res, ylim=c(-8,8))  ### MA Plot
    })
    output$volcano <- renderPlot({
        res_all2 <- data.frame(res_all()) %>% mutate(threshold = padj < input$padj_th) %>% arrange(padj)
        top_genes <- res_all2[1:input$gene_n,]
        ggplot(res_all2) +  ### Volcano Plot
            geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
            xlab("log2 fold change") + 
            ylab("-log10 adjusted p-value") +
            geom_text(data=top_genes, aes(x = log2FoldChange, y = -log10(padj),label=symbol), hjust = 0, nudge_x = 0.1) +
            theme(legend.position = "none", 
                  plot.title = element_text(size = rel(1.5), hjust = 0.5), 
                  axis.title = element_text(size = rel(1.25)))
    })
    output$h_heatmap <- renderPlot({
        res_sig <- res_all() %>% filter(padj < input$padj_th) %>% arrange(padj)
        sig_norm_counts <- dge_normalized[res_sig$ensgene, ]
        heat_colors <- brewer.pal(n = 6, name = input$brewer_pal)
        pheatmap(sig_norm_counts,   ### Hierarchical Heatmap
                 color = heat_colors, 
                 cluster_rows = T, 
                 show_rownames = F,
                 annotation = select(metadata, condition), 
                 scale = "row")
    })
    output$table <- DT::renderDT({
            res_all() %>% filter(padj < input$padj_th) %>% arrange(padj)  ### DE Table
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
