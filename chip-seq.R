## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Instalar chipseeker y paquete de anotación de Arabidopsis thaliana
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")

library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28


## Leer fichero de picos
LEC1.peaks <- readPeakFile(peakfile = "LEC1.narrowPeak",header=FALSE)

## Definir la región que se considera promotor entorno al TSS
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)

## Anotación de los picos
                              tssRegion=c(-1000, 1000),
                              TxDb=txdb)

plotAnnoPie(LEC1.peakAnno)
plotAnnoBar(LEC1.peakAnno)
plotDistToTSS(LEC1.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")
upsetplot(LEC1.peakAnno)

plotPeakProf2(peak = LEC1.peaks, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)

## Convertir la anotación a data frame
LEC1.annotation <- as.data.frame(LEC1.peakAnno)
head(LEC1.annotation)

target.genes <- LEC1.annotation$geneId[LEC1.annotation$annotation == "Promoter"]

write(x = target.genes,file = "LEC1_target_genes.txt")

## Enriquecimiento funcional. 
library(clusterProfiler)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("enrichplot")

library(org.At.tair.db)
library(enrichplot)

LEC1.enrich.go <- enrichGO(gene = target.genes,
                           OrgDb         = org.At.tair.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable      = FALSE,
                           keyType = "TAIR")

barplot(LEC1.enrich.go,showCategory = 10)
dotplot(LEC1.enrich.go,showCategory = 10)
emapplot(pairwise_termsim (LEC1.enrich.go),showCategory = 15, cex_label_category=0.5)
cnetplot(LEC1.enrich.go,showCategory = 15, cex_label_category=0.5, cex_label_gene=0.5)

LEC1.enrich.kegg <- enrichKEGG(gene  = target.genes,
                               organism = "ath",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05)
df.LEC1.enrich.kegg <- as.data.frame(LEC1.enrich.kegg)
head(df.LEC1.enrich.kegg)

## ChIPpeakAnno es un paquete de R de Bioconductor que implementa 
## análisis de los resultados del procesamiento de los datos de ChIP-seq
## alternativo y complementario al presentado anteriormente. 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")
