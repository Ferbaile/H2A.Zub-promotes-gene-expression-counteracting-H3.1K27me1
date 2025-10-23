library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

################ WT 48HAI H3K27me3 IP1 #################
wt.me3.IP1.peaks <- readPeakFile(peakfile = 
                                   "seeds_me1/seed_48HAI_me1_IP1_peaks.narrowPeak", header=FALSE)
wt.me3.IP1.peakAnno <- annotatePeak(peak = wt.me3.IP1.peaks, 
                                    addFlankGeneInfo = T, flankDistance = 1,
                                    tssRegion=c(-500, 500), TxDb=txdb)

## Function to make pie charts
plotAnnoPie(wt.me3.IP1.peakAnno)

## Data frame conversion of the peak file generated and target genes extraction 

wt.me3.IP1.annotation <- as.data.frame(wt.me3.IP1.peakAnno)
target.genes <- unique(c(wt.me3.IP1.annotation[grep("Exon", 
                                                    wt.me3.IP1.annotation$annotation), "geneId"],
                         wt.me3.IP1.annotation[grep("Promoter", 
                                                    wt.me3.IP1.annotation$annotation), "geneId"]))

write(x = target.genes.wt.me3.IP1,file = "wt/wt_me3_IP1_target_genes.txt")

## Annotation of genes covered by broad peaks.
write(wt.me3.IP1.annotation[wt.me3.IP1.annotation$width > 2000, ]
      $flank_geneIds, file = "wt_me3_overlap_genes_IP1.txt")

wt.me3.multiple.overlap.IP1 <- read.table(file = 
                                            "wt_me3_overlap_genes_IP1.txt")[[1]]
wt.me3.multiple.overlap.IP1 <- unique(wt.me3.multiple.overlap.IP1)

wt.me3.IP1.genes.complete <- union(wt.me3.multiple.overlap.IP1,target.genes)

## After obtaining target genes for IP2, we intersect both lists and obtain 
## the marked genes in both replicates that we used in this study. 
wt.me3.total.genes <- intersect(wt.me3.IP1.genes.complete, 
                                wt.me3.IP2.genes.complete)



########## venn diagram function ########
pairwise.venn.diagram <- function(set1, set2, set1.name,set2.name)
{
  grid.newpage()
  draw.pairwise.venn(area1 = length(set1),
                     area2 = length(set2),
                     cross.area = length(intersect(set1,set2)),
                     category = c(set1.name, set2.name), fontface = c(2,2,2),
                     fontfamily = rep("sans", 3), alpha = c(0.4,0.4),
                     fill=c("blue2","blue2"), cat.pos = c(-150,150), 
                     cat.dist = rep(0.05, 2), cat.cex= rep(1.2,2),
                     cat.fontfamily = rep("sans", 2),
                     lwd= rep(4, 2), cex = rep(1.8,3) )
  
}



                                
                                
                                
                                