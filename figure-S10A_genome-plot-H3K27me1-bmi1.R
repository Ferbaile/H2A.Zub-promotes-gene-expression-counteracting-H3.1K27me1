library(limma)
library(dplyr)
library(tidyverse)
library(rtracklayer)

######### DIFFERENTIAL ANALYSIS ########

 ## Matrix from multiBigwigSummary (deepTools)

bw.me1 <- read.table("data/multiBigwigSummary_H3K27me1_1000bp.tsv", header = T)
colnames(bw.me1) <- c("chr" ,"start","end","wt_IP1_me1","bmi1abc_IP1_me1", "bmi1abc_IP2_me1","wt_IP2_me1","seed_48HAI_me1_IP1", "seed_48HAI_me1_IP2")
head(bw.me1)
log.cpm <- log2(bw.me1[,c(4,7,5,6,8,9)]+1)


## Experimental design for differential analysis with LIMMA
experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3)))
colnames(experimental.design) <- c("wt_me1","bmi1abc_me1", "wt_48h_me1")

linear.fit <- lmFit(log.cpm, experimental.design)


## To do the contrasts between genotypes, we have to specify the names of
## both conditions separated by -

contrast.matrix <- makeContrasts(bmi1abc_me1-wt_me1, bmi1abc_me1-wt_48h_me1,
                                 wt_48h_me1-wt_me1, 
                                 levels = c("wt_me1","bmi1abc_me1", "wt_48h_me1"))

## contrasts.fit and eBayes functions are used for computing the FC and 
## p-vaules for each comparison. 

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## With topTable we obtain the p-value and FC. Then, a granges object is 
## created to make the genome plot. 

bmi1.wt <- topTable(contrast.results, number = nrow(log.cpm), coef = 1, sort.by = "none")

bmi1.wt.gr <- GRanges(seqnames = bw.me1$chr, ranges = IRanges(bw.me1$start, bw.me1$end ),strand = NULL, bmi1.wt ) 


## Simplifying the Granges and making it a data.frame
bmi1.wt.me1 <- as.data.frame(bmi1.wt.gr)
bmi1.wt.ord <- rbind(bmi1.wt.me1[bmi1.wt.me1$seqnames == 1, ], bmi1.wt.me1[bmi1.wt.me1$seqnames == 2, ], bmi1.wt.me1[bmi1.wt.me1$seqnames == 3, ], 
                     bmi1.wt.me1[bmi1.wt.me1$seqnames == 4, ], bmi1.wt.me1[bmi1.wt.me1$seqnames == 5, ])

bmi1.wt.ord.simple <- bmi1.wt.ord[ ,c(1,2,3,6,9)]
colnames(bmi1.wt.ord.simple) <- c("chrom","start", "end", "Mval", "pval")
bmi1.wt.ord.simple <- bmi1.wt.ord.simple %>% 
  mutate(chrom = str_c("Chr", chrom))
head(bmi1.wt.ord.simple)


## Function to plot chr as facets in ggplot

facet_chr <- function(chr = chr, rows = NULL, scales = "free_x", 
                      space = "free_x", breaks = 10e6, switch = NULL) 
{
  require(scales)
  list(facet_grid(rows = vars({{ rows }}), cols = vars({{ chr }}), 
                  scales = scales, space = space, switch = switch),
       scale_x_continuous(expand = c(0, 0), limits = c(0, NA), 
                          breaks = scales::breaks_width(breaks), 
                          labels = scales::label_number(scale_cut =
                                                          cut_short_scale())),
       xlab("Genomic coordinates (bp)"))
}

## Data and Function to create heterochromatin ribbon
centromeres <- tibble(chr = str_c("Chr",1:5), position = c(15086045, 
                                                           3607929, 13587786, 3956021, 11725024))
heterochromatin <- read_tsv("data/heterochomatin.bed", 
                            col_names = c("chr","start","end"), col_types = cols())

TAIR10_chrs <- str_c("Chr", 1:5)

geom_heterochromatin <- function(cent = centromeres, hchrom = 
                                   heterochromatin, chr_col = chr, chrs = TAIR10_chrs) 
{
  cent <- cent %>% filter(chr %in% chrs) %>% 
    dplyr::rename({{ chr_col }} := chr)
  hchrom <- hchrom %>% filter(chr %in% chrs) %>% 
    dplyr::rename({{ chr_col }} := chr)
  list(geom_vline(data = cent, aes(xintercept = position), 
                  linetype = "dashed", color = "black"),
       geom_rect(data = hchrom, aes(x = NULL, y = NULL, xmin = start,
                                    xmax = end, ymin = -Inf, ymax = Inf), 
                 color = NA, fill = "black", alpha = 0.1))
}

## Extraction of differential methylated peaks (DMP) using a p-value > 0.01
##and a fold-change bigger than 1.5

dmps <- filter(bmi1.wt.ord.simple, pval < 0.01, abs(Mval) > log2(1.5), 
               chrom %in% TAIR10_chrs)
n_dmps_down <- filter(dmps, Mval < 0) %>% nrow()
n_dmps_up   <- filter(dmps, Mval > 0) %>% nrow()

## Drawing the genome plot

ggplot(bmi1.wt.ord.simple %>% filter(chrom %in% TAIR10_chrs),
       aes(x = start, y = Mval)) + facet_chr(chr = chrom) +
  geom_heterochromatin(chr_col = chrom, chrs = TAIR10_chrs) +
  geom_hline(yintercept = c(0), color = "grey10", linetype = "dashed") +
  geom_point(size = 0.1, shape = 1, alpha = 0.5, color = "grey50") +
  geom_point(data = dmps, color = "#1DC2AF", size = 0.85, shape = 16) +
  coord_cartesian(ylim = c(-2,2)) +
  labs(title = "bmi1_me1 vs WT_7D_me1", y = "M value (log2FC)", 
       subtitle = str_c("green points have a p-value < 0.01 & |Mval| > 
log2(1.5) - (gain:", n_dmps_up, " loss:", n_dmps_down, "), total = ", 
                        nrow(bmi1.wt.ord.simple), " windows"))

