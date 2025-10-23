library(ballgown)
library(BiocGenerics)

## Transcriptomic data manipulation with ballgown package.
pheno.data <- read.csv("pheno_data.csv") # Table with samples and conditions
bg.data <- ballgown (dataDir = ".", samplePattern ="sample", 
                     pData=pheno.data)

gene.expression <- gexpr(bg.data)
## Name the columns of gene.expression with colnames() with sample names

## Upper quantiles normalization
upper.quantiles <- vector(mode="numeric" ,length=ncol(gene.expression))

for(i in 1:ncol(gene.expression))
{
  upper.quantiles[i] <- quantile(gene.expression[,i],probs=0.75)
}

## Column division by its superior quantile
for(i in 1:ncol(gene.expression))
{
  gene.expression[,i] <- gene.expression[,i] / upper.quantiles[i]
}

########## DIFFERENTIAL ANALYSIS. LIMMA ##########

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4)))
colnames(experimental.design) <- c("wt_0h", "wt_48h", "WT_7D", "abc_7D")

linear.fit <- lmFit(log.gene.expression, experimental.design)
contrast.matrix <- makeContrasts(abc_7D-wt_0h, abc_7D-wt_48h, abc_7D-WT_7D,
                                 levels=c("wt_0h", "wt_48h", "WT_7D", "abc_7D"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

wt0h.abc7D <- topTable(contrast.results, number = nrow(gene.expression), 
                       coef = 1, sort.by = "logFC") #Note that is sorted by LogFC

############## bmi1 vs. WT_0h ################
fc.wt0h.abc7D <- wt0h.abc7D$logFC
pv.wt0h.abc7D <- wt0h.abc7D$adj.P.Val
genes.ids.1 <- rownames(wt0h.abc7D)
names(fc.wt0h.abc7D) <- genes.ids.1

activated.genes.wt0h.abc7D <- genes.ids.1[fc.wt0h.abc7D > log2(fc) & 
                                            pv.wt0h.abc7D < 0.05]
repressed.genes.wt0h.abc7D <- genes.ids.1[fc.wt0h.abc7D < - log2(fc) & 
                                            pv.wt0h.abc7D < 0.05]



