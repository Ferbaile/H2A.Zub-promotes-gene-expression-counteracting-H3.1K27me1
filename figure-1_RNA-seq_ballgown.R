library(ballgown)
library(BiocGenerics)
library(limma)

## Transcriptomic data manipulation with ballgown package.
pheno.data <- read.csv("pheno_data.csv") # Table with samples and conditions
bg.data <- ballgown (dataDir = ".", samplePattern ="sample", 
                     pData=pheno.data)

gene.expression <- gexpr(bg.data)

gene.expression <- read.table(row.names = 1, "~/Library/CloudStorage/OneDrive-UniversityofWarwick/Collaboration-Sevilla-Warwick/2025-10-Genome-Biology-publication/scripts_for_github/data/FPKM_imbibed_seeds_7D_abc_clf_noplastids.tsv", 
                              header = T)

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



##### Generation of a Mean expression. matrix  ###########
## Calculamos la matrix de expresión media. 
head(gene.expression)

wt.0h <- (gene.expression[,1] + gene.expression[,2]+ gene.expression[,3])/3
wt.48h <- (gene.expression[,4] + gene.expression[,5] +  gene.expression[,6])/3
wt.7D <- (gene.expression[,7] + gene.expression[,8] +  gene.expression[,9])/3
abc.7D <- (gene.expression[,10] + gene.expression[,11])/2


mean.expression <- matrix(c(wt.0h,wt.48h, wt.7D, abc.7D),ncol=4)
colnames(mean.expression) <- c("wt_0h", "wt_48h", "WT_7D", "abc_7D")
rownames(mean.expression) <- rownames(gene.expression)
head(mean.expression)


########## DIFFERENTIAL ANALYSIS. LIMMA ##########
log.gene.expression <- log2(gene.expression[,1:11] + 1)
head(log.gene.expression)

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
fc <- 1.5
fc.wt0h.abc7D <- wt0h.abc7D$logFC
pv.wt0h.abc7D <- wt0h.abc7D$adj.P.Val
genes.ids.1 <- rownames(wt0h.abc7D)
names(fc.wt0h.abc7D) <- genes.ids.1

activated.genes.wt0h.abc7D <- genes.ids.1[fc.wt0h.abc7D > log2(fc) & 
                                            pv.wt0h.abc7D < 0.05]
repressed.genes.wt0h.abc7D <- genes.ids.1[fc.wt0h.abc7D < - log2(fc) & 
                                            pv.wt0h.abc7D < 0.05]



############  BOXPLOT CODE  ###############

boxplot(log2(mean.expression[activated.genes.wt0h.abc7D, ] +1) , cex.lab=1.2, cex.axis = 1.1, lwd=3, col= rainbow(20), 
        pch=1, cex = 0.7, main = "Upregulated bmi1abc",outline=F, las = 2, cex.main = 2,ylab="log2(FPKM + 1)", boxwex = 0.65)



