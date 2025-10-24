library(clusterProfiler)
library(org.At.tair.db)
library(rrvgo)

########### GO TERM ENRICHMENT WITH CLUSTERPROFILER ##########

## Using enrichGO function to compute enrichments and generate GO table. 

atha.universe <- unique(select(org.At.tair.db,columns = c("GO"), 
                               keys=keys(org.At.tair.db,keytype = "TAIR"))[["TAIR"]])

enrich.go <- enrichGO(gene = target.genes , ## Gene list input
                      universe = atha.universe, OrgDb = org.At.tair.db,
                      ont = "BP", pAdjustMethod = "BH",pvalueCutoff  = 0.05,
                      qvalueCutoff = 0.05, readable = F, keyType = "TAIR")

enrich.go.result <- as.data.frame(enrich.go)
write.table(enrich.go.result, file = "GOterms/GO_table.txt", sep = "\t", 
            row.names = F,quote = F)

## Now with the simplify method, we can remove redundant terms

simply.go <- simplify(enrich.go, cutoff=0.7, by="p.adjust", select_fun=min, 
                      measure="Wang", semData = NULL)
write.table(simply.go, file = "GOterms/GO_table_SIMPLY.txt", sep = "\t", 
            row.names = F, quote = F)

########   REVIGO   ########
## Reducing the redundancy of GO using rrvgo package

simMatrix <- calculateSimMatrix(enrich.go.result$`GO ID`,
                                orgdb="org.At.tair.db", ont="BP", method="Rel")

## Group terms based on similarity
scores <- setNames(-log10(enrich.go.result$`q-value`), 
                   enrich.go.result$`GO ID`)
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7,
                                orgdb="org.At.tair.db")

## Use reducedTerms matrix to plot GO terms enrichment. Treemaps are space-filling visualization of hierarchical structures

treemapPlot(reducedTerms, title="H3K27me1-bmi1abc genes 5987", 
            size = "score") ## choose score to represent -log10(q-val)



############# HEATMAP GENERATION WITH GO TERMS ####################

# Split the simply.go table that simplify provides to create the annotation with 
# down and upregulated genes in bmi1abc vs. WT7D
## First, for downregulated genes,

simply.go <- as.data.frame(simply.go)
head(simply.go)

df.d.simpl <- cbind.data.frame(intersect(key.terms.d,simply.go$ID),
                               simply.go[simply.go$ID %in% intersect(key.terms.d,simply.go$ID),"geneID" ])
colnames(df.d.simpl) <- c("GO", "genes")
head(df.d.simpl)

setdiff(key.terms.d,simply.go$ID)

df.d.simpl_separado <- df.d.simpl %>%
  mutate(genes = strsplit(genes, "/")) %>%  # Dividir la columna genes
  unnest(genes)  # Expandir en filas

go.down.bmi1.simpl <- as.data.frame(df.d.simpl_separado)
dim(go.down.bmi1.simpl)


## The same for genes upreg in bmi1abc
simply.go <- as.data.frame(simply.go)
head(simply.go)

df.u.simpl <- cbind.data.frame(intersect(key.terms.u,simply.go$ID),
                               simply.go[simply.go$ID %in% intersect(key.terms.u,simply.go$ID),"geneID" ])
colnames(df.u.simpl) <- c("GO", "genes")
head(df.u.simpl)

setdiff(key.terms.u,simply.go$ID)

df.u.simpl_separado <- df.u.simpl %>%
  mutate(genes = strsplit(genes, "/")) %>%  # Divide gene column
  unnest(genes)  # expand in rows

go.u.bmi1.simpl <- as.data.frame(df.u.simpl_separado)
dim(go.u.bmi1.simpl)


## Lets include replicates in the heatmap
go.up.bmi1.expr <- cbind(mean.expression[go.u.bmi1.simpl$genes, ], go.u.bmi1.simpl$GO)
go.up.bmi1.expr <- na.exclude(go.up.bmi1.expr)

go.down.bmi1.expr <- cbind(mean.expression[go.down.bmi1.simpl$genes, ], go.down.bmi1.simpl$GO)
go.down.bmi1.expr <- na.exclude(go.down.bmi1.expr)
colnames(go.up.bmi1.expr)[5] <- "GO_terms"
colnames(go.down.bmi1.expr)[5] <- "GO_terms"

go.bmi1.complete.simply <- rbind(go.up.bmi1.expr,go.down.bmi1.expr)

gene.expression <- read.table(file = "/Volumes/FerBaile/Calonje_ubuntu/Documentos/imbibed_REF6_paper2023/rnaseq_abc_imbibed/FPKM_imbibed_seeds_7D_abc_clf_noplastids.txt")

genes.names <- substr(rownames(go.bmi1.complete.simply), start = 1, stop = 9) 
gene.expression[genes.names, 1:11]


# Create a colour palette for GO terms
library(circlize)
library(ComplexHeatmap)

go_terms <- unique(unlist(strsplit(as.character(go.bmi1.complete$GO_terms), ",")))

# Create annotation for GO terms
anotacion_GO <- as.data.frame(go.bmi1.complete$GO_terms)

rownames(anotacion_GO) <- rownames(go.bmi1.complete)
colnames(anotacion_GO) <- "GO_terms"

# Give colours to the annotation object
unique(go.bmi1.complete$GO_terms)


key.terms.u <- c("GO:0019915", "GO:0010876", "GO:0044242" ,"GO:0006869","GO:0010344",
                 "GO:0010431","GO:0022611","GO:0009845","GO:0010262",
                 "GO:0051093","GO:0045165","GO:2000026")

key.terms.d <- c("GO:0015979","GO:0046148","GO:0010109",
                 "GO:0009639","GO:0009733","GO:0007623","GO:0009637",
                 "GO:0010087","GO:0006833","GO:0052545","GO:0010232","GO:0010345",
                 "GO:0009657","GO:0032544","GO:0010119","GO:0010118","GO:0002376","GO:0010375",
                 "GO:0009704","GO:0010016")

go.terms <- c(key.terms.u, key.terms.d)

ann_colors = list(
  GO_terms = c("GO:0019915" = "#F6F892", "GO:0044242" = "#FAFFD0", "GO:0010431" = "#49FF00", "GO:0022611" = "#AAFFAC" , "GO:0009845" = "forestgreen" , "GO:2000026" = "#CD762F",
               "GO:0010262" = "#7DBB7E","GO:0051093" = "#FF5700","GO:0045165" = "#F5A800","GO:0010876" = "#FFDB00","GO:0006869" = "#F8FF00","GO:0010344" = "#D4CD43",
               "GO:0015979" = "#FFD6D6","GO:0009657" = "#1E00C8","GO:0046148" = "#FF7575","GO:0010109" = "#FD4447","GO:0032544" = "#00E9FF","GO:0009639" = "#FF00D0","GO:0009733" = "#D91560",
               "GO:0010119" = "#0066FF","GO:0010118" = "#8ED4FF","GO:0007623" = "#FF89E2","GO:0009637" = "#FFBCFC","GO:0002376" = "#24A0C2","GO:0010087" = "#7500FF","GO:0006833" = "#A000FF",
               "GO:0010375" = "#75B8FF","GO:0052545" = "#CAA1E4","GO:0009704" = "#272883","GO:0010232" = "#865497","GO:0010345" = "#CC00FF","GO:0010016" = "#BCCFFF")
  
)



library(RColorBrewer)

## DRAWING THE HM
pheatmap(na.exclude(log2(gene.expression[genes.names, 1:11]+1)), 
         scale = "row" , cluster_rows = F , col = rev(brewer.pal(n = 10, name = "RdYlBu")), 
         show_rownames = F, main = "DEGs bmi1abc GOs", cluster_cols = F , fontsize_row = 3,  
         cutree_rows = NA, treeheight_row = 10, border_color = NA, breaks = c(-1.5, 1.5))
