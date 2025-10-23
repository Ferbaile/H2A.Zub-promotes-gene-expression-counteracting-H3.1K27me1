#################### METAGENEPLOTS ##################

average.signal.metagene <- function(bigwig.file, gene.set,txdb)
{
  ## Extract genes data 
  genes.data <- genes(txdb)
  gene.set.data <- subset(genes.data, gene_id %in% gene.set)
  
  ## TSS data extraction 
  gene.set.tss <- resize(gene.set.data, width=1, fix='start')
  nrow(as.data.frame(gene.set.tss))
  
  ## Centering the TSS upstream and downstream environment to 2kb
  around.genes.tss <- gene.set.tss
  start(around.genes.tss) <- start(gene.set.tss) - 2000
  end(around.genes.tss) <- end(gene.set.tss) + 2000 
  feature.recentered <- around.genes.tss
  
  ## Import the CPM value from .bw files
  cvglists <- sapply(bigwig.file, import, 
                     format="BigWig", 
                     which=feature.recentered, 
                     as="RleList")
  
  ## TSS Signal extraction. (This may be slow)
  sig <- featureAlignedSignal(cvglists, around.genes.tss, 
                              upstream=2000, downstream=2000,n.tile=50) 
  sig.tss <- sig[[1]][rowMeans(sig[[1]]) < 50 ,1:25]
  mean.sig.tss <- colMeans(sig.tss,na.rm = TRUE)
  
  ## TES data extraction 
  gene.set.tes <- resize(gene.set.data, width=1, fix='end')
  nrow(as.data.frame(gene.set.tes))
  
  ## Centering the TSS upstream and downstream environment to 2kb
  around.genes.tes <- gene.set.tes
  start(around.genes.tes) <- start(gene.set.tes) - 2000
  end(around.genes.tes) <- end(gene.set.tes) + 2000 
  feature.recentered <- around.genes.tes
  
  ## Import the CPM value from .bw files
  cvglists <- sapply(bigwig.file, import, 
                     format="BigWig", 
                     which=feature.recentered, 
                     as="RleList")
  
  ## TSS Signal extraction. (This may be slow)
  sig <- featureAlignedSignal(cvglists, around.genes.tes, 
                              upstream=2000, downstream=2000,n.tile=50) 
  sig.tes <- sig[[1]][rowMeans(sig[[1]]) < 50 ,26:50]
  mean.sig.tes <- colMeans(sig.tes,na.rm = TRUE)

  ## Gene body data extraction
  gene.body.signal <- matrix(0,nrow=length(gene.set.data$gene_id),ncol=100)
  
  for(i in 1:length(gene.set.data$gene_id))
  {
    ## Printing current gene and extracting the features for the next two 
    ## genes
    print(i)
    current.gene.data <- gene.set.data[c(i,i),]
    
    ## Getting widths of the features
    widths <- width(current.gene.data)
    
    ## Getting the signal for the current gene data
    cvglists <- sapply(bigwig.file, import, 
                       format="BigWig", 
                       which=current.gene.data, 
                       as="RleList")
    
    
    if(widths[1] > 100)
    {
      sig <- featureAlignedSignal(cvglists, current.gene.data, n.tile=100)
      gene.body.signal[i,] <- sig[[1]][1,]
    }
  }
  
  gene.body.signal <- gene.body.signal[rowMeans(gene.body.signal) < 50,]
  mean.sig.gene.body <- colMeans(gene.body.signal)
  
  ## Merging the signal from the TSS, the gene body and TES
  metagene.signal <- c(mean.sig.tss,mean.sig.gene.body,mean.sig.tes)
  
  metagene.complete.signal <- list(sig.tss, gene.body.signal, sig.tes)
  names(metagene.complete.signal) <- c("TSS", "gene body", "TES")
  
  output <- list(metagene.signal, metagene.complete.signal)
  names(output) <- c("mean signal", "complete signal")
  return(output)
}

## Usage of function with an example 
h3k27me1.me1.wt <- average.signal.metagene(bigwig.file= "h3k27me1_wt.bw", 
                                           gene.set = H3k27me1.genes,
                                           txdb=txdb)

