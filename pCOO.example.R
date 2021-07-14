require(RNOmni)
require(data.table)
require(randomForest)
require(fmsb)
require(Matrix)
################################
# functions
get.ave <- function(sc, cluster) {
  dat <- as.data.table(t(as.matrix(sc)))
  dat$cluster <- cluster
  dat.ave <- setDT(dat)[, lapply(.SD, mean, na.rm = TRUE), by = cluster]
  cl <- dat.ave$cluster
  dat.ave <- as.matrix(dat.ave[, -1])
  rownames(dat.ave) <- cl
  dat.ave <- t(dat.ave)
  return(dat.ave)
}

getHVGs <- function(sc, block = NULL, n=2000) {
  if (is.null(block)) {
    modgene <- scran::modelGeneVar(sc)
  } else {
    modgene <- scran::modelGeneVar(sc, block = block)
  }
  g=scran::getTopHVGs(modgene, n= n)
  return(g)
 
}

radarplot <- function(comb, groups, title, cols=NULL, ma=NULL, mi=NULL) {
  d=t(comb[,groups, drop=F])
  #dd=as.numeric(comb)
  #ma=ceiling(100*max(dd))/100
  #mi=floor(100*min(dd))/100
  if (is.null(ma)) {ma=0.4}
  if (is.null(mi)) {mi=0}
  
  d=rbind(rep(ma, ncol(d)), rep(mi, ncol(d)), d)
  rownames(d)[1:2]=c('max','min')
  
  if (is.null(cols)) {
    cols=scales::hue_pal()(length(groups))
    names(cols)=groups
  }
  
  s=seq(from = mi, to = ma, by = 0.1)
  print(fmsb::radarchart(as.data.frame(d), title=title, pcol= cols[groups],
                         pfcol=adjustcolor(cols[groups], alpha.f = 0.2) , 
                         plwd=1 , plty=1, cglcol="black", seg= length(s)-1,
                         axistype = 1, # label axis
                         axislabcol = 1, #pty = 32,
                         caxislabels = s, centerzero = TRUE
  )
  )
  legend('topright',legend = rownames(d)[-c(1:2)], col=cols[groups], lty=1, bty = 'n')
}

###########################################
# read single cell data matrix
nsc <- readRDS('../data/nsc.raw.rds')
nsc <- do.call(cbind, nsc)

# read annotation file
nsc.anno <- read.csv('../data/normal_anno.csv', stringsAsFactors = F, row.names = 1)
# keep epithelial cells only
nsc.anno <- nsc.anno[grepl('PT|IC|DL|tAL|TAL|DCT|CNT|PC', nsc.anno$label), ]
nsc <- nsc[, rownames(nsc.anno)]

# normalize by total UMI and log transform
nsc <- t(t(nsc)/colSums(nsc))
nsc <- log2(5000*(nsc)+1)

# read bulk RNA-seq data, in TPM
bulk.tpm <- readRDS('../data/tcga.rcc.tpm.rds')
bulk.meta <- bulk.tpm$meta
bulk.tpm <- bulk.tpm$tpm

# keep common genes in both bulk and sc
comm.genes <- intersect(rownames(nsc), rownames(bulk.tpm))
nsc <- nsc[comm.genes, ]
bulk.tpm <- bulk.tpm[comm.genes, ]

# get cluster average expression to filter lowly expressed genes
nsc.ave <- get.ave(nsc, nsc.anno$label)
flt <- apply(nsc.ave, 1, function(x) any(x>0.2))
kept.genes <- rownames(nsc.ave[flt, ])

# inverse normal transformation for both bulk and sc data
# nsc.norm <- apply(nsc[kept.genes,], 2, rntransform) # GenABEL was removed from cran depository
# tpm.norm <- apply(tpm[kept.genes,], 2, rntransform)
nsc.norm <- apply(nsc[kept.genes,], 2, RankNorm)
tpm.norm <- apply(bulk.tpm[kept.genes,], 2, RankNorm)
rownames(nsc.norm) <- kept.genes
rownames(tpm.norm) <- kept.genes

# get highly variable genes
hvgs  <- getHVGs(nsc[kept.genes,], block = nsc.anno$sample)

# the sample sizes for each cell type are not balanced. so randomly select cells for large population 
# run this multiple times and calculate average
sel.i=lapply(unique(nsc.anno$label), function(x) {
  i=which(nsc.anno$label==x)
  set.seed(1)
  sample(i, min(200, length(i)))
})
sel.i=unlist(sel.i)
rf.mod <- randomForest(t(nsc.norm[hvgs, sel.i]), as.factor(nsc.anno$label[sel.i]), 
                  classwt = table(nsc.anno$label[sel.i])/200)
pred <- predict(rf.mod, t(tpm.norm[hvgs, ]), type='prob')
pred.ave <- get.ave(t(pred), bulk.meta$group)

radarplot(comb = pred.ave, groups = colnames(pred.ave), title = 'TCGA RCCs')
