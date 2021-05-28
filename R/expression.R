writeEdgeR <- function(df, out){
  fdr <- p.adjust(df$PValue,method="BH")
  write.table(
    x = data.frame(names=rownames(df), 
                   df, 
                   fdr = fdr,
                   log10p = -log10(df$PValue),
                   log10fdr = -log10(fdr)), 
    file = out, 
    quote = FALSE, 
    row.names = FALSE, 
    sep = "\t")
}



edgeRead <- function(f, row.names="index"){
  x <- read.delim(f, row.names=row.names)
  # dim(x)
  # dim(x[complete.cases(x),])
  x <- x[complete.cases(x),]
  return(x)
}

deseqRead <- function(f, row.names="index"){
  return(round(edgeRead(f,row.names)))
}


writeDEseq <- function(df, out){
  write.table(
    x = data.frame(names=rownames(df),df), 
    file = out, 
    quote = FALSE, 
    row.names = FALSE, 
    sep = "\t")
}




edgeR_prepare <- function(f, row.names="index", min.count = 1, min.total.count=10, min.prop=0.25, names, levels_order){
  x <- edgeRead(f,row.names)
  
  # group <- factor(c(rep("FA",2),rep("FE",3),rep("FL",2), rep("FP",2), rep("MA",2),rep("ME",3),rep("ML",2),rep("MP",2)), 
  #                 levels = c("FE","FL","FP","FA","ME","ML","MP","MA"))
  group <- factor(names, 
                  levels = levels_order)
  
  y <- DGEList(counts=x,group=group)
  keep <- filterByExpr(y, min.count = min.count, min.total.count=min.total.count, min.prop=min.prop)#; summary(keep)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y); 
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y,design, robust=TRUE)
  
  return(list(y=y, design=design, fit=fit))
}