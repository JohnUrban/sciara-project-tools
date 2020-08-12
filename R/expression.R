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

