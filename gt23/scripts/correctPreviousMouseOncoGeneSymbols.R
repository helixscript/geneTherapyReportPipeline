f <- read.table('mouseOncoGeneSymbolLookup', sep = '\t', header = TRUE)
o <- gt23::mouseOncoGenesList

write(date(), file = 'log', append = FALSE)
r <- unlist(lapply(o, function(x){
  if(toupper(x) %in% f$Input){
    i <- match(toupper(x), toupper(f$Input))
    write(paste(x, '-> ', f[i,]$Symbol), append = TRUE, file = 'log')
    return(f[i,]$Symbol)
  } else {
    return(x)
  }
}))

rm(f, o)