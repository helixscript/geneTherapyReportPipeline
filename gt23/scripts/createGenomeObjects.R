library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm9)
library(tidyverse)
set.seed(46)

createRandomFragments <- function(refGenome, n = 1, fragWidth = 1){
  
  min_seqs <- grep("_", seqnames(refGenome), fixed=TRUE, invert=TRUE, value=TRUE)
  refGenome@user_seqnames <- setNames(min_seqs, min_seqs)
  refGenome@seqinfo <- refGenome@seqinfo[min_seqs]
  
  fragsPerChromosome <-  ceiling(n / length(refGenome))
  g <- unlist(GRangesList(lapply(1:length(refGenome), function(x){
    starts <- sample(1:(length(refGenome[[x]]) - fragWidth), fragsPerChromosome, replace = TRUE)
    GRanges(seqnames =  names(refGenome)[x], strand = '*', ranges = IRanges(start=starts, end = (starts + fragWidth - 1)))
  })))
  
  p <- paste0(seqnames(g), ':', start(g), '-', end(g))
  s <- as.character(BSgenome::getSeq(refGenome, g))
  i <- grepl('N', s)
  p <- p[!i]
  s <- s[!i]
  names(s) <- p
  s
}

hg38.randomFragments.1000000.100 <- sample(createRandomFragments(refGenome = BSgenome.Hsapiens.UCSC.hg38, n = 1500000, fragWidth = 100), 1000000)
mm9.randomFragments.1000000.100  <- sample(createRandomFragments(refGenome = BSgenome.Mmusculus.UCSC.mm9, n = 1500000, fragWidth = 100), 1000000)

save(hg38.randomFragments.1000000.100, file='../data/hg38.randomFragments.1000000.100.RData', compress = TRUE, compression_level = 9)
save(mm9.randomFragments.1000000.100, file='../data/mm9randomFragments.1000000.100.RData', compress = TRUE, compression_level = 9)



# Manualy edit gbff files and create duplicate records for genes with joined ranges.

d <- paste0(readLines('annotations/geneReferences/Oct2019/GCF_001458475.1_KKKWG1_genomic.gbff'), collapse = '\n')
g <- unlist(str_split(d, '  gene  '))
d <- GenomicRanges::makeGRangesFromDataFrame(bind_rows(lapply(g[2:length(g)], function(x){
       range <- unlist(str_match_all(x, '\\d+'))
       strand <- ifelse(grepl('complement\\(<?\\d+\\.\\.\\d+\\)', x), '-', '+')
       
       gene <- str_match(x, 'gene=\"([^\\"]+)"')[,2]
       product <-  str_match(x, 'product=\"([^\\"]+)"')[,2]
       protein_id <- str_match(x, 'protein_id=\"([^\\"]+)"')[,2]
      # geneID <- str_match(x, 'db_xref=\"([^\\"]+)"')[,2]
       locusTag <- str_match(x, 'locus_tag=\"([^\\"]+)"')[,2]
  
       name <- gene
       if(is.na(name)) name <- product
       if(is.na(name)) name <- protein_id
       if(is.na(name)) name <- geneID
       if(is.na(name)) name <- locusTag
       if(is.na(name)) name <- 'Unknown'
       
       if(as.integer(range[1]) > as.integer(range[2])) return(tibble())
       tibble(seqnames = 'chr1', start = as.integer(range[1]), end = as.integer(range[2]), strand = strand, name = name, name2 = locusTag)
    })), keep.extra.columns = TRUE)
  
    
pros <- data.frame(d) %>%
        mutate(regStart = ifelse(strand == '+', start-51, end+1),
               regEnd   = ifelse(strand == '+', start-1, end+51)) %>%
        mutate(start = regStart, end = regEnd, name = paste(name, 'PRO'), name2 = paste(name2, 'PRO')) %>%
        select(-width, -regStart, -regEnd) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Remove promotors which overlap with coding regions.
o <- findOverlaps(pros, d)
pros <- pros[-unique(o@from)]

d <- c(d, pros)

Kingella_kingae_KKKWG1.refSeqGenesGRanges <- d
Kingella_kingae_KKKWG1.refSeqGenesGRanges.exons <- d

save(Kingella_kingae_KKKWG1.refSeqGenesGRanges, file = paste0('../data/Kingella_kingae_KKKWG1.refSeqGenesGRanges.RData'), compress = TRUE, compression_level = 9)
save(Kingella_kingae_KKKWG1.refSeqGenesGRanges.exons, file = paste0('../data/Kingella_kingae_KKKWG1.refSeqGenesGRanges.exons.RData'), compress = TRUE, compression_level = 9)


