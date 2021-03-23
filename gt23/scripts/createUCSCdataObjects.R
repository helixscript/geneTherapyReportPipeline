library(GenomicRanges)
library(tidyverse)
options(stringsAsFactors = FALSE)

humanGeneFilter <- function(d){
  humanGenes <- gsub('\\.\\d+$', '', gt23::hg38.refSeqGenesGRanges$name)
  d$name <- gsub('\\.\\d+$', '', d$name)
  subset(d, toupper(name) %in% toupper(humanGenes))
}

createRefSeqObjects <- function(genomeLabel, file, humanGeneFilter = FALSE){
  d <- read.table(file, sep = '\t', header = FALSE, quote = '')

  names(d) <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart',
                'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2',
                'cdsStartStat', 'cdsEndStat', 'exonFrames')

  if(humanGeneFilter) d <- humanGeneFilter(d)

  g <- makeGRangesFromDataFrame(d, seqnames.field = 'chrom', start.field = 'txStart',
                                   end.field = 'txEnd', strand.field = 'strand',
                                   keep.extra.columns = TRUE)

  e <- group_by(d, 1:nrow(d)) %>%
       unnest(exonStarts = strsplit(exonStarts, ","), exonEnds = strsplit(exonEnds, ",")) %>%
       mutate(name = paste('exon', 1:n())) %>%
       ungroup() %>%
       makeGRangesFromDataFrame(seqnames.field = 'chrom', start.field = 'exonStarts',
                                end.field = 'exonEnds', strand.field = 'strand',
                                keep.extra.columns = TRUE)

  assign(paste0(genomeLabel, '.refSeqGenesGRanges'), g)
  assign(paste0(genomeLabel, '.refSeqGenesGRanges.exons'), e)

  save(list = paste0(genomeLabel, '.refSeqGenesGRanges'), file = paste0('../data/', genomeLabel, '.refSeqGenesGRanges.RData'), compress = TRUE, compression_level = 9)
  save(list = paste0(genomeLabel, '.refSeqGenesGRanges.exons'), file = paste0('../data/', genomeLabel, '.refSeqGenesGRanges.exons.RData'), compress = TRUE, compression_level = 9)
}

# Create special empty vector object for instances when no refGenomes are to be used.
none <- vector(mode = "character", length = 0)
save(none,  file = paste0('../data/none.RData'))

createRefSeqObjects('sacCer3', 'annotations/geneReferences/March_2021/sacCer3.refSeq.curated.txt.gz')
createRefSeqObjects('hg38',    'annotations/geneReferences/June2019/hg38.refSeq.curated.txt.gz')
createRefSeqObjects('hg18',    'annotations/geneReferences/June2019/hg18.refGene.txt.gz')
createRefSeqObjects('mm9',     'annotations/geneReferences/June2019/mm9.refGene.txt.gz')
createRefSeqObjects('susScr3', 'annotations/geneReferences/June2019/susScr3.refGene.txt.gz')
createRefSeqObjects('macFas5', 'annotations/geneReferences/June2019/macFas5.refGene.txt.gz')
createRefSeqObjects('canFam3', 'annotations/geneReferences/June2019/canFam3.refGene.txt.gz')
createRefSeqObjects('canFam3.humanXeno', 'annotations/geneReferences/June2019/canFam3.xenoRefGene.txt.gz', humanGeneFilter = TRUE)
createRefSeqObjects('macFas5.humanXeno', 'annotations/geneReferences/June2019/macFas5.xenoRefGene.txt.gz', humanGeneFilter = TRUE)

# Create human oreganno gene regulation GRanges object.
d <- read.table('annotations/regulation/June2019/hg38.oreganno.txt.gz', header = FALSE)
hg38.oreganno <- GenomicRanges::makeGRangesFromDataFrame(d, seqnames.field = 'V2', start.field = 'V3',
                                                         end.field = 'V4', strand.field = 'V6', 
                                                         starts.in.df.are.0based = TRUE)
hg38.oreganno$id <- d$V7
save(hg38.oreganno, file='../data/hg38.oreganno.RData', compress = TRUE, compression_level = 9)


# Create human CCDS GRanges object.
d <- read.table('annotations/ccds/June2019/hg38.ccdsGene.txt.gz', header = FALSE); 
hg38.ccds <- GenomicRanges::makeGRangesFromDataFrame(dplyr::bind_rows(lapply(1:nrow(d), function(x){
  r <- d[x,]
  data.frame('name'     = r$V2,
             'seqnames' = r$V3,
             'strand'   = r$V4,
             'start'    = unlist(strsplit(r$V10, ',')),
             'end'      = unlist(strsplit(r$V11, ',')))
})), keep.extra.columns = TRUE)
save(hg38.ccds, file='../data/hg38.ccds.RData', compress = TRUE, compression_level = 9)


