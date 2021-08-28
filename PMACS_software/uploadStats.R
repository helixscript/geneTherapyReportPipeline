library(RMySQL)
options(useFancyQuotes = FALSE, stringsAsFactors = FALSE)

#outputDir <- '210302_MN01490_0008_A000H3FHWW_hg38'

args <- commandArgs(trailingOnly=TRUE)
outputDir <- args[1]


if(! file.exists(file.path(outputDir, 'stats.txt'))) 
  stop(paste0('Error -- ', file.path(outputDir, 'stats.txt'), ' not found.'))

o <- outputDir
o <- unlist(strsplit(o, '/'))
finalDir <- o[length(o)]

if(! grepl('_', finalDir)) stop('Error -- output dir does not contain an underscore.')

logFile <- paste0('output/', finalDir, '.uploadStats.log')
write(date(), file = logFile, append = FALSE)

o <- unlist(strsplit(finalDir, '_'))
refGenome <- o[length(o)]

miseqid <- paste0(o[1:(length(o)-1)], collapse = '_')

o <- read.table(file.path(outputDir, 'stats.txt'), sep = '\t', header = TRUE)

dbConn  <- dbConnect(MySQL(), group='intsites_miseq')

for( x in (1:nrow(o))) {
  x <- o[x,]
  x[is.na(x)] <- 'NULL'
  write.table(data.frame(x), file = 't', sep = ',')
  comm <- paste0(c("insert into intSiteCallerStats (sampleName, barcoded, p_qTrimmed, l_qTrimmed, LTRed, linkered, ltredlinkered, lenTrim, vTrimed, uniqL, uniqP, uniqP30,  lLen, pLen, readsAligning, readsWithGoodAlgnmts, numProperlyPairedAlignments, chimeras, numAllSingleReads, numAllSingleSonicLengths, numUniqueSites, multihitReads, multihitSonicLengths, multihitClusters, totalSonicLengths, totalEvents, linkerSequence, bcSeq, gender, primer, ltrBit, largeLTRFrag, vectorSeq, refGenome, miseqid) values (", 
                   paste(sQuote(c(x[1:(length(x)-2)], miseqid)), collapse = ', '), ')'), collapse = '')

  comm <- gsub("'NULL'", "NULL", comm)

  message(comm, '\n\n')
  
  out <- tryCatch(
    {
       a <- DBI::dbSendQuery(dbConn, comm)
       return(a)
    },
    error=function(cond) {
       return(cond)
    }
   )
}

for(x in (1:nrow(o))) {
  x <- o[x,]
  r <- DBI::dbGetQuery(dbConn, paste0('select * from intSiteCallerStats where sampleName = "', x[1], '"'))
  if(nrow(r) == 0) write(paste0('Error - could not insert results for sample ', x[1]), file = logFile, append = TRUE)
}
