library(dplyr)
library(ggplot2)
library(ShortRead)
library(scales)
library(rmarkdown)
library(RMySQL)
library(gt23)
options(stringsAsFactors = FALSE)


nBarCodes <- 30
LTRwidth  <- 14
readSamplePlotN <- 10000
I1.readSamplePlotWidth <- 12
R1.readSamplePlotWidth <- 100
R2.readSamplePlotWidth <- 100
intSiteDB.group <- 'intsites_miseq'

runDir <- '/media/sequencing/Illumina/210302_MN01490_0008_A000H3FHWW'
o <- unlist(strsplit(runDir, '/'))
runID <- o[length(o)]


sampleSheet <- tibble()
LTRsampleTable <- tibble()

o <- readLines(paste0(runDir, '/SampleSheet.csv'))
if(any(grepl('\\[metaData\\]', o))){
  o <- o[(grep('\\[metaData\\]', o)+1):length(o)]
  sampleSheet <- read.table(textConnection(o), sep = ',', header = TRUE)
}

LTRsampleTable <- data.frame(table(paste0(sampleSheet$primer, sampleSheet$ltrBit)))
names(LTRsampleTable) <- c('LTR', 'samples')




I1 <- readFastq(file.path(runDir, 'Data/Intensities/BaseCalls/Undetermined_S0_L001_I1_001.fastq.gz'))
R1 <- readFastq(file.path(runDir, 'Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz'))
R2 <- readFastq(file.path(runDir, 'Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz'))

barCodeTbl <- data.frame(sort(table(paste(as.character(I1@sread))), decreasing = TRUE))[1:nBarCodes,]
names(barCodeTbl) <- c('I1 barcodes', 'counts')

linkerCodeTbl <- data.frame(sort(table(paste(as.character(subseq(R1@sread, 1, 20)))), decreasing = TRUE))[1:nBarCodes,]
names(linkerCodeTbl) <- c('R1 1-20', 'counts')

LTRtbl <- data.frame(sort(table(paste(as.character(subseq(R2@sread, 1, LTRwidth)))), decreasing = TRUE))[1:nBarCodes,]
names(LTRtbl) <- c(paste0('R2 1-', LTRwidth), 'counts')

counts <- bind_cols(barCodeTbl, linkerCodeTbl, LTRtbl)
names(counts) <-  c('I1 barcodes', 'counts', 'R1 1-20', 'counts', paste0('R2 1-', LTRwidth), 'counts')

counts$`I1 barcodes` <- as.character(counts$`I1 barcodes`)
counts$`R1 1-20`     <- as.character(counts$`R1 1-20`)
counts$`R2 1-10`     <- as.character(counts$`R2 1-10`)

readSamplePlot <- function(reads, n){
  ds <- as.character(sample(unique(reads), n))
  dp <- lapply(strsplit(sort(ds), ''), function(x){ tibble(base = x, n = 1:length(x)) })
  dp <- bind_rows(mapply(function(x, n){ x$read <- n; x}, dp, 1:length(dp), SIMPLIFY = FALSE))
  dp$base <- factor(dp$base, levels = c('A', 'T', 'C', 'G', 'N'))

  ggplot(dp, aes(n, read, fill = base)) + theme_bw() + geom_tile() +
    scale_fill_manual(values =  c('red', 'green', 'blue', 'gold', 'gray50')) +
    scale_x_continuous(limits = c(0, width(reads[1])), expand = c(0, 0)) +
    scale_y_continuous(label=comma, limits = c(0, n), expand = c(0, 0)) +
    labs(x = 'Position', y = 'Reads') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

I1plot <- readSamplePlot(subseq(I1@sread, 1, I1.readSamplePlotWidth), readSamplePlotN)
R1plot <- readSamplePlot(subseq(R1@sread, 1, R1.readSamplePlotWidth), readSamplePlotN)
R2plot <- readSamplePlot(subseq(R2@sread, 1, R2.readSamplePlotWidth), readSamplePlotN)

rm(R1, R2, I1)


invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
runStats    <- dbGetQuery(dbConn, paste0('select * from intSiteCallerStats where miseqid ="', runID,'"'))





getDBgenomicFragments <- function(samples, intSiteDB.group){
  options(useFancyQuotes = FALSE)
  
  dbConn <- DBI::dbConnect(RMySQL::MySQL(), group=intSiteDB.group)
  
  intSiteSamples <- DBI::dbGetQuery(dbConn, 'select * from samples')
  intSiteSamples$GTSP <- gsub('\\-\\d+$', '', intSiteSamples$sampleName)
  
  sampleTable <-  base::subset(intSiteSamples, sampleName %in% samples)
  
  sampleIDs <- sampleTable$sampleID
  
  if(length(sampleIDs) == 0) return(GenomicRanges::GRanges())
  
  replicateQuery <- paste('samples.sampleID in (', paste0(sampleIDs, collapse = ','), ')')
  
  q <- sprintf("select position, chr, strand, breakpoint, count, refGenome,
               sampleName from sites left join samples on
               sites.sampleID = samples.sampleID
               left join pcrbreakpoints on
               pcrbreakpoints.siteID = sites.siteID
               where (%s)", replicateQuery)
  
  sites <- DBI::dbGetQuery(dbConn, q)
  
  DBI::dbDisconnect(dbConn)
  
  if(nrow(sites) == 0) return(GRanges())
  
  sites$sampleName <- paste0(sites$refGenome, '-', sites$sampleName)
  
  
  sites$GTSP      <- as.character(sub('\\-\\d+$', '', sites$sampleName))
  sites$patient   <- 'PositiveControl'
  sites$timePoint <- 'd0'
  sites$cellType  <- 'control'
  
  sites <- dplyr::bind_cols(sites, expandTimePoints(sites$timePoint))
  
  intSites <-
    GenomicRanges::GRanges(seqnames = S4Vectors::Rle(sites$chr),
                           ranges   = IRanges::IRanges(start = pmin(sites$position, sites$breakpoint),
                                                       end   = pmax(sites$position, sites$breakpoint)),
                           strand   = S4Vectors::Rle(sites$strand))
  
  sites <- dplyr::rename(sites, reads = count)
  GenomicRanges::mcols(intSites) <- data.frame(dplyr::select(sites, refGenome, reads, patient, sampleName,
                                                             GTSP, cellType, timePoint, timePointDays, timePointMonths))
  intSites
}


controlSites <- tibble()

if(nrow(runStats) != 0){
  runStats <- bind_rows(lapply(split(runStats, runStats$sampleName), function(x){
    o <- select(x, sampleName, barcoded, LTRed, linkered, ltredlinkered)
    a <- tidyr::spread(select(x, numUniqueSites, refGenome), refGenome, numUniqueSites)
    names(a) <- paste(names(a), 'sites')
    bind_cols(o[1,], a)
  }))
  
  posControls <- sampleSheet[grepl('Positive', sampleSheet$alias, ignore.case = TRUE),]$alias
    
  controlSites <- getDBgenomicFragments(posControls, 'intsites_miseq') %>%
  stdIntSiteFragments() %>%
  calcSampleAbunds()
  
  if(length(controlSites) > 0){
    controlSites <- data.frame(controlSites) %>% mutate(posid = paste0(seqnames, strand, start)) %>% select(sampleName, posid, reads, estAbund)
  } 
  
}


rmarkdown::render('seqRunAnalysis.Rmd',
                  output_file = 'test.pdf',
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = runID))



