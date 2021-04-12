sequencingArchiveDir <- '/media/sequencing/Illumina/'
intSiteDB.group      <- 'intsites_miseq'
sampleDB.group       <- 'specimen_management'
reportOutputDir      <- '/media/lorax/data/export/projects'
Rscript              <- '/opt/R-3.4.4-20180823/lib64/R/bin/Rscript'
softwareDir          <- '/media/lorax/data/software/geneTherapyReports'

seqRunAnalysesToKeep <- 25
nBarCodes <- 30
LTRwidth  <- 14
readSamplePlotN <- 10000
I1.readSamplePlotWidth <- 12
R1.readSamplePlotWidth <- 70
R2.readSamplePlotWidth <- 70


if(file.exists(file.path(reportOutputDir, 'LOCK2'))) q()
system(paste('touch', file.path(reportOutputDir, 'LOCK2')))

options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(RMySQL)
library(stringr)
library(knitr)
library(dplyr)
library(parallel)
library(lubridate)
library(readr)
library(ggplot2)
library(ShortRead)
library(scales)
library(rmarkdown)
library(gt23)
source(file.path(softwareDir, 'lib.R'))

# Create trial directories.
if(! dir.exists(file.path(reportOutputDir, 'trials'))) dir.create(file.path(reportOutputDir, 'trials'))
if(! dir.exists(file.path(reportOutputDir, 'trials', 'seqRuns'))) dir.create(file.path(reportOutputDir, 'trials', 'seqRuns'))


# Select sequencing runs to analyze.
d <- tibble(dir = list.dirs(sequencingArchiveDir, recursive = FALSE, full.names = TRUE),
            run = sapply(dir, le, '/'),
            date = unlist(lapply(dir, function(x){
              o <- unlist(strsplit(x, '/'))
              unlist(strsplit(o[length(o)], '_'))[1]
            }))) %>% 
  filter(grepl('^\\d\\d\\d\\d\\d\\d_M', run)) %>% 
  top_n(seqRunAnalysesToKeep, wt = date) %>%
  arrange(desc(date))

dbConn  <- dbConnect(MySQL(), group = sampleDB.group)
runsToInclude <- unname(unlist(dbGetQuery(dbConn, "select miseqid from seqRunAnalysis where updateReport = '1'")))

invisible(lapply(1:nrow(d), function(x){
  x <- d[x,]
  
  runID <- le(x$dir, '/')
  if(file.exists(file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'trial.rds')) & ! runID %in% runsToInclude) return('report exists')
  message('\n\n', runID, '\n\n')
  
  sampleSheet <- tibble()
  LTRsampleTable <- tibble()
  controlSites <- tibble()
  runStats <- tibble()
  R1plot <- ggplot()
  R2plot <- ggplot()
  
  if(file.exists(file.path(x$dir, 'SampleSheet.csv'))){
    o <- readLines(file.path(x$dir, 'SampleSheet.csv'))
    
    if(any(grepl('\\[metaData\\]', o))){
      o <- o[(grep('\\[metaData\\]', o)+1):length(o)]
      sampleSheet <- read.table(textConnection(o), sep = ',', header = TRUE)
      
      if('primer' %in% names(sampleSheet) & 'ltrBit' %in% names(sampleSheet))
      {
        LTRsampleTable <- data.frame(table(paste0(sampleSheet$primer, sampleSheet$ltrBit)))
        names(LTRsampleTable) <- c('LTR', 'samples')
        LTRsampleTable$LTR <- as.character(LTRsampleTable$LTR )
        LTRwidth <- max(nchar(LTRsampleTable$LTR))
      }
    }
  }
  
  
  LTRtbl <- tibble()
  linkerCodeTbl <- tibble()
  barCodeTbl <- tibble()
  
  
  R1.files <- list.files(file.path(x$dir, 'Data/Intensities/BaseCalls/'), pattern = '_R1', full.names = TRUE)
  if(length(R1.files) > 0){
    R1 <- Reduce('append', lapply(R1.files, function(file){
      strm <- FastqStreamer(file)
      o <- DNAStringSet()
      repeat {
        fasta <- yield(strm)@sread
        if (length(fasta) == 0) break
        o <- Reduce('append', list(o, fasta[width(fasta) >= R1.readSamplePlotWidth]))
      }
      return(o)
    }))
    
    if(length(R1) > 0){
      linkerCodeTbl <- data.frame(sort(table(paste(as.character(subseq(R1, 1, 20)))), decreasing = TRUE))[1:nBarCodes,]
      names(linkerCodeTbl) <- c('R1 1-20', 'R1 counts')
      linkerCodeTbl[,1] <- as.character(linkerCodeTbl[,1])
      R1plot <- readSamplePlot(subseq(R1, 1, R1.readSamplePlotWidth), readSamplePlotN)
      rm(R1)
    } 
  }
  
  R2.files <- list.files(file.path(x$dir, 'Data/Intensities/BaseCalls/'), pattern = '_R2', full.names = TRUE)
  if(length(R2.files) > 0){
    R2 <- Reduce('append', lapply(R2.files, function(file){
      strm <- FastqStreamer(file)
      o <- DNAStringSet()
      repeat {
        fasta <- yield(strm)@sread
        if (length(fasta) == 0) break
        o <- Reduce('append', list(o, fasta[width(fasta) >= R2.readSamplePlotWidth]))
      }
      return(o)
    }))
    
    if(length(R2) > 0){
      LTRtbl <- data.frame(sort(table(paste(as.character(subseq(R2, 1, LTRwidth)))), decreasing = TRUE))[1:nBarCodes,]
      names(LTRtbl) <- c(paste0('R2 1-', LTRwidth), 'R2 counts')
      LTRtbl[,1] <- as.character(LTRtbl[,1])
      R2plot <- readSamplePlot(subseq(R2, 1, R2.readSamplePlotWidth), readSamplePlotN)
      rm(R2)
    } 
  }
  
  I1.files <- list.files(file.path(x$dir, 'Data/Intensities/BaseCalls/'), pattern = '_I1', full.names = TRUE)
  if(length(I1.files) > 0){
    I1 <- Reduce('append', lapply(I1.files, function(file){
      strm <- FastqStreamer(file)
      o <- DNAStringSet()
      repeat {
        fasta <- yield(strm)@sread
        if (length(fasta) == 0) break
        o <- Reduce('append', list(o, fasta))
      }
      return(o)
    }))
    
    if(length(I1) > 0){
      barCodeTbl <- data.frame(sort(table(paste(as.character(I1))), decreasing = TRUE))[1:nBarCodes,]
      names(barCodeTbl) <- c('I1 barcodes', 'I1 counts')
      barCodeTbl[,1] <- as.character(barCodeTbl[,1])
      rm(I1)
    } 
  }
  
  counts <- bind_cols(barCodeTbl, linkerCodeTbl, LTRtbl)
  
  if(nrow(LTRsampleTable) != 0){
    invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
    dbConn  <- dbConnect(MySQL(), group=intSiteDB.group)
    runStats <- dbGetQuery(dbConn, paste0('select * from intSiteCallerStats where miseqid ="', runID,'"'))
    
    if(nrow(runStats) != 0){
      runStats <- bind_rows(lapply(split(runStats, runStats$sampleName), function(x){
        o <- select(x, sampleName, barcoded, LTRed, linkered, ltredlinkered)
        a <- tidyr::spread(select(x, numUniqueSites, refGenome), refGenome, numUniqueSites)
        names(a) <- paste(names(a), 'sites')
        bind_cols(o[1,], a)
      }))
      
      posControls <- sampleSheet[grepl('Positive', sampleSheet$alias, ignore.case = TRUE),]$alias
      
      if(length(posControls) > 0){  
        controlSites <- getDBgenomicFragments(posControls, 'intsites_miseq') %>%
          stdIntSiteFragments() %>%
          calcSampleAbunds()
        
        if(length(controlSites) > 0){
          controlSites <- data.frame(controlSites) %>% mutate(posid = paste0(seqnames, strand, start)) %>% select(sampleName, posid, reads, estAbund)
        }
      }
    }
  }
  
  if(! dir.exists(file.path(reportOutputDir, 'trials', 'seqRuns', runID))) dir.create(file.path(reportOutputDir, 'trials', 'seqRuns', runID))
  
  file.remove(file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'report.pdf'))
  
  rmarkdown::render(file.path(softwareDir, 'seqRunAnalysis.Rmd'),
                    output_file = file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'report.pdf'),
                    params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                  'title' = runID))
  
  if(file.exists(file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'report.pdf'))){
    invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
    dbConn  <- dbConnect(MySQL(), group=sampleDB.group)
    DBI::dbSendQuery(dbConn,paste0('insert into seqRunAnalysis (miseqid, updateReport) values ("', runID, '", "0") ON DUPLICATE KEY UPDATE updateReport="0"'))
  }         
             
  gc()
  
  saveRDS(list(trial = le(x$dir, '/'), lastSampleSeqDate = as.character(ymd(x$date))), 
          file = file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'trial.rds'))
}))

file.remove(file.path(reportOutputDir, 'LOCK2'))