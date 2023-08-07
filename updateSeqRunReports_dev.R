#!/opt/R-3.4.4-20180823/lib64/R/bin/Rscript

sequencingArchiveDir <- '/media/sequencing/Illumina/'
intSiteDB.group      <- 'intsites_miseq'
sampleDB.group       <- 'specimen_management'
reportOutputDir      <- '/media/lorax/data/export/projects'
Rscript              <- '/opt/R-3.4.4-20180823/lib64/R/bin/Rscript'
softwareDir          <- '/media/lorax/data/software/geneTherapyReports'
logFile              <- file.path(reportOutputDir, 'log')


if(file.exists(file.path(reportOutputDir, 'PMACS.lock'))) stop('PMACS.lock file found -- aborting.')

seqRunAnalysesToKeep <- 75
nBarCodes <- 30
LTRwidth  <- 14
readSamplePlotN <- 10000
I1.readSamplePlotWidth <- 12
R1.readSamplePlotWidth <- 70
R2.readSamplePlotWidth <- 70


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

write(paste0(date(), ' -- Starting seq run report creator'), file = logFile, append = TRUE)

# Create trial directories.
if(! dir.exists(file.path(reportOutputDir, 'trials'))) dir.create(file.path(reportOutputDir, 'trials'))
if(! dir.exists(file.path(reportOutputDir, 'trials', 'seqRuns'))) dir.create(file.path(reportOutputDir, 'trials', 'seqRuns'))


# Select sequencing runs to analyze.
# NexSeq runs typically do not have the M after the leading date.
d <- tibble(dir = list.dirs(sequencingArchiveDir, recursive = FALSE, full.names = TRUE),
            run = sapply(dir, le, '/'),
            date = unlist(lapply(dir, function(x){
              o <- unlist(strsplit(x, '/'))
              unlist(strsplit(o[length(o)], '_'))[1]
            }))) %>% 
  #filter(grepl('^\\d\\d\\d\\d\\d\\d_M', run)) %>% 
  filter(grepl('^\\d\\d\\d\\d\\d\\d', run)) %>% 
  top_n(seqRunAnalysesToKeep, wt = date) %>%
  arrange(desc(date))

dbConn  <- dbConnect(MySQL(), group = sampleDB.group)
runsToInclude <- unname(unlist(dbGetQuery(dbConn, "select miseqid from seqRunAnalysis where updateReport = '1'")))

cluster <- makeCluster(15)
clusterExport(cluster, c('R1.readSamplePlotWidth', 'R2.readSamplePlotWidth', 'I1.readSamplePlotWidth', 'logFile'))

write(paste0('Processing ', nrow(d), ' sequencing runs.'), file = logFile, append = TRUE)
write.table(d, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = logFile, append = TRUE)

invisible(lapply(1:nrow(d), function(x){
  x <- d[x,]
  
  browser()
  
  runID <- le(x$dir, '/')
  write(paste0('Starting ', runID), file = logFile, append = TRUE)
  
  
  # Skip this run if the final output is present and it is not in the runsToIncude vector from the seqRunAnalysis 
  # database table which would overide.
  
  if(file.exists(file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'trial.rds')) & ! runID %in% runsToInclude){
    write('Skipping this run because output exists and not in override.', file = logFile, append = TRUE)  
    return('report exists')
  }

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
  R2.files <- list.files(file.path(x$dir, 'Data/Intensities/BaseCalls/'), pattern = '_R2', full.names = TRUE)
  
  
  R1.size <- as.numeric(stringr::str_extract(gdata::humanReadable(sum(unlist(lapply(R1.files, file.size))), standard = 'Unix', units = 'G'), '[\\d\\.]+'))
  R2.size <- as.numeric(stringr::str_extract(gdata::humanReadable(sum(unlist(lapply(R2.files, file.size))), standard = 'Unix', units = 'G'), '[\\d\\.]+'))

  if(R1.size > 5 | R2.size > 5){
  write('Skipping this run because the data size is too large.', file = logFile, append = TRUE)  
    return(NA)
  }
  
  
  if(length(R1.files) > 0){
    
    write('Streaming in R1 files...', file = logFile, append = TRUE)
    
    R1 <- unlist(parLapply(cluster, R1.files, function(file){
            library(ShortRead)
            r <- readFastq(file)
            r <- r[width(r) >= R1.readSamplePlotWidth]
            as.character(narrow(r, 1, R1.readSamplePlotWidth)@sread)
          }))
    
    if(length(R1) > 0){
      linkerCodeTbl <- data.frame(sort(table(paste(substr(R1, 1, 20))), decreasing = TRUE))[1:nBarCodes,]
      names(linkerCodeTbl) <- c('R1 1-20', 'R1 counts')
      linkerCodeTbl[,1] <- as.character(linkerCodeTbl[,1])
      R1plot <- readSamplePlot(R1, readSamplePlotN)
      rm(R1)
      gc()
    } 
  }
  
  
  if(length(R2.files) > 0){
    write('Streaming in R2 files...', file = logFile, append = TRUE)
    
    R2 <- unlist(parLapply(cluster, R2.files, function(file){
      library(ShortRead)
      r <- readFastq(file)
      r <- r[width(r) >= R2.readSamplePlotWidth]
      as.character(narrow(r, 1, R2.readSamplePlotWidth)@sread)
    }))
    
    if(length(R2) > 0){
      #LTRtbl <- data.frame(sort(table(paste(as.character(subseq(R2, 1, LTRwidth)))), decreasing = TRUE))[1:nBarCodes,]
      LTRtbl <- data.frame(sort(table(paste(substr(R2, 1, LTRwidth))), decreasing = TRUE))[1:nBarCodes,]
      names(LTRtbl) <- c(paste0('R2 1-', LTRwidth), 'R2 counts')
      LTRtbl[,1] <- as.character(LTRtbl[,1])
      R2plot <- readSamplePlot(R2, readSamplePlotN)
      rm(R2)
      gc()
    } 
  }
  
  
  I1.files <- list.files(file.path(x$dir, 'Data/Intensities/BaseCalls/'), pattern = '_I1', full.names = TRUE)
  if(length(I1.files) > 0){
    write('Streaming in I1 files...', file = logFile, append = TRUE)
    
    I1 <- unlist(parLapply(cluster, I1.files, function(file){
      library(ShortRead)
      r <- readFastq(file)
      as.character(r@sread)
      ### r <- r[width(r) >= I1.readSamplePlotWidth]
      ### as.character(narrow(r, 1, I1.readSamplePlotWidth)@sread)
    }))
    
    
    if(length(I1) > 0){
      write('Processing I1 data.', file = logFile, append = TRUE)
      
      #barCodeTbl <- data.frame(sort(table(paste(as.character(I1))), decreasing = TRUE))[1:nBarCodes,]
      barCodeTbl <- data.frame(sort(table(I1), decreasing = TRUE))[1:nBarCodes,]
      
      names(barCodeTbl) <- c('I1 barcodes', 'I1 counts')
      barCodeTbl[,1] <- as.character(barCodeTbl[,1])
      rm(I1)
      gc()
    } 
  }
  
  counts <- bind_cols(barCodeTbl, linkerCodeTbl, LTRtbl)
  if(nrow(counts) == 0) return()
  
  if(nrow(LTRsampleTable) != 0){
    write('Processing LTR sample table.', file = logFile, append = TRUE)
    invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
    dbConn  <- dbConnect(MySQL(), group=intSiteDB.group)
    runStats <- dbGetQuery(dbConn, paste0('select * from intSiteCallerStats where miseqid ="', runID,'"'))
    
    if(nrow(runStats) != 0){
      write('Processing run stats.', file = logFile, append = TRUE)
      runStats <- bind_rows(lapply(split(runStats, runStats$sampleName), function(x){
        o <- select(x, sampleName, barcoded, LTRed, linkered, ltredlinkered)
        a <- tidyr::spread(select(x, numUniqueSites, refGenome), refGenome, numUniqueSites)
        names(a) <- paste(names(a), 'sites')
        bind_cols(o[1,], a)
      }))
      
      posControls <- sampleSheet[grepl('Positive', sampleSheet$alias, ignore.case = TRUE),]$alias
      
      if(length(posControls) > 0){
        write('Processing positive controls.', file = logFile, append = TRUE)
        controlSites <- getDBgenomicFragments(posControls, 'intsites_miseq') 
        if(length(controlSites) > 0) controlSites <- stdIntSiteFragments(controlSites)
        if(length(controlSites) > 0) controlSites <- calcSampleAbunds(controlSites)
        
        if(length(controlSites) > 0){
          controlSites <- data.frame(controlSites) %>% mutate(posid = paste0(seqnames, strand, start)) %>% select(sampleName, posid, reads, estAbund)
        } else {
          controlSites <- data.frame()
        }
      }
    }
  }
  
  if(! dir.exists(file.path(reportOutputDir, 'trials', 'seqRuns', runID))) dir.create(file.path(reportOutputDir, 'trials', 'seqRuns', runID))
  
  file.remove(file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'report.pdf'))
  
  write('Knitting report.', file = logFile, append = TRUE)
  
  rmarkdown::render(file.path(softwareDir, 'seqRunAnalysis.Rmd'),
                    output_file = file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'report.pdf'),
                    params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                  'title' = runID))
  
  if(file.exists(file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'report.pdf'))){
    write('Updating database.', file = logFile, append = TRUE)
    invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
    dbConn  <- dbConnect(MySQL(), group=sampleDB.group)
    q <- paste0('insert into seqRunAnalysis (miseqid, updateReport) values ("', runID, '", "0") ON DUPLICATE KEY UPDATE updateReport="0"')
    message(q)
    DBI::dbSendQuery(dbConn, q)
  }         
             
  gc()
  
  write(paste0('Completed ', runID, '\n\n'), file = logFile, append = TRUE)
  
  saveRDS(list(trial = le(x$dir, '/'), lastSampleSeqDate = as.character(ymd(x$date))), 
          file = file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'trial.rds'))
}))

write(paste0(date(), ' -- updateSeqRunReports done.'), file = logFile, append = TRUE)
