sequencingArchiveDir <- '/media/sequencing/Illumina/'
sampleDB.group       <- 'specimen_management'
intSiteDB.group      <- 'intsites_miseq'
reportOutputDir      <- '/media/lorax/data/export/projects'
Rscript              <- '/opt/R-3.4.4-20180823/lib64/R/bin/Rscript'
softwareDir          <- '/media/lorax/data/software/geneTherapyReports'
CPUs                 <- 10
legacyData           <- '/media/lorax/data/software/geneTherapyReports/geneTherapyData/legacyData.rds'

seqRunAnalysesToKeep <- 100
nBarCodes <- 30
LTRwidth  <- 14
readSamplePlotN <- 10000
I1.readSamplePlotWidth <- 12
R1.readSamplePlotWidth <- 100
R2.readSamplePlotWidth <- 100


if(file.exists(file.path(reportOutputDir, 'LOCK'))) q()
system(paste('touch', file.path(reportOutputDir, 'LOCK')))

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
source('lib.R')

log <- file.path(reportOutputDir, 'reportCreation.log')
write(date(), file = log, append = FALSE)




# Create trial directories.
if(! dir.exists(file.path(reportOutputDir, 'trials'))) dir.create(file.path(reportOutputDir, 'trials'))
if(! dir.exists(file.path(reportOutputDir, 'trials', 'seqRuns'))) dir.create(file.path(reportOutputDir, 'trials', 'seqRuns'))



# Select sequencing runs to analyze.
d <- tibble(dir = list.dirs(sequencingArchiveDir, recursive = FALSE, full.names = TRUE),
                            date = unlist(lapply(dir, function(x){
                              o <- unlist(strsplit(x, '/'))
                              unlist(strsplit(o[length(o)], '_'))[1]
                            }))) %>% 
     filter(grepl('^\\d\\d\\d\\d', date)) %>% 
     top_n(seqRunAnalysesToKeep, wt = date) %>%
     arrange(desc(date))


lapply(1:nrow(d), function(x){
  x <- d[x,]
  
  runID <- le(x$dir, '/')
  sampleSheet <- tibble()
  LTRsampleTable <- tibble()
  controlSites <- tibble()
  runStats <- tibble()
  
  if(file.exists(file.path(x$dir, 'SampleSheet.csv'))){
    o <- readLines(file.path(x$dir, 'SampleSheet.csv'))
    
    if(any(grepl('\\[metaData\\]', o))){
      o <- o[(grep('\\[metaData\\]', o)+1):length(o)]
      sampleSheet <- read.table(textConnection(o), sep = ',', header = TRUE)
      
      if('primer' %in% names(sampleSheet) & 'ltrBit' %in% names(sampleSheet))
      {
        LTRsampleTable <- data.frame(table(paste0(sampleSheet$primer, sampleSheet$ltrBit)))
        names(LTRsampleTable) <- c('LTR', 'samples')
        LTRwidth <- max(nchar(LTRsampleTable$LTR))
        
      }
    }
  }
  
  o <- tryCatch({
         I1 <- readFastq(file.path(x$dir, 'Data/Intensities/BaseCalls/Undetermined_S0_L001_I1_001.fastq.gz'))
         R1 <- readFastq(file.path(x$dir, 'Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz'))
         R2 <- readFastq(file.path(x$dir, 'Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz'))
       }, error = function(e) {
          return(NA)
       })
  
  barCodeTbl <- data.frame(sort(table(paste(as.character(I1@sread))), decreasing = TRUE))[1:nBarCodes,]
  names(barCodeTbl) <- c('I1 barcodes', 'counts')
  
  linkerCodeTbl <- data.frame(sort(table(paste(as.character(subseq(R1@sread, 1, 20)))), decreasing = TRUE))[1:nBarCodes,]
  names(linkerCodeTbl) <- c('R1 1-20', 'counts')
  
  LTRtbl <- data.frame(sort(table(paste(as.character(subseq(R2@sread, 1, LTRwidth)))), decreasing = TRUE))[1:nBarCodes,]
  names(LTRtbl) <- c(paste0('R2 1-', LTRwidth), 'counts')
  
  counts <- bind_cols(barCodeTbl, linkerCodeTbl, LTRtbl)
  names(counts) <-  c('I1 barcodes', 'counts', 'R1 1-20', 'counts', paste0('R2 1-', LTRwidth), 'counts')
  
  counts[,1] <- as.character(counts[,1])
  counts[,3] <- as.character(counts[,3])
  counts[,5] <- as.character(counts[,5])

  I1.readSamplePlotWidth <- as.integer(names(sort(table(width(I1)), decreasing = TRUE)[1]))
  
  o <- tryCatch({
    ### I1plot <- readSamplePlot(subseq(I1@sread, 1, I1.readSamplePlotWidth), readSamplePlotN)
    R1plot <- readSamplePlot(subseq(R1@sread, 1, R1.readSamplePlotWidth), readSamplePlotN)
    R2plot <- readSamplePlot(subseq(R2@sread, 1, R2.readSamplePlotWidth), readSamplePlotN)
  }, 
  error = function(e) {
    return(NA)
  })
  
  rm(R1, R2, I1)
  
  if(nrow(LTRsampleTable) != 0){
    invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
    dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
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
  
  browser()
  
  rmarkdown::render('seqRunAnalysis.Rmd',
                    output_file = file.path(reportOutputDir, 'trials', 'seqRuns', runID, 'report.pdf'),
                    params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                  'title' = runID))
  })  
  
  
  
  
  
  






legacyData <- readRDS(legacyData)

# Create a list of all GTSPs that passed through the INSPIIRED pipeline.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples    <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
intSitesamples    <- unique(gsub('\\-\\d+$', '', intSitesamples))
intSitesamplesTbl <- dbGetQuery(dbConn, 'select * from samples where sampleName like "%GTSP%"')

# Read in sample data.
dbConn  <- dbConnect(MySQL(), group = sampleDB.group)
samples <- dbGetQuery(dbConn, paste0('select * from gtsp')) %>% filter(SpecimenAccNum %in% intSitesamples)


# Create a table of samples, seq dates, trials and patients.
intSitesamplesTbl$date <- ymd(unlist(lapply(strsplit(intSitesamplesTbl$miseqid, '_'), '[', 1)))
intSitesamplesTbl$SpecimenAccNum <- sub('\\-\\d+$', '', intSitesamplesTbl$sampleName)
intSitesamplesTbl <- left_join(intSitesamplesTbl, select(samples, SpecimenAccNum, Patient, Trial), by = 'SpecimenAccNum')




invisible(sapply(unique(samples$Trial), function(x) dir.create(file.path(reportOutputDir, 'trials', x))))




# Idenitfy Trials and patient reports which are no longer defined in the databse and delete them in the output directory.
invisible(lapply(list.files(reportOutputDir, recursive = TRUE, pattern = '.pdf$'), function(x){
  o <- unlist(strsplit(x, '/'))
  
  file <- paste0(reportOutputDir, '/', paste0(o[1:(length(o)-1)], collapse = '/'), '/', 'sampleData.csv')
  ### message(file)
  
  if(file.exists(file)){
    sd <- readr::read_csv(file, col_types = cols(), trim_ws = FALSE) %>% 
          filter(! SpecimenAccNum %in% legacyData$GTSP) %>%  # Remove legacy data from report table since it would not be in db.
          arrange(SpecimenAccNum)
    sd2 <- subset(samples, SpecimenAccNum %in% sd$SpecimenAccNum) %>% select(SpecimenAccNum, Trial, CellType, Timepoint) %>% arrange(SpecimenAccNum) 

    delete <- FALSE
    if (nrow(sd) != nrow(sd2)){
       delete <- TRUE
    } else {
      if(! all(sd == sd2)){
        delete <- TRUE
      }
    }
    
    if(delete){
      #browser()
      message('\nDeleting patient directiory due to db mismatch: ', paste0(o[1:(length(o)-1)], collapse = '/'))
      unlink(paste0(reportOutputDir, '/', paste0(o[1:(length(o)-1)], collapse = '/')), recursive = TRUE)
    }
  }
  
  # Trial / patient does not exist in databae -- delete patient directory
  if(nrow(subset(samples, Trial == o[length(o)-2] & Patient == o[length(o)-1])) == 0){ 
    message('Deleting patient directiory ', paste0(o[1:(length(o)-1)], collapse = '/'))
    unlink(paste0(reportOutputDir, '/', paste0(o[1:(length(o)-1)], collapse = '/')), recursive = TRUE)
  }
  
  # Trial does not exist in databae -- delete trial directory.
  if(nrow(subset(samples, Trial == o[length(o)-2])) == 0){ 
    message('Deleting trial directory ', paste0(o[1:(length(o)-2)], collapse = '/'))
    unlink(paste0(reportOutputDir, '/', paste0(o[1:(length(o)-2)], collapse = '/')), recursive = TRUE)
  }
}))




# Create a table of reports to create / update.
r <- bind_rows(lapply(split(intSitesamplesTbl, paste(intSitesamplesTbl$Trial, intSitesamplesTbl$Patient)), function(x){
       file <-  paste0(file.path(reportOutputDir, 'trials', x$Trial[1], x$Patient[1], x$Patient[1]), '.pdf')
       tibble(trial = x$Trial[1], patient = x$Patient[1], 
              mtime =  file.info(file)$mtime,
              update = ifelse(is.na(mtime), TRUE, file.info(file)$mtime < max(x$date)),
              reportPath = file)
    })) %>% filter(update == TRUE)

r <- r[! is.na(r$patient),]
r <- r[! is.na(r$trial),]
o <- tibble()

if(nrow(r) > 0){
  r$n <- ntile(1:nrow(r), CPUs)
  cluster <- makeCluster(CPUs)
  clusterExport(cluster, c('Rscript', 'log', 'sampleDB.group', 'intSiteDB.group', 'reportOutputDir', 'softwareDir'))
  
  o <- bind_rows(parLapply(cluster, split(r, r$n), function(x){
  #o <- bind_rows(lapply(split(r, r$n), function(x){
         library(dplyr)
         bind_rows(lapply(split(x, paste(x$patient, x$trial)), function(x2){
      
           command <- paste0(Rscript,  ' ', softwareDir, '/geneTherapySubjectReport/report.R ',
                             '--specimenDB  ', sampleDB.group, ' --intSiteDB ', intSiteDB.group, ' ',
                             '--patient "', x2$patient, '" --trial "', x2$trial, '" ',
                             '--outputDir ', file.path(reportOutputDir, 'trials', x2$trial), ' ',
                             '--reportFile geneTherapySubjectReport/report.Rmd --legacyData geneTherapyData/legacyData.rds')
           
           write(c(command, '\n'), file = log, append = TRUE)
    
           # Remove previous result.
           # Report maker software creates patient diretory in output folder.
           unlink(file.path(reportOutputDir, 'trials', x2$trial, x2$patient), recursive = TRUE)
           
           o <- tryCatch({
                 system(command)
                 'Success'
                }, error = function(e) {
                  'Failed'
               })
           
           if(file.exists(paste0(file.path(reportOutputDir, 'trials', x2$trial, x2$patient, x2$patient), '.pdf'))){
              result <- 'success'
              system(paste0('touch ', file.path(reportOutputDir, 'trials', x2$trial, 'updateData')))
           } else {
             ### browser()
             result <- 'failed'
           }
           
           tibble(trial = x2$trial, patient = x2$patient, result = result)
    }))
  }))
  
  write.table(o, sep = '\t', file = log, quote = FALSE, row.names = FALSE, append = TRUE)
}

# Create trial level pages.
# fairly quick -- fine to regenerate all reports with each script run.
createTrialReports <- lapply(list.dirs(file.path(reportOutputDir, 'trials'), recursive = FALSE), function(x){
  tryCatch({
    o <-unlist(strsplit(x, '/'))
    rmarkdown::render('trialReport.Rmd',
                      output_file = file.path(x, 'index.html'),
                      params = list('title' = gsub('_', ' ', o[length(o)]),
                                    'trial' = o[length(o)],
                                    'reportPath' = x))
    return('Success')
  }, error = function(e) {
    return('Failed')
  })
})

createTrialReports <- data.frame(result = unlist(createTrialReports))
createTrialReports$trial <- list.dirs(file.path(reportOutputDir, 'trials'), recursive = FALSE)



if(! dir.exists(file.path(reportOutputDir, 'dashboard'))) dir.create(file.path(reportOutputDir, 'dashboard'))

# Create overview page.
rmarkdown::render('dashboard.Rmd', output_file = file.path(reportOutputDir, 'dashboard', 'index.html'),  
                  params = list('reportOutputDir' = reportOutputDir,
                                'sequencingArchiveDir' = sequencingArchiveDir))

file.remove(file.path(reportOutputDir, 'LOCK'))