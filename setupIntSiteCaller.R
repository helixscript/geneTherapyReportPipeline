# Define remote hosts and paths.
remoteHost                <- 'everett\\@microb120.med.upenn.edu'
remoteHostDataDir         <- '/media/sequencing/Illumina'
intSiteCaller.hoursToWait <- 5

options(stringsAsFactors = FALSE)
library(RMySQL)
write(date(), file = 'log')
wd <- getwd()

# Retrieve sequencing run id from command line.
args = commandArgs(trailingOnly=TRUE)
### args <- '210302_MN01490_0008_A000H3FHWW'
if(length(args) == 0) stop('Error - no sequencing run id was provided.')
seqRun <- args[1]

system('rm -rf intSiteCaller_vectors')
system('git clone https://github.com/BushmanLab/intSiteCaller_vectors')

# Setup and run intSiteCaller.
source('setupIntSiteCaller.lib.R')

s <- metaData()

dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSiteSamples <- dbGetQuery(dbConn, 'select sampleName, refGenome from samples')
intSiteSamples <- paste0(intSiteSamples$sampleName, '/', intSiteSamples$refGenome)

intSiteRuns <- dbGetQuery(dbConn, 'select miseqid, refGenome from samples')
intSiteRuns <- paste0(intSiteRuns$miseqid, '/', intSiteRuns$refGenome)


for(refGenome in unique(s$refGenome)){
  
  s$id <- paste0(s$alias, '/', refGenome)

  if(paste0(seqRun, '/', refGenome) %in% intSiteRuns){
    write(paste0(seqRun, '/', refGenome, 'run is already in the database.'), file = 'log', append = TRUE)
    stop('Run duplication error -- see log.')
  }
  
  duplicateSamples <- s$id[s$id %in% intSiteSamples]
  if(length(duplicateSamples)){
    write(paste0('Error - the following samples are already in the intSite database: ',
                 paste0(duplicateSamples, collapse=', ')), file = 'log', append = TRUE)
    stop('Sample duplication error -- see log.')
  }
  
  if(! file.exists(paste0('output/', seqRun, '_', refGenome, '/error.txt'))) runIntSiteCaller(paste0('output/', seqRun, '_', refGenome), refGenome, s)
  
  # Wait for intSiteCaller to finish, check every 5 minutes (300 seconds).
  done <- FALSE
  timeElapsed <- 0
  write(paste0('Waiting for intSiteCaller to finish ', seqRun, ' refGenome:', refGenome), file = 'log', append = TRUE)
  
  for (w in rep(300, (intSiteCaller.hoursToWait * 12))){
    Sys.sleep(w)
    timeElapsed <- timeElapsed + w
    write(paste0(timeElapsed/60, ' minutes passed waiting for intSiteCaller to finish. Looking for ',
                 'output/', seqRun, '_', refGenome, '/error.txt'), file = 'log', append = TRUE)
    if(file.exists(paste0('output/', seqRun, '_', refGenome, '/error.txt'))){
      done <- TRUE
      break
    }
  }
  
  # Kill all the intSiteCaller jobs if it has not completed within the aloted time.
  if(done == FALSE){
    write(paste0('intSiteCaller did not finish after ', intSiteCaller.hoursToWait, ' hours. Stopping all processes.'), file = 'log', append = TRUE)
    sapply(system(paste0('bjobs -w | grep ', seqRun, ' | cut -c1-8'), intern = TRUE), function(x) system(paste0('bkill ', x)))
    stop('Time elapsed error -- see log')
  } else {
    write('intSiteCaller completed.', file = 'log', append = TRUE)
    setwd(paste0('output/', seqRun, '_', refGenome))
    system('Rscript ../../intSiteCaller/check_stats.R > stats.txt')
    setwd(wd)
    
    write('Uploading intSites to database...', file = 'log', append = TRUE)
    
    system(paste0('Rscript intSiteUploader/intSiteUploader.R output/', seqRun, '_', refGenome, ' > ',
                  'output/', seqRun, '_', refGenome, '.uploadLog 2>&1'))
       
    # Recreate the intSiteSamples list which should now include the uploaded samples.
    dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
    intSiteSamples <- dbGetQuery(dbConn, 'select sampleName, refGenome from samples')
    intSiteSamples <- paste0(intSiteSamples$sampleName, '/', intSiteSamples$refGenome)
    
    if(! all(s$id %in% intSiteSamples)){
      msg <- paste0('These samples did not upload to the database: ', paste0(s$id[! s$id %in% intSiteSamples], collapse = ', '))
      write(msg, file = 'log', append = TRUE)
      stop(msg)
    } else {
      write('All samples uploaded to intSite database.', file = 'log', append = TRUE)
    }

    
    # Upload stats to database.
    system(paste0('Rscript ./uploadStats.R output/', seqRun, '_', refGenome))
  }
}

