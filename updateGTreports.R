sequencingArchiveDir <- '/media/sequencing/Illumina/'
sampleDB.group       <- 'specimen_management'
intSiteDB.group      <- 'intsites_miseq'
reportOutputDir      <- '/media/lorax/data/export/projects'
Rscript              <- '/opt/R-3.4.4-20180823/lib64/R/bin/Rscript'
softwareDir          <- '/media/lorax/data/software/geneTherapyReports'
CPUs                 <- 10
legacyData           <- '/media/lorax/data/software/geneTherapyReports/geneTherapyData/legacyData.rds'

if(file.exists(file.path(reportOutputDir, 'LOCK'))) q()
system(paste('touch', file.path(reportOutputDir, 'LOCK')))

log <- file.path(reportOutputDir, 'log')
write(date(), file = log)

# Update sequencing run reports.
# Writes log to {reportOutputDir}/log
system(file.path(softwareDir, 'updateSeqRunReports.R'))

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
  
  # Skip seq run trials.
  if(grepl('seqRuns', x)) return(NA)

  o <- unlist(strsplit(x, '/'))
  trial <- o[length(o)-2]
  patient <-  o[length(o)-1]
  
    
  # Retrieve snap show of sample data from the time the report was last created.
  file <- paste0(reportOutputDir, '/', paste0(o[1:(length(o)-1)], collapse = '/'), '/', 'sampleData.csv')
  
  
  if(file.exists(file)){
    sd <- readr::read_csv(file, col_types = cols(), trim_ws = FALSE) %>% 
          filter(! SpecimenAccNum %in% legacyData$GTSP) %>%  # Remove legacy data from report table since it would not be in db.
          arrange(SpecimenAccNum)
    sd2 <- subset(samples, SpecimenAccNum %in% sd$SpecimenAccNum) %>% select(SpecimenAccNum, Trial, CellType, Timepoint) %>% arrange(SpecimenAccNum) 

    availableIntSiteSamples <- unique(subset(intSitesamplesTbl, Patient == patient & Trial == trial)$SpecimenAccNum)
    
    # Rebuild report if different intSiteData data is available. 
    delete <- FALSE
    
    if(length(availableIntSiteSamples) != length(sd$SpecimenAccNum)){
      delete <- TRUE 
    } else {
      if(! all(sort(availableIntSiteSamples) == sort(sd$SpecimenAccNum))) delete <- TRUE
    }
    
    if (nrow(sd) != nrow(sd2)){
       delete <- TRUE
    } else {
      if(length(sd) == length(sd2)){
        if(! all(sd == sd2)){
          delete <- TRUE
        }
      }
    }
    
    if(delete){
      message('\nDeleting patient directiory due to db mismatch: ', paste0(o[1:(length(o)-1)], collapse = '/'))
      unlink(paste0(reportOutputDir, '/', paste0(o[1:(length(o)-1)], collapse = '/')), recursive = TRUE)
    }
  }
  
  # Trial / patient does not exist in database -- delete patient directory
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
  clusterExport(cluster, c('Rscript', 'log', 'sampleDB.group', 'intSiteDB.group', 'reportOutputDir', 'softwareDir', 'intSitesamplesTbl', 'sampleDB.group'))
  
  o <- bind_rows(parLapply(cluster, split(r, r$n), function(x){
  #o <- bind_rows(lapply(split(r, r$n), function(x){
         library(dplyr)
         library(RMySQL)
         bind_rows(lapply(split(x, paste(x$patient, x$trial)), function(x2){
      
           command <- paste0(Rscript,  ' ', softwareDir, '/geneTherapySubjectReport/report.R ',
                             '--specimenDB  ', sampleDB.group, ' --intSiteDB ', intSiteDB.group, ' ',
                             '--patient "', x2$patient, '" --trial "', x2$trial, '" ',
                             '--outputDir ', file.path(reportOutputDir, 'trials', x2$trial), ' ',
                             '--noSiteReportFile ', file.path(softwareDir, 'geneTherapySubjectReport', 'report_noSites.Rmd'), ' ',
                             '--reportFile ', file.path(softwareDir, 'geneTherapySubjectReport', 'report.Rmd'), ' ',
                             '--legacyData ', file.path(softwareDir,  'geneTherapyData', 'legacyData.rds'))
           
           message(command)
           write(c(command, '\n'), file = log, append = TRUE)
  
      
           invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
           dbConn  <- dbConnect(MySQL(), group=sampleDB.group)
           
           
           # Create a flag in the specimen database to update run analyses.
           sapply(unique(subset(intSitesamplesTbl, Patient == x2$patient & Trial == x2$trial)$miseqid), function(run){
             comm <- paste0('insert into seqRunAnalysis (miseqid, updateReport) values ("', run, '", "1") ON DUPLICATE KEY UPDATE updateReport="1"')
             DBI::dbSendQuery(dbConn, comm)
           })
           
           dbDisconnect(dbConn)
           
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
             result <- 'failed'
           }
           
           tibble(trial = x2$trial, patient = x2$patient, result = result)
    }))
  }))
  
  write.table(o, sep = '\t', file = log, quote = FALSE, row.names = FALSE, append = TRUE)
}


# Update data exchange archives.
system(file.path(softwareDir, 'updateDataExchanges.R'))


# Create trial level pages.
# fairly quick -- fine to regenerate all reports with each script run.

dirs <- list.dirs(file.path(reportOutputDir, 'trials'), recursive = FALSE, full.names = TRUE)

for(d in dirs){
  o <- unlist(strsplit(d, '/'))
  message(o[length(o)], '\n\n')
  rmarkdown::render(file.path(softwareDir, 'trialReport.Rmd'),
                    output_file = file.path(d, 'index.html'),
                    params = list('title' = gsub('_', ' ', o[length(o)]),
                                  'trial' = o[length(o)],
                                  'reportPath' = d))
  gc()
}


#createTrialReports <- data.frame(result = unlist(createTrialReports))
#createTrialReports$trial <- list.dirs(file.path(reportOutputDir, 'trials'), recursive = FALSE)


if(! dir.exists(file.path(reportOutputDir, 'dashboard'))) dir.create(file.path(reportOutputDir, 'dashboard'))

# Create overview page.
rmarkdown::render(file.path(softwareDir, 'dashboard.Rmd'), output_file = file.path(reportOutputDir, 'dashboard', 'index.html'),  
                  params = list('reportOutputDir' = reportOutputDir,
                                'sequencingArchiveDir' = sequencingArchiveDir))

file.remove(file.path(reportOutputDir, 'LOCK'))