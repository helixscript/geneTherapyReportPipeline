library(dplyr)
library(readr)
library(RMySQL)

dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

projectDir <- '/media/lorax/data/export/projects'
outputDir  <- '/media/lorax/data/export/projects/exchange'
d <- data.frame(file = system(paste0('find ', projectDir, ' -name intSites.csv '), intern = TRUE))

d$trial   <- unlist(lapply(strsplit(as.character(d$file), '/'), function(x) x[length(x)-2]))
d$subject <- unlist(lapply(strsplit(as.character(d$file), '/'), function(x) x[length(x)-1]))

invisible(lapply(split(d, d$trial), function(x){

  f <- unlist(strsplit(as.character(x[1,]$file), '/'))
  trialName <- f[length(f) - 2]
  trialDir <- paste0(f[1:(length(f)-2)], collapse = '/')
  
  if(dir.exists(file.path(outputDir, trialName))) unlink(file.path(outputDir, trialName), recursive = TRUE)
  if(file.exists(paste0(file.path(outputDir, trialName), '.zip'))) file.remove(paste0(file.path(outputDir, trialName), '.zip'))
  
  dir.create(file.path(outputDir, trialName))
  dir.create(file.path(outputDir, trialName, 'reports'))
  reports <- system(paste0('find ', trialDir, ' -name *.pdf '), intern = TRUE)
  
  invisible(lapply(reports, function(r) system(paste0('cp ', r, ' ', file.path(outputDir, trialName, 'reports')))))
  
  r <- bind_rows(lapply(as.character(x$file), function(x2){
         x2 <- read_csv(x2)
         x2$patient <- as.character(x2$patient)
         x2$cellType <- as.character(x2$cellType)
         x2$dataSource <- as.character(x2$dataSource)
         x2$timePoint <- as.character(x2$timePoint)
         x2
       })) %>% left_join(select(samples, SpecimenAccNum, SamplePatientCode), by = c('GTSP' = 'SpecimenAccNum'))
  
  o <- select(r, seqnames, start, strand, refGenome, reads, patient, SamplePatientCode, GTSP, cellType,
              timePoint, estAbund, relAbund, nearestFeature,	inFeature, nearestFeatureStrand, inFeatureExon,	
              nearestFeatureDist,	nearestOncoFeature,	nearestOncoFeatureDist)
   
  names(o) <- c('chromosome',	'position',	'strand','refGenome', 'reads', 'subject', ',externalSampleID',
                'internalSampleID', 'cellType',	'timePoint',	'estAbund',	'relAbund',	'nearestFeature',	'inFeature',	
                'nearestFeatureStrand',	'inFeatureExon',	'nearestFeatureDist',	'nearestOncoFeature', 'nearestOncoFeatureDist')
  
  write_tsv(o, file.path(outputDir, trialName, 'intSites.tsv'), col_names = TRUE)
  
  system(paste0('zip -r ', paste0(file.path(outputDir, trialName), '.zip'), ' ', file.path(outputDir, trialName)))
}))