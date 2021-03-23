runIntSiteCaller <- function(outputDir, refGenome, s){
  
  if(dir.exists(outputDir)) stop('Error - analysis directory already exists.')
  dir.create(outputDir)
  dir.create(paste0(outputDir, '/Data'))
  
  # Test ssh keys and data presence.
  r <- system(paste('ssh ', remoteHost, ' ls ', remoteHostDataDir, ' > /dev/null 2>&1'), intern = FALSE)
  if(r == 255) stop('Could not connect to remote host.')
  if(r == 2) stop('Remote data directory not present.')
  
  r <- system(paste('ssh ', remoteHost, ' ls ', 
                    paste0(remoteHostDataDir, '/', seqRun, '/Data/Intensities/BaseCalls/Undetermined* > /dev/null 2>&1')), intern= FALSE)
  if(r == 2) stop('Remote data files not present. No Undetermined FASTQ files found.')
  
  r <- system(paste('ssh ', remoteHost, ' ls ', paste0(remoteHostDataDir, '/', seqRun, '/Data/Intensities/BaseCalls/Undetermined*')), intern = TRUE)
  
  # Force metaData to a single reference genome.
  s$refGenome <- refGenome
  
  # Write table to analysis directory.
  write.table(s, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE, file = paste0(outputDir, '/sampleInfo.tsv'))
  
  # Transfer vector files.
  if(! all(unique(s$vectorSeq) %in% list.files('intSiteCaller_vectors'))) stop('All the vector files in the sample info file were not found in the vector file directory.')
  
  invisible(sapply(unique(s$vectorSeq), function(x){
    system(paste0('cp  intSiteCaller_vectors/', x, '  ', outputDir, '/'))
  }))
  
  # Retrieve sequencing files from server.
  invisible(sapply(r, function(x){
    message('Retrieving ', x)
    system(paste0('scp ', remoteHost, ':', x, ' ', paste0(outputDir, '/Data/')))
  }))
  
  # Check for the expected number of transfered files.
  if(! length(list.files(paste0(outputDir, '/Data/'))) == 3) stop('Error - the expected number of FASTQ files were not retrieved.')
  
  # Update INSPIIRED control template file and copy to the analysis directory.
  o <- readLines('INSPIIRED.yml')
  write(sub('runIDplaceHolder', seqRun, o), file = paste0(outputDir, '/INSPIIRED.yml'))
  
  # Change into the analysis directory and start intSiteCaller.
  setwd(outputDir)
  
  i <- unlist(strsplit(seqRun, '[-_]'))
  i <- paste0(i[length(i)], '_', refGenome)
  
  message('Starting intSiteCaller with run id: ', i)
  system(paste0('Rscript ../../intSiteCaller/intSiteCaller.R -j ', i))
  setwd(wd)
  return(TRUE)
}


metaData <- function(){
  o <- system(paste0('ssh ', remoteHost, ' "cat ', remoteHostDataDir, '/', seqRun, '/SampleSheet.csv"'), intern = TRUE)
  o <- unlist(strsplit(o, '[\r\n]+'))
  if(! any(grepl('\\[metaData\\]', o))) stop('Error - could not retrieve [metaData] block from SampleSheet.')
  o <- o[(grep('\\[metaData\\]', o)+1):length(o)]
  s <- read.table(textConnection(o), sep = ',', header = TRUE)

  # Override provided gender with male in order to analyze Y chromosome.
  s$gender <- 'm'

  # Test for missing vector file names.
  # Assume all file names are at least 3 characters long.
  if(! all(nchar(s$vectorSeq) >= 3)) stop('Error -- there appears to be missing vector file entries.')

  # Update genomes to the latest supported genomes.
  s$refGenome <- ifelse(grepl('^hg', s$refGenome, ignore.case = TRUE), 'hg38', s$refGenome)
  s$refGenome <- ifelse(grepl('^mm', s$refGenome, ignore.case = TRUE), 'mm9',  s$refGenome)
  s$refGenome <- ifelse(grepl('^sacCer', s$refGenome, ignore.case = TRUE), 'sacCer3',  s$refGenome)
  
  # Switch reference genome to Yeast for positive controls
  i <- grep('PositiveControl', s$alias, ignore.case = TRUE)
  if(length(i) > 0) s[i,]$refGenome <- 'sacCer3'
  
  s
}
