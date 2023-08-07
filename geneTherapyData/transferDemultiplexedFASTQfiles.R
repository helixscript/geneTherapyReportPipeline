
localRunDirs  <- '~/test'
outputDir     <- '/media/lorax/data/BushmanGeneTherapy/demultiplexedINSPIIREDsamples'
remoteRunDirs <- system("ssh everett@microb120.med.upenn.edu ls /media/lorax/data/BushmanGeneTherapy/demultiplexedINSPIIREDsamples", intern = TRUE)
remoteFiles   <- system("ssh everett@microb120.med.upenn.edu find /media/lorax/data/BushmanGeneTherapy/demultiplexedINSPIIREDsamples -name *.fastq.gz", intern = TRUE)
remoteFiles   <- unlist(lapply(strsplit(remoteFiles, '/'), function(x) x[length(x)]))

f <- list.files(localRunDirs, pattern = '*.fastq.gz$', recursive = TRUE)
f <- f[! grepl('sacCer3', f)]
f <- f[grepl('GTSP\\d+\\-\\d+_', f)]

invisible(sapply(f, function(x){
  o <- unlist(strsplit(x, '/'))
  a <- unlist(strsplit(o[1], '_'))
  o[1] <- paste0(a[1:(length(a)-1)], collapse = '_')
  
  if(! o[1] %in% remoteRunDirs) system(paste0('ssh everett@microb120.med.upenn.edu mkdir ', file.path(outputDir, o[1])))
  
  if(! o[4] %in% remoteFiles) system(paste0('scp ', localRunDirs, '/', x, ' everett@microb120.med.upenn.edu:', file.path(outputDir, o[1]), '/'))
}))

