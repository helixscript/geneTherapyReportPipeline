
o <- readLines('/media/sequencing/Illumina/220603_MN01490_0081_A000H3WT25/SampleSheet.csv')


o[which(grepl('\\[metaData\\]', o))+1:(length(o)-(which(grepl('\\[metaData\\]', o))))]