#' Pretty print numbers with commas.
#' 
#' @param number Floating point of integer.
#' 
#' @export
ppNum <- function (n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)


#' Add position id, ie. chrX+231231 to GRange object.
#'
#' @param GRange Grange object.
#'
#' @export
addPositionID <- function (gr) {
  gr$posid <- paste0(seqnames(gr), strand(gr), start(flank(gr, -1, start = T)))
  gr
}




#' Create a integration site location plot.
#'
#' @param gr GRange object.
#' @param chromosomeLengths List with chromosome names which stores chromosomal lengths. 
#' @param alpha Alpha value of integration tick marks.
#' @param siteColor Color name or hex code for integration tick marks. 
#'
#' @export
intSiteDistributionPlot<- function (gr, chromosomeLengths, alpha = 0.8, siteColor = "black") 
{
  library(GenomicRanges)
  library(ggplot2)
  library(gtools)
  library(dplyr)
  d <- GenomicRanges::as.data.frame(gr)[, c("start", "seqnames")]
  d <- suppressWarnings(bind_rows(d, bind_rows(lapply(names(chromosomeLengths)[!names(chromosomeLengths) %in% 
                                                                                 unique(as.character(d$seqnames))], function(x) {
                                                                                   data.frame(start = 1, seqnames = x)
                                                                                 }))))
  d$seqnames <- factor(d$seqnames, levels = mixedsort(names(chromosomeLengths)))
  d <- lapply(split(d, d$seqnames), function(x) {
    lines <- data.frame(x = rep(x$start, each = 2), y = rep(c(0, 
                                                              1), nrow(x)), g = rep(1:nrow(x), each = 2), seqnames = x$seqnames[1])
    box <- data.frame(boxYmin = 0, boxYmax = 1, boxXmin = 1, 
                      boxXmax = chromosomeLengths[[as.character(x$seqnames[1])]], 
                      seqnames = x$seqnames[1])
    list(lines = lines, box = box)
  })
  sites <- do.call(rbind, lapply(d, "[[", 1))
  boxes <- do.call(rbind, lapply(d, "[[", 2))
  ggplot() + theme_bw() + geom_line(data = sites, alpha = alpha, 
                                    color = siteColor, aes(x, y, group = g)) + geom_rect(data = boxes, 
                                                                                         color = "black", alpha = 0, mapping = aes(xmin = boxXmin, 
                                                                                                                                   xmax = boxXmax, ymin = boxYmin, ymax = boxYmax)) + 
    facet_grid(seqnames ~ ., switch = "y") + scale_x_continuous(expand = c(0, 
                                                                           0)) + labs(x = "Genomic position", y = "") + theme(axis.text.y = element_blank(), 
                                                                                                                              axis.ticks.y = element_blank(), panel.grid.major = element_blank(), 
                                                                                                                              panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                                                                                              panel.border = element_blank(), strip.text.y = element_text(size = 12, 
                                                                                                                                                                                          angle = 180), strip.background = element_blank())
}


#' Abbreviate numbers with scientific suffixes.  
#'
#' @param number Floating point or integer.
#' @param rounding Boolean, round provided number if TRUE.
#' @param digits Number of significant digits to round number to if rounding set to TRUE.
#'
#' @export
numShortHand <- function (number, rounding = F, digits = ifelse(rounding, NA, 2)) 
{
  # https://stackoverflow.com/questions/11340444/is-there-an-r-function-to-format-number-using-unit-prefix/29932218
  
  lut <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06, 
           0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21, 
           1e+24, 1e+27)
  pre <- c("y", "z", "a", "f", "p", "n", "u", "m", "", "k", 
           "M", "G", "T", "P", "E", "Z", "Y", NA)
  ix <- findInterval(number, lut)
  if (ix > 0 && ix < length(lut) && lut[ix] != 1) {
    if (rounding == T && !is.numeric(digits)) {
      sistring <- paste(round(number/lut[ix]), pre[ix])
    }
    else if (rounding == T || is.numeric(digits)) {
      sistring <- paste(signif(number/lut[ix], digits), 
                        pre[ix])
    }
    else {
      sistring <- paste(number/lut[ix], pre[ix])
    }
  }
  else {
    sistring <- as.character(number)
  }
  sistring <- gsub("\\s", "", sistring)
  return(sistring)
}



#' Standardize genomic fragments from intSite experiments using the gintools package.
#'
#' @param frags GenomicRange object containing genomic fragments where fragments with positive strands 
#' are assumed to have their integration positions on the left and break points on the right (opposite assumed from negative strand).
#'
#' @param CPUs Number of CPUs to use for break point standardizations.
#' @param countsCol Name of meta data column containing the number of times fragments were observed. 
#'
#' @importFrom magrittr '%>%'
#' @export
stdIntSiteFragments <- function(frags, CPUs = 10, countsCol = 'reads'){ 
  # Setup parallelization.
  cluster <- parallel::makeCluster(CPUs)
  parallel::clusterExport(cl = cluster, envir = environment(), varlist = c('countsCol'))

  
  # Standardize break points by sample.
  frags <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, as.list(split(frags, frags$sampleName)), function(x){
    gintools::refine_breakpoints(x, counts.col = countsCol)
  })))
  
  
  # Standardize fragment start positions.
  frags <- gintools::standardize_sites(frags)


  # Stop the cluster.
   parallel::stopCluster(cluster)

  # Create a postion id now that fragment boundaries have been standardized.
  frags$posid <- paste0(GenomicRanges::seqnames(frags), GenomicRanges::strand(frags), GenomicRanges::start(GenomicRanges::flank(frags, -1, start = T)))


  # Now that the fragment positions have been adjusted, merge ranges and re-tally the counts.
  # Here we sort by sampleName which is the replicate level sample id so that the same fragment found in two or more
  # replicates are represeneted more than once.
  dplyr::group_by(data.frame(frags), sampleName, start, end, strand) %>%
  dplyr::mutate(reads = sum(reads)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}



#' Collapse technical replicates and estimate abundances by number of unique fragments associated with integration positions.
#'
#' @param f GRange object of standardized genomic fragments.
#' @return GRange object of integration positions including estimated abundances.
#' @importFrom magrittr '%>%'
#' @export
collapseReplicatesCalcAbunds <- function(f){
  # Conversion of GRange object to data frame creates width column.
  # Here we sum replicate specific fragment widths to estimate abundances.
  f <- data.frame(f)
  f$start       <- ifelse(as.character(f$strand) == '+', f$start, f$end)
  f$end         <- f$start
  f$sampleWidth <- paste0(f$sampleName, '/', f$width)
  dplyr::group_by(f, GTSP, posid) %>%
  dplyr::mutate(reads = sum(reads)) %>%
  dplyr::mutate(estAbund = dplyr::n_distinct(sampleWidth)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-sampleName, -sampleWidth) %>%
  dplyr::group_by(GTSP) %>%
  dplyr::mutate(relAbund = (estAbund / sum(estAbund))*100) %>%
  dplyr::ungroup() %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}


#' Estimate abundances by number of unique fragments associated with integration positions for each sample.
#'
#' @param f GRange object of standardized genomic fragments.
#' @return GRange object of integration positions including estimated abundances.
#' @importFrom magrittr '%>%'
#' @export
calcSampleAbunds <- function (f) 
{
  f <- data.frame(f)
  f$start <- ifelse(as.character(f$strand) == "+", f$start,  f$end)
  f$end <- f$start
  f$sampleWidth <- paste0(f$sampleName, "/", f$width)
  
  dplyr::group_by(f, sampleName, posid) %>% 
    dplyr::mutate(reads = sum(reads)) %>% 
    dplyr::mutate(estAbund = dplyr::n_distinct(sampleWidth)) %>% 
    dplyr::slice(1) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-sampleWidth) %>% 
    dplyr::group_by(sampleName) %>% 
    dplyr::mutate(relAbund = (estAbund/sum(estAbund)) * 100) %>% 
    dplyr::ungroup() %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}


  
 


#' Create a list of default file mappings.
#' 
#' An empty character vector object named 'none' is available for instances where 
#' gene lists are not available.
#'
#' @return List of default gt23 file mappings.
#'
#' @export
defaultGenomeFileMappings <- function(){
   list('hg38'    = list('genes'             = 'hg38.refSeqGenesGRanges', 
                         'exons'             = 'hg38.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'hg38.oncoGeneList',
                         'lymphomaGenesList' = 'hg38.lymphomaGenesList'),
        'hg18'    = list('genes'             = 'hg18.refSeqGenesGRanges', 
                         'exons'             = 'hg18.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'hg38.oncoGeneList',
                        'lymphomaGenesList'  = 'hg38.lymphomaGenesList'),
        'mm9'     = list('genes'             = 'mm9.refSeqGenesGRanges', 
                         'exons'             = 'mm9.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'mm9.oncoGeneList',
                         'lymphomaGenesList' = 'none'),
        'susScr3' = list('genes'             = 'susScr3.refSeqGenesGRanges', 
                         'exons'             = 'susScr3.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'none',
                         'lymphomaGenesList' = 'none'),
        'sacCer3' = list('genes'             = 'sacCer3.refSeqGenesGRanges', 
                         'exons'             = 'sacCer3.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'none',
                         'lymphomaGenesList' = 'none'),
        'canFam3' = list('genes'             = 'canFam3.humanXeno.refSeqGenesGRanges', 
                         'exons'             = 'canFam3.humanXeno.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'hg38.oncoGeneList',
                         'lymphomaGenesList' = 'hg38.lymphomaGenesList'),
        'macFas5' = list('genes'             = 'macFas5.humanXeno.refSeqGenesGRanges', 
                         'exons'             = 'macFas5.humanXeno.refSeqGenesGRanges.exons',
                         'oncoGeneList'      = 'hg38.oncoGeneList',
                         'lymphomaGenesList' = 'hg38.lymphomaGenesList'),
        'Kingella_kingae_KKKWG1' = list('genes'             = 'Kingella_kingae_KKKWG1.refSeqGenesGRanges', 
                                        'exons'             = 'Kingella_kingae_KKKWG1.refSeqGenesGRanges.exons',
                                        'oncoGeneList'      = 'none',
                                        'lymphomaGenesList' = 'none'))
        
 }
 
 
#' Validate that each named data object in genome file mapping list is present in the gt23 package.
#'
#' The special object 'none' is allowed which is an empty character vector object.
#'
#' @return Error if one or more named data objects are not present.
#'
#' @export
 validateGenomeFileMap <- function(x){
   availableObjects <- data(package='gt23')$results[,3]
   
   invisible(lapply(x, function(i){
     i <- unlist(i)
     i <- i[! is.na(i)]
     if(any(! i %in% availableObjects)) 
       stop(paste0(paste0(i[! i %in% availableObjects], collapse = ', '), ' could not be found in the package database. ',
                   'Availabe data objects include: ', paste0(availableObjects, collapse = ', ')))
   }))
 }
 
 
 #' Calculate nearest genomic and onocogene features. 
 #'
 #' @return GRange object updated with genomic feature annotations.
 #'
 #' @export
 annotateIntSites <- function(sites, CPUs = 20, genomeFileMap = NULL, oncoGeneList = NULL,
                              oncoGeneListSide = 'either', lymphomaGenesListSide = 'either'){
   
   if(is.null(genomeFileMap)) genomeFileMap <- gt23::defaultGenomeFileMappings()
   gt23::validateGenomeFileMap(genomeFileMap)
   
   # Here we convert GRanges to data frames after annoations are added and then back to a single 
   # GRange object so that different metadata cols can be concatenated where nonexistinting columns will be NA.
   
   GenomicRanges::makeGRangesFromDataFrame(dplyr::bind_rows(lapply(split(sites, paste(sites$patient, sites$refGenome)), function(x){
     
     message(paste('Starting', x$patient[1], '/', x$refGenome[1]))
     
     # Use the refGenome id to retrieve the appropriate gene boundary data.
     if(! x$refGenome[1] %in% names(genomeFileMap)) stop(paste('Error -- ', x$refGenome[1], ' is not defined in the genomeFileMap list.'))
                                                         
     genome_refSeq      <- eval(parse(text = paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$genes)))
     genome_refSeqExons <- eval(parse(text = paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$exons)))
     
     if(is.null(oncoGeneList)){
        oncoGeneList <- eval(parse(text = paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$oncoGeneList)))
     } else {
       message('Using custom oncogene list with ', dplyr::n_distinct(oncoGeneList), ' gene names.')
     }
     
     lymphomaGenesList  <- eval(parse(text = paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$lymphomaGenesList)))
     
     # Gene list check.
     message(paste0('oncoGeneList: ',      paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$oncoGeneList), ' - ', length(oncoGeneList), '  genes\n',
                    'lymphomaGenesList: ', paste0('gt23::', genomeFileMap[[x$refGenome[1]]]$lymphomaGenesList), ' - ', length(lymphomaGenesList), '  genes\n'))
             
     
     
     # Setup parallelization.
     cluster <- parallel::makeCluster(CPUs)
     parallel::clusterExport(cl=cluster, envir = environment(), varlist = c('CPUs', 'genome_refSeq', 'genome_refSeqExons', 'oncoGeneList', 'lymphomaGenesList'))
    
     # Create a splitting vector for parallelization.
     x$s <- dplyr::ntile(seq_along(x), CPUs)
     x$s <- dplyr::ntile(seq_along(x), CPUs)
     names(x) <- NULL
     
     x <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, split(x, x$s), function(x2){
     ###x <- unlist(GenomicRanges::GRangesList(lapply(split(x, x$s), function(x2){
          library(GenomicRanges)
          library(gt23)
          library(dplyr)
       
          # Create an order column to ensure that the incoming ranges keep the same order.
          x2$n <- 1:length(x2)
       
          # Nearest gene boundary
          x2 <- gt23::nearestGenomicFeature(x2, subject = genome_refSeq, subject.exons = genome_refSeqExons)
          x2 <- x2[order(x2$n)]
       
         # Nearest oncogene
         if(length(oncoGeneList) > 0){
           o <- gt23::nearestGenomicFeature(x2, subject = genome_refSeq, subject.exons = genome_refSeqExons, 
                                            geneList = oncoGeneList, subjectSide = oncoGeneListSide)
           o <- o[order(o$n)]
           
           d <- data.frame(GenomicRanges::mcols(x2))
           d <- d[order(d$n),]
           stopifnot(all(o$n == d$n))
           d$nearestOncoFeature       <- o$nearestFeature
           d$nearestOncoFeatureDist   <- o$nearestFeatureDist
           d$nearestOncoFeatureStrand <- o$nearestFeatureStrand
           GenomicRanges::mcols(x2) <- d
         }
     
         # Nearest lymphoma gene
         if(length(lymphomaGenesList) > 0){
           o <- gt23::nearestGenomicFeature(x2, subject = genome_refSeq, subject.exons = genome_refSeqExons, 
                                            geneList = lymphomaGenesList, subjectSide = lymphomaGenesListSide)
           o <- o[order(o$n)]
       
           d <- data.frame(GenomicRanges::mcols(x2))
           d <- d[order(d$n),]
           stopifnot(all(o$n == d$n))
           d$nearestlymphomaFeature       <- o$nearestFeature
           d$nearestlymphomaFeatureDist   <- o$nearestFeatureDist
           d$nearestlymphomaFeatureStrand <- o$nearestFeatureStrand
           GenomicRanges::mcols(x2) <- d
         }
       
         x2$n <- NULL
         x2
     })))
     
     x$s      <- NULL
     names(x) <- NULL
     parallel::stopCluster(cluster)
     data.frame(x)
   })), keep.extra.columns = TRUE)
 }
 
 
 
 
 

