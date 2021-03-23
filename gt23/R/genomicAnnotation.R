#' Annotate genomic ranges
#'
#' Create a dataframe of sequencing run ids associated with 1 or more patient identifiers.
#'
#' @param query GRange object containing ranges from a genome supported by the gt23 package.
#' @param subject gt23 genome specific gene boundary object.
#' @param subject.exons gt23 genome specific gene exon object.
#' @param subjectSide Which subject boundary should be tested ('start', 'end', 'midpoint', 'either').
#' @param geneList Case insensitive list of gene names to filter nearest features against.
#'
#' @return An updated query GRange object.
#'
#' @export
nearestGenomicFeature <- function(query, subject = NULL, subject.exons = NULL, subjectSide = 'either', geneList = NULL){
  
  if(is.null(subject))       stop('subject parameter can not be NULL.')
  if(is.null(subject.exons)) stop('subject.exons parameter can not be NULL.')
  if(! is.null(geneList))    subject <- GenomicRanges::subset(subject, toupper(name2) %in% toupper(geneList))
  
  message('Starting nearestGenomicFeature() -- subject contains ', length(ranges), ', subject.exon contains ', length(subject.exons), ' ranges.')
  
  # If subjectSide is not set to either, collapse the subject ranges to single positions
  if (subjectSide %in% c("start", "end", "midpoint")) {
    message('subjectSide parameter: ', subjectSide)
    options(warn=-1)
    if (subjectSide == "start") subject <- GenomicRanges::flank(subject, width = -1)
    if (subjectSide == "end")   subject <- GenomicRanges::flank(subject, width = -1, start = FALSE)
    if (subjectSide == "midpoint") ranges(subject) <- IRanges(mid(ranges(subject)), width = 1)
    options(warn=0)
  }
  
  options(stringsAsFactors = FALSE)
  
  query.df  <- GenomicRanges::as.data.frame(query)
  subject.df <- GenomicRanges::as.data.frame(subject)
  
  query.df$strand <- as.character(query.df$strand)
  subject.df$strand <- as.character(subject.df$strand)
  
  subject.df$name2 <- as.character(subject.df$name2)
  
  
  subject.exons.df <- GenomicRanges::as.data.frame(subject.exons)
  query.df$inFeature            <- FALSE
  query.df$nearestFeature       <- 'None.found'
  query.df$nearestFeatureStrand <- 'None.found'
  query.df$inFeatureExon        <- FALSE
  query.df$inFeatureSameOrt     <- FALSE
  query.df$nearestFeatureStart  <- Inf
  query.df$nearestFeatureEnd    <- Inf
  query.df$nearestFeatureDist   <- Inf
  
  o  <- suppressWarnings(GenomicRanges::nearest(query, subject, select='all', ignore.strand=TRUE))
  
  if(length(o) > 0){
    
    createCol <- function(a, b, n){
      paste0(unique(cbind(a, b))[,n], collapse=',')
    }
    
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
         dplyr::mutate(gene   = createCol(subject.df[subjectHits,]$name2, subject.df[subjectHits,][['strand']], 1),
                       strand = createCol(subject.df[subjectHits,]$name2, subject.df[subjectHits,][['strand']], 2),
                       hitStart = min(subject.df[subjectHits,][['start']]),
                       hitEnd   = max(subject.df[subjectHits,][['end']])) %>% 
         dplyr::ungroup() %>%
         dplyr::select(queryHits, gene, strand, hitStart, hitEnd) %>% 
         dplyr::distinct() %>% 
         data.frame()
    
    if(nrow(a) > 0){
      query.df[a$queryHits,]$nearestFeature       <- a$gene
      query.df[a$queryHits,]$nearestFeatureStrand <- a$strand
      query.df[a$queryHits,]$nearestFeatureStart  <- a$hitStart
      query.df[a$queryHits,]$nearestFeatureEnd    <- a$hitEnd
    }
  }
  
  o <- suppressWarnings(GenomicRanges::findOverlaps(query, subject, select='all', ignore.strand=TRUE, type='any'))
  
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = paste(unique(subject.df[subjectHits,]$name2), collapse=',')) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    if(nrow(a) > 0) query.df[a$queryHits,]$inFeature <- TRUE
  }
  
  o <- suppressWarnings(GenomicRanges::distanceToNearest(query,  subject, select='all', ignore.strand=TRUE))
  
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::top_n(-1, distance) %>%
      dplyr::ungroup() %>%
      dplyr::select(queryHits, distance) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    
    if(nrow(a) > 0) query.df[a$queryHits,]$nearestFeatureDist <- a$distance
  }
  
  
  # JKE
  query.df$nearestFeatureBoundary <- ifelse(abs(query.df$start - query.df$nearestFeatureStart) > 
                                              abs(query.df$start - query.df$nearestFeatureEnd),   
                                            query.df$nearestFeatureEnd,  
                                            query.df$nearestFeatureStart)
  
  query.df$nearestFeatureDist <- query.df$nearestFeatureDist * sign(query.df$start - query.df$nearestFeatureBoundary)
  query.df$nearestFeatureDist <- ifelse(query.df$nearestFeatureStrand=='+', query.df$nearestFeatureDist, query.df$nearestFeatureDist * -1)
  
  query.df$nearestFeatureStart    <- NULL
  query.df$nearestFeatureEnd      <- NULL
  query.df$nearestFeatureBoundary <- NULL
  
  # In exon THIS ???
  o <- suppressWarnings(GenomicRanges::findOverlaps(query, subject.exons, select='all', ignore.strand=TRUE, type='any'))
 
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = paste(unique(subject.exons.df[subjectHits,]$name2), collapse=',')) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    #query.df[a$queryHits,]$inFeatureExon  <- a$gene
    if(nrow(a) > 0) query.df[a$queryHits,]$inFeatureExon  <- TRUE
  }
  
  # In TU ort
  # There may be cases where a site overlaps two or more features which have the sampe orientation.
  # ie. +,+,+ and we want to reduce these down to a single unique sign for comparison.

  query.df$nearestFeatureStrand2 <- unlist(lapply(strsplit(query.df$nearestFeatureStrand, ','), function(x){ paste(unique(x), collapse=',')}))
  i <- grepl(',', query.df$nearestFeatureStrand2)
  if(any(i)) query.df[i,]$nearestFeatureStrand2 <- '*'
  query.df$inFeatureSameOrt <- query.df$strand == query.df$nearestFeatureStrand2
  query.df$nearestFeatureStrand2 <- NULL
  
  GenomicRanges::makeGRangesFromDataFrame(query.df, keep.extra.columns = TRUE)
}




#' Annotate genomic ranges
#'
#' Create a dataframe of sequencing run ids associated with 1 or more patient identifiers.
#'
#' @param d Data fram containing genomic ranges to be added to UCSC track.
#' @param abundCuts Cut points for estimated abundance (estAbund) values. 
#' @param posColors Color codes for binned abundances (positive integrations).
#' @param negColors Color codes for binned abundances (positive integrations).
#' @param title Track title.
#' @param outputFile Track output file.
#' @param visibility Track default visibility (0 - hide, 1 - dense, 2 - full, 3 - pack, and 4 - squish).
#' @param position Deafult track position.
#' @param padSite Number of NTs to pad sites with for increased visibility.
#' @param siteLabel Text to appear next to sites, ie. 'Patient X, chr12+1052325'.
#' 
#' @return Nothing.
#'
#' @export
createIntUCSCTrack <- function(d, abundCuts = c(5,10,50), 
                                  posColors = c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB"),
                                  negColors = c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000"),
                                  title = 'intSites', outputFile = 'track.ucsc', visibility = 1, 
                                  position = 'chr7:127471196-127495720', padSite = 0,
                                  siteLabel = NA){

  # Check function inputs.
  if(length(posColors) != length(negColors)) 
    stop('The pos and neg color vectors are not the same length.')
  
  if(length(abundCuts) != length(posColors) - 1) 
    stop('The number of aundance cut offs must be one less than the number of provided colors.')
  
  if(! all(c('start', 'end', 'strand', 'seqnames', 'estAbund') %in% names(d))) 
    stop("The expected column names 'start', 'end', 'strand', 'seqnames', 'estAbund' were not found.") 
  
  if(is.na(siteLabel) | ! siteLabel %in% names(d)) 
    stop('The siteLabel parameter is not defined or can not be found in your data.')
  
  
  # Cut the abundance data. Abundance bins will be used to look up color codes.
  # We flank the provided cut break points with 0 and Inf in order to bin all values outside of breaks.
  cuts <- cut(d$estAbund, breaks = c(0, abundCuts, Inf), labels = FALSE)
  
  
  # Convert Hex color codes to RGB color codes. 
  # col2rgb() returns a matrix, here we collapse the columns into comma delimited strings.
  #   grDevices::col2rgb(posColors)
  #         [,1] [,2] [,3] [,4]
  #   red    140  103   66   29
  #   green  157  104   52    0
  #   blue   255  227  199  171
  
  posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
  negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')
  
  
  # Create data fields needed for track table.
  d$score <- 0
  d$color <- ifelse(d$strand == '+', posColors[cuts], negColors[cuts])
  
  # Pad the site n NTs to increase visibility.
  if(padSite > 0){
    d$start <- floor(d$start - padSite/2)
    d$end   <- ceiling(d$end + padSite/2)
  }
  
  # Define track header.
  trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s\nbrowser position %s",
                       title, title, visibility, position)
  
  # Write out track table.
  write(trackHead, file = outputFile, append = FALSE)
  write.table(d[, c('seqnames', 'start', 'end', siteLabel, 'score', 'strand', 'start', 'end', 'color')], 
              sep = '\t', col.names = FALSE, row.names = FALSE, file = outputFile, append = TRUE, quote = FALSE)
}

