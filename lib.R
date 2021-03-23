le <- function(x, s){
  o <- unlist(strsplit(x, s))
  o[length(o)]
}


readSamplePlot <- function(reads, n){
  ds <- as.character(sample(unique(reads), n))
  dp <- lapply(strsplit(sort(ds), ''), function(x){ tibble(base = x, n = 1:length(x)) })
  dp <- bind_rows(mapply(function(x, n){ x$read <- n; x}, dp, 1:length(dp), SIMPLIFY = FALSE))
  dp$base <- factor(dp$base, levels = c('A', 'T', 'C', 'G', 'N'))
  
  ggplot(dp, aes(n, read, fill = base)) + theme_bw() + geom_tile() +
    scale_fill_manual(values =  c('red', 'green', 'blue', 'gold', 'gray50')) +
    scale_x_continuous(limits = c(0, width(reads[1])), expand = c(0, 0)) +
    scale_y_continuous(label=comma, limits = c(0, n), expand = c(0, 0)) +
    labs(x = 'Position', y = 'Reads') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}



getDBgenomicFragments <- function(samples, intSiteDB.group){
  options(useFancyQuotes = FALSE)
  
  dbConn <- DBI::dbConnect(RMySQL::MySQL(), group=intSiteDB.group)
  
  intSiteSamples <- DBI::dbGetQuery(dbConn, 'select * from samples')
  intSiteSamples$GTSP <- gsub('\\-\\d+$', '', intSiteSamples$sampleName)
  
  sampleTable <-  base::subset(intSiteSamples, sampleName %in% samples)
  
  sampleIDs <- sampleTable$sampleID
  
  if(length(sampleIDs) == 0) return(GenomicRanges::GRanges())
  
  replicateQuery <- paste('samples.sampleID in (', paste0(sampleIDs, collapse = ','), ')')
  
  q <- sprintf("select position, chr, strand, breakpoint, count, refGenome,
               sampleName from sites left join samples on
               sites.sampleID = samples.sampleID
               left join pcrbreakpoints on
               pcrbreakpoints.siteID = sites.siteID
               where (%s)", replicateQuery)
    
    sites <- DBI::dbGetQuery(dbConn, q)
    
    DBI::dbDisconnect(dbConn)
    
    if(nrow(sites) == 0) return(GRanges())
    
    sites$sampleName <- paste0(sites$refGenome, '-', sites$sampleName)
    
    
    sites$GTSP      <- as.character(sub('\\-\\d+$', '', sites$sampleName))
    sites$patient   <- 'PositiveControl'
    sites$timePoint <- 'd0'
    sites$cellType  <- 'control'
    
    sites <- dplyr::bind_cols(sites, expandTimePoints(sites$timePoint))
    
    intSites <-
      GenomicRanges::GRanges(seqnames = S4Vectors::Rle(sites$chr),
                             ranges   = IRanges::IRanges(start = pmin(sites$position, sites$breakpoint),
                                                         end   = pmax(sites$position, sites$breakpoint)),
                             strand   = S4Vectors::Rle(sites$strand))
    
    sites <- dplyr::rename(sites, reads = count)
    GenomicRanges::mcols(intSites) <- data.frame(dplyr::select(sites, refGenome, reads, patient, sampleName,
                                                               GTSP, cellType, timePoint, timePointDays, timePointMonths))
    intSites
  }