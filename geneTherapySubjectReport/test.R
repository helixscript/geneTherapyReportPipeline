library(RMySQL)
library(dplyr)
options(useFancyQuotes = FALSE)

dbGroup <- 'intsites_miseq'

sampleIDs <- c('GTSP5030', 'GTSP5061', 'GTSP5068')


sql <- paste0("select samples.sampleName, samples.refGenome, multihitpositions.multihitID, ",
              "multihitlengths.length from multihitlengths left join multihitpositions on ",
              "multihitpositions.multihitID = multihitlengths.multihitID left join samples on ",
              "samples.sampleID = multihitpositions.sampleID where sampleName like ", 
              paste0(sQuote(paste0(sampleIDs, '-%')), collapse = ' or sampleName like '))

dbConn  <- dbConnect(MySQL(), group = args$intSiteDB)
sites.multi <- unique(dbGetQuery(dbConn, sql))
dbDisconnect(dbConn)

if( nrow(sites.multi) > 0 ){
  
  replicateLevelTotalAbunds <- 
    group_by(data.frame(intSites.std.reps), sampleName, posid) %>%
    mutate(estAbund = n_distinct(width)) %>%
    ungroup() %>%
    group_by(sampleName) %>%
    summarise(totalReplicateAbund = sum(estAbund)) %>%
    ungroup()
  
  sites.multi$GTSP <- str_extract(sites.multi$sampleName, 'GTSP\\d+') 
  
  sites.multi <- left_join(sites.multi, replicateLevelTotalAbunds, by = 'sampleName') %>%
    group_by(sampleName, multihitID) %>%
    mutate(estAbund = length(unique(length)),
           relAbund = ((estAbund / sum(totalReplicateAbund[1], estAbund))*100)) %>%
    summarise(totalReplicateCells = sum(totalReplicateAbund[1], estAbund), relAbund = relAbund[1], GTSP = GTSP[1]) %>%
    ungroup() %>%
    left_join(select(sampleData, CellType, Timepoint, SpecimenAccNum), c("GTSP" = "SpecimenAccNum")) %>%
    filter(relAbund >= 20) %>%
    select(-GTSP) %>%
    arrange(desc(relAbund)) %>%
    mutate(relAbund = sprintf("%.1f%%", relAbund))
  
  names(sites.multi) <- c('Replicate', 'Multihit id', 'Total cells in replicate', 'Multihit relative Abundance', 'Celltype', 'Timepoint')
}
