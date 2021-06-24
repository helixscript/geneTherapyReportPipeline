library(RMySQL)

intSiteDB.group  <- 'intsites_miseq'


invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group=intSiteDB.group)

runID <- '210611_MN01490_0021_A000H3FK2C'


commands <- c("delete from multihitlengths where multihitID in (select multihitID from multihitpositions where sampleID in (select sampleID from samples where miseqid='xxx'))",
              "delete from multihitpositions where sampleID in (select sampleID from samples where miseqid='xxx')",
              "delete from pcrbreakpoints where siteID in (select siteID from sites  where sampleID in (select sampleID from samples where miseqid='xxx'))",
              "delete from sites where sampleID in (select sampleID from samples where miseqid='xxx')",
              "delete from samples where miseqid='xxx'",
              "delete from intsitecallerstats where miseqid='xxx'")

commands <- gsub('xxx', runID, commands)

invisible(sapply(commands, function(comm){
  DBI::dbSendQuery(dbConn, comm)
}))


