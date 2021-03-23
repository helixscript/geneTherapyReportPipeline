# gt23

Basic workflow from sample list:  

```
library(dplyr)
samples  <- c('GTSP2683', 'GTSP2684', 'GTSP2685', 'GTSP2686')
intSites <- gt23::getDBgenomicFragments(samples, 'specimen_management', 'intsites_miseq') %>%
            gt23::stdIntSiteFragments() %>%
            gt23::collapseReplicatesCalcAbunds() %>%
            gt23::annotateIntSites()
```