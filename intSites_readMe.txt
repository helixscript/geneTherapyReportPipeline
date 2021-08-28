Each row in the intSite.tsv spreadsheet describes an integration site in the genome (refGenome). 
The first two fields (chromosome and position) define the position of an integration site while 
the strand field defines in the vector integrated in the forward (+) or reverse (-) orientation. 
The (reads) field reports the number of sequencing reads which supporting an integration site.

An estimated abundance (estAbund) is provided for each site which is a lower end estimate of the 
number of cells harboring an integration. The relative abundance (relAbund) is the estimated 
proportion of cells within a sample harobring an integration. 

The (nearestFeature) field is the name of nearest annotated gene or RNA where the (nearestFeatureDist) 
field contains the number of NT between the integration site and the nearest feature boundary. 
The (inFeature) field is set to TRUE when the integration is within a gene's transcription unit while 
the (inFeatureExon) field is set to TRUE if the integration site is within an annotated gene exon.

The (nearestOncoFeature) and (nearestOncoFeatureDist) fields are similiar to the (nearestFeature) and 
(nearestFeatureDist) except they refer to the nearest oncogene to the integration annoated for 
the reference genome.

