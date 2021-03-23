# intSiteDistributionPlot <- function (gr, chromosomeLengths, alpha = 0.8, siteColor = 'black')
# {
#   library(GenomicRanges)
#   library(ggplot2)
#   library(gtools)
#   library(dplyr)
#   
#   d <- GenomicRanges::as.data.frame(gr)[, c("start", "seqnames")]
#   
#   # Add mock rows for chromosomes missing from the data.
#   d <- suppressWarnings(bind_rows(d, bind_rows(lapply(names(chromosomeLengths)[!names(chromosomeLengths) %in% unique(as.character(d$seqnames))], 
#                                                       function(x){
#                                                         data.frame(start=1, seqnames=x)
#                                                       }))))
#   
#   
#   d$seqnames <- factor(d$seqnames, levels = mixedsort(names(chromosomeLengths)))
#   
#   d <- lapply(split(d, d$seqnames), function(x) {
#     lines <- data.frame(x=rep(x$start, each=2), y=rep(c(0,1), nrow(x)), g=rep(1:nrow(x), each=2), seqnames=x$seqnames[1])
#     box <- data.frame(boxYmin = 0, boxYmax = 1, boxXmin = 1,
#                       boxXmax = chromosomeLengths[[as.character(x$seqnames[1])]], seqnames=x$seqnames[1])
#     list(lines = lines, box = box)
#   })
#   
#   sites <- do.call(rbind, lapply(d, "[[", 1))
#   boxes <- do.call(rbind, lapply(d, "[[", 2))
#   
#   ggplot() +
#     theme_bw() +
#     geom_line(data=sites, alpha=alpha, color=siteColor, aes(x, y, group=g)) +
#     geom_rect(data  = boxes,
#               color = "black",
#               alpha = 0,
#               mapping = aes(xmin = boxXmin, xmax = boxXmax, ymin = boxYmin, ymax = boxYmax)) +
#     facet_grid(seqnames~., switch='y') +
#     scale_x_continuous(expand = c(0, 0)) +
#     labs(x='Genomic position', y='') +
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.background = element_blank(),
#           panel.border=element_blank(),
#           strip.text.y = element_text(size = 12, angle = 180),
#           strip.background = element_blank())
# }
# 
# 
# 
# 
# 
# 
# 
# # Function for plotting the data obtained from genomicHeatmap2dataframe().
# genericHeatmap <- function(d, tileColors=colorRampPalette(c('blue4', 'white', 'gold3'))(21)){
#   
#   heatmap_dims <- function(p) {
#     .x <- as.character(p$mapping$x)
#     .y <- as.character(p$mapping$y)
#     ncols <- length(unique(p$data[[.x]]))
#     nrows <- length(unique(p$data[[.y]]))
#     return(list(ncols=ncols, nrows=nrows))
#   }
#   
#   make_square <- function(p, fudge=1) {
#     dims <- heatmap_dims(p)
#     p + ggplot2::theme(aspect.ratio = (dims$nrows/dims$ncols)*fudge)
#   }
#   
#   p <- ggplot(d, aes(sample, test, fill=score)) +
#     theme_bw() +
#     labs(x='', y='') +
#     geom_tile(color='black') +
#     scale_x_discrete(position = "top") +
#     scale_fill_gradientn(colors=tileColors,
#                          breaks=seq(from=0, to=1, by=0.2),
#                          labels=seq(from=0, to=1, by=0.2),
#                          limits=c(0,1),
#                          name='Color Key') +
#     theme(axis.text.x=element_text(angle=90, hjust=0),
#           #legend.position="bottom",
#           plot.background = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank())
#   
#   make_square(p)
# }