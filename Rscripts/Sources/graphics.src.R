## Graphic utils
# Author: Binyamin Knisbacher
library(ggplot2)
library(ggdendro)
library(ggpubr)
library(ggthemes)

### ggplot2 custom themes ###
#Note: to work, must be before other calls to theme()
no_bg_theme = theme(panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_text(colour="black", size=10),
                    axis.text.y = element_text(colour="black", size=10),
                    axis.ticks = element_line(colour="black"))

### Palettes ###
diff_colors_hexa = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#000000",
"#800000", "#008000", "#000080", "#808000", "#800080", "#008080", "#808080",
"#C00000", "#00C000", "#0000C0", "#C0C000", "#C000C0", "#00C0C0", "#C0C0C0",
"#400000", "#004000", "#000040", "#404000", "#400040", "#004040", "#404040",
"#200000", "#002000", "#000020", "#202000", "#200020", "#002020", "#202020",
"#600000", "#006000", "#000060", "#606000", "#600060", "#006060", "#606060",
"#A00000", "#00A000", "#0000A0", "#A0A000", "#A000A0", "#00A0A0", "#A0A0A0",
"#E00000", "#00E000", "#0000E0", "#E0E000", "#E000E0", "#00E0E0", "#E0E0E0")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#924900", "#490092", "#920000", "#24ff24", "#000000", "#db6d00", "#006ddb")

### Functions ###

#create multiplot from list of plots
#source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#Distance fuctions for heatmap (For dist= argument in heatmap.2)
dist.pear <- function(x) as.dist(1-cor(t(x))) #for heatap.2's distfun parameter
dist.spear <- function(x) as.dist(1-cor(t(x), method="spearman")) #for heatmap.2's distfun parameter


generateHeatmap <- function(mytable, heat_out_pdf=NULL, cellValMat=NULL, width=6, height=6){
  .e <- environment()

  if(is.null(cellValMat)){
    p = heatmap.2(as.matrix(mytable),col=rev(heat.colors(256)), dendrogram='none', Rowv=F, Colv=F,
                  na.color="black", scale="none", density.info="none",
                  trace="none", cexCol=0.8, cexRow=0.8, margins=c(5,8), adjCol=c(NA,0.5), lhei = c(2,6))
  } else {
    p = heatmap.2(as.matrix(mytable),col=rev(heat.colors(256)), dendrogram='none', Rowv=F, Colv=F,
                  na.color="black", scale="none", density.info="none", cellnote=cellValMat, notecol="black",
                  trace="none", cexCol=0.8, cexRow=0.8, margins=c(5,8), adjCol=c(NA,0.5), lhei = c(2,6))
  }

  # p = heatmap(mytable, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
  if(!is.null(heat_out_pdf)){
    pdf(heat_out_pdf, width=width, height=height)
    p
    dev.off()
  } else {
    print(p)
  }
  return(p)
}
