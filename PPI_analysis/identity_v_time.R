rm(list = ls())
require(ggplot2)
require(grid)
require(latticeExtra)

this.chain = "C"

last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('Wild Type' = "#000000", 'Non-Bound' = mycols[1], 'Low Stability' = mycols[4])

get.data <- function(this.folder, which.chain) {
  start <- this.folder
  dirs <- list.files(start)
  
  ab.count <- c()
  stab <- c()
  ev.binding <- c()
  an.binding <- c()
  name <- c()
  chain <- c()
  replicate <- c()
  identity <- c()
  count <- 0
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    final.letters <- sapply(dat$name, last.letter)
    dat <- dat[final.letters == which.chain, ] 
    final.letters <- final.letters[final.letters == which.chain]
    
    count <- count + 1
    
    identity <- append(identity, dat$identity)
    ab.count <- append(ab.count, dat$count);
    stab <- append(stab, dat$stability);
    ev.binding <- append(ev.binding, dat$evolved_interaction);
    an.binding <- append(an.binding, dat$ancestral_interaction);
    name <- append(name, dat$name)
    chain <- append(chain, final.letters)
    replicate <- append(replicate, rep(as.character(count), length(final.letters)))
  }
  
  tmp.data <- data.frame(names=name, identity=identity, ab.count=ab.count, ev.binding=ev.binding, an.binding=an.binding, stab=stab, 
                         chain=chain, replicate=replicate)
  
  return(tmp.data)
  
}

identity.plot <- function(df, name) {
  graphics.off()
  the_pointsize=18
  theme_set(theme_bw(base_size=the_pointsize))
  old_theme <- theme_update(panel.border=element_blank(),
                            axis.line=element_line(),
                            panel.grid.minor=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.background=element_blank(),
                            panel.border=element_blank(),
                            axis.line=element_line())
  
  g <- ggplot(df, aes(x=x, y=y, fill=set, color=set)) + 
    geom_point(alpha = 0.2, size=1.5) + 
    scale_colour_manual(values=cbbPalette) + 
    scale_fill_manual(values=cbbPalette) +
    stat_smooth(se = FALSE, lwd=1.5) 
  g <- g + theme(strip.background=element_blank())
  g <- g + ylab('Sequence Identity (%)')
  g <- g + xlab('Time (Mutations Attempted)')
  g <- g + scale_x_continuous(breaks=seq(0, 1000, 100), limits=c(0, 1000))
  g <- g + scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1))
  g <- g + theme(panel.border=element_blank(), axis.line=element_line())
  g <- g + theme(axis.title.x = element_text(size=24, vjust=-1))
  g <- g + theme(axis.text.x = element_text(size=24))
  g <- g + theme(axis.title.y = element_text(size=24, vjust=2))
  g <- g + theme(axis.text.y = element_text(size=24))
  g <- g + theme(axis.line = element_line(colour = 'black', size = 1))
  g <- g + theme(axis.ticks = element_line(colour = 'black', size = 1))
  g <- g + theme(plot.margin=unit(c(1.5, 1.5, 1.5, 1.5), "lines"))
  g <- g + theme(axis.ticks.margin = unit(0.25, "cm"))
  g <- g + theme(legend.position = c(0.75, 0.9),
                 legend.title=element_blank(), 
                 legend.key = element_blank(), 
                 legend.text=element_text(size=22),
                 legend.key.size = unit(1, "cm"))
  
  ggsave(g, file=paste('~/Sandbox/complex_divergence_simul/figures/identity_v_time_', name, '.pdf', sep=''), width=12, height=10)
  return(g)
}

WT <- get.data('~/Sandbox/complex_divergence_simul/data/WT_data/', this.chain)
UnS <- get.data('~/Sandbox/complex_divergence_simul/data/UnS_data/', this.chain)
UnB <- get.data('~/Sandbox/complex_divergence_simul/data/UnB_data/', this.chain)

plot.data <- data.frame(x=c(WT$ab.count, UnB$ab.count, UnS$ab.count),
                        y=c(WT$identity, UnB$identity, UnS$identity),
                        set=c(rep('Wild Type', length(WT$identity)), 
                              rep('Non-Bound', length(UnB$identity)), 
                              rep('Low Stability', length(UnS$identity))
                              )
                        )

##Put in levels so the legend is in the right order
plot.data$set <- factor(plot.data$set, levels = c('Wild Type', 'Low Stability', 'Non-Bound'))

identity.plot(plot.data, this.chain)
