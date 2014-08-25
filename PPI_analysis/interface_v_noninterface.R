rm(list = ls())
require(ggplot2)
require(grid)
require(latticeExtra)

this.chain = "C"

last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('Wild Type' = "#000000", 'Non-Bound' = mycols[1], 'Low Stability' = mycols[4])

last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}


get.data <- function(this.folder, which.chain) {
  start <- this.folder
  dirs <- list.files(start)
  
  survival.divergence <- c()
  status.divergence <- c()
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    final.letters <- sapply(dat$name, last.letter)
    dat <- dat[final.letters == which.chain, ] 
    
    divergence <- 1 - dat$identity
    survived <- dat$ancestral_interaction <= survival.value
    cutoff.surv <- max(divergence[survived])
    cutoff.divergence <- min(divergence[!survived & divergence >= cutoff.surv])
    cutoff.divergence[is.infinite(cutoff.divergence) | is.na(cutoff.divergence)] <- max(divergence)
    
    survival.divergence <- append(survival.divergence, cutoff.divergence)
    status.divergence <- append(status.divergence, as.numeric(!cutoff.divergence == max(divergence)))
  }
  
  tmp.survival.data <- data.frame(survival.divergence=survival.divergence,
                                  status.divergence=status.divergence
  )
  
  return(tmp.survival.data)
}

interface.plot <- function(df) {
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
  
  g <- ggplot(df, aes(x=int.frac, fill=id, color=id)) + 
    geom_density(alpha=0.2) + 
    scale_colour_manual(values=cbbPalette) + 
    scale_fill_manual(values=cbbPalette)
  g <- g + theme(strip.background=element_blank())
  g <- g + ylab('Probability Density')
  g <- g + xlab('Odds of accepting Non-interface to Interface Mutations')
  g <- g + scale_x_continuous(breaks=seq(0, 7, 1), limits=c(0, 7))
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
  
  ggsave(g, file=paste('~/Sandbox/complex_divergence_simul/figures/interface_odds.pdf', sep=''), width=12, height=10)
  return(g)
}

##Get all of the data
survival.data.WT <- get.data('~/Sandbox/complex_divergence_simul/data/WT_data/', this.chain)
survival.data.UnB <- get.data('~/Sandbox/complex_divergence_simul/data/UnB_data/', this.chain)
survival.data.UnS <- get.data('~/Sandbox/complex_divergence_simul/data/UnS_data/', this.chain)

plot.data <- data.frame(id=c(rep('Wild Type', length(WT$interface)), 
                             rep('Non-Bound', length(UnB$interface)), 
                             rep('Low Stability', length(UnS$interface))), 
                        int.frac=c(WT$non_interface/WT$interface, 
                                   UnB$non_interface/UnB$interface, 
                                   UnS$non_interface/UnS$interface))

##Put in levels so the legend is in the right order
plot.data$id <- factor(plot.data$id, levels = c('Wild Type', 'Low Stability', 'Non-Bound'))


interface.plot(plot.data)