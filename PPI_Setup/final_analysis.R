rm(list = ls())
survival.value = -7

last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

get.data <- function(this.folder) {
  start <- this.folder
  dirs <- list.files(start)
  
  ab.count <- c()
  stab <- c()
  ev.binding <- c()
  an.binding <- c()
  name <- c()
  chain <- c()
  replicate <- c()
  survival <- c()
  survival.replicate <- c()
  count <- 0
  
  for(i in dirs) {
    dat<-read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    final.letters <- sapply(dat$name, last.letter)
    count <- count + 1
    
    ab.count <- append(ab.count, dat$count);
    stab <- append(stab, dat$stability);
    ev.binding <- append(ev.binding, dat$evolved_interaction);
    an.binding <- append(an.binding, dat$ancestral_interaction);
    name <- append(name, dat$name)
    chain <- append(chain, final.letters)
    replicate <- append(replicate, rep(as.character(count), length(final.letters)))
    survival.replicate <- append(survival.replicate, rep(count, 1000))
    for(i in 1:1000){ if(i > min(which(dat$ancestral_interaction > survival.value))) survival <- append(survival, 0) 
                      else survival <- append(survival, 1)}
    
  }
  
  tmp.data <- data.frame(names=name, ab.count=ab.count, ev.binding=ev.binding, an.binding=an.binding, stab=stab, 
                         chain=chain, replicate=replicate, set=rep('LS', length(replicate)))
  
  tmp.survival.data <- data.frame(x=1:1000, y=rowMeans(as.data.frame(split(survival, survival.replicate))))
  
  return(tmp.survival.data)
  
}

survival.data.B1 <- get.data('B1_ancest_data/')
survival.data.B2 <- get.data('B2_ancest_data/')
survival.data.B3 <- get.data('B3_ancest_data/')
survival.data.B4 <- get.data('B4_ancest_data/')
survival.data.UnB <- get.data('UnB_data/')

survival.data <- data.frame(x=rep(1:1000, 5), y=c(survival.data.B1$y, survival.data.B2$y, survival.data.B3$y, survival.data.B4$y, survival.data.UnB$y), 
                            replicate=c(rep('b1',1000), rep('b2', 1000), rep('b3', 1000), rep('b4', 1000), rep('UnB', 1000)))


survival.lines <- function(df) {
  require(ggplot2)
  require(grid)
  
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
  
  g <- ggplot(df, aes(x=x, y=y, color=replicate)) + geom_line(size=2.5)
  g <- g + theme(strip.background=element_blank())
  g <- g + ylab('Ancestral Binding')
  g <- g + xlab('Time (Mutations Attempted)')
  g <- g + scale_x_continuous(breaks=seq(0, 1000, 25), limits=c(0, 150))
  #  g <- g + scale_y_continuous(breaks=seq(500, 1000, 100), limits=c(550, 900))
  g <- g + theme(panel.border=element_blank(), axis.line=element_line())
  g <- g + theme(axis.title.x = element_text(size=24, vjust=-1))
  g <- g + theme(axis.text.x = element_text(size=24))
  g <- g + theme(axis.title.y = element_text(size=24, vjust=2))
  g <- g + theme(axis.text.y = element_text(size=24))
  g <- g + theme(axis.line = element_line(colour = 'black', size = 1))
  g <- g + theme(axis.ticks = element_line(colour = 'black', size = 1))
  g <- g + theme(plot.margin=unit(c(1.5, 1.5, 1.5, 1.5), "lines"))
  g <- g + theme(axis.ticks.margin = unit(0.25, "cm"))
  g <- g + theme(legend.position = "none")
  
  ggsave(g, file='~/Desktop/survival.pdf', width=10, height=10)
  return(g)
}

survival.lines(survival.data)