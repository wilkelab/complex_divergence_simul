require(splines)
require(quantreg)

rm(list = ls())
this.chain = "A"
alpha = 0.05/2

last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

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
  count <- 0
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    final.letters <- sapply(dat$name, last.letter)
    dat <- dat[final.letters == which.chain, ] 
    final.letters <- final.letters[final.letters == which.chain]
    
    count <- count + 1
    
    ab.count <- append(ab.count, dat$count);
    stab <- append(stab, dat$stability);
    ev.binding <- append(ev.binding, dat$evolved_interaction);
    an.binding <- append(an.binding, dat$ancestral_interaction);
    name <- append(name, dat$name)
    chain <- append(chain, final.letters)
    replicate <- append(replicate, rep(as.character(count), length(final.letters)))
  }
  
  tmp.data <- data.frame(names=name, ab.count=ab.count, ev.binding=ev.binding, an.binding=an.binding, stab=stab, 
                         chain=chain, replicate=replicate, set=rep('WT', length(replicate)))
  
  return(tmp.data)
  
}

survival.data.WT <- get.data('~/Desktop/WT_data/', this.chain)

plot.data <- data.frame(x=rep(survival.data.WT$ab.count, 2), 
                        y=c(survival.data.WT$ev.binding, survival.data.WT$an.binding), 
                        id=c(rep('Evolved', dim(survival.data.WT)[1]), rep('Ancestral', dim(survival.data.WT)[1]))
)

plot.data <- plot.data[order(plot.data$x), ]

plot.data.ancestral <- plot.data[plot.data$id == 'Ancestral', ]

degree.freedom <- 3

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.ancestral)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.ancestral)
  y.fit.ancestral <- X %*% fit$coef
  plot.data.ancestral <- cbind(plot.data.ancestral, y.fit.ancestral)
}

plot.data.evolved <- plot.data[plot.data$id == 'Evolved', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.evolved)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.evolved)
  y.fit.evolved <- X %*% fit$coef
  plot.data.evolved <- cbind(plot.data.evolved, y.fit.evolved)
}

plot.data <- data.frame(x=c(plot.data.ancestral$x, plot.data.evolved$x),
                        y=c(plot.data.ancestral$y, plot.data.evolved$y),
                        ysmooth=c(plot.data.ancestral[, 5], plot.data.evolved[, 5]),
                        ymin=c(plot.data.ancestral[, 4], plot.data.evolved[, 4]),
                        ymax=c(plot.data.ancestral[, 6], plot.data.evolved[, 6]),
                        id=c(rep('Ancestral', length(plot.data.ancestral$x)), rep('Evolved', length(plot.data.evolved$x))))

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
  
  g <- ggplot(df, aes(x=x, y=y, color=id, fill=id)) + geom_point(size=1.4) + geom_line(aes(y=ysmooth), size=1.4) + geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.2, color=NA)
  g <- g + theme(strip.background=element_blank())
  g <- g + ylab('Binding Energy')
  g <- g + xlab('Time (Mutations Attempted)')
  g <- g + scale_x_continuous(breaks=seq(0, 1000, 100), limits=c(0, 1000))
  g <- g + scale_y_continuous(breaks=seq(-16., 0., 2.), limits=c(-16., 0.))
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
  
  ggsave(g, file=paste('~/Desktop/binding_v_time_', this.chain, '.pdf', sep=''), width=10, height=10)
  return(g)
}

survival.lines(plot.data)