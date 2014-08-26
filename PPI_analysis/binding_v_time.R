require(splines)
require(quantreg)
require(latticeExtra)

rm(list = ls())
this.chain = "C"
alpha = 0.318/2

mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)

cbbPalette <- c('WT-Ancestral' = "#000000", 'NonB-Ancestral' = mycols[1], 'LowS-Ancestral' = mycols[4])

##Determine the last letter because that's where the chain is encoded
last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

##Function that grabs the data from the files in the specified folder
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

##Get all of the data
survival.data.WT <- get.data('~/Sandbox/marcotte/complex_divergence_simul/data/WT_data/', this.chain)
survival.data.UnB <- get.data('~/Sandbox/marcotte/complex_divergence_simul/data/UnB_data/', this.chain)
survival.data.UnS <- get.data('~/Sandbox/marcotte/complex_divergence_simul/data/UnS_data/', this.chain)

##Put the data in a usable format for later
plot.data <- data.frame(x=c(rep(survival.data.WT$ab.count, 2), 
                            rep(survival.data.UnB$ab.count, 2), 
                            rep(survival.data.UnS$ab.count, 2)), 
                        y=c(survival.data.WT$ev.binding, 
                            survival.data.WT$an.binding, 
                            survival.data.UnB$ev.binding, 
                            survival.data.UnB$an.binding, 
                            survival.data.UnS$ev.binding, survival.data.UnS$an.binding), 
                        id=c(rep('WT-Evolved', dim(survival.data.WT)[1]), 
                             rep('WT-Ancestral', dim(survival.data.WT)[1]), 
                             rep('NonB-Evolved', dim(survival.data.UnB)[1]), 
                             rep('NonB-Ancestral', dim(survival.data.UnB)[1]), 
                             rep('LowS-Evolved', dim(survival.data.UnS)[1]), 
                             rep('LowS-Ancestral', dim(survival.data.UnS)[1]))
)

plot.data <- plot.data[order(plot.data$x), ]

##Start of WT
plot.data.ancestral.wt <- plot.data[plot.data$id == 'WT-Ancestral', ]

degree.freedom <- 4

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.ancestral.wt)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.ancestral.wt)
  y.fit.ancestral <- X %*% fit$coef
  plot.data.ancestral.wt <- cbind(plot.data.ancestral.wt, y.fit.ancestral)
}

plot.data.evolved.wt <- plot.data[plot.data$id == 'WT-Evolved', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.evolved.wt)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.evolved.wt)
  y.fit.evolved <- X %*% fit$coef
  plot.data.evolved.wt <- cbind(plot.data.evolved.wt, y.fit.evolved)
}

##Start of Non-Bound
plot.data.ancestral.nb <- plot.data[plot.data$id == 'NonB-Ancestral', ]

degree.freedom <- 4

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.ancestral.nb)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.ancestral.nb)
  y.fit.ancestral <- X %*% fit$coef
  plot.data.ancestral.nb <- cbind(plot.data.ancestral.nb, y.fit.ancestral)
}

plot.data.evolved.nb <- plot.data[plot.data$id == 'NonB-Evolved', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.evolved.nb)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.evolved.nb)
  y.fit.evolved <- X %*% fit$coef
  plot.data.evolved.nb <- cbind(plot.data.evolved.nb, y.fit.evolved)
}

##Start of Low Stability
plot.data.ancestral.ls <- plot.data[plot.data$id == 'LowS-Ancestral', ]

degree.freedom <- 4

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.ancestral.ls)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.ancestral.ls)
  y.fit.ancestral <- X %*% fit$coef
  plot.data.ancestral.ls <- cbind(plot.data.ancestral.ls, y.fit.ancestral)
}

plot.data.evolved.ls <- plot.data[plot.data$id == 'LowS-Evolved', ]

X <- model.matrix(y ~ bs(x, df=degree.freedom), data=plot.data.evolved.ls)

for (tau in c(alpha, 0.5, 1-alpha)) {
  fit <- rq(y ~ bs(x, df=degree.freedom), tau=tau, data=plot.data.evolved.ls)
  y.fit.evolved <- X %*% fit$coef
  plot.data.evolved.ls <- cbind(plot.data.evolved.ls, y.fit.evolved)
}

##Assemble data together for plotting
plot.data <- data.frame(x=c(plot.data.ancestral.wt$x, 
                            plot.data.ancestral.nb$x,
                            plot.data.ancestral.ls$x),
                        y=c(plot.data.ancestral.wt$y, 
                            plot.data.ancestral.nb$y,
                            plot.data.ancestral.ls$y),
                        ysmooth=c(plot.data.ancestral.wt[, 5], 
                                  plot.data.ancestral.nb[, 5],
                                  plot.data.ancestral.ls[, 5]),
                        ymin=c(plot.data.ancestral.wt[, 4], 
                               plot.data.ancestral.nb[, 4],
                               plot.data.ancestral.ls[, 4]),
                        ymax=c(plot.data.ancestral.wt[, 6], 
                               plot.data.ancestral.nb[, 6],
                               plot.data.ancestral.ls[, 6]),
                        id=c(rep('WT-Ancestral', length(plot.data.ancestral.wt$x)), 
                             rep('NonB-Ancestral', length(plot.data.ancestral.nb$x)),
                             rep('LowS-Ancestral', length(plot.data.ancestral.ls$x))))

##Put in levels so the legend is in the right order
plot.data$id <- factor(plot.data$id, levels = c('WT-Ancestral', 'LowS-Ancestral', 'NonB-Ancestral'))

##Plot function 
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
  
  g <- ggplot(df, aes(x=x, y=y, color=id, fill=id)) + 
    #geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.2, color=NA)
    geom_point(alpha=0.2, size=1.5) +
    geom_line(aes(y=ysmooth), size=2) + 
    scale_colour_manual(values=cbbPalette) + 
    scale_fill_manual(values=cbbPalette)
  g <- g + theme(strip.background=element_blank())
  g <- g + ylab('Binding Energy')
  g <- g + xlab('Time (Mutations Attempted)')
  g <- g + scale_x_continuous(breaks=seq(0, 1000, 100), limits=c(0, 1000))
  g <- g + scale_y_continuous(breaks=seq(-16., 16., 4.), limits=c(-16., 12.))
  g <- g + theme(panel.border=element_blank(), axis.line=element_line())
  g <- g + theme(axis.title.x = element_text(size=32, vjust=-0.25))
  g <- g + theme(axis.text.x = element_text(size=24, vjust=1.3))
  g <- g + theme(axis.title.y = element_text(size=32, vjust=-0.5))
  g <- g + theme(axis.text.y = element_text(size=24, hjust=1))
  g <- g + theme(axis.line = element_line(colour = 'black', size = 1.5))
  g <- g + theme(axis.ticks = element_line(colour = 'black', size = 1.5))
  g <- g + theme(plot.margin=unit(c(0.25, 0.25, 0.75, -0.5), "lines"))
  g <- g + theme(axis.ticks.margin = unit(0.25, "cm"))
  g <- g + theme(legend.position = c(0.25, 0.9),
                 legend.title=element_blank(), 
                 legend.key = element_blank(), 
                 legend.text=element_text(size=24),
                 legend.key.size = unit(1, "cm"))
  
  ggsave(g, file=paste('~/Sandbox/marcotte/complex_divergence_simul/figures/binding_v_time_', this.chain, '.pdf', sep=''), width=10, height=10)
  return(g)
}

survival.lines(plot.data)