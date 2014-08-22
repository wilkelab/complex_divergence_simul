require(survival)
require(latticeExtra)

rm(list = ls())

mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)
mycols <- c("#000000", mycols[1], mycols[4])


survival.value = -7.5
this.chain = "A"
treatment = 'WT'

last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

get.data <- function(this.folder, which.chain) {
  start <- this.folder
  dirs <- list.files(start)
  
  survival.divergence <- c()
  status.divergence <- c()
  survival.divergence.interface <- c()
  status.divergence.interface <- c()
  survival.divergence.noninterface <- c()
  status.divergence.noninterface <- c()
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    final.letters <- sapply(dat$name, last.letter)
    dat <- dat[final.letters == which.chain, ] 
    
    divergence <- 1 - dat$identity
    survived <- which(dat$ancestral_interaction <= survival.value)
    cutoff.divergence <- max(divergence[survived])
    cutoff.divergence[is.infinite(cutoff.divergence) | is.na(cutoff.divergence)] <- max(divergence)
    survival.divergence <- append(survival.divergence, cutoff.divergence)
    status.divergence <- append(status.divergence, as.numeric(!cutoff.divergence == max(divergence)))
    
    divergence.interface <- 1 - dat$interface_identity
    survived <- which(dat$ancestral_interaction <= survival.value)
    cutoff.divergence.interface <- max(divergence.interface[survived])
    cutoff.divergence.interface[is.infinite(cutoff.divergence.interface) | is.na(cutoff.divergence.interface)] <- max(divergence.interface)
    survival.divergence.interface <- append(survival.divergence.interface, cutoff.divergence.interface)
    status.divergence.interface <- append(status.divergence.interface, as.numeric(!cutoff.divergence.interface == max(divergence.interface)))
    
    divergence.noninterface <- 1 - dat$non_interface_identity
    survived <- which(dat$ancestral_interaction <= survival.value)
    cutoff.divergence.noninterface <- max(divergence.noninterface[survived])
    cutoff.divergence.noninterface[is.infinite(cutoff.divergence.noninterface) | is.na(cutoff.divergence.noninterface)] <- max(divergence.noninterface)
    survival.divergence.noninterface <- append(survival.divergence.noninterface, cutoff.divergence.noninterface)
    status.divergence.noninterface <- append(status.divergence.noninterface, as.numeric(!cutoff.divergence.noninterface == max(divergence.noninterface)))
  }
  
  tmp.survival.data <- data.frame(survival.divergence=survival.divergence,
                                  status.divergence=status.divergence,
                                  survival.divergence.interface=survival.divergence.interface,
                                  status.divergence.interface=status.divergence.interface,
                                  survival.divergence.noninterface=survival.divergence.noninterface,
                                  status.divergence.noninterface=status.divergence.noninterface
  )
  
  return(tmp.survival.data)
}

##Get all of the data
survival.data.WT <- get.data(paste('~/Sandbox/complex_divergence_simul/data/', treatment, '_data/', sep=''), this.chain)

survival.data <- data.frame(time=c(survival.data.WT$survival.divergence, 
                                   survival.data.WT$survival.divergence.interface,
                                   survival.data.WT$survival.divergence.noninterface), 
                            status=c(survival.data.WT$status.divergence, 
                                     survival.data.WT$status.divergence.interface,
                                     survival.data.WT$status.divergence.noninterface), 
                            replicate=c(rep('Total Identity', length(survival.data.WT$status.divergence)), 
                                        rep('Interface Identity', length(survival.data.WT$status.divergence.interface)),
                                        rep('Non-Interface Identity', length(survival.data.WT$status.divergence.noninterface)))
)

fit = survfit(Surv(time,status)~replicate, data=survival.data)
print(survdiff(Surv(time,status)~replicate, data=survival.data))

pdf(paste('~/Sandbox/complex_divergence_simul/figures/survival_v_divergence_sites_', treatment, '_', this.chain, '.pdf', sep=''), height=11, width=12)
par(mar=c(5,5,1,2)+0.1)
par(mgp=c(3, 1, 0))
par(family = 'Helvetica')

plot(fit, xlab="Divergence (1 - Identity)", 
     ylab="Survival Probability", 
     col=c(mycols[3], mycols[2], mycols[1]), 
     cex.lab=2,
     mark=19,
     axes=F,
     xlim=c(0,1),
     lwd=2.5
)

axis( 1, 
      cex.axis=2,
      lwd=2)
axis( 2, 
      cex.axis=2,
      lwd=2)

legend(0, 0.15, c('Total', 'Interface', 'Non-Interface'), col=c(mycols[1], mycols[3], mycols[2]), lty=1, cex=2, lwd=2.5, bty = "n")

dev.off()