require(survival)
require(latticeExtra)

rm(list = ls())

mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)
mycols <- c("#000000",mycols[1], mycols[4])


survival.value = -7.5
this.chain = "C"

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
    
    cutoff.divergence <- 1 - min(dat$identity[which(dat$ancestral_interaction < survival.value)])
    cutoff.divergence[is.infinite(cutoff.divergence) | is.na(cutoff.divergence)] <- max(1 - dat$identity)
    
    survival.divergence <- append(survival.divergence, cutoff.divergence)
    status.divergence <- append(status.divergence, as.numeric(!cutoff.divergence == max(1 - dat$identity)))
  }
  
  tmp.survival.data <- data.frame(survival.divergence=survival.divergence,
                                  status.divergence=status.divergence
  )
  
  return(tmp.survival.data)
}

##Get all of the data
survival.data.WT <- get.data('~/Sandbox/complex_divergence_simul/data/WT_data/', this.chain)
survival.data.UnB <- get.data('~/Sandbox/complex_divergence_simul/data/UnB_data/', this.chain)
survival.data.UnS <- get.data('~/Sandbox/complex_divergence_simul/data/UnS_data/', this.chain)

print(mean(survival.data.UnB$survival.divergence))

survival.data <- data.frame(time=c(survival.data.WT$survival.divergence, 
                                   survival.data.UnB$survival.divergence,
                                   survival.data.UnS$survival.divergence), 
                            status=c(survival.data.WT$status.divergence, 
                                     survival.data.UnB$status.divergence,
                                     survival.data.UnS$status.divergence), 
                            replicate=c(rep('WT-Survival', length(survival.data.WT$status.divergence)), 
                                        rep('NonB-Survival', length(survival.data.UnB$status.divergence)),
                                        rep('LowS-Survival', length(survival.data.UnS$status.divergence)))
)

fit = survfit(Surv(time,status)~replicate, data=survival.data)
print(survdiff(Surv(time,status)~replicate, data=survival.data))

pdf(paste('~/Sandbox/complex_divergence_simul/figures/survival_v_divergence_max_', this.chain, '.pdf', sep=''), height=11, width=12)
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

legend(0, 0.15, c('Wild Type', 'Low Stability', 'Non-Bound'), col=c(mycols[1], mycols[3], mycols[2]), lty=1, cex=2, lwd=2.5, bty = "n")

dev.off()