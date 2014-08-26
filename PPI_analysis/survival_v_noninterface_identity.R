require(survival)
require(latticeExtra)

rm(list = ls())

mycols <- dput(ggplot2like(n = 5, h.start = 0, l = 65)$superpose.line$col)
mycols <- c("#000000",mycols[1], mycols[4])


survival.value = -7.5
this.chain = "A"

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
    
    divergence <- 1 - dat$non_interface_identity
    survived <- dat$ancestral_interaction <= survival.value
    cutoff.divergence <- max(divergence[survived])
    cutoff.divergence[is.infinite(cutoff.divergence) | is.na(cutoff.divergence)] <- max(divergence)
    
    survival.divergence <- append(survival.divergence, cutoff.divergence)
    status.divergence <- append(status.divergence, as.numeric(!cutoff.divergence == max(divergence)))
  }
  
  tmp.survival.data <- data.frame(survival.divergence=survival.divergence,
                                  status.divergence=status.divergence
  )
  
  return(tmp.survival.data)
}

##Get all of the data
survival.data.WT <- get.data('~/Sandbox/marcotte/complex_divergence_simul/data/WT_data/', this.chain)
survival.data.UnB <- get.data('~/Sandbox/marcotte/complex_divergence_simul/data/UnB_data/', this.chain)
survival.data.UnS <- get.data('~/Sandbox/marcotte/complex_divergence_simul/data/UnS_data/', this.chain)

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

pdf(paste('~/Sandbox/marcotte/complex_divergence_simul/figures/survival_v_noninterface_divergence_', this.chain, '.pdf', sep=''), height=11, width=12)
par(mar=c(4,6,0,2))
par(mgp=c(3, 1.2, 0))

plot(fit,
     col=c(mycols[3], mycols[2], mycols[1]), 
     cex.lab=2,
     mark=19,
     axes=F,
     xlim=c(0,1.0),
     lwd=6,
     cex=2.5,
     font=2
)

axis( 1, 
      cex.axis=2.25,
      lwd=5,
      at=seq(0, 1.0, by=0.2),
      labels=seq(0, 100, by=20),
      line=-1)
axis( 2, 
      cex.axis=2.25,
      lwd=5,
      at=seq(0, 1.0, by=0.2),
      labels=seq(0, 100, by=20),
      las=1)

mtext("Divergence (100 - % Amino Acid Identity)", side=1, line = 2.5, cex=3, font=2)
mtext("% Binding Ancestor", side=2, line = 3.75, cex=3, font=2)

legend(0, 0.15, c('Wild Type', 'Low Stability', 'Non-Bound'), col=c(mycols[1], mycols[3], mycols[2]), lty=1, cex=2, lwd=10, bty = "n")

dev.off()