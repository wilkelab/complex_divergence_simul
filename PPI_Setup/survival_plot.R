require(survival)

rm(list = ls())
survival.value = -7
this.chain = "A"

last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

get.data <- function(this.folder, which.chain) {
  start <- this.folder
  dirs <- list.files(start)
  
  survival.count <- c()
  survival.divergence <- c()
  status.divergence <- c()
  status.count <- c()
  
  for(i in dirs) {
    dat <- read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
    final.letters <- sapply(dat$name, last.letter)
    dat <- dat[final.letters == which.chain, ] 
    
    cutoff.count <- min(dat$count[which(dat$ancestral_interaction > survival.value)])
    cutoff.divergence <- 1 - min(dat$identity[which(dat$ancestral_interaction > survival.value)])
    
    cutoff.count[is.infinite(cutoff.count) | is.na(cutoff.count)] <- max(dat$count)
    cutoff.divergence[is.infinite(cutoff.divergence) | is.na(cutoff.divergence)] <- max(1 - dat$identity)
    
    survival.divergence <- append(survival.divergence, cutoff.divergence)
    status.divergence <- append(status.divergence, as.numeric(!cutoff.divergence == max(1 - dat$identity)))
    
    survival.count <- append(survival.count, cutoff.count)
    status.count <- append(status.count, as.numeric(!cutoff.count == max(dat$count)))
  }
  
  tmp.survival.data <- data.frame(survival.count=survival.count,
                                  survival.divergence=survival.divergence,
                                  status.count=status.count, 
                                  status.divergence=status.divergence
  )
  
  return(tmp.survival.data)
}

survival.data.WT <- get.data('~/Desktop/WT_data/', this.chain)
survival.data.UnB <- get.data('~/Desktop/UnS_data/', this.chain)

survival.data <- data.frame(time=c(survival.data.WT$survival.count, survival.data.UnB$survival.count), 
                            status=c(survival.data.WT$status.count, survival.data.UnB$status.count), 
                            replicate=c(rep('WT', length(survival.data.WT$status.count)), rep('UnS', length(survival.data.UnB$status.count)))
)

#survival.data <- data.frame(time=c(survival.data.UnB$survival.count), 
#                            status=c(survival.data.UnB$status.count), 
#                            replicate=c(rep('UnB', length(survival.data.UnB$status.count)))
#)

fit = survfit(Surv(time,status)~replicate, data=survival.data)
print(survdiff(Surv(time,status)~replicate, data=survival.data))

pdf(paste('~/Desktop/survival_', this.chain, '.pdf', sep=''), height=11, width=12)
par(mar=c(5,5,1,2)+0.1)
par(mgp=c(3, 1, 0))
par(family = 'Helvetica')

plot(fit, xlab="Time (Mutations Attempted)", 
     ylab="Survival Probability", 
     col=c('red', 'blue'), 
     cex.lab=1.75,
     mark=19,
     axes=F,
     xlim=c(0,1000),
     lwd=1.5
)

axis( 1, 
      cex.axis=1.75,
      at = seq(0, 1000, 100),
      lwd=1.5)
axis( 2, 
      cex.axis=1.75,
      #at = seq(0, 1, 0.1),
      lwd=1.5)

legend( 0.1, 0.1, c('UnB', 'WT'), col=c('red', 'blue'), lty=1, cex=1.4, lwd=2)

dev.off()