rm(list = ls())

last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

start <- 'B4_ancest_data/'
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
  for(i in 1:1000){ if(i > min(which(dat$ancestral_interaction > 0))) survival <- append(survival, 1) else survival <- append(survival, 0)}

}

all.data.B4 <- data.frame(names=name, ab.count=ab.count, ev.binding=ev.binding, an.binding=an.binding, stab=stab, chain=chain, replicate=replicate, set=rep('LS', length(replicate)))

survival.data.B4 <- data.frame(x=1:1000, y=rowMeans(as.data.frame(split(survival, survival.replicate))))