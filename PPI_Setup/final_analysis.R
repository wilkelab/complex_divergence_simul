last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

start <- 'B2_ancest_data/'
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
  dat<-read.table(paste(start, i, sep=''), sep='\t', header=T, stringsAsFactors=F);
  final.letters <- sapply(dat$name, last.letter)
  count <- count + 1

  ab.count <- append(ab.count, dat$count);
  stab <- append(stab, dat$stability);
  ev.binding <- append(ev.binding, dat$evolved_interaction);
  an.binding <- append(an.binding, dat$ancestral_interaction);
  name <- append(name, dat$name)
  chain <- append(chain, final.letters)
  replicate <- append(replicate, rep(factor(count), length(final.letters)))
}

all.data.B2 <- data.frame(names=name, ab.count=ab.count, ev.binding=ev.binding, an.binding=an.binding, stab=stab, chain=chain, replicate=replicate, set=rep('HS', length(replicate)))
