last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

start <- 'both/'
dirs <- list.files(start)

ab.counts <- c()
stab <- c()
binding <- c()
name <- c()
chain <- c()

for(i in dirs) {
  dat<-read.table(paste(start, i, "/ancestral_comparisons.txt", sep=''), sep='\t', header=T, stringsAsFactors=F);
  final.letters <- sapply(dat$name, last.letter)

  ab.counts <- append(ab.counts, dat$count);
  stab <- append(stab, dat$stability);
  binding <- append(binding, dat$interaction);
  name <- append(name, dat$name)
  chain <- append(chain, final.letters)
}

all.data <- data.frame(names=name, ab.counts=ab.counts, binding=binding, stab=stab, chain=chain)
