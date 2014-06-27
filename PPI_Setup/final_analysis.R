last.letter <- function(this.string) {tmp.length <- nchar(this.string); substring(this.string, tmp.length, tmp.length)}

dat <- read.table('ancestral_comparisons.txt', sep='\t', header=T, , stringsAsFactors=F)

final.letters <- sapply(dat$name, last.letter)

as <- final.letters == "A"
cs <- final.letters == "C"

print(cor.test(dat[as,]$count, dat[as,]$interaction))
print(cor.test(dat[cs,]$count, dat[cs,]$interaction))
