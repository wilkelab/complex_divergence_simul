require(stringr)

dat <- read.table('ancestral_comparisons.txt', sep='\t', header=T)

as <- str_sub(levels(dat$name), start=-1) == "A"
cs <- str_sub(levels(dat$name), start=-1) == "C"

print(cor.test(dat[as,]$count, dat[as,]$interaction))
print(cor.test(dat[cs,]$count, dat[cs,]$interaction))
