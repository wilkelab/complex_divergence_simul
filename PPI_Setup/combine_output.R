start <- 'both/'
dirs <- list.files(start)

ab.counts <- c()
rel.counts <- c()
stab1 <- c()
stab2 <- c()
binding <- c()
name <- c()

for(i in dirs) {
  dat<-read.table(paste(start, i, "/add_counts.txt", sep=''), sep='\t', header=T, stringsAsFactors=F); 
  ab.counts <- append(ab.counts, dat$absolute_count); 
  rel.counts <- append(rel.counts, dat$relative_count);
  stab1 <- append(stab1, dat$stab1); 
  stab2 <- append(stab2, dat$stab2); 
  binding<-append(binding, dat$binding); 
  name <- append(name, dat$mutant)
  print(dim(dat))
}

print(length(name))
print(length(ab.counts))
print(length(rel.counts))
print(length(binding))
print(length(stab1))
print(length(stab2))

all.data.evolved <- data.frame(names=name, ab.counts=ab.counts, rel.counts=rel.counts, binding=binding, stab1=stab1, stab2=stab2)
