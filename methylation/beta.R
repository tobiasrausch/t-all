library(minfi)

path = "../../fastq/methylome/"
targets = read.table("../../fastq/methylome/samples.tsv", header=T)
print(summary(targets))
targets$Basename = file.path(path, targets$Basename)

rgset = read.metharray.exp(targets = targets, verbose=T)
mset = preprocessIllumina(rgset)
mset = mapToGenome(mset)
#plot(as.matrix(getQC(mset))
betaval = data.frame(getBeta(mset, type = "Illumina"))
betaval$ILMNID = rownames(betaval)
write.table(betaval, file="methylation.beta.tsv", sep="\t", quote=FALSE, row.names=FALSE)
