library(DESeq2)
library(vsn)
library(ggplot2)

# Rscript diffRNA.R gene.count sample.info [all|1|2]

# Load count matrix and sample info
args=commandArgs(trailingOnly=TRUE)
x = read.table(args[1], comment.char='$', quote='$', sep="\t", header=T)
s = read.table(args[2], header=T)
atype = args[3]
rownames(x) = x$gene
rownames(s) = s$id
s$inirel = factor(s$inirel, levels=c("INI", "REL"))
s$reltype = factor(s$reltype, levels=c(1, 2))
geneinfo = x[,1:5]
x = x[,6:ncol(x)]

# Filter unexpressed genes
print("Filter unexpressed genes (not absolutely necessary)")
print(nrow(x))
x = x[rowSums(x) >= 10,]
print(nrow(x))

# Subset if necessary
if (atype == "1") {
 s = s[s$reltype == 1,]
 x = x[,colnames(x) %in% s$id]
} else if (atype == "2") {
 s = s[s$reltype == 2,]
 x = x[,colnames(x) %in% s$id]
}

# Build DESeq object
dds = DESeqDataSetFromMatrix(countData = x, colData = s, design = ~ patient + inirel)
print("DESeq2 version")
print(metadata(dds)[["version"]])

# Variance stabilizing transformation
vsd = vst(dds, blind=F)

# Plot PCA
pcaData = plotPCA(vsd, intgroup=c("fusion", "reltype", "patient", "inirel"), returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))
png(paste0("pca.rna.type_", atype, ".png"), width=800, height=800)
p=ggplot(data=pcaData, aes(x = PC1, y = PC2, color=reltype, fill=patient, shape=inirel)) + geom_point(size=3, stroke=1.5)
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))
p=p+scale_fill_brewer(palette="Paired")
p=p+scale_shape_manual(labels=c("Initial", "Relapse"), values=c(24, 25))
p=p+scale_color_manual(labels=c("Type1", "Type2"), values=c("black", "grey"))
p=p+guides(fill = guide_legend(override.aes = list(shape = 21)))
p=p+guides(color = guide_legend(override.aes = list(shape = 21)))
p=p+facet_wrap(~fusion)
p=p+labs(color="T-ALL Type", fill="Patient", shape="Time")
p
dev.off()

# Differential expression
dds = DESeq(dds)
print(resultsNames(dds))  # List the coefficients

# Results
res = results(dds, alpha=0.05)
print(summary(res))

# QC
resLFC = lfcShrink(dds, coef="inirel_REL_vs_INI", type="apeglm")
png(paste0("ma.rna.type_", atype, ".png"))
plotMA(resLFC, ylim=c(-2,2))
dev.off()

# Significant results
resSig = subset(res, padj<0.05)
resSig = resSig[order(resSig$pvalue),]

# Write results
resSig = as.data.frame(resSig)
resSig$gene = rownames(resSig)
resSig = merge(geneinfo, resSig, by=c("gene"))
write.table(resSig, file=paste0("res.sig.rna.type_", atype, ".tsv"), col.names=T, row.names=F, quote=F, sep="\t")
