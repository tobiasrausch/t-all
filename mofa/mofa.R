library(DESeq2)
library(vsn)
library(MOFA2)

args=commandArgs(trailingOnly=TRUE)

ntop = 2500

# Samples
s = read.table(args[1], header=T)
colnames(s)[1] = "sample"
rownames(s) = s$sample
s$inirel = factor(s$inirel, levels=c("INI", "REL"))
s$reltype = factor(s$reltype, levels=c(1, 2))
s$patient = factor(s$patient)

# RNA 
rna = read.table(args[2], comment.char='$', quote='$', sep="\t", header=T)
rownames(rna) = rna$gene
geneinfo = rna[,1:5]
rna = rna[,rownames(s)]
dds = DESeqDataSetFromMatrix(countData = rna, colData = s, design = ~ patient + inirel)
vsd = vst(dds, blind=F)
rv = rowVars(assay(vsd))
selectIDX = order(rv, decreasing=T)[seq_len(min(ntop, length(rv)))]
rna = assay(vsd)[selectIDX,]
print(head(rna))

# ATAC
atac = read.table(args[3], comment.char='$', quote='$', sep="\t", header=T)
colnames(atac)[4] = "peak"
rownames(atac) = atac$peak
peakinfo = atac[,1:4]
atac = atac[,rownames(s)]
dds = DESeqDataSetFromMatrix(countData = atac, colData = s, design = ~ patient + inirel)
vsd = vst(dds, blind=F)
rv = rowVars(assay(vsd))
selectIDX = order(rv, decreasing=T)[seq_len(min(ntop, length(rv)))]
atac = assay(vsd)[selectIDX,]
print(head(atac))

# MOFA object
multiomic = list()
multiomic$atac = as.matrix(atac)
multiomic$rna = as.matrix(rna)
print(lapply(multiomic,dim))
MOFAobject = create_mofa(multiomic)

# Overview
png("data.overview.png")
plot_data_overview(MOFAobject)
dev.off()

# Options
data_opts = get_default_data_options(MOFAobject)
model_opts = get_default_model_options(MOFAobject)
model_opts$num_factors = 9
train_opts = get_default_training_options(MOFAobject)
train_opts$convergence_mode = "slow"
train_opts$seed = 42

# Train the model
MOFAobject = prepare_mofa(MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAobject = run_mofa(MOFAobject)
samples_metadata(MOFAobject) = s
print(MOFAobject)

# Factor correlation
png("factor.correlation.png")
plot_factor_cor(MOFAobject)
dev.off()

# Variance explained by data modality
png("variance.explained.png")
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
dev.off()

# Variance explained by factor
png("variance.by.factor.png")
plot_variance_explained(MOFAobject, max_r2=15)
dev.off()

# Covariates
png("covariates.png")
correlate_factors_with_covariates(MOFAobject, covariates = c("patient", "fusion", "sex", "reltype", "inirel"), plot="log_pval")
dev.off()

# Factor values by inirel
png("inirel.png")
plot_factor(MOFAobject, factors = 1, color_by = "inirel", add_violin = TRUE,  dodge = TRUE)
dev.off()

# Plot gene weights
png("factor3.rna.png")
plot_weights(MOFAobject, view = "rna", factor = 3, nfeatures = 10)
dev.off()

# Heatmap
png("heatmap.factor1.png")
plot_data_heatmap(MOFAobject, view = "rna", factor = 1, features = 25, denoise = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()
