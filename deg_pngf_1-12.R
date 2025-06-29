# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE70757", GSEMatrix = TRUE, AnnotGPL = FALSE)
if (length(gset) > 1L) {
  idx <- grep("GPL10787",
    fixed = TRUE,
    attr(gset, "names")
  )
} else {
  idx <- 1L
}
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000XXXXXXXXXXXX111XXX"
sml <- strsplit(gsms, split = "", fixed = TRUE)[[1L]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[, sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5L] > 100L) ||
  (qx[6L] - qx[1L] > 50L && qx[2L] > 0L)
if (LogC) {
  ex[which(ex <= 0L)] <- NaN
  exprs(gset) <- log2(ex)
}

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("CTRL", "CASE"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~ group + 0L, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design) # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[2L], groups[1L], sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)

tT <- subset(tT, select = c("ID", "adj.P.Val", "P.Value", "t", "logFC", "GENE", "GENE_SYMBOL"))
write.table(tT, file = file.path("extdata", "pngf_1m-12m.csv"), row.names = FALSE, sep = ",")

writable <- subset(tT, select = c("GENE_SYMBOL", "logFC", "P.Value"))
writable <- writable[writable[["GENE_SYMBOL"]] != "", ]
names(writable) <- c("gene_name", "logFC", "P.Value")
writable <- writable[order(writable$gene_name, -writable$logFC), ]
writable <- writable[!duplicated(writable$gene_name), ]
write.table(writable, file = file.path("extdata", "pngf_1m-12m_clean.csv"), row.names = FALSE, sep = ",")
