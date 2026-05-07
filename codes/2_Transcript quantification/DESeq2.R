library(DESeq2)
mycounts <- read.table("Ppup.CA.SGD.txt", header = T, row.names = "gene_id",sep="\t")
mycounts<-as.data.frame(lapply(mycounts,as.integer),row.names=rownames(mycounts))
head(mycounts)
condition <- data.frame(condition = factor(rep(c("control", 'treat'), each = 3), levels = c('control', 'treat')))
colData <- data.frame(row.names = colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'treat', 'control'))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'Ppup.CA.SGD.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)

res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.05),'sig'] <- 'none'

res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'Ppup.CA.SGD.DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)

res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

write.table(res1_up, file = 'Ppup.CA.SGD.DESeq2.up.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'Ppup.CA.SGD.DESeq2.DESeq2.down.txt', sep = '\t', col.names = NA, quote = FALSE)
