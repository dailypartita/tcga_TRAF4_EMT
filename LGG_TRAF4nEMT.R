setwd("C:/Users/yangkaixin1/Desktop/tcga_coad/github_version")
library(survival)
library(limma)
library(gplots)
library(corrplot)

# input counts, and filter genes whose expression is zero in more than half the samples
# 读取COAD的RNAseq数据，去掉超过一半样本中不表达的基因
rna = as.matrix(read.table("LGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header=T, row.names=1, sep="\t"))
rna = rna[-1,]
rownames(rna) = sapply(rownames(rna), function(x) unlist(strsplit(x,"\\|"))[[1]])
rna = type.convert(rna)
x = t(apply(rna,1,as.numeric))
r = as.numeric(apply(x,1,function(i) sum(i == 0)))
remove = which(r > dim(rna)[2]*0.5)
rna = rna[-remove,]

# 获取肿瘤和对照样本的colindex，提取肿瘤组RNA表达数据
n_index <- which(substr(colnames(rna),14,14) == "1")
t_index <- which(substr(colnames(rna),14,14) == "0")
tumor <- rna[,t_index]
z_rna <- matrix(nrow=nrow(tumor), ncol=ncol(tumor))
colnames(z_rna) <- colnames(rna[,t_index])
rownames(z_rna) <- rownames(rna[,t_index])

# 计算所有基因的zscore
for(i in 1:nrow(tumor)){
  for(j in 1:ncol(tumor)){
    z_rna[i,j] <- (tumor[i, j] - mean(tumor[i,])) / sd(tumor[i,])
  }
}

# 提取EMT相关基因的zscore
emt_gene = c("TRAF4", "CDH1", "CDH2", "VIM", "FN1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "SNAI3", "TWIST1", "TWIST2")
emt_zrna = lapply(emt_gene, function(gene){
  gene_zrna = z_rna[which(rownames(z_rna) == gene),]
  return(gene_zrna)
})
emt_zrna = data.frame(setNames(emt_zrna, emt_gene))
emt_zrna = emt_zrna[order(emt_zrna$TRAF4),]

# 绘制热图
# x11()
data = data.matrix(t(emt_zrna))
data = ifelse(data > 4, 4, data)
data = ifelse(data < -4, -4, data)
pdf(file = "LGG_TRAF4_EMT_zscore_heatmap.pdf", width = 10, height = 10)
heatmap.2(data,
          Colv = NULL,
          Rowv = NULL,
          trace="none",
          density="none",
          col=bluered(40),
          breaks = seq(-1,1,0.05), 
          labCol = NA,
          dendrogram = "none" # remove dendrogram and clustering
          )
dev.off()

# 绘制基因表达相关矩阵
emt_zrna_cor = cor(emt_zrna)
pdf(file = "LGG_TRAF4_EMT_zscore_cor.pdf", width = 10, height = 10)
corrplot(emt_zrna_cor)
dev.off()

# read the Clinical file, in this case i transposed it to keep the clinical feature title as column name
# 读取临床数据，输出TRAF4 & EMT的zscore和临床数据
# clinical <- read.csv("COAD.merged_only_clinical_clin_format.txt",sep = "\t")
# clinical = t(clinical)
# write.csv(emt_zrna, "COAD_RNAseq_TRAF4nEMT_zscore.csv")
# write.csv(clinical,"COAD_clinical.csv")
