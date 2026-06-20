install.packages("WebGestaltR")
# 1. CARREGAR AS BIBLIOTECAS
library(tximport)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(WebGestaltR)
library(dplyr)

# 2. CONFIGURAR METADADOS E CAMINHOS
samples <- read.table("samples.tsv", header=TRUE, sep="\t")
tx2gene <- read.csv("reference/tx2gene.csv")

# onde estão os arquivos quant.sf
files <- file.path("salmon_results", "results", "salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

# 3. IMPORTAR DADOS DO SALMON
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# 4. CRIAR OBJETO DESEQ2 E FILTRAR RUÍDOS
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 5. RODAR ANÁLISE ESTATÍSTICA
dds <- DESeq(dds)

# 6. EXTRAIR RESULTADOS DAS COMPARAÇÕES
res_MEKi_vs_ctrl  <- results(dds, contrast=c("condition", "MEKi", "ctrl"))
res_FLT3i_vs_ctrl <- results(dds, contrast=c("condition", "FLT3i", "ctrl"))
res_MEKi_vs_FLT3i <- results(dds, contrast=c("condition", "MEKi", "FLT3i"))

# Criar pasta para salvar resultados caso não exista
if(!dir.exists("results")) dir.create("results")
if(!dir.exists("results/deseq2")) dir.create("results/deseq2", recursive = TRUE)

# Salvar tabelas em CSV
write.csv(as.data.frame(res_MEKi_vs_ctrl),  "results/deseq2/MEKi_vs_ctrl.csv")
write.csv(as.data.frame(res_FLT3i_vs_ctrl), "results/deseq2/FLT3i_vs_ctrl.csv")
write.csv(as.data.frame(res_MEKi_vs_FLT3i), "results/deseq2/MEKi_vs_FLT3i.csv")


# GERAÇÃO DOS GRÁFICOS
vsd <- vst(dds, blind=FALSE)

# GRAFICO 1: PCA
png("results/deseq2/PCA_GSE161218.png", width=800, height=600)
p_pca <- plotPCA(vsd, intgroup="condition") + theme_bw() + geom_point(size=4)
print(p_pca)
dev.off()

# GRAFICO 2: VOLCANO PLOT
df_volcano <- as.data.frame(res_MEKi_vs_ctrl)
png("results/deseq2/Volcano_MEKi_vs_ctrl.png", width=800, height=600)
plot(df_volcano$log2FoldChange, -log10(df_volcano$pvalue), pch=20, col=rgb(0,0,0,0.3),
     xlab="log2 Fold Change", ylab="-log10(p-value)", main="Volcano Plot: MEKi vs Control")
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-1, 1), col="blue", lty=2)
dev.off()

# GRÁFICO 3: MA PLOT
png("results/deseq2/MAPlot_MEKi_vs_ctrl.png", width=800, height=600)
plotMA(res_MEKi_vs_ctrl, ylim=c(-5,5), main="MA Plot: MEKi vs Control")
dev.off()

# GRÁFICO 4: HEATMAP DOS TOP 20 GENES
top20_indices <- head(order(res_MEKi_vs_ctrl$pvalue, decreasing=FALSE), 20)
matriz_heatmap <- assay(vsd)[top20_indices, ]
png("results/deseq2/Heatmap_Top20_Expressao.png", width=900, height=700)
pheatmap(matriz_heatmap, cluster_rows=TRUE, cluster_cols=TRUE, scale="row",
         annotation_col=data.frame(Tratamento=samples$condition, row.names=samples$sample),
         main="Top 20 Genes com Maior Diferença de Expressão")
dev.off()


# ENRIQUECIMENTO COM WEBGESTALT

# 1. Preparar a tabela para o WebGestalt (usando MEKi vs ctrl)
df_dados <- as.data.frame(res_MEKi_vs_ctrl)

# LIMPEZA: Remover o ponto final e a versão dos IDs Ensembl (ENSG000001.12 vira ENSG000001)
df_dados$gene_id <- sub("\\..*", "", rownames(df_dados))

# 2. Calcular o Score de Rank
df_dados$score <- -log10(df_dados$pvalue) * sign(df_dados$log2FoldChange)

# 3. Remover linhas com valores vazios (NA)
lista_gsea <- na.omit(df_dados[, c("gene_id", "score")])

# Criar a pasta para o WebGestalt se não existir
if(!dir.exists("results/webgestalt_kegg")) dir.create("results/webgestalt_kegg", recursive = TRUE)

# 4. Rodar o WebGestalt oficial por código (GSEA contra o KEGG)
WebGestaltR(
  enrichMethod     = "GSEA",
  organism         = "hsapiens",
  enrichDatabase   = "pathway_KEGG",
  interestGene     = lista_gsea,
  interestGeneType = "ensembl_gene_id",
  outputDirectory  = "results/webgestalt_kegg"
)


