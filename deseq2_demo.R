library(tximport)
library(rtracklayer)
library(DESeq2)

quant_paths = c(
  "mi4_1" = "data/Mi4_1/quant.sf",
  "mi4_2" = "data/Mi4_2/quant.sf",
  "pm4_1" = "data/Pm4_1/quant.sf",
  "pm4_2" = "data/Pm4_2/quant.sf"
)

gene_model_path = "/share/data/dmel-all-r6.52.gtf.gz"
gene_model = import(gene_model_path)

gene_model = subset(
  gene_model,
  type == "mRNA"
)

tx_annotation = data.frame(
  TXNAME= gene_model$transcript_id,
  GENEID = gene_model$gene_symbol
)

salmon_result = tximport(
  files = quant_paths,
  type = "salmon",
  tx2gene = tx_annotation
)

## Tell DESeq2 what to compare by setting metadata
metadata = data.frame(
  row.names = colnames(salmon_result$abundance),
  celltype = c("Mi4", "Mi4", "Pm4", "Pm4")
)

degobj = DESeqDataSetFromTximport(salmon_result,
                                   colData = metadata,
                                   design = ~ celltype)

## Filter lowly expressed genes
not_low_in_at_least_n = 2
tokeep = rowSums(counts(degobj) >= 40) >= not_low_in_at_least_n
degobj = degobj[tokeep,]

degs = DESeq(degobj)
degtbl = results(degs, alpha = 0.01)

plotMA(degtbl)

library(pheatmap)
sig_tpm = salmon_result$abundance[row.names(degtbl), ]
pheatmap(sig_tpm)




