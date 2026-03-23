args <- commandArgs(trailingOnly = TRUE)
path = args[1]
thresholdFC = as.numeric(args[2])
thresholdpval= as.numeric(args[3])

library(DESeq2)
library(EnhancedVolcano)

result <- paste(path,"/results/",sep="")
setwd(result)
table <- paste(result,"/host_mapped_23-24nt_smallrna.tsv",sep="")

main_rna_tab <- read.delim(table,row.names = NULL)

pvalstr=gsub("\\.","_",as.character(thresholdpval))  

row.names(main_rna_tab) <- main_rna_tab$Sequence_Id
main_rna_tab$Sequence <- NULL
main_rna_tab$Sequence_Length <- NULL
main_rna_tab$Sequence_Id <- NULL


# ---- Define fixed control group ----
control_prefix <- "Mock"
control_cols <- grep(paste0("^", control_prefix), colnames(main_rna_tab), value = TRUE)

# ---- Detect all other group prefixes dynamically ----
all_prefixes <- gsub("_[0-9]+$", "", colnames(main_rna_tab))  # remove trailing digits
all_groups <- unique(all_prefixes)
comparison_groups <- setdiff(all_groups, control_prefix)

# ---- Run comparisons: Mock vs each group ----
for (grp in comparison_groups) {
  exp_cols <- grep(paste0("^", grp), colnames(main_rna_tab), value = TRUE)
  if (length(exp_cols) == 0) {
    warning(paste("No columns found for group:", grp))
    next
  }
  cat("Comparison:", control_prefix, "vs", grp, "\n")
  # Combine columns and create colData
  count_matrix <- main_rna_tab[ , c(control_cols, exp_cols)]
  group_labels <- c(rep("control", length(control_cols)), rep("infected", length(exp_cols)))
  colData <- data.frame(condition = factor(group_labels))
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~condition)
  dds <- DESeq(dds,parallel=TRUE)
  res <- results(dds,parallel=TRUE)
  res <- res[order(res$log2FoldChange), ]
  df <- as.data.frame(res)
  
  # Volcano plot
  plotname<-paste0("host_mapped_23-24nt_smallrna_Volcano_", grp, "_vs_Mock_adjpval",pvalstr,"_FC", thresholdFC,".png")
  ggsave(plotname,EnhancedVolcano(df,
                  lab = rownames(df),
                  x = 'log2FoldChange',
                  y = 'padj',
                  pCutoff = thresholdpval,
                  FCcutoff = thresholdFC))
  
  # Save significant DE miRNAs
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  df <- df[df$padj <thresholdpval,]
  df <- df[abs(df$log2FoldChange) >= thresholdFC, ]
  df$Sequence_Id <- rownames(df)
  df$regulation <- with(df,ifelse(log2FoldChange > 0,'Up',"Down"))
  output_file <- paste0("host_mapped_23-24nt_smallrna_DESeq2_", grp, "_vs_Mock_adjpval",pvalstr,"_FC", thresholdFC, ".csv")
  write.csv(df, output_file, row.names = FALSE)
  
  # writing DEseq result per condtion & regulation
  # output_file_up <- paste0("DESeq2_", grp, "_vs_Mock_pval",pvalstr,"_FC", thresholdFC, "_UP.csv")
  # output_file_dw <- paste0("DESeq2_", grp, "_vs_Mock_pval",pvalstr,"_FC", thresholdFC, "_DOWN.csv")
  # write.csv(subset_df_up, output_file_up, row.names = FALSE)
  # write.csv(subset_df_dw, output_file_dw, row.names = FALSE)
  
  cat("Finished comparison:", control_prefix, "vs", grp, "\n")
}
