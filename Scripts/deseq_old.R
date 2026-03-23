args <- commandArgs(trailingOnly = TRUE)
path = args[1]
thresholdFC = as.numeric(args[2])
pval= as.numeric(args[3])

library(DESeq2)
library(EnhancedVolcano)

result <- paste(path,"/results/",sep="")
setwd(result)
target <- paste(path,"/targeting/",sep="")
table <- paste(result,"miRtable.tsv",sep="")

main_rna_tab <- read.delim(table,row.names = NULL)

pvalstr=gsub("\\.","_",as.character(pval))

#aggregate mirna by names
rna_tab <- main_rna_tab
rna_tab$Sequence <- NULL
agg_df <- aggregate(rna_tab[-1], by=list(rna_tab$miR_Name), FUN=sum)
row.names(agg_df) <- agg_df$Group.1
agg_df$Group.1 <- NULL 

# ---- Define fixed control group ----
control_prefix <- "Mock"
control_cols <- grep(paste0("^", control_prefix), colnames(agg_df), value = TRUE)

# ---- Detect all other group prefixes dynamically ----
all_prefixes <- gsub("_[0-9]+$", "", colnames(agg_df))  # remove trailing digits
all_groups <- unique(all_prefixes)
comparison_groups <- setdiff(all_groups, control_prefix)

# ---- Run comparisons: Mock vs each group ----
for (grp in comparison_groups) {
  exp_cols <- grep(paste0("^", grp), colnames(agg_df), value = TRUE)
  if (length(exp_cols) == 0) {
    warning(paste("No columns found for group:", grp))
    next
  }
  
  # Combine columns and create colData
  count_matrix <- agg_df[ , c(control_cols, exp_cols)]
  group_labels <- c(rep("control", length(control_cols)), rep("infected", length(exp_cols)))
  colData <- data.frame(condition = factor(group_labels))
  # preparing table to filter for target analysis later in the pipeline
  rna_tab_exp <- main_rna_tab[,c(control_cols, exp_cols,"Sequence","miR_Name")]
  rna_tab_exp<-rna_tab_exp[rowSums(rna_tab_exp[,c(control_cols, exp_cols)])>0,]
  rna_tab_exp <- rna_tab_exp[,c("Sequence","miR_Name")]
  
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$log2FoldChange), ]
  df <- as.data.frame(res)
  
  # Volcano plot
  plotname<-paste0("Volcano_", grp, "_vs_Mock_pval",pvalstr,"_FC", thresholdFC,".png")
  ggsave(plotname,EnhancedVolcano(df,
                  lab = rownames(df),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  pCutoff = pval,
                  FCcutoff = thresholdFC))
  
  # Save significant DE miRNAs
  df <- df[!is.na(df$pvalue) & !is.na(df$log2FoldChange), ]
  df <- df[df$pvalue <pval,]
  subset_df_up <- df[df$log2FoldChange >= thresholdFC, ]
  subset_df_dw <- df[df$log2FoldChange <= -thresholdFC, ]
  subset_df_up$miRNA <- rownames(subset_df_up)
  subset_df_dw$miRNA <- rownames(subset_df_dw)

  # filtering and writing table for the target analysis next
  DE_mir_tab_exp_up <- rna_tab_exp[rna_tab_exp$miR_Name %in% subset_df_up$miRNA, ]
  DE_mir_tab_exp_dw <- rna_tab_exp[rna_tab_exp$miR_Name %in% subset_df_dw$miRNA, ]
  output_tab_up <- paste0(target, grp, "_UP.tsv")
  output_tab_dw <- paste0(target, grp, "_DOWN.tsv")
  write.table(DE_mir_tab_exp_up , output_tab_up,col.names=FALSE, row.names = FALSE,quote=FALSE,sep="\t")
  write.table(DE_mir_tab_exp_dw, output_tab_dw,col.names=FALSE, row.names = FALSE,quote=FALSE,sep="\t")
  
  # writing DEseq result per condtion & regulation
  output_file_up <- paste0("DESeq2_", grp, "_vs_Mock_pval",pvalstr,"_FC", thresholdFC, "_UP.csv")
  output_file_dw <- paste0("DESeq2_", grp, "_vs_Mock_pval",pvalstr,"_FC", thresholdFC, "_DOWN.csv")
  write.csv(subset_df_up, output_file_up, row.names = FALSE)
  write.csv(subset_df_dw, output_file_dw, row.names = FALSE)
  
  cat("Finished comparison:", control_prefix, "vs", grp, "\n")
}
