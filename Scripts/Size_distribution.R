args <- commandArgs(trailingOnly = TRUE)
distribution_file=args[1]
outfile=file.path(dirname(distribution_file), "Percent_size _of_sRNA_relative_mapped_reads_in_each_replicate.svg")

packages <- c("readr", "ggplot2", "reshape2")
options(repos = c(CRAN = "https://cloud.r-project.org"))
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

distribution_percent <- read_delim(distribution_file)

distribution_percent$size<- as.factor(distribution_percent$size)

data_plt <- melt(distribution_percent,id.vars="size", variable.name = "replicates")
data_plt$condition <- gsub("_[0-9]", "", data_plt$replicates)
data_plt$size <- gsub("$", "-nt", data_plt$size)

ggplot(data_plt,aes(x=size,y=value,fill=condition,group=replicates))+
  geom_col(position = "dodge",colour = "black",linewidth = 0.1)+scale_fill_manual(values=c('#ff9f4e','#ff2a0e','#1f77b4'))+
  labs(y="% of mapped reads",x="Read size",fill="Condition")
ggsave(outfile,width=10,height=7)


