normalized_counts <- read.table("GSE255581.txt", header = TRUE, sep ="\t")

genename <- c("COQ2","GCLM","LPCAT3","MAP1LC3C","PCBP1","PCBP2","ACSL3","ACSL4","DPP4","FTH1","FTL","GPX4","HMOX1","NCOA4","KEAP1","SLC7A11","SCD","SCD5")

library(tidyr)

heatmap_in2<-normalized_counts%>%
  filter(gene %in% genename)%>%
  pivot_longer(-1,names_to = "Sample",values_to = "Exp")%>%
  group_by(gene)%>%
  mutate(Exp=(Exp-mean(Exp))/sd(Exp))%>%
  ungroup()

tidyHeatmap::heatmap(.data = heatmap_in2,
                     .row = gene,
                     .column=Sample,
                     .value = Exp,
                     cluster_columns = FALSE,
                     show_row_names=TRUE,
                     palette_value = c("blue", "white", "red"))