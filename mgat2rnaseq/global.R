library(shiny)
library(tidyverse)
library(shinythemes)
library(DT)
#library(plotly)


library(tidyverse)
#library(DESeq2)
#library(org.Mm.eg.db)
#library(KEGGREST)
#library(biomaRt)
library(gghalves)
library(ggpubr)
select <- dplyr::select


##import files
remove_list = c("IKOHL_1")

metadata = read_csv("metadata.csv")%>%
  column_to_rownames("Sample_ID")%>%
  mutate(genotype = fct_relevel(genotype, c("WT","FF","KO","IKO")))%>%
  mutate(bile_acid = fct_relevel(bile_acid, c("low","high")))%>%
  filter(!rownames(.) %in% remove_list)%>%
  mutate(geno = case_when(genotype %in% c("FF","IKO") ~ "Intestine",
                          genotype %in% c("WT","KO") ~ "Global"))
# 
# 
# countdata = read_delim("featurecount3.txt", skip = 1)%>%
#   column_to_rownames("Geneid")%>%
#   select(-(1:5))%>%
#   rownames_to_column("rowname")%>%
#   mutate(rowname = gsub("\\..*","",rowname))%>%
#   column_to_rownames("rowname")
# 
# colnames(countdata) <- gsub("\\./aligned/+|\\.bam", "", colnames(countdata))
# 
# countdata = select(countdata, -remove_list)


##annotation files
ann_df <- read_csv("ann.csv")
# ann_dict <- setNames(ann_df$external_gene_name, ann_df$ensembl_gene_id)
# 
# 
# ann_dict <- ann_dict[!is.na(ann_dict)]




##normalize size factors of countdata
# dds <-  DESeqDataSetFromMatrix(countData = countdata,
#                                colData = metadata,
#                                design = ~b_genotype + diet + diet:b_genotype)
# dds <- estimateSizeFactors(dds)
# countdata_n <-  counts(dds, normalized = TRUE)
# 
# count_named = countdata_n%>%
#   as.data.frame()%>%
#   rownames_to_column("ensembl_gene_id")%>%
#   inner_join(ann_df, by = "ensembl_gene_id")%>%
#   select(ensembl_gene_id, external_gene_name, everything())
count_named <- read_csv("counttable_n.csv")

metadata_2 <- metadata%>%
  rownames_to_column("key")


plot_boxes <- function(gene){
  gene_gg <- count_named%>%
    filter(str_detect(str_to_lower(external_gene_name), str_to_lower(gene)))%>%
    select(-external_gene_name, -ensembl_gene_id, -description,-entrezgene_id)%>%
    gather()%>%
    left_join(metadata_2)%>%
    mutate(value = as.numeric(value))%>%
    mutate(b_genotype = fct_relevel(b_genotype, c("WTL","KOL","KOH")))
  
  
  
  comparisons = list(c("WTL","KOH"), c("WTL", "KOL"), c("KOL", "KOH"))
  
  ggplot(gene_gg, aes(x = b_genotype, y = value, color = b_genotype))+
    geom_half_violin(side = "l", aes(fill = b_genotype))+
    geom_boxplot(position = position_nudge(x = 0.22), width=0.2)+
    geom_half_point(side = "r")+
    facet_grid(geno~diet)+
    scale_fill_manual(values=c("#4E79A7","#59A14F","#F28E2B"),
                      labels = c("WT Low BA", "M2KO Low BA", "M2KO High BA"))+
    scale_color_manual(values=c("#4E79A7","#59A14F","#F28E2B"))+
    guides(color = "none")+
    stat_compare_means(method = "kruskal.test")+
    ggtitle("Gene Expression")

}




