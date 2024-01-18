library(shiny)
library(tidyverse)
library(shinythemes)
library(DT)
#library(plotly)
library(ggplot2)

library(tidyverse)
#library(DESeq2)
#library(org.Mm.eg.db)
#library(KEGGREST)
#library(biomaRt)
library(gghalves)
library(ggpubr)
#library(png)
library(ComplexHeatmap)
select <- dplyr::select


msig_df <- read_csv("msig_df.csv")

##import files
remove_list <-  c("IKOHL_1")
hm_group_order <- c("WTC","KOLC","KOHC",
                    "FFC","IKOLC","IKOHC",
                    "WTL","KOLL","KOHL",
                    "FFL","IKOLL","IKOHL",
                    "WTH","KOLH","KOHH",
                    "FFH","IKOLH","IKOHH")
# pathways list
pathwaylist <- c("HALLMARK_ADIPOGENESIS", "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_ANDROGEN_RESPONSE", 
                 "HALLMARK_ANGIOGENESIS", "HALLMARK_APICAL_JUNCTION", "HALLMARK_APICAL_SURFACE", 
                 "HALLMARK_APOPTOSIS", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", 
                 "HALLMARK_COAGULATION", "HALLMARK_COMPLEMENT", "HALLMARK_DNA_REPAIR", 
                 "HALLMARK_E2F_TARGETS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                 "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE", 
                 "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_G2M_CHECKPOINT", 
                 "HALLMARK_GLYCOLYSIS", "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_HEME_METABOLISM", 
                 "HALLMARK_HYPOXIA", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
                 "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                 "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_KRAS_SIGNALING_DN", 
                 "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MTORC1_SIGNALING", 
                 "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_MYOGENESIS", 
                 "HALLMARK_NOTCH_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                 "HALLMARK_P53_PATHWAY", "HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_PEROXISOME", 
                 "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_PROTEIN_SECRETION", 
                 "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_SPERMATOGENESIS", 
                 "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                 "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_UV_RESPONSE_DN", 
                 "HALLMARK_UV_RESPONSE_UP", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", 
                 "HALLMARK_XENOBIOTIC_METABOLISM")
metadata <-  read_csv("metadata.csv")%>%
  column_to_rownames("Sample_ID")%>%
  mutate(genotype = fct_relevel(genotype, c("WT","FF","KO","IKO")))%>%
  mutate(bile_acid = fct_relevel(bile_acid, c("low","high")))%>%
  filter(!rownames(.) %in% remove_list)%>%
  mutate(geno = case_when(genotype %in% c("FF","IKO") ~ "Intestine",
                          genotype %in% c("WT","KO") ~ "Global"))


metadata_3 <- metadata%>%
  rownames_to_column("Sample_ID")%>%
  select(-Sample_ID, -sample_number, -b_genotype, -e_genotype)%>%
  distinct(group, .keep_all = TRUE)%>%
  #arrange(match(group, hm_group_order))%>%
  column_to_rownames("group")

metadata_3 <- metadata_3[hm_group_order, , drop = FALSE]

dds_vst2 <- read.csv("dds_matrix.csv",row.names = 1)
dds_vst_cname = colnames(dds_vst2)
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
ann_dict <- setNames(ann_df$external_gene_name, ann_df$ensembl_gene_id)


ann_dict <- ann_dict[!is.na(ann_dict)]


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







entrezid2ensembl_dict = setNames(ann_df$ensembl_gene_id,ann_df$entrezgene_id)



# plot box function
plot_boxes <- function(gene){
  gene_gg <- count_named%>%
    filter(str_detect(external_gene_name, regex(paste0("\\b", gene, "\\b"), ignore_case = TRUE)))%>%
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

#plot heatmap function
generate_heatmap <- function(dds_vst3) {
  # Create Heatmap Annotation
  column_ha <- HeatmapAnnotation(
    Diet = metadata_3$diet,
    Bile_acid = metadata_3$bile_acid,
    Genotype = metadata_3$genotype,
    col = list(
      Bile_acid = c("low" = "#5778a4", "high" = "#e49444"),
      Diet = c("low_fat" = "#F7D96A", "chow" = "#964B00", "high_fat" = "#499894"),
      Genotype = c("WT" = "#FF9D9A", "KO" = "#E15759", "FF" = "#D4A6C8", "IKO" = "#B07AA1")
    )
  )
  
  # Generate the Heatmap
  hm <- Heatmap(transpose_dataframe(dds_vst3),
                cluster_columns = FALSE,
                top_annotation = column_ha,
                column_names_side = "top",
                column_names_rot = 45)
  
  # Draw the Heatmap
  draw(hm)
}

#process matrix for heatmap
process_dds_vst <- function(path_gene_df_list) {
  dds_vst3 <- dds_vst2 %>%
    rownames_to_column("gene") %>%
    filter(gene %in% path_gene_df_list) %>%
    mutate(gene = ann_dict[gene]) %>%
    column_to_rownames("gene") %>%
    apply(1, scale) %>%
    as.data.frame()
  
  rownames(dds_vst3) <- dds_vst_cname
  
  return(dds_vst3)
}


#make pathway list
generate_path_gene_list <- function(pathway) {
  path_gene_df_list <- msig_df %>%
    filter(gs_name == pathway) %>%
    mutate(ensembl_gene_id = entrezid2ensembl_dict[entrez_gene]) %>%
    pull(ensembl_gene_id)
  
  return(path_gene_df_list)
}


#transpose dataframe
transpose_dataframe <- function(df) {
  # Store original column and row names
  col_names <- colnames(df)
  row_names <- rownames(df)
  
  # Transpose the dataframe to a matrix, then back to a dataframe
  df_transposed <- as.data.frame(t(as.matrix(df)))
  
  # Reassign the original column and row names
  colnames(df_transposed) <- row_names
  rownames(df_transposed) <- col_names
  
  return(df_transposed)
}