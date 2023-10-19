get_type <- function(variable){ 
  sz <- as.integer(length(variable)) #length of your variable 
  tof <- typeof(variable)            #typeof your variable 
  cls <- class(variable)             #class of your variable 
  isc <- is.character(variable)      #what is.character() has to say about it.  
  d <- dim(variable)                 #dimensions of your variable 
  isv <- is.vector(variable) 
  if (is.matrix(variable)){  
    d <- dim(t(variable))             #dimensions of your matrix
  }    
  #observations ----> datatype 
  if (sz>=1 && tof == "logical" && cls == "logical" && isv == TRUE){ return("vector of logical") } 
  if (sz>=1 && tof == "integer" && cls == "integer" ){ return("vector of integer") } 
  if (sz==1 && tof == "double"  && cls == "Date" ){ return("Date") } 
  if (sz>=1 && tof == "raw"     && cls == "raw" ){ return("vector of raw") } 
  if (sz>=1 && tof == "double"  && cls == "numeric" ){ return("vector of double") } 
  if (sz>=1 && tof == "double"  && cls == "array" ){ return("vector of array of double") } 
  if (sz>=1 && tof == "character"  && cls == "array" ){ return("vector of array of character") } 
  if (sz>=0 && tof == "list"       && cls == "data.frame" ){ return("data.frame") } 
  if (sz>=1 && isc == TRUE         && isv == TRUE){ return("vector of character") } 
  if (sz>=1 && tof == "complex"    && cls == "complex" ){ return("vector of complex") } 
  if (sz==0 && tof == "NULL"       && cls == "NULL" ){ return("NULL") } 
  if (sz>=0 && tof == "integer"    && cls == "factor" ){ return("factor") } 
  if (sz>=1 && tof == "double"     && cls == "numeric" && isv == TRUE){ return("vector of double") } 
  if (sz>=1 && tof == "double"     && cls == "matrix"){ return("matrix of double") } 
  if (sz>=1 && tof == "character"  && cls == "matrix"){ return("matrix of character") } 
  if (sz>=1 && tof == "list"       && cls == "list" && isv == TRUE){ return("vector of list") } 
  if (sz>=1 && tof == "closure"    && cls == "function" && isv == FALSE){ return("closure/function") } 
  return("custom container") 
} 
assert <- function(a, b){ 
  if (a == b){ 
    cat("P") 
  } 
  else{ 
    cat("\nFAIL!!!  Sniff test:\n") 
    sz <- as.integer(length(variable))   #length of your variable 
    tof <- typeof(variable)              #typeof your variable 
    cls <- class(variable)               #class of your variable 
    isc <- is.character(variable)        #what is.character() has to say about it. 
    d <- dim(variable)                   #dimensions of your variable 
    isv <- is.vector(variable) 
    if (is.matrix(variable)){  
      d <- dim(t(variable))                   #dimensions of your variable 
    } 
    if (!is.function(variable)){ 
      print(paste("value: '", variable, "'")) 
    } 
    print(paste("get_type said: '", a, "'")) 
    print(paste("supposed to be: '", b, "'")) 
    
    cat("\nYour pointer to memory has properties:\n")  
    print(paste("sz: '", sz, "'")) 
    print(paste("tof: '", tof, "'")) 
    print(paste("cls: '", cls, "'")) 
    print(paste("d: '", d, "'")) 
    print(paste("isc: '", isc, "'")) 
    print(paste("isv: '", isv, "'")) 
    quit() 
  } 
}

ddsres_to_df <- function(dds, contrast){
  results(dds,contrast)%>%
    as.data.frame()%>%
    na.omit()%>%
    rownames_to_column("ensembl_gene_id")%>%
    arrange(pvalue)%>%
    inner_join(ann_df, by = "ensembl_gene_id")%>%
    na.omit()%>%
    distinct()
}

get_dds_diets <- function(dds){
  res_h <- ddsres_to_df(dds, list(c("bile_acid_high_vs_low","bile_acidhigh.diethigh_fat")))%>%
    select(entrezgene_id, stat)%>%
    arrange(desc(stat))%>%
    #filter(abs(stat)>2)%>%
    deframe()
  
  res_l <- ddsres_to_df(dds, list(c("bile_acid_high_vs_low","bile_acidhigh.dietlow_fat")))%>%
    select(entrezgene_id, stat)%>%
    arrange(desc(stat))%>%
    #filter(abs(stat)>2)%>%
    deframe()
  
  res_c <- ddsres_to_df(dds, list(c("bile_acid_high_vs_low")))%>%
    select(entrezgene_id, stat)%>%
    arrange(desc(stat))%>%
    
    
    #filter(abs(stat)>2)%>%
    deframe()
  
  
  
  
  
  res_all <- list(c = res_c, hf = res_h, lf =res_l)
  
  
  return(res_all)
}

pathview_3 <- function(gene_data = kegg_data, 
                       pathways = c("00062"),
                       out.suffix = "pv", 
                       save = "KOHKOL"){
  current_wd <- getwd()
  save_path <- paste0(current_wd,"/", save)
  if (!file.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)  
  }
  setwd(save_path)
  pathview(gene.data = gene_data,
           pathway.id = pathways,
           species = "mmu",
           out.suffix = out.suffix,
           gene.idtype = "entrez",
           limit = list(gene = 5),
           bins = list(gene = 10))
  setwd(current_wd)
  return(save_path)
}

split_df <- function(name_list, suffix, count = countdata, meta = metadata){
  
  countdata_out = count%>%
    select(colnames(.)[colnames(.) %in% name_list])
  
  metadata_out <- meta%>%
    filter(rownames(.) %in% name_list)
  
  countdata_new <- paste0("countdata_", suffix)
  metadata_new <- paste0("metadata_", suffix)
  
  assign(countdata_new, countdata_out, envir = .GlobalEnv)
  assign(metadata_new, metadata_out, envir = .GlobalEnv)
  
  return(print(c(dim(countdata_out),dim(metadata_out))))
}

do_hmgsea <- function(genelist){
  gsea_res <- GSEA(genelist, TERM2GENE = msig_df,pvalueCutoff = 1)
  result <- setReadable(gsea_res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")%>%
    filter(ID %in% pathway_list)%>%
    mutate(ID = fct_relevel(ID, 
                            
                            "HALLMARK_MYC_TARGETS_V1",
                            "HALLMARK_MYC_TARGETS_V2",
                            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                            "HALLMARK_MITOTIC_SPINDLE",
                            "HALLMARK_MTORC1_SIGNALING",
                            "HALLMARK_P53_PATHWAY"
    ))%>%
    
    return(result)
}

heatplot_gsea <- function(gsearesult,genelist,name){
  

    heatplot(gsearesult, foldChange = genelist, label_format = 100)+
    scale_y_discrete(position = "right")+
    coord_flip()+
    ggtitle(name)+
    theme_classic()+
    theme(axis.text.y = element_text(size  = 3),
          axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     size = 8),)
  
  ggsave(paste0(name,".png"), dpi = 600, width = 8, height = 20)    
}


do_dds <- function(count, meta, design, n_cutoff = 1){
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData = meta,
                                design = design) 
  dds <-  dds[rowSums(counts(dds)) > n_cutoff]
  dds <- DESeq(dds)
  
  return(dds)
}

do_vst <- function(count, meta, design, n_cutoff = 1){
  dds <-  DESeqDataSetFromMatrix(countData = count,
                                 colData = meta,
                                 design = design)
  dds <- dds[rowSums(counts(dds))  > n_cutoff]
  dds <- varianceStabilizingTransformation(dds, blind = TRUE)
  #normalized_counts <- counts(dds, normalized =TRUE)
  return(assay(dds))
}

do_dds_norm <- function(count, meta, design, n_cutoff = 1){
  dds <-  DESeqDataSetFromMatrix(countData = count,
                                 colData = meta,
                                 design = design)
  
  dds <- dds[rowSums(counts(dds))  > n_cutoff]
  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized =TRUE)
  return(normalized_counts)
}


heatplot_vst <- function(gsearesult,genelist,name){
  
  
  heatplot(gsearesult, foldChange = genelist, label_format = 100)+
    scale_y_discrete(position = "right")+
    coord_flip()+
    ggtitle(name)+
    theme_classic()+
    theme(axis.text.y = element_text(size  = 3),
          axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     size = 8),)
  
  ggsave(paste0(name,".png"), dpi = 600, width = 8, height = 20)    
}

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
