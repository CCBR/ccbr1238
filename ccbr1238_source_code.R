################################################################################
# formatting
################################################################################
## shorten any list of names
SHORTEN_NAMES<-function(list_in,to_remove){
  shortened_names=sub(to_remove,"",list_in)
  return(shortened_names)
}

################################################################################
# initial QC
################################################################################
## input the counts matrix
FORMAT_COUNTS_MATRIX<-function(contrast_id,sampleinfo){
  # set extension
  extensions=c(paste0("__",dedup_status,"__",norm_type_cutandrun,".bed"))
  
  # set variables
  rawcountsmatrix=paste0(car_subpath,contrast_id,extensions,"/",
                         contrast_id,extensions,"_",method,"countsmatrix.txt")
  # prep counts
  rawcounts = read.csv(rawcountsmatrix,
                       header = TRUE,sep="\t",
                       comment.char = "#", 
                       strip.white = TRUE,
                       check.names = FALSE,
                       colClasses = "character")
  rawcounts = as.data.frame(rawcounts)
  rawcounts %>% column_to_rownames("peakID") -> rawcounts
  
  # convert character to numeric to integer
  x = matrix(as.numeric(as.matrix(rawcounts)),ncol=ncol(rawcounts))
  x = matrix(mapply(x,FUN=as.integer),ncol=ncol(rawcounts))
  x = as.data.frame(x)
  colnames(x) = colnames(rawcounts)
  rownames(x) = rownames(rawcounts)
  rawcounts = x
  
  return(rawcounts)
}

## generate RLA and PCOA plots
GENERATE_RLE_PLOT<-function(input_data,input_samples,input_title,stage_in){
  #input_data=filtered;input_samples=sample_list;input_title="Fig. Before Normalization"
  
  #set colors
  x=subset(groups_df,replicate %in% input_samples)$group
  colors <- brewer.pal(length(unique(x)), "Set2")
  x=as.factor(x)#set pheno data
  
  #Plot results
  par(mfrow=c(1, 2), oma=c(3, 2, 0, 0)+0.1)
  if (stage_in=="before_norm"){
    plotRLE(as.matrix(input_data), outline=FALSE, ylim=c(-.5, .5), col=colors[x],las=2, cex.axis = .8)
    plotPCA(as.matrix(input_data), col=colors[x], cex=.8)
  } else if (stage_in=="upper_quant"){
    plotRLE(input_data, outline=FALSE, ylim=c(-.5, .5), col=colors[x],las=2, cex.axis = .8)
    plotPCA(input_data, col=colors[x], cex=.8)
  } else if (stage_in=="DESEQ"){
    plotRLE(counts(input_data, normalize=TRUE), outline=FALSE, ylim=c(-.5, .5), col=colors[x],las=2, cex.axis = .8)
    plotPCA(counts(input_data,normalize=TRUE), col=colors[x], cex=.8)
  }
  mtext(input_title, side=2, outer=TRUE, adj=0)  
}

## generate RMS boxplots
GENERATE_BOXPLOTS<-function(sampleinfo,filtered){
  #set lib
  #determine lib reduction factor
  if (mean(colSums(filtered))>10000000){
    lib_factor=1e6
  } else if (mean(colSums(filtered))>1000000){
    lib_factor=1e5
  } else if (mean(colSums(filtered))>100000){
    lib_factor=1e4
  } else if (mean(colSums(filtered))>10000){
    lib_factor=1e3
  } else if (mean(colSums(filtered))>1000){
    lib_factor=1e2
  } else if (mean(colSums(filtered))>100){
    lib_factor=1e1
  } else {
    lib_factor=1e1
  }
  print(paste0("the lib",lib_factor," ",mean(colSums(filtered))))
  
  # filter
  sampleinfo=sampleinfo[sampleinfo$replicate==colnames(filtered),]
  sampleinfo$library_size=colSums(filtered)/lib_factor
  sampleinfodf = as.data.frame(sampleinfo)
  rownames(sampleinfo) = sampleinfo$sampleid
  pander(sampleinfodf,style="rmarkdown")
  
  # melt data
  rawcounts_logcpm = log2(cpm(filtered))
  cpm_melt=reshape2::melt(rawcounts_logcpm)
  colnames(cpm_melt)=c("geneID","sampleid","log2cpm")
  
  # print boxplots
  p = ggplot(cpm_melt,aes(x=sampleid,y=log2cpm)) + 
    geom_boxplot(fill=as.factor(as.numeric(as.factor(sampleinfo$group))+1)) +
    theme_classic() +
    coord_flip()
  print(p)
}

## FORMAT_COUNTS_MATRIX
## GENERATE_RLE_PLOT
## GENERATE_BOXPLOTS
MAIN_PREP_QC_CORE<-function(groups_in,stage_in){
  # set sample_list
  sample_list=subset(groups_df,group%in%groups_in)$replicate
  
  # read in rawcounts
  fpath=paste0(pipeliner_dir,"DEG_ALL/RawCountFile_RSEM_genes.txt")
  raw_counts=read.csv(fpath,sep="\t")
  
  colnames(raw_counts)=SHORTEN_NAMES(colnames(raw_counts),"_expected_count")
  raw_counts=raw_counts[,c("symbol",sample_list)]

  out_raw_counts=separate(raw_counts,col="symbol",into=c("ENSEMBL","SYMBOL"),sep="[|]")
  fpath=paste0(output_dir,"replicate_prenormalized_",csid,".csv")
  write.csv(out_raw_counts,fpath)
  
  # filter
  ## CPM is calcualted as "how many counts would I get for a gene if the sample 
  #had a library size of 1M".
  raw_counts=column_to_rownames(raw_counts)
  filtered=ceiling(raw_counts)
  cpm_counts=edgeR::cpm(as.matrix(filtered))
  log_cpm_counts=log2(cpm_counts)
  
  # filter by threshold
  sample_count_threshold=round(ncol(raw_counts)*sample_consensus_threshold+.5)-1
  print(paste0("Filtering for genes with min reads in at least ",sample_count_threshold," samples."))
  keep=rowSums(cpm_counts>0.5)>sample_count_threshold

  filtered=filtered[keep,]
  head(filtered)
  
  if (stage_in=="before_norm"){
    # generate RLE plot
    GENERATE_RLE_PLOT(filtered,sample_list,"Fig. Before Normalization",stage_in)
    
    # generate boxplots
    sampleinfo=subset(groups_df,group%in%group_list)
    rownames(sampleinfo)=NULL
    GENERATE_BOXPLOTS(sampleinfo,filtered)
  } else if (stage_in %in% c("upper_quant","DESEQ")){
    # run upper quant norm
    x=subset(groups_df,replicate %in% sample_list)$group
    x=as.factor(x)
    set <- newSeqExpressionSet(as.matrix(filtered),
                               phenoData = data.frame(x, row.names=colnames(filtered)))
    
    if (stage_in=="upper_quant"){
      set_u <- betweenLaneNormalization(set, which="upper")
      GENERATE_RLE_PLOT(set_u,sample_list,"Fig. UpperQuant Normalization",stage_in)  
    } else { # run DESEQ
      dds <- DESeqDataSetFromMatrix(countData = counts(set),
                                    colData = pData(set),
                                    design = ~ x)
      dds <- DESeq(dds)
      GENERATE_RLE_PLOT(dds,sample_list,"Fig. DESEq2 Normalization",stage_in)
      
      # normalize count matrix
      Normalized_counts_matrix<-as.data.frame(counts(dds, norm=TRUE))
      nrow(Normalized_counts_matrix)
      
      # add flags for passing
      sample_count_threshold=round(length(sample_list)*sample_consensus_threshold+.5)-1
      Normalized_counts_matrix$sample_threshold=rowSums(Normalized_counts_matrix[,sample_list] > (read_minimum_threshold-1))
      Normalized_counts_matrix$sample_threshold_flag=ifelse(Normalized_counts_matrix$sample_threshold > sample_consensus_threshold,"Y","N")
      head(Normalized_counts_matrix)
      
      fpath=paste0(output_dir,"replicate_normalized_",csid,".csv")
      write.csv(Normalized_counts_matrix,fpath)
    }
  }
}

################################################################################
# primary differential
################################################################################
# reads diff results from pipeliner
PRIMARY_DIFFERENTIAL<-function(contrast_in){
  # id contrast
  print(paste0("** Processing ", contrast_in," **"))
  
  # read in DEG
  contrast_hypen=gsub("_vs_","-",contrast_in)
  fpath=paste0(pipeliner_dir,"DEG_",contrast_in,"_0.5_0.5/",
               "DESeq2_DEG_",contrast_in,"_all_genes.txt")
  deg_df=read.csv(fpath,sep="\t")
  
  # split EID and gene, remove duplicate gene/symbol name
  deg_df = deg_df %>% separate(col=ensid_gene,
                               into=c("ensembleID","SYMBOL"),
                               sep="[|]")
  deg_df=deg_df[,c(1:2,4:ncol(deg_df))]
  
  # determine sig
  deg_df$significance="N"
  deg_df$significance[abs(deg_df$log2fc)>log2_cutoff & deg_df$fdr<fdr_cutoff]="Y"
  
  # print out result
  print(paste0("--the number of significant genes found: ",nrow(subset(deg_df,significance=="Y"))))
  
  # write out DEG
  fpath=paste0(output_dir,"DEG_",contrast_in,".csv")
  write.csv(deg_df,fpath,row.names=FALSE)
}

COMPARATIVE_VENN_DIAGRAM<-function(compare_list_in,subtitle_in){
  # for each of the contrasts, read in the DEG file and generate sig lists
  for (contrast_id in compare_list_in){
    sample_name=gsub("_0pt5mM-SCR_0pt5mM","",contrast_id)
    sample_name=gsub("_4mM-SCR_0pt5mM","4",sample_name)
    print(paste0("** Processing ", sample_name," **"))
    
    # read DEG from primary differential
    fpath=paste0(output_dir,"DEG_",contrast_id,".csv")
    deg_df=read.csv(fpath)
    
    # create gene list for sig values only
    tmp_list=subset(deg_df,significance=="Y")$SYMBOL
    
    # filter NA
    tmp_list=tmp_list[!is.na(tmp_list)]
    
    # assign to SH name list
    assign(paste0(sample_name,"_siglist"),tmp_list)
  }
  # create a sample list
  sample_list=gsub("_0pt5mM-SCR_0pt5mM","",compare_list_in)
  sample_list=gsub("_4mM-SCR_0pt5mM","4",sample_list)
  sample_sig_list=paste0(sample_list,"_siglist")
  
  # List of genes
  x <- list(A = get(sample_sig_list[1]), B = get(sample_sig_list[2]), C=get(sample_sig_list[3]))
  
  # Venn diagram with custom category names
  p = ggVennDiagram(x, color = 1, lwd = 0.7,
                    category.names = c(sample_list[1], sample_list[2],sample_list[3])) + 
    scale_fill_gradient(low = "red", high = "blue")
  full_title=paste0("Significantly Differentiated Genes:\n",subtitle_in)
  pf = p + ggtitle(full_title)
  
  # save and print venn diagram
  print(pf)
  fpath=paste0(output_dir,"venndiagram_",subtitle,"_genes.pdf")
  ggsave(fpath,pf)
  
  # find intersections
  union_all=union(get(sample_sig_list[1]),get(sample_sig_list[2]))
  union_all=union(union_all,get(sample_sig_list[3]))
  length(union_all)
  
  in_all=intersect(get(sample_sig_list[1]),get(sample_sig_list[2]))
  in_all=intersect(in_all,get(sample_sig_list[3]))
  length(in_all)
  
  cross_12=setdiff(intersect(get(sample_sig_list[1]),get(sample_sig_list[2])),in_all)
  length(cross_12)
  
  cross_13=setdiff(intersect(get(sample_sig_list[1]),get(sample_sig_list[3])),in_all)
  length(cross_13)
  
  cross_23=setdiff(intersect(get(sample_sig_list[2]),get(sample_sig_list[3])),in_all)
  length(cross_23)
  
  only_1=setdiff(get(sample_sig_list[1]),in_all)
  only_1=setdiff(only_1,cross_12)
  only_1=setdiff(only_1,cross_13)
  length(only_1)
  
  only_2=setdiff(get(sample_sig_list[2]),in_all)
  only_2=setdiff(only_2,cross_12)
  only_2=setdiff(only_2,cross_23)
  length(only_2)
  
  only_3=setdiff(get(sample_sig_list[3]),in_all)
  only_3=setdiff(only_3,cross_13)
  only_3=setdiff(only_3,cross_23)
  length(only_3)
  
  # create df
  venn_df=data.frame(union_all)
  venn_df$status=""
  venn_df$status[venn_df$union_all %in% in_all]="sig_all"
  venn_df$status[venn_df$union_all %in% cross_12]=paste0("sig_",sample_list[1],"_and_",sample_list[2])
  venn_df$status[venn_df$union_all %in% cross_13]=paste0("sig_",sample_list[1],"_and_",sample_list[3])
  venn_df$status[venn_df$union_all %in% cross_23]=paste0("sig_",sample_list[2],"_and_",sample_list[3])
  venn_df$status[venn_df$union_all %in% only_1]=paste0("sig_",sample_list[1])
  venn_df$status[venn_df$union_all %in% only_2]=paste0("sig_",sample_list[2])
  venn_df$status[venn_df$union_all %in% only_3]=paste0("sig_",sample_list[3])
  head(venn_df)
  
  # write out df
  fpath=paste0(output_dir,"venn_diagram_significance_",subtitle,".csv")
  write.csv(venn_df,fpath)
}

DUAL_VENN_DIAGRAM<-function(compare_list_in,subtitle_in){
  # for each of the contrasts, read in the DEG file and generate sig lists
  merged_df=data.frame()
  for (contrast_id in compare_list_in){
    sample_name=gsub("_0pt5mM-SCR_0pt5mM","",contrast_id)
    print(paste0("** Processing ", sample_name," **"))
    
    # read DEG from primary differential
    fpath=paste0(output_dir,"DEG_",contrast_id,".csv")
    deg_df=read.csv(fpath)
    
    # create gene list for sig values only
    sub_df=subset(deg_df,significance=="Y")
    tmp_list=sub_df$SYMBOL
    
    # filter NA
    tmp_list=tmp_list[!is.na(tmp_list)]
    
    # assign to SH name list
    assign(paste0(sample_name,"_siglist"),tmp_list)
    
    sub_df
    merged_df=rbind(merged_df,sub_df)
  }
  
  # create a sample list
  sample_list=gsub("_0pt5mM-SCR_0pt5mM","",compare_list_in)
  sample_sig_list=paste0(sample_list,"_siglist")
  
  # List of genes
  x <- list(A = get(sample_sig_list[1]), B = get(sample_sig_list[2]))
  
  # Venn diagram with custom category names
  p = ggVennDiagram(x, color = 1, lwd = 0.7,
                    category.names = c(sample_list[1], sample_list[2])) + 
    scale_fill_gradient(low = "red", high = "blue")
  full_title=paste0("Significantly Differentiated Genes:\n",subtitle_in)
  pf = p + ggtitle(full_title)
  
  # save and print venn diagram
  print(pf)
  fpath=paste0(output_dir,"venndiagram_",subtitle,"_genes.pdf")
  ggsave(fpath,pf)
  
  # find intersections
  union_all=union(get(sample_sig_list[1]),get(sample_sig_list[2]))
  length(union_all)
  
  in_all=intersect(get(sample_sig_list[1]),get(sample_sig_list[2]))
  length(in_all)
  
  cross_12=setdiff(intersect(get(sample_sig_list[1]),get(sample_sig_list[2])),in_all)
  length(cross_12)
  
  only_1=setdiff(get(sample_sig_list[1]),in_all)
  only_1=setdiff(only_1,cross_12)
  length(only_1)
  
  only_2=setdiff(get(sample_sig_list[2]),in_all)
  only_2=setdiff(only_2,cross_12)
  length(only_2)
  
  # create dual sig df
  
  # create df
  venn_df=data.frame(union_all)
  venn_df$status=""
  venn_df$status[venn_df$union_all %in% in_all]="sig_all"
  venn_df$status[venn_df$union_all %in% cross_12]=paste0("sig_",sample_list[1],"_and_",sample_list[2])
  venn_df$status[venn_df$union_all %in% only_1]=paste0("sig_",sample_list[1])
  venn_df$status[venn_df$union_all %in% only_2]=paste0("sig_",sample_list[2])
  head(venn_df)
  
  # write out df
  fpath=paste0(output_dir,"venn_diagram_significance_",subtitle,".csv")
  write.csv(venn_df,fpath)
}

RUN_DESEQ2<-function(counts_in,sample_list){
  # counts_in=filtered;
  
  # create groups df
  sampleinfo=subset(groups_df,replicate %in% sample_list)[,c("replicate","group")]
  head(sampleinfo)
  colnames(sampleinfo)=c("sampleid","group")
  sampleinfo$group=as.factor(sampleinfo$group)
  
  # run DESEQ
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_in),
                                colData = sampleinfo,
                                design = ~ group)
  dds <- DESeq(dds)
  
  return(dds)
}

GROUP_DIFFERENTIAL<-function(group1_in,group2_in,control_in){
  
  #group1_in="SH1_0pt5mM"; group2_in="SH4_0pt5mM"; control_in="SCR_0pt5mM"
  
  # set contrasts
  contrast_1=paste0(group1_in,"-",control_in)
  contrast_2=paste0(group2_in,"-",control_in)
  group1_group2=paste0(group1_in,"_vs_",group2_in)
  sample_list=subset(groups_df,group %in% c(group1_in,group2_in,control_in))$replicate
  
  # read in counts matrix
  fpath=paste0(pipeliner_dir,"DEG_ALL/RawCountFile_RSEM_genes.txt")
  rawcounts=read.csv(fpath,sep="\t")
  colnames(rawcounts)=SHORTEN_NAMES(colnames(rawcounts),"_expected_count")
  rownames(rawcounts)=rawcounts$symbol
  head(rawcounts)
  
  # subset for groups
  rawcounts=rawcounts[,c(sample_list)]
  
  # filter
  filtered=ceiling(rawcounts)
  cpm_counts=edgeR::cpm(as.matrix(filtered))
  log_cpm_counts=log2(cpm_counts)
  keep=rowSums(cpm_counts>0.5)>2
  filtered=filtered[keep,]
  head(filtered)
  
  # run DESEQ
  dds=RUN_DESEQ2(filtered,sample_list)
  
  # pull out contrasts
  # https://support.bioconductor.org/p/101494/
  res <- as.data.frame(results(dds, contrast=c("group",group1_in,group2_in)))
  res$merged_name=rownames(res)
  
  # split EID and gene, remove duplicate gene/symbol name
  deg_df = res %>% separate(col=merged_name,
                                            into=c("ensembleID","SYMBOL"),
                                            sep="[|]")
  
  # determine sig
  deg_df$significance="N"
  deg_df$significance[abs(deg_df$log2FoldChange)>log2_cutoff & deg_df$padj<fdr_cutoff]="Y"
  
  # calculated gsea ranking score
  deg_df$gsea_ranking_score=-log10(deg_df$pvalue)*sign(deg_df$log2FoldChange)
  
  # print out result
  print(paste0("--the number of significant genes found: ",nrow(subset(deg_df,significance=="Y"))))
  
  # rename cols
  colnames(deg_df)=gsub("log2FoldChange","log2fc",colnames(deg_df))
  
  # output df
  fpath=paste0(output_dir,"DEG_",group1_group2,".csv")
  write.csv(deg_df,fpath)
}

################################################################################
# GSEA
################################################################################
# load necessary files for mus musculus
LOAD_KEGG_MMU_FILES<-function(){
  #https://cran.r-project.org/web/packages/pathfindR/vignettes/intro_vignette.html
  #https://cran.r-project.org/web/packages/pathfindR/vignettes/obtain_data.html
  ## Save both as RDS/SIF files for later use
  RDS_genes=paste0("~/../../Volumes/ccbr1238/scripts/mmu_kegg_genes.RDS")
  RDS_desc=paste0("~/../../Volumes/ccbr1238/scripts/mmu_kegg_descriptions.RDS")
  
  # check both RDS files have been created; if not create them
  if (!file.exists(RDS_genes) || !file.exists(RDS_genes)){
    print("Creating RDS files")
    gsets_list <- get_gene_sets_list(source = "KEGG",org_code = "mmu")
    
    mmu_kegg_genes <- gsets_list$gene_sets
    mmu_kegg_descriptions <- gsets_list$descriptions
    
    saveRDS(mmu_kegg_genes, RDS_genes)
    saveRDS(mmu_kegg_descriptions, RDS_desc)
  }
  
  # load the files
  if (file.exists(RDS_genes) & file.exists(RDS_genes) & file.exists(SIF_PIN)){
    print("Loading RDS files")
    mmu_kegg_genes <- readRDS(RDS_genes)
    mmu_kegg_descriptions <- readRDS(RDS_desc)
  } else{
    print("Error creating RDS files")
    exit()
  }
}

# set the annotation dbs
DB_LOOKUP<-function(t2g,species_in){
  # t2g=db_id; species_in=species_in
  # generate gene lists for C1
  # generate gene lists for C2 with subtypes biocarta, kegg, reactome, wiki
  # generate gene lists for C5 with subtypes MF, BP, CC
  # generate gene lists for Hallmark
  if (t2g=="C1"){
    db_out=msigdbr(species = species_in, category = "C1") %>% 
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C2:BIOCARTA"){
    db_out=msigdbr(species = species_in, category = "C2", subcategory = "BIOCARTA") %>% 
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C2:KEGG"){
    db_out=msigdbr(species = species_in, category = "C2", subcategory = "KEGG") %>% 
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C2:REACTOME"){
    db_out=msigdbr(species = species_in, category = "C2", subcategory = "REACTOME") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C2:WIKIPATHWAYS"){
    db_out=msigdbr(species = species_in, category = "C2", subcategory = "WIKIPATHWAYS") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C5:MF"){
    db_out=msigdbr(species = species_in,  category = "C5", subcategory = "GO:MF") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C5:BP"){
    db_out=msigdbr(species = species_in,  category = "C5", subcategory = "GO:BP") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C5:CC"){
    db_out=msigdbr(species = species_in,  category = "C5", subcategory = "GO:CC") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="H"){
    db_out=msigdbr(species = species_in, category = "H") %>% 
      dplyr::select(gs_name,ensembl_gene)  
  } else{
    print("DB does not exist. Please review")
  }
  
  return(db_out)
}

# GSEA needs ENTREZids, this converts ENSEMBL to ENTREZ
CAPTURE_ENTREZids<-function(input_df){
  #input_df=deg_df_sub
  
  # search for ENTREZ by ENSEMBL
  gene_df_e <- bitr(input_df[,c("ENSEMBL")], fromType = "ENSEMBL",
                    toType = c("ENTREZID"),
                    OrgDb = get(species_db))
  
  # if any genes did not map, try by SYMBOL
  missing_df=subset(input_df,ENSEMBL %ni% gene_df_e$ENSEMBL)
  if (nrow(missing_df)>0){
    gene_df_s <- bitr(input_df[,c("SYMBOL")], fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = get(species_db))
    
    # merge dfs together
    gene_output=merge(gene_df_e,gene_df_s,all=TRUE,by="ENTREZID")
    
  }
  
  # remove duplicated ENSEMBL ID's that have a gene symbol
  dup_list=gene_output[duplicated(gene_output$ENSEMBL),]
  filter_list=dup_list[is.na(dup_list$SYMBOL),]$ENTREZID
  gene_output2=subset(gene_output, ENTREZID %ni% filter_list)
  
  # fill in missing gene symbols
  filter_list=gene_output2[is.na(gene_output2$SYMBOL),]
  for (i in rownames(filter_list)){
    eid=filter_list[i,"ENSEMBL"]
    #pull symbol from deg
    gene_output2[i,"SYMBOL"]=unique(subset(input_df,ENSEMBL==eid)$SYMBOL)
  }
  
  # fill in missing ENSEMBL
  filter_list=gene_output2[is.na(gene_output2$ENSEMBL),]
  for (i in rownames(filter_list)){
    sym=filter_list[i,"SYMBOL"]
    #pull symbol from deg
    gene_output2[i,"ENSEMBL"]=subset(input_df,SYMBOL==sym)$ENSEMBL[1]
  }
  
  #remove dups
  output_df=gene_output2[!duplicated(gene_output2[,c("ENSEMBL")]),]
  
  # rename cols to match deg
  colnames(output_df)=c("ENTREZID","ENSEMBL","SYMBOL")
  return(output_df)
}

CAPTURE_KEGGids<-function(input_df){
  #input_df=deg_anno_df
  
  # search for ENTREZ by ENSEMBL
  gene_df_e <- bitr_kegg(geneID = input_df$ENTREZID, 
                         fromType='ncbi-geneid', toType='kegg', organism=species_short)
  
  #remove dups
  output_df=gene_df_e[!duplicated(gene_df_e[,c("kegg")]),]
  
  # rename cols to match deg
  colnames(output_df)=c("ENTREZID","KEGG")
  return(output_df)
}

CAPTURE_ENSEMBLids<-function(input_df){
  # change result to list 
  entrez_list=""
  for (rowid in rownames(input_df)){
    entrez_list=append(entrez_list,
                       strsplit(input_df[rowid,"core_enrichment"],"[/]")[[1]])
  }
  entrez_list
  
  # deduplicate list
  entrez_list=entrez_list[!duplicated(entrez_list)]
  
  # create df
  gene_df=data.frame(entrez_list)
  colnames(gene_df)=c("ENTREZID")
  head(gene_df)
  
  # search for ENTREZ by ENSEMBL
  gene_df_e <- bitr(gene_df[,c("ENTREZID")], fromType = "ENTREZID",
                    toType = c("ENSEMBL"),
                    OrgDb = get(species_db))
  
  # sub the original ENTREZIDs with ENSEMBLids
  input_df$core_enrichment=gsub("/",";",input_df$core_enrichment)
  for (rowid in rownames(gene_df_e)){
    ezid=gene_df_e[rowid,"ENTREZID"]
    eid=gene_df_e[rowid,"ENSEMBL"]
    
    input_df$core_enrichment=gsub(ezid,eid,input_df$core_enrichment)
  }
  return(input_df)
}

CAPTURE_GENEids<-function(input_df){
  # input_df=t2g_df
  
  # create list of E'IDS
  eid_list=input_df$EIDS_included
  
  # collapse list, take unique values
  tmp_collapsed_list=paste0(unlist(eid_list),collapse=";")
  tmp_collapsed_list=as.list(unique(strsplit(tmp_collapsed_list, ";"))[[1]])
  length(tmp_collapsed_list)
  tmp_unique_list=tmp_collapsed_list[!duplicated(tmp_collapsed_list)]
  length(tmp_unique_list)  
  
  # create df
  tmp_df=as.data.frame(do.call(rbind,tmp_unique_list))
  
  # convert to SYMBOL
  df_annotated=bitr(tmp_df[,c("V1")], fromType = "ENSEMBL",
                    toType = c("SYMBOL"),
                    OrgDb = get(species_db))
  
  # replace EIDs with SYMBOL
  for (rid in rownames(df_annotated)){
    eid_list=gsub(df_annotated[rid,"ENSEMBL"],df_annotated[rid,"SYMBOL"],eid_list)
  }
  input_df$EIDS_included=eid_list
  
  return(input_df)
}

# run GSEA analysis
RUN_GSEA_ANALYSIS<-function(gsea_genelist,t2g,anno_db){
  # gsea_genelist=ranked_list; t2g=db_id; 
  # run GSEA
  if (t2g=="C2:KEGG"){
    if (species_in=="Mus musculus"){
      result=GSEA(geneList=gsea_genelist,
                  pvalueCutoff = padj_cutoff,
                  eps=0,
                  pAdjustMethod = "BH",
                  TERM2GENE = anno_db)
    } else{
      #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
      result=gseKEGG(geneList=gsea_genelist,
                     pvalueCutoff=padj_cutoff,
                     eps=0,
                     pAdjustMethod="BH", 
                     organism=species_short,
                     verbose=FALSE,
                     keyType="ncbi-geneid")
    }
  } else if (t2g=="C2:REACTOME"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
    result=gsePathway(gene=gsea_genelist, 
                      pvalueCutoff = padj_cutoff,
                      eps=0,
                      pAdjustMethod = "BH", 
                      verbose = FALSE)
    
  } else if (t2g=="C2:WIKIPATHWAYS"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/wikipathways-analysis.html
    result=gseWP(gene=gsea_genelist, 
                 pvalueCutoff = padj_cutoff,
                 eps=0,
                 pAdjustMethod = "BH",
                 organism=species_in)
    
  } else if ((t2g=="C5:MF") | (t2g=="C5:BP") | (t2g=="C5:CC")){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
    ont_id=strsplit(t2g,":")[[1]][2]
    result=gseGO(geneList=gsea_genelist,
                 pvalueCutof = padj_cutoff,
                 eps=0,
                 pAdjustMethod = "BH",
                 OrgDb= get(species_db),
                 ont=ont_id,
                 verbose= FALSE)
    
  } else if ((t2g=="C1") | (t2g=="C2:BIOCARTA") | (t2g=="H")){
    pulled_db=DB_LOOKUP(t2g)
    result=GSEA(geneList=gsea_genelist,
                pvalueCutoff = padj_cutoff,
                eps=0,
                pAdjustMethod = "BH",
                TERM2GENE = pulled_db)
  } else {
    print (paste0(t2g, ": DB selected is not valid"))
  }
  
  return(result)
}

# create plots
CREATE_FGSEA_PLOTS<-function(t2g,result_in,ranked_list,msigdbr_list){
  # run pathway analysis
  topPathwaysUp <- result_in[ES > 0][head(order(pval), n=5), pathway]
  topPathwaysDown <- result_in[ES < 0][head(order(pval), n=5), pathway]
  
  # generate plots, DT's for up and down
  path_name=str_wrap(gsub("_"," ",topPathwaysUp[1]), width = 30)
  p1 = plotEnrichment(msigdbr_list[[topPathwaysUp[1]]], ranked_list) + 
    labs(title=paste0("Up-regulated\n",path_name),cex=.5)
  path_name=str_wrap(gsub("_"," ",topPathwaysDown[1]), width = 30)
  p2 = plotEnrichment(msigdbr_list[[topPathwaysDown[1]]], ranked_list) + 
    labs(title=paste0("Down-regulated\n",path_name))
  title1=text_grob(paste0("Top pathways for ", t2g), size = 15, face = "bold")
  grid.arrange(
    p1,p2,
    top=title1,
    nrow=1
  )
}

CREATE_GSEA_PLOTS<-function(t2g,result_in,contrast_id){
  # t2g=db_id; result_in=gseaRes
  # create dot plots for all DB, ridgeplots for specific DB's
  if(nrow(result_in)==0){
    pf = ggparagraph( paste0("\n\n\n No Sig Results for GSEA:",t2g,
                             "\n-",contrast_id), 
                      size = 20, face = "bold")
  } else{
    p1 = dotplot(result_in,
                 title=paste0(contrast_id,"\nGSEA:",t2g),
                 font.size = 6, showCategory=2, split=".sign",orderBy="p.adjust") +
      facet_grid(.~.sign)
    p2 = ridgeplot(result_in, label_format = 30, showCategory = 4, orderBy="p.adjust") +
      labs(x = "Enrichment distribution for top 5 pathways") + 
      theme(text = element_text(size=6),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=5.5))
    pcol <- cowplot::plot_grid(
      p1 + theme(legend.position="none"),
      p2 + theme(legend.position="none"),
      nrow = 2
    )
    legend<-get_legend(p1)
    pf=cowplot::plot_grid(pcol,legend,rel_widths = c(3, .4))
  }
  print(pf)
  
  fpath=paste0(output_dir,"gsea_plots_",t2g,"_",contrast_id,".png")
  ggsave(fpath,pf)
}

# create output dt that summarized pathways
OUTPUT_GSEA_FGSEA_DF<-function(result_in,db_id,analysis_type){
  # result_in=fgseaRes; analysis_type="fgsea"
  # result_in=gseaRes; analysis_type="gsea"
  t2g_df=as.data.frame(result_in)
  
  # if there are pathways identified, analyze
  if (nrow(t2g_df)>0){
    # separate cols for gsea
    if (analysis_type=="gsea"){
      t2g_df=t2g_df %>% 
        separate(leading_edge,
                 c("percent_included","percent_list","percent_signal"),
                 sep=", ")
      # remove value= in cols
      t2g_df$percent_included=sub("tags=","",t2g_df$percent_included)
      t2g_df$percent_list=sub("list=","",t2g_df$percent_list)
      t2g_df$percent_signal=sub("signal=","",t2g_df$percent_signal)
      
      # convert the ENTREZids with ENSEMBLids
      if (db_id!="C2:KEGG"){
        t2g_df=CAPTURE_ENSEMBLids(t2g_df)
      }
      
      # determine length
      t2g_df$N_total=round(as.numeric(t2g_df$setSize)/(as.numeric(gsub("%","",
                                                                       t2g_df$percent_included))/100))
      colnames(t2g_df)
      
      # collapse EIDs
      for (rowid in rownames(t2g_df)){
        t2g_df[rowid,"core_enrichment"]=paste(t2g_df[rowid,"core_enrichment"][[1]],collapse = ';')
      }
      t2g_df$core_enrichment=gsub("\\/",";",t2g_df$core_enrichment)
      
      # subset cols
      select_cols=c("Description","p.adjust","enrichmentScore","NES",
                    "setSize","N_total","core_enrichment")
      t2g_df_sub=t2g_df[,select_cols]
      
      # rename cols
      select_cols=c("pathway","padj","ES","NES",
                    "N_included","N_total","EIDS_included")
      colnames(t2g_df_sub)=select_cols
      head(t2g_df_sub)
      
      # set values to round
      round_cols=c("padj","ES","NES")
      
    } else{
      # calculate number of genes in set
      t2g_df$set_total=sapply(t2g_df$leadingEdge,length)
      
      # collapse the EIDS
      for (rowid in rownames(t2g_df)){
        t2g_df[rowid,"EIDS_included"]=paste(unlist(t2g_df[rowid,"leadingEdge"]),
                                            collapse = ';')
      }
      
      # subset cols
      select_cols=c("pathway","padj","ES","NES","set_total","size","EIDS_included")
      t2g_df_sub=t2g_df[,select_cols]
      
      # rename cols
      select_cols=c("pathway","padj","ES","NES","N_included","N_total","EIDS_included")
      colnames(t2g_df_sub)=select_cols
      
      # set values to round
      round_cols=c("padj","ES","NES")
    }
    
    #round cols
    for (colid in round_cols){
      t2g_df_sub[,colid]=signif(t2g_df_sub[,colid], digits=3)
    }
    
    # replace EID's with gene symbols
    t2g_df_sub=CAPTURE_GENEids(t2g_df_sub)
    
    # set other params
    t2g_df_sub$annotation_db=db_id
    t2g_df_sub$significance="N"
    t2g_df_sub$significance[t2g_df_sub$padj<padj_cutoff]="Y"
    
    # calculate percent genes included
    t2g_df_sub$N_total=as.numeric(t2g_df_sub$N_total)
    t2g_df_sub$percent_included=(t2g_df_sub$N_included/t2g_df_sub$N_total)*100
    t2g_df_sub$percent_included=signif(t2g_df_sub$percent_included, digits=2)
    head(t2g_df_sub)
    
    # sort by p.adjust
    output_df=t2g_df_sub[order(t2g_df_sub$padj),]
    
    # reorder cols
    select_cols=c("significance","annotation_db","pathway","padj",
                  "ES","NES","N_included","N_total","percent_included",
                  "EIDS_included")
    t2g_df_sub=t2g_df_sub[,select_cols]
  } else{
    t2g_df_sub="No Significant Genes"
  }
  return(t2g_df_sub)
}

# main function
MAIN_FGSEA_GSEA_ANALYSIS<-function(contrast_id,db_list){
  # contrast_id="SH1_0pt5mM-SCR_0pt5mM"; db_id=db_list[4]
  print(paste0("--working on ", contrast_id))
  
  # read in deg df
  fpath=paste0(output_dir,"DEG_",contrast_id,".csv")
  deg_df=read.csv(fpath)
  
  # change inf to highest or lowest values
  # https://www.biostars.org/p/276783/
  deg_df$gsea_ranking_score[deg_df$gsea_ranking_score=="-Inf"]="toolow"
  deg_df$gsea_ranking_score[deg_df$gsea_ranking_score=="Inf"]="toohigh"
  
  gsea_num_list=as.numeric(deg_df$gsea_ranking_score)
  gsea_num_list=gsea_num_list[!is.na(gsea_num_list)]
  
  deg_df$gsea_ranking_score[deg_df$gsea_ranking_score=="toolow"]=min(gsea_num_list)
  deg_df$gsea_ranking_score[deg_df$gsea_ranking_score=="toohigh"]=max(gsea_num_list)
  deg_df$gsea_ranking_score=as.numeric(deg_df$gsea_ranking_score)
  range(deg_df$gsea_ranking_score)
  
  # prep df
  deg_df_sub=deg_df %>% 
    separate("ensembleID",c("ENSEMBL","ID"),sep="[.]")
  
  ## add entrezIDS for GSEA, KEGGids for GSEA
  gene_ref_db=CAPTURE_ENTREZids(deg_df_sub)
  deg_anno_df=merge(gene_ref_db,deg_df_sub,by="SYMBOL")
  gene_ref_db=CAPTURE_KEGGids(deg_anno_df)
  deg_anno_df=merge(gene_ref_db,deg_anno_df,by="ENTREZID")
  head(deg_anno_df)
  
  # sort df by score
  deg_df_sort=deg_anno_df[order(deg_anno_df$log2fc,decreasing=TRUE),]
  colnames(deg_df_sort)=gsub("[.]y","",colnames(deg_df_sort))
  head(deg_df_sort)
  
  # create ranked list
  ranked_list=deg_df_sort$log2fc
  
  # pull db and generate pathways
  merged_fgsea_df=data.frame()
  merged_gsea_df=data.frame()
  
  # db_id=db_list[1]
  for (db_id in db_list){
    print(paste0("----processing ", db_id))
    anno_db=DB_LOOKUP(db_id,species_in)
    msigdbr_list=split(anno_db$ensembl_gene,anno_db$gs_name)
    
    # run fgsea, plot, and merge df
    # https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
    # https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
    names(ranked_list)=as.character(deg_df_sort$ENSEMBL)
    fgseaRes <- fgsea(msigdbr_list, ranked_list,minSize  = minSize_gene_set)
    CREATE_FGSEA_PLOTS(db_id,fgseaRes,ranked_list,msigdbr_list)
    merged_fgsea_df=rbind(merged_fgsea_df,
                          OUTPUT_GSEA_FGSEA_DF(fgseaRes,db_id,analysis_type="fgsea"))
    head(merged_fgsea_df)
    
    # run gsea, plot and merge df
    if (db_id=="C2:KEGG"){
      names(ranked_list)=as.character(deg_df_sort$ENSEMBL)
    } else{
      names(ranked_list)=as.character(deg_df_sort$ENTREZID)
    }
    gseaRes <- RUN_GSEA_ANALYSIS(ranked_list,db_id,anno_db)
    CREATE_GSEA_PLOTS(db_id,gseaRes,contrast_id)
    merged_gsea_df=rbind(merged_gsea_df,
                         OUTPUT_GSEA_FGSEA_DF(gseaRes,db_id,analysis_type="gsea"))
    head(merged_gsea_df)
    
    # breakpoint for the select_pathway function versus standard analysis
    #https://stephenturner.github.io/deseq-to-fgsea/
  }
  
  # sort and output top 10 pathways
  output_fgsea_df=merged_fgsea_df[order(merged_fgsea_df$padj),][1:10,]
  output_gsea_df=merged_gsea_df[order(merged_gsea_df$padj),][1:10,]
  
  # clean output dt
  merged_fgsea_df$pathway=gsub("GOMF_||GOBP_||GOCC_||KEGG_","",merged_fgsea_df$pathway)
  merged_gsea_df$pathway=gsub("GOMF_||GOBP_||GOCC_||KEGG_","",merged_gsea_df$pathway)
  
  # save top pathway info, print DT
  caption_title=paste0("Top 10 pathways across all annotation databases for FGSEA ", contrast_id)
  DT::datatable(output_fgsea_df, extensions = 'Responsive', 
                caption=htmltools::tags$caption(paste0(caption_title),
                                                style="color:gray; font-size: 18px" ),
                rownames=F)
  fpath=paste0(output_dir,"GO_KEGG_pathways_fgsea_",contrast_id,".csv")
  write.table(merged_fgsea_df,fpath,sep=",",row.names = FALSE)
  
  caption_title=paste0("Top 10 pathways across all annotation databases for GSEA ", contrast_id)
  DT::datatable(output_gsea_df, extensions = 'Responsive', 
                caption=htmltools::tags$caption(paste0(caption_title),
                                                style="color:gray; font-size: 18px" ),
                rownames=F)
  fpath=paste0(output_dir,"GO_KEGG_pathways_gsea_",contrast_id,".csv")
  write.table(merged_gsea_df,fpath,sep=",",row.names = FALSE)
  
}

######################################################################
# functions heatmaps
######################################################################
# Overwrites the pheatmap defaults
DRAW_COLNAMES_45 <- function (coln, gaps, ...) {
  "Overwrites body of pheatmap:::draw_colnames, customizing it my liking"
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 25, gp = gpar(...))
  return(res)
}

SAVE_PHEATMAP_PDF <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# creates heatmap for multiple replicates
PLOT_HEAT_MAP<-function(df_in,show_names="ON",title_in="",cluster_by_rows="ON",fpath=""){
  #df_in=counts_matrix_complete;show_names="OFF";title_in="";cluster_by_rows="ON";fpath=fpath
  
  ####################
  # formatting
  #####################
  # Overwrite pheatmaps default draw_colnames with new version
  assignInNamespace(x="draw_colnames", value="DRAW_COLNAMES_45",ns=asNamespace("pheatmap")) 
  
  # Heatmap Color Gradients 
  paletteLength <- 1000
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  ####################
  # metadata
  ####################
  # Creating Dataframe to map samplenames to groups
  meta = groups_df
  groups <- data.frame(as.factor(meta$group))
  colnames(groups) <- "Groups"
  rownames(groups) <- meta$sampleid
  
  # Creating Group Column Annotation Colors
  columnColors <- c("lightpink","lightblue","orange","purple","red","green","darkblue","brown")
  names(columnColors) <- unique(groups$Groups)
  anno_colors <- list(Groups = columnColors)
  
  # set title
  if(title_in==""){
    title_in=paste0("Genes Included (N=",nrow(df_in),")")
  }
  ####################
  # function
  ####################
  if (show_names=="OFF" && cluster_by_rows=="ON"){
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 5, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = FALSE)
  } else if (show_names=="ON" && cluster_by_rows=="ON") {
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 5, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = TRUE)
  } else if (show_names=="OFF" && cluster_by_rows=="OFF") {
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 5, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,cluster_rows=F,annotation_colors = anno_colors, 
               show_rownames = FALSE)
  }
  
  # set a fpath if it's missing
  if (fpath==""){
    fpath=paste0(img_dir,"replicate_heatmap_",contrast_id,".pdf")
  }
  
  SAVE_PHEATMAP_PDF(p, fpath)
}

# main function to prep and run plot_heat_map
MAIN_REPLICATE_HEATMAPS<-function(sample_list,gene_list,scale_flag,name_flag){
  # read in counts matrix
  fpath=paste0(output_dir,"replicate_normalized_",csid,".csv")
  counts_matrix=read.csv(fpath,sep=",")
  counts_matrix=separate(counts_matrix,col="X",into=c("ENSEMBL","SYMBOL"),sep="[|]")
  nrow(counts_matrix)
  
  # filter for genes meeting threshold
  counts_matrix_filtered=subset(counts_matrix,sample_threshold_flag=="Y")
  nrow(counts_matrix_filtered)
  
  # filter for gene_list provided
  print("**Subsetting based on gene_list provided")
  counts_matrix_subset=subset(counts_matrix_filtered,SYMBOL %in% gene_list)
  missing_genes=gene_list[gene_list%ni%unique(counts_matrix_subset$SYMBOL)]
  if(length(missing_genes)>0){print(paste0("The following genes are missing:",missing_genes))}
  print("The following genes are included:")
  print(unique(counts_matrix_subset$SYMBOL))
  
  # pull rowname
  rownames(counts_matrix_subset)=make.unique(counts_matrix_subset$SYMBOL)
  head(counts_matrix_subset)
  
  # Subet for samples
  print("**Subsetting based on sample_list provided")
  counts_matrix_subset=counts_matrix_subset[,sample_list]
  
  # scale, if needed
  if (scale_flag=="ON"){
    print("**Performing scaling")
    
    # ztransform df
    counts_matrix_complete=t(scale(t(counts_matrix_subset)))
    
    # fix any nan or inf
    counts_matrix_complete[is.nan(counts_matrix_complete)] <- 0
    counts_matrix_complete[counts_matrix_complete=="Inf"] <- 0
    range(counts_matrix_complete)
    
    # set fpath
    fpath=paste0(img_dir,"replicate_heatmap_withscale_")
    
  } else{
    print("**No scaling will be performed")
    counts_matrix_complete=counts_matrix_subset
    
    # set fpath
    fpath=paste0(img_dir,"replicate_heatmap_withoutscale_")
  }
  
  # shorten col names
  colnames(counts_matrix_complete)=gsub("pt",".",colnames(counts_matrix_complete))
  
  print("**Generating heatmaps")
  # write out file
  fpath_f=paste0(fpath,"heatmap.csv")
  write.csv(counts_matrix_complete,fpath_f)
  
  # generate heatmap
  fpath=paste0(fpath,"heatmap.pdf")
  PLOT_HEAT_MAP(df_in=counts_matrix_complete,
                show_names=name_flag,
                title_in="",
                cluster_by_rows="ON",
                fpath=fpath)
}

######################################################################
# IPA Analysis
######################################################################
CREATE_ZSCORE_PLOTS<-function(top_df,mol_num,sample_id,exp_type){
  #top_df=ipa_genes_specific; mol_num=nrow(top_df)
  # create plot
  #https://stackoverflow.com/questions/40469757/change-the-shape-of-legend-key-for-geom-bar-in-ggplot2
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
  "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9",
  "#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  subPalette=cbPalette[length(top_df$Molecule.Type)]
  scale_factor <-2
  
  #create plot without legend
  p1 = top_df %>%
    mutate(Upstream.Regulator = fct_reorder(Upstream.Regulator, Activation.z.score)) %>%
    ggplot(aes(x = Upstream.Regulator))+
    geom_bar(aes(y = Activation.z.score,fill=Molecule.Type),colour="red", stat="identity")+
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor,color=""),
               shape=20, cex=8)+ # create dummy points so pvalue legend symbol appears
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor),
               shape=20, cex=8,color="blue")+
    scale_color_manual(name = "log10(pvalue) x \nActivation Z Score",
                       values = c(" " = "blue")) + # add legened for pvalue point
    scale_fill_brewer(direction = -1, palette = "Pastel2")+ # change palette for bars
    coord_flip()+ # flip x and y axis
    scale_y_continuous(sec.axis = sec_axis(~ .*1, #create a sec "x" axis for the log10(pvalue data)
                   labels = number_format(scale=scale_factor), 
                   name="log10(pvalue) x sign of Activation Z Score"))+
    xlab("Upstream Regulator")+
    ylab("Activation Z Score")+
    theme_bw() + # remove grey background
    theme(axis.line.x.top = element_line(color = "blue"),
          axis.ticks.x.top= element_line(color = "blue"), # change top x axis color to match pvalues
          axis.text.x.top = element_text(color = "blue"),
          axis.title.x.top=element_text(colour="blue"))+
    theme(axis.line.x.bottom = element_line(color = "red"),
          axis.ticks.x.bottom= element_line(color = "red"), # change bottom x axis color to match z scores
          axis.text.x.bottom = element_text(color = "red"),
          axis.title.x.bottom=element_text(colour="red"))+
    ggtitle(paste0(exp_type,": ",sample_id)) + # add title
    theme(plot.title = element_text(size=22, face="bold")) + #  increase title
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    labs(color="log10(pvalue) x \nActivation Z Score",
         fill="Activation Z Score")+ # change the order of the legend so log10 is on top
    guides(fill  = guide_legend(order = 2),color = guide_legend(order = 1))+
    theme(legend.position = "none") # remove legend
  
  #legend with FILL only
  p2 = top_df %>%
    mutate(Upstream.Regulator = fct_reorder(Upstream.Regulator, Activation.z.score)) %>%
    ggplot(aes(x = Upstream.Regulator))+
    geom_bar(aes(y = Activation.z.score,fill=Molecule.Type),colour="red", stat="identity")+
    scale_fill_brewer(direction = -1, palette = "Pastel2")+ # change palette for bars
    coord_flip()+ # flip x and y axis
    scale_y_continuous(sec.axis = sec_axis(~ .*1, #create a secondary "x" axis for the log10(pvalue data)
                   labels = number_format(scale=scale_factor), 
                   name="log10(pvalue) x sign of Activation Z Score"))+
    xlab("Upstream Regulator")+
    ylab("Activation Z Score")+
    theme_bw() + # remove grey background
    theme( axis.line.x.top = element_line(color = "blue"),
           axis.ticks.x.top= element_line(color = "blue"), # change top x axis color to match pvalues
           axis.text.x.top = element_text(color = "blue"),
           axis.title.x.top=element_text(colour="blue"))+
    theme( axis.line.x.bottom = element_line(color = "red"),
           axis.ticks.x.bottom= element_line(color = "red"), # change bottom x axis color to match z scores
           axis.text.x.bottom = element_text(color = "red"),
           axis.title.x.bottom=element_text(colour="red"))+
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    theme(legend.direction = "vertical", legend.box = "horizontal")+
    theme(legend.title = element_text(size=12,color="red",face="bold"))
  
  # legend for COLOR only
  p3 = top_df %>%
    mutate(Upstream.Regulator = fct_reorder(Upstream.Regulator, Activation.z.score)) %>%
    ggplot(aes(x = Upstream.Regulator))+
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor,color=""),
               shape=20, cex=8)+ # create dummy points so pvalue legend symbol appears
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor),
               shape=20, cex=8,color="blue")+
    scale_color_manual(name = "log10(pvalue) x \nsign of Activation Z Score",
                       values = c(" " = "blue")) + # add legened for pvalue point
    coord_flip()+ # flip x and y axis
    scale_y_continuous(sec.axis = sec_axis(~ .*1, #create a secondary "x" axis for the log10(pvalue data)
                   labels = number_format(scale=scale_factor), 
                   name="log10(pvalue) x sign of Activation Z Score"))+
    theme_bw() + # remove grey background
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    labs(color="log10(pvalue) x \nActivation Z Score",
         fill="Activation Z Score")+ # change the order of the legend so log10 is on top
    theme(legend.direction = "vertical", legend.box = "horizontal")+
    theme(legend.title = element_text(size=12,color="blue",face="bold"))
  
  # get two legends
  leg2=get_legend(p2)
  leg3=get_legend(p3)
  
  # combine legends
  leg32 = plot_grid(leg3, leg2,
                    nrow = 2)
  
  # create final plot
  final_p = plot_grid(p1,
                      leg32,
                      nrow = 1,
                      align = "h",
                      axis = "t",
                      rel_widths = c(1, 0.3))
  
  
  # save and print
  ext=gsub(" ","_",exp_type)
  fpath=paste0(img_dir,sample_id,"_",ext,"_molecules.png")
  ggsave(fpath,final_p)
  print(final_p)
}

MAIN_ZSCORE_PLOTS<-function(sample_id,gene_list,exp_type){
  # read in IPA file
  ## File generated from UPSTREAM analysis
  ## Select and download Casual Networks
  fpath=paste0(ipa_dir,sample_id,"_upstream.txt")
  ipa_df=read.csv(fpath,sep="\t",skip=2)
  head(ipa_df)
  
  # remove values without z scores, pvalues
  ipa_df_filt=ipa_df[complete.cases(ipa_df$Activation.z.score),]
  ipa_df_filt=ipa_df_filt[complete.cases(ipa_df_filt$p.value.of.overlap),]
  
  # subset for specific gene list
  ipa_genes_specific=subset(ipa_df_filt,Upstream.Regulator %in% gene_list)
  missing_genes=gene_list[gene_list%ni%unique(ipa_genes_specific$Upstream.Regulator)]
  if(length(missing_genes)>0){print(paste0("The following genes were not found as a Master regulator: ",
                                           missing_genes))}
  
  # add calc cols
  ipa_genes_specific$signA=ipa_genes_specific$Activation.z.score<0
  ipa_genes_specific$signA=gsub("TRUE",1,ipa_genes_specific$signA)
  ipa_genes_specific$signA=gsub("FALSE",-1,ipa_genes_specific$signA)
  ipa_genes_specific$logpvalue=log(ipa_genes_specific$p.value.of.overlap,
                                   100)*as.numeric(ipa_genes_specific$signA)
  
  # create plot
  CREATE_ZSCORE_PLOTS(ipa_genes_specific,nrow(ipa_genes_specific),sample_id,exp_type)
  
  #create DT
  select_cols=c("Molecule.Type","Upstream.Regulator","Activation.z.score","p.value.of.overlap","logpvalue")
  out_df=ipa_genes_specific[,select_cols]
  out_df$logpvalue=signif(out_df$logpvalue, 4)
  colnames(out_df)=c("Molecule.Type","Upstream.Regulator","Activation.z.score","pvalue","log10.pvalue")  
  DT::datatable(out_df,caption = paste0(exp_type,": ",sample_id))
}

######################################################################
# FGSEA/GSEA Followup
########################################################################
PLOT_FGSEA_GSEA_PATHWAYS<-function(sample_list,pathway_list,listid){
  
  # pull fgsea and gsea results
  go_df=data.frame()
  analysis_list=c("gsea","fgsea")
  for (aid in analysis_list){
    for (sid in sample_list){
      fpath=paste0(GO_dir,"GO_KEGG_pathways_",aid,"_",sid,"-SCR_0pt5mM.csv")
      tmp_df=read.csv(fpath)
      
      # add sampleid
      tmp_df$sample_id=gsub("pt",".",sid)
      
      # convert fsea path names to gsea path names
      tmp_df$pathway=gsub("_"," ",tolower(tmp_df$pathway))
      
      if(nrow(go_df)==0){
        go_df=tmp_df
      } else{
        go_df=rbind(go_df,tmp_df)
      }
    }  
  }
  head(go_df)
  
  # filter for pathway list
  go_filt=subset(go_df,pathway %in% pathway_list)[,c("pathway","NES","sample_id")]
  head(go_filt)
  
  # sort
  go_filt = go_filt %>% arrange(desc(pathway))
  go_filt$pathway <- factor(go_filt$pathway, levels = unique(go_filt$pathway))
  head(go_filt)
  
  # check pathways
  missing_paths=pathway_list[pathway_list%ni%unique(go_filt$pathway)]
  if (length(missing_paths)>0){print(paste0("There are missing pathways: ", missing_paths))}
  
  # plot
  p = ggplot(data = go_filt, aes(x = pathway, y = NES, fill = sample_id)) +
    geom_bar(stat = "identity", position = position_dodge(),width=0.5)  +
    labs(x = "Pathway", y = "NES") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_text(face="bold", size = 12),
          axis.title.y = element_text(face="bold", size = 12),
          legend.title = element_text(face="bold", size = 10)) +
    theme_bw() +
    theme(aspect.ratio = 1.5)+  
    coord_flip()
  p
  
  # save plots
  sample_ext=gsub("_0.5mM","",paste0(unique(go_filt$sample_id),collapse="_"))
  fpath=paste0(img_dir,"pathway_",sample_ext,"_",listid,".png")
  ggsave(fpath,p)
  print(p)
}

######################################################################
# functions secondary pathway analysis
######################################################################
# creates heatmap formatted df with only log2fc values
create_heatmap_df<-function(df_in,extra_filter=""){
  df_out=dplyr::select(df_in,contains("log"))
  
  # cleanup col names
  output_list=sub('--log2FoldChange', '',colnames(df_out))
  cntrl=strsplit(output_list,"_vs_")
  for (ct in cntrl){
    output_list=sub(ct[2], '',output_list)
    output_list=sub('_vs_', '',output_list)
  }
  
  # use extra filter
  if (extra_filter != ""){
    output_list=sub(extra_filter, '',output_list)
  }
  
  colnames(df_out)=output_list
  return(df_out)
}

# creates df of log2fc and pvalues
create_output_df<-function(df_in,n_in,extra_filter=""){
  output_list=colnames(df_in)
  
  # split treat from control
  cntrl=strsplit(output_list,"_vs_")
  for (ct in cntrl){
    metric=strsplit(ct[2],"--")[[1]][2]
    output_list=sub(ct[2], paste0("_",metric),output_list)
    output_list=sub('_vs_', '',output_list)
  }
  
  # cleanup col names
  output_list=sub('log2FoldChange', 'log2FC',output_list)
  
  # use extra filter
  extra_filter="_without_IFNb"
  if (extra_filter != ""){
    output_list=sub(extra_filter, '',output_list)
  }
  
  colnames(df_in)=output_list
  
  # round df
  for (colid in colnames(df_in)){
    df_in[,colid]=signif(df_in[,colid], digits=3)
  }
  
  # if number of pathways is less than n_in, adjust
  if (nrow(df_in)<n_in){
    n_in=nrow(df_in)
  }
  
  # print out table
  caption_title=paste0("Expression values for ",n_in," genes")
  p = DT::datatable(df_in, extensions = 'Responsive', 
                    caption=htmltools::tags$caption(paste0(caption_title) ,
                                                    style="color:gray; font-size: 18px" ))
  return(p)
}

# creates heatmap
generate_heat_map<-function(df_in,show_names="ON",title_in="",cluster_by_rows="ON",fpath=""){
  
  ####################
  # formatting
  #####################
  # Overwrite pheatmaps default draw_colnames with new version
  assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap")) 
  
  # Heatmap Color Gradients 
  paletteLength <- 1000
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  ####################
  # metadata
  ####################
  # Creating Dataframe to map samplenames to groups
  meta = groups_df
  groups <- data.frame(as.factor(meta$group))
  colnames(groups) <- "Groups"
  rownames(groups) <- meta$sampleid
  
  # Creating Group Column Annotation Colors
  columnColors <- c("lightpink","lightblue","orange","purple","red","green")
  names(columnColors) <- unique(groups$Groups)
  anno_colors <- list(Groups = columnColors)
  
  # set title
  if(title_in==""){
    title_in=paste0("Significant Genes (N=",nrow(df_in),")")
  }
  ####################
  # function
  ####################
  if (show_names=="OFF" && cluster_by_rows=="ON"){
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 7, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = FALSE)
  } else if (show_names=="ON" && cluster_by_rows=="ON") {
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 5, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = TRUE)
  } else if (show_names=="OFF" && cluster_by_rows=="OFF") {
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 7, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,cluster_rows=F,annotation_colors = anno_colors, 
               show_rownames = FALSE)
  }
  if (fpath==""){
    fpath=paste0(img_dir,"heatmap_",contrast_id,".pdf")
  }
  save_pheatmap_pdf(p, fpath)
}

# creates heatmap for multiple replicates
generate_replicate_heatmaps<-function(contrast_id,peak_type,scale,gene_list_subset="",sample_subset=""){
  # split contrast
  contrast_1=strsplit(contrast_id,"_vs_")[[1]][1]
  contrast_2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  # define extension
  contrast_extension=paste0(contrast_id,"__",dedup_status,"__",peak_type,"_peaks.bed")
  
  # read in counts matrix
  fpath=paste0(output_car_dir,"DESEQ_norm_counts_",contrast_id,".csv")
  counts_matrix=read.csv(fpath)
  head(counts_matrix)
  colnames(counts_matrix)=c("peakID",colnames(counts_matrix)[2:ncol(counts_matrix)])
  
  # replace X at the beginning of cols
  colnames(counts_matrix)=gsub("X","",colnames(counts_matrix))
  head(counts_matrix)
  
  # subset for samples, if needed
  if (length(sample_subset)>1){
    sample_list=sample_subset
    counts_matrix=counts_matrix[,c("peakID",sample_list)]
  } else{
    sample_list=subset(groups_df,group%in% c(contrast_1,contrast_2))$sampleid
  }
  
  # add annotation information
  fpath=paste0(output_car_dir,"peak_annotation_",contrast_id,".csv")
  peak_df=read.csv(fpath)[,c("peakID","SYMBOL","shortAnno")]
  head(peak_df)
  counts_matrix_anno=merge(counts_matrix,peak_df,by="peakID")
  head(counts_matrix_anno)
  
  # pull rowname
  rownames(counts_matrix_anno)=counts_matrix_anno$peakID
  counts_matrix=counts_matrix_anno[,2:ncol(counts_matrix_anno)]
  head(counts_matrix_anno)
  
  # read in sig gene list
  fpath=paste0(output_dir,"deeptools_sig_",contrast_id,".csv")
  deeptools_gene_list=read.csv(fpath)$x
  head(deeptools_gene_list)
  
  # subset for gene list
  counts_matrix_subset=subset(counts_matrix_anno,SYMBOL %in% deeptools_gene_list)
  
  # subset for non-distal annotated peaks
  counts_matrix_subset=subset(counts_matrix_subset,shortAnno!="Distal")
  
  # if needed, subset for immune gene list
  if (length(gene_list_subset)>1){
    print("Subsetting based on gene_list provided")
    counts_matrix_subset=subset(counts_matrix_subset,SYMBOL %in% gene_list_subset)
    
    print(paste0("Total number of peaks: ", nrow(counts_matrix_subset)))
    print(paste0("Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
    
    # keep only one instance per symbol
    counts_matrix_subset=counts_matrix_subset[!duplicated(counts_matrix_subset$SYMBOL), ]
    print(paste0("Total number of peaks after subsetting only unique genes: ", nrow(counts_matrix_subset)))
    
    # rename rows for printing
    rownames(counts_matrix_subset)=make.unique(paste0(counts_matrix_subset$SYMBOL," (",counts_matrix_subset$shortAnno,")"))
  } else{
    print(paste0("Total number of peaks: ", nrow(counts_matrix_subset)))
    print(paste0("Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
  }
  
  # scale if necessary
  if (scale=="ON"){
    print("performing scaling")
    # ztransform df
    counts_matrix_complete=t(scale(t(counts_matrix_subset[,sample_list])))
    (head(counts_matrix_complete))
  } else{
    print("no scaling will be performed")
    (head(counts_matrix_subset[,sample_list]))
    counts_matrix_complete=counts_matrix_subset[,sample_list]
  }
  
  # define output name, run heatmap
  if (scale=="ON"){
    fpath=paste0(img_dir,"heatmap_withscale_")
  } else{
    fpath=paste0(img_dir,"heatmap_withoutscale_")
  }
  if (length(gene_list_subset)>1){
    fpath=paste0(fpath,"subset_",gene_list_name,"_",contrast_id,".pdf")
    print(fpath)
    
    # generate heatmap
    generate_heat_map(counts_matrix_complete,
                      show_names="ON",
                      title_in="",
                      cluster_by_rows="ON",
                      fpath=fpath)
  } else{
    fpath=paste0(fpath,contrast_id,".pdf")
    
    # generate heatmap
    generate_heat_map(counts_matrix_complete,
                      show_names="OFF",
                      title_in="",
                      cluster_by_rows="ON",
                      fpath=fpath)
  }
}

# creates heatmap for multiple replicates in RNASeq data
generate_replicate_heatmaps_rna<-function(contrast_id,scale_flag,gene_list_subset="",sample_subset=""){
  # split contrast
  contrast_1=strsplit(contrast_id,"_vs_")[[1]][1]
  contrast_2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  # read in counts matrix
  counts_matrix=read.csv(paste0(input_dir,"DEG_ALL/RawCountFile_RSEM_genes_filtered.txt"),sep="\t")
  head(counts_matrix)
  
  # split EID and SYMBOL
  counts_matrix=separate(counts_matrix,col="symbol",into=c("EID","SYMBOL"),sep="[|]")
  
  # keep one gene
  counts_matrix=counts_matrix %>% distinct(SYMBOL,.keep_all=TRUE)
  head(counts_matrix)
  
  # subset for samples, if needed
  if (length(sample_subset)>1){
    sample_list=sample_subset
    counts_matrix=counts_matrix[,c("SYMBOL",sample_list)]
    print(head(counts_matrix))
  } else{
    sample_list=subset(groups_df,group%in% c(contrast_1,contrast_2))$sampleid
  }
  
  # pull rowname
  rownames(counts_matrix)=make.unique(counts_matrix$SYMBOL)
  head(counts_matrix)
  
  # read in sig gene list
  deeptools_gene_list=read.csv(car_deeptools)$x
  head(deeptools_gene_list)
  
  # subset for gene list
  counts_matrix_subset=subset(counts_matrix,SYMBOL %in% deeptools_gene_list)
  
  # if needed, subset for immune gene list
  if (length(gene_list_subset)>1){
    print("Subsetting based on gene_list provided")
    counts_matrix_subset=subset(counts_matrix_subset,SYMBOL %in% gene_list_subset)
    
    # keep only one instance per symbol
    #counts_matrix_subset=counts_matrix_subset[!duplicated(counts_matrix_subset$SYMBOL), ]
    
    # rename rows for printing
    rownames(counts_matrix_subset)=make.unique(paste0(counts_matrix_subset$SYMBOL))
  }
  
  # scale if necessary
  if (scale_flag=="ON"){
    print("performing scaling")
    # ztransform df
    x=t(counts_matrix_subset[,sample_list])
    counts_matrix_complete=t(scale(x))
    print(head(counts_matrix_complete))
  } else{
    print("no scaling will be performed")
    (head(counts_matrix_subset[,sample_list]))
    counts_matrix_complete=counts_matrix_subset[,sample_list]
  }
  
  # remove na values
  counts_matrix_complete=counts_matrix_complete[complete.cases(counts_matrix_complete), ]
  
  #output counts matrix
  if (length(gene_list_subset)>1){
    fpath=paste0(output_rna_dir,"heatmap_normcounts_",contrast_id,".txt")
  } else{
    fpath=paste0(output_rna_dir,"heatmap_normcounts_subset_",contrast_id,".txt")
  }
  write.csv(counts_matrix_complete,fpath)
  
  # define output name, run heatmap
  if (scale_flag=="ON"){
    fpath=paste0(img_dir,"heatmap_withscale_")
  } else{
    fpath=paste0(img_dir,"heatmap_withoutscale_")
  }
  if (length(gene_list_subset)>1){
    fpath=paste0(fpath,"subset_",gene_list_name,"_",contrast_id,".pdf")
    
    # generate heatmap
    generate_heat_map(counts_matrix_complete,
                      show_names="ON",
                      title_in="",
                      cluster_by_rows="ON",
                      fpath=fpath)
  } else{
    fpath=paste0(fpath,contrast_id,".pdf")
    
    # generate heatmap
    generate_heat_map(counts_matrix_complete,
                      show_names="OFF",
                      title_in="",
                      cluster_by_rows="ON",
                      fpath=fpath)
  }
}

#create heatmap and DT
generate_heat_map_select<-function(select_deg,contras){
  # prep df for heatmap generation
  rownames(select_deg)=select_deg$ensid_gene
  log_name=paste0(contras[1],"_vs_",contras[2],"--log2FoldChange")
  p_name=paste0(contras[1],"_vs_",contras[2],"--padj")
  select_deg=select_deg[,c("log2fc","fdr")]
  colnames(select_deg)=c(log_name,p_name)
  
  # create heatmap df of sig genes in pathway list
  sig_df=subset(select_deg, get(log_name) > log2fc_cutoff | get(log_name) < -log2fc_cutoff)
  sig_df=subset(sig_df, get(p_name) < padj_cutoff )
  heat_df=create_heatmap_df(sig_df,"")
  
  # run functions
  generate_heat_map(heat_df)
  p = create_output_df(select_deg,nrow(select_deg),"")
  return(p)
}

# create gseaplot
gsea_plus_plots_select<-function(deg,t2g,path_id,contras){
  
  # create GSEA genelist
  gsea_genelist=PREP_GSEA_GL(deg,t2g=t2g)
  
  # run GSEA, save plots
  result=gsea_plus_plot(gl=gsea_genelist,t2g=t2g,contrast_in=contras,select_flag="ON")
  
  #get rowname of pathway
  result_df=as.data.frame(result)
  rownames(result_df) <- NULL
  row_id=as.numeric(rownames(subset(result_df,ID==path_id))[[1]])
  path_desc=subset(result_df,ID==path_id)$Description[[1]]
  
  # create plot
  p1=gseaplot(result, by = "all", title = paste0(path_id,": ",path_desc),
              geneSetID =row_id)
  pf=p1&theme(text = element_text(size=8),
              axis.text.x = element_text(size=8),
              axis.text.y = element_text(size=8),
              axis.text =element_text(size=8),
              axis.title=element_text(size=8,face="bold")) 
  print(pf)
}

main_selectpath_function<-function(cntrl_in,treat_in,type_in,t2g,path_id){
  t2g="C1"
  type_in="GSEA"
  path_id="chr18q21"
  cntrl_in=cntrl
  treat_in=treatment
  
  # set contrast
  contras=c(treat_in,cntrl_in)
  
  # read in datatable created during GSEA/ORA plotting, get list of gene names in pathway selected
  ttl_abbrev=sub(" ","_",sub(":","_",t2g))
  fpath=paste0(output_rna_dir,type_in,"_",contras[1],"-",contras[2],"_table_",ttl_abbrev,".txt")
  path_df=read.csv(fpath,sep="\t")
  
  #check pathway exists
  if (nrow(subset(path_df,ID==path_id))==0){
    stop(paste0("The selected pathway (",path_id,") does not exist in the annotation database (",t2g,
                "). Please select a valid combination"))
  }
  
  # create gene list from pathway
  genes_in_pathway=strsplit(subset(path_df,ID==path_id)$core_enrichment,"/")[[1]]
  
  # read in created deg
  fpath=paste0(output_rna_dir,"DESeq2_",contras[1],"-",contras[2],"_DEG_allgenes.txt")
  deg=read.csv(fpath,header=TRUE,sep="\t")
  
  # subset deg for genes,
  #convert ENTREZID if necessary
  special_dbs=c("C1","C2:BIOCARTA","H")
  if (t2g %ni% special_dbs){
    gene_ref_db = CAPTURE_ENTREZids(deg)
    genes_in_pathway=subset(gene_ref_db,ENTREZID %in% genes_in_pathway)$gene
  } 
  select_deg=subset(deg, gene %in% genes_in_pathway)
  
  # create heatmaps
  p = generate_heat_map_select(select_deg,contras)
  
  # create gseaPlot
  gsea_plus_plots_select(deg,t2g,path_id,contras)
  
  # return DT
  return(p)
}