###########################################################################
## 
###########################################################################
handle_packages<-function(list.of.packages){
  for (pkg in list.of.packages){
    if ( !(pkg %in% rownames(installed.packages()))){
      print(paste0("Installing: ", pkg))
      tryCatch(
        {
          install.packages(pkg,local = FALSE,ask=FALSE,update=FALSE)
        }, error={
          BiocManager::install(pkg,ask=FALSE,update=FALSE)
        }
      )
    }
    print(paste0("Loading: ",pkg))
    invisible(lapply(pkg, library, character.only = TRUE))
  }
}

shorten_names<-function(input_df){
  colnames(input_df)=gsub(paste0("_",unique(rna_meta$harvestID)[1],"_S[0-9]*"),"",colnames(input_df))
  colnames(input_df)=gsub(paste0("_",unique(rna_meta$harvestID)[2],"_S[0-9]*"),"",colnames(input_df))
  colnames(input_df)=gsub(paste0("_",unique(rna_meta$harvestID)[3],"_S[0-9]*"),"",colnames(input_df))
  return(input_df)
}

###########################################################################
## COUNTS
###########################################################################
counts_dat_to_matrix <- function(counts_tbl) {
  gene_id <- GeneName <- NULL
  counts_dat <- counts_tbl %>%
    # deseq2 requires integer counts
    dplyr::mutate(dplyr::across(
      dplyr::where(is.numeric),
      \(x) as.integer(round(x, 0))
    )) %>%
    as.data.frame()
  row.names(counts_dat) <- counts_dat %>% dplyr::pull("gene_id")
  # convert counts tibble to matrix
  counts_mat <- counts_dat %>%
    dplyr::select(-c(gene_id, GeneName)) %>%
    as.matrix()
  return(counts_mat)
}

run_deseq2 <- function(counts,sampleinfo) {
  DESeqDataSetFromMatrix(countData = as.matrix(counts),
                         colData = sampleinfo[,c("sampleid","group")],
                         design = ~ group)
  return(renee_ds)
}

filter_low_counts <- function(counts_dat, min_counts = 0) {
  gene_id <- count <- count_sum <- NULL
  
  # find genes passing filtering
  genes_above_threshold <- counts_dat %>%
    tidyr::pivot_longer(!c("gene_id", "GeneName"),
                        names_to = "sample_id", values_to = "count"
    ) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(count_sum = sum(count)) %>%
    dplyr::filter(count_sum >= min_counts) %>%
    dplyr::pull(gene_id)
  
  return( counts_dat %>%
            dplyr::filter(gene_id %in% (genes_above_threshold)))
}

###########################################################################
## PCA
###########################################################################
PLOT_RLE_PCA<-function(input_data,input_title,inputcol){
  tmp_df=shorten_names(input_data)
  
  #Plot results
  par(mfrow=c(1, 2), oma=c(3, 2, 0, 0)+0.1)
  EDASeq::plotRLE(as.matrix(tmp_df), outline=FALSE, ylim=c(-.5, .5), col=colors[inputcol],las=2, cex.axis = .8)
  EDASeq::plotPCA(as.matrix(tmp_df), col=colors[inputcol],las=2, cex.axis = .8)
  
  mtext(input_title, side=2, outer=TRUE, adj=0)  
}  

PLOT_PCA<-function(dds,samplemeta,exclusion_list="",contrast_id,flag="raw",select_col=""){
  try(if (exclusion_list==""){exclusion_list=samplemeta$filename})
  
  if (flag == "raw"){
    # default for nsub is 1000, if there are less than 1000 rows this will error
    if (nrow(dds) < 1000){
      rld <- vst(dds,nsub=nrow(dds))
    } else{
      rld <- vst(dds)
    }
    assayrld = as.data.frame(assay(rld))
    
    plottitle=paste0("Raw Data\n",contrast_id)
  } else{
    assayrld = as.data.frame(counts(dds, normalize=TRUE))
    #assayrld=assayrld[,colnames(assayrld) %in%exclusion_list]
    plottitle=paste0("Normalized Data\n",contrast_id)
  }
  
  # analysis of variance
  assayrld$row_variance = rowVars(as.matrix(assayrld))
  assayrld = arrange(assayrld,desc(row_variance))
  zero_variance_rows=assayrld$row_variance<1e-5
  assayrld$row_variance = NULL
  assayrld = assayrld[!zero_variance_rows,]
  if (nrow(assayrld) > 500){
    assayrld=assayrld[1:500,]
  }
  
  #plot PCA
  pca=prcomp(t(assayrld),scale. = T)
  m.pc1 = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
  m.pc2 = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
  xlab=paste0("PC1(",m.pc1,"%)")
  ylab=paste0("PC2(",m.pc2,"%)")
  p = ggplot(pca$x,aes(x=PC1,y=PC2,label=rownames(pca$x)))+
    geom_point(col=as.factor(as.numeric(as.factor(samplemeta[,select_col]))+1))+
    xlab(xlab)+ylab(ylab)+
    ggtitle(plottitle)+
    geom_text_repel(max.overlaps = 20,size=2)+
    theme_light()
  print(p)
  
  # save final image
  fpath=paste0(img_dir,"pcoa_",gsub("-","",contrast_id),".png")
  ggsave(fpath,p)
}
###########################################################################
## DEG
###########################################################################
db_lookup<-function(t2g){
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

#set genelist for GSEA
deg2geneList<-function(input_df,t2g){
  # create genelist
  gsea_genelist=input_df$log2FoldChange
  
  if ((t2g=="H") | (t2g=="C1") | (t2g=="C2:BIOCARTA")) {
    names(gsea_genelist)=as.character(input_df$ENSEMBL)
  } else{
    input_df=capture_entrezids(input_df)
    names(gsea_genelist)=as.character(input_df$ENTREZID)
  }
  gsea_genelist=sort(gsea_genelist,decreasing=TRUE)
  return(gsea_genelist)
}

# create df of ENSEMBL, ENTREZ and SYMBOL
## different annotations are required for different db's
capture_entrezids<-function(input_df){
  
  # search for ENTREZ by ENSEMBL
  colnames(input_df)=gsub("GeneName","SYMBOL",colnames(input_df))
  gene_df_e <- clusterProfiler::bitr(input_df[,c("ENSEMBL")], 
                                     fromType = "ENSEMBL",
                                     toType = c("ENTREZID"),
                                     OrgDb="org.Mm.eg.db")
  
  # if any genes did not map, try by SYMBOL
  missing_df=subset(input_df,ENSEMBL %ni% gene_df_e$ENSEMBL)
  if (nrow(missing_df)>0){
    gene_df_s <- clusterProfiler::bitr(input_df[,c("SYMBOL")], fromType = "SYMBOL",
                                       toType = c("ENTREZID"),
                                       OrgDb = "org.Mm.eg.db")
    
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
  colnames(output_df)=c("ENTREZID","ENSEMBL","gene")
  return(output_df)
}

## Main function
RUN_CLUSTERPROFILER<-function(t2g){
  print(paste0("--",t2g))
  
  # run GSEA
  result=data.frame()
  if (t2g=="C2:KEGG"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
    result=clusterProfiler::gseKEGG(geneList=gl,
                                    pvalueCutoff=padj_cutoff,
                                    eps=0,
                                    pAdjustMethod="BH", 
                                    organism=species_short,
                                    verbose=FALSE)
    
  } else if (t2g=="C2:REACTOME"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
    result=ReactomePA::gsePathway(gene=gl, 
                                  pvalueCutoff = padj_cutoff,
                                  eps=0,
                                  pAdjustMethod = "BH", 
                                  verbose = FALSE)
    
  } else if (t2g=="C2:WIKIPATHWAYS"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/wikipathways-analysis.html
    result=clusterProfiler::gseWP(gene=gl, 
                                  pvalueCutoff = padj_cutoff,
                                  eps=0,
                                  pAdjustMethod = "BH",
                                  organism=species_in)
    
  } else if ((t2g=="C5:MF") | (t2g=="C5:BP") | (t2g=="C5:CC")){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
    ont_id=strsplit(t2g,":")[[1]][2]
    result=clusterProfiler::gseGO(geneList=gl,
                                  pvalueCutof = padj_cutoff,
                                  eps=0,
                                  pAdjustMethod = "BH",
                                  OrgDb= get(org_db),
                                  ont=ont_id,
                                  verbose= FALSE)
    
  } else if ((t2g=="C1") | (t2g=="C2:BIOCARTA") | (t2g=="H")){
    pulled_db=db_lookup(t2g)
    result=clusterProfiler::GSEA(geneList=gl,
                                 pvalueCutoff = padj_cutoff,
                                 eps=0,
                                 pAdjustMethod = "BH",
                                 TERM2GENE = pulled_db)
  } else {
    print (paste0(t2g, ": DB selected is not valid"))
  }
  
  # create dot plots, ridgeplots for specific DB's
  if(nrow(result)==0){
    pf = ggpubr ::ggparagraph( paste0("\n\n\n No Sig Results for ",t2g), 
                               size = 20, face = "bold")
  } else{
    resultdf=as.data.frame(result)
    fpath=sub(":","_",paste0(output_dir,"GSEA_",t2g,".txt"))
    write.table(resultdf,file=fpath,quote=FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
    
    p1 = dotplot(result,
                 title=paste0("\nGSEA:",t2g),
                 font.size = 6, showCategory=2, split=".sign",orderBy="p.adjust") +
      facet_grid(.~.sign)
    p2 = ridgeplot(result, label_format = 30, showCategory = 4, orderBy="p.adjust") +
      labs(x = "Enrichment distribution for top 5 pathways") + 
      theme(text = element_text(size=6),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=5.5))
    pcol <- cowplot::plot_grid(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               nrow = 2)
    legend<-get_legend(p1)
    pf=cowplot::plot_grid(pcol,legend,rel_widths = c(3, .4))
    print(pf)
    
    fpath=sub(":","_",paste0(output_dir,"GSEA_",t2g,".png"))
    ggsave(fpath,pf)
  }
}  
  
###########################################################################
## IPA
###########################################################################
create_top_list<-function(sub_df,mol_num){
  # find top values for each molecule type
  for (mol_id in mol_list){
    # subset for molecule
    sub_df2=subset(sub_df,Molecule.Type==mol_id)
    
    # sort by p value
    sorted_df=sub_df2[order(sub_df2$p.value.of.overlap),]
    
    # add top value to df
    if (exists("top_df")){
      top_df=full_join(top_df,sorted_df[mol_num,])
    } else{
      top_df=sorted_df[mol_num,]
    }
  }
  # add log10 of pvalue
  top_df$signA=top_df$Activation.z.score<0
  top_df$signA=gsub("TRUE",1,top_df$signA)
  top_df$signA=gsub("FALSE",-1,top_df$signA)
  top_df$logpvalue=log(top_df$p.value.of.overlap,100)*as.numeric(top_df$signA)
  return(top_df)
}

create_zscore_pvalue_plot<-function(top_df,mol_num){
  # create plot
  #https://stackoverflow.com/questions/40469757/change-the-shape-of-legend-key-for-geom-bar-in-ggplot2
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                 "#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  subPalette=cbPalette[length(top_df$Molecule.Type)]
  scale_factor <-2
  
  # clear NA's
  top_df=top_df[complete.cases(top_df$Activation.z.score),]
  
  #create plot without legend
  p1 = top_df %>%
    mutate(Upstream.Regulator = fct_reorder(Upstream.Regulator, Activation.z.score)) %>%
    ggplot(aes(x = Upstream.Regulator))+
    geom_bar(aes(y = Activation.z.score,fill=Molecule.Type),colour="red", stat="identity")+
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor,color=""),shape=20, cex=8)+ # create dummy points so pvalue legend symbol appears
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor),shape=20, cex=8,color="blue")+
    scale_color_manual(name = "log10(pvalue) x \nActivation Z Score",values = c(" " = "blue")) + # add legened for pvalue point
    scale_fill_brewer(direction = -1, palette = "Pastel2")+ # change palette for bars
    coord_flip()+ # flip x and y axis
    scale_y_continuous(sec.axis = sec_axis(~ .*1, #create a secondary "x" axis for the log10(pvalue data)
                                           labels = number_format(scale=scale_factor), 
                                           name="log10(pvalue) x sign of Activation Z Score"))+
    xlab("Upstream Regulator")+
    ylab("Activation Z Score")+
    theme_bw() + # remove grey background
    theme( axis.line.x.top = element_line(color = "blue"),axis.ticks.x.top= element_line(color = "blue"), # change top x axis color to match pvalues
           axis.text.x.top = element_text(color = "blue"),axis.title.x.top=element_text(colour="blue"))+
    theme( axis.line.x.bottom = element_line(color = "red"),axis.ticks.x.bottom= element_line(color = "red"), # change bottom x axis color to match z scores
           axis.text.x.bottom = element_text(color = "red"),axis.title.x.bottom=element_text(colour="red"))+
    ggtitle(paste0("Top ", max(mol_num), " Upstream Regulators by Molecule Type")) + # add title
    theme(plot.title = element_text(size=22, face="bold")) + #  increase title
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    labs(color="log10(pvalue) x \nActivation Z Score",fill="Activation Z Score")+ # change the order of the legend so log10 is on top
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
    theme( axis.line.x.top = element_line(color = "blue"),axis.ticks.x.top= element_line(color = "blue"), # change top x axis color to match pvalues
           axis.text.x.top = element_text(color = "blue"),axis.title.x.top=element_text(colour="blue"))+
    theme( axis.line.x.bottom = element_line(color = "red"),axis.ticks.x.bottom= element_line(color = "red"), # change bottom x axis color to match z scores
           axis.text.x.bottom = element_text(color = "red"),axis.title.x.bottom=element_text(colour="red"))+
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    theme(legend.direction = "vertical", legend.box = "horizontal")+
    theme(legend.title = element_text(size=12,color="red",face="bold"))
  
  # legend for COLOR only
  p3 = top_df %>%
    mutate(Upstream.Regulator = fct_reorder(Upstream.Regulator, Activation.z.score)) %>%
    ggplot(aes(x = Upstream.Regulator))+
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor,color=""),shape=20, cex=8)+ # create dummy points so pvalue legend symbol appears
    geom_point(aes(y = logpvalue/scale_factor,size=logpvalue/scale_factor),shape=20, cex=8,color="blue")+
    scale_color_manual(name = "log10(pvalue) x \nsign of Activation Z Score",values = c(" " = "blue")) + # add legened for pvalue point
    coord_flip()+ # flip x and y axis
    scale_y_continuous(sec.axis = sec_axis(~ .*1, #create a secondary "x" axis for the log10(pvalue data)
                                           labels = number_format(scale=scale_factor), 
                                           name="log10(pvalue) x sign of Activation Z Score"))+
    theme_bw() + # remove grey background
    theme(legend.key=element_rect(fill=NA)) +# remove grey from pvalue background
    labs(color="log10(pvalue) x \nActivation Z Score",fill="Activation Z Score")+ # change the order of the legend so log10 is on top
    theme(legend.direction = "vertical", legend.box = "horizontal")+
    theme(legend.title = element_text(size=12,color="blue",face="bold"))
  
  # get two legends
  leg2=get_legend(p2)
  leg3=get_legend(p3)
  
  # combine legends
  leg32 = plot_grid(leg3, leg2,
                    nrow = 2)
  
  final_p = plot_grid(p1,
                      leg32,
                      nrow = 1,
                      align = "h",
                      axis = "t",
                      rel_widths = c(1, 0.3))
  
  
  # save and print
  fpath=paste0(ipa_dir,"top_",max(mol_num),"_molecules.png")
  ggsave(fpath,final_p)
  print(final_p)
}

