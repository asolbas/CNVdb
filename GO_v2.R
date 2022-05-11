if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db") #human databse
BiocManager::install("clusterProfiler")
BiocManager::install("rrvgo")
BiocManager::install("GO.db")
install.packages("data.table")
# Install and load data.table
library("data.table")
library("org.Hs.eg.db") #human databse
library("GO.db") #db with GO_ID-GO_term relations
library("rrvgo") ##Bubbleplot representation
library("clusterProfiler")
library("ggplot2")
theme_set(theme_bw())

#IMPORT DATA FILES -------------------------------------------------------------
#Data bases
DUP_db<- read.table(
  file = 'DB_DUP_PC.tsv',
  sep = '\t', header = TRUE, stringsAsFactors = FALSE,check.names=FALSE)
#add a unique id for each cnv
DUP_db<-cbind(ID=seq.int(nrow(DUP_db)),DUP_db)

DEL_db<- read.table(
  file = 'DB_DEL_PC.tsv',
  sep = '\t', header = TRUE, stringsAsFactors = FALSE,check.names=FALSE)
#add unique id for each cnv
DEL_db<-cbind(ID=seq.int(nrow(DEL_db)),DEL_db)

# GENE ANNOTATIONS--------------------------------------------------------------

#Sophia annotation file
sophia <- read.table(
  file = 'sophia_clinical_exome_ces_annotated.bed',
  sep = '\t', header = FALSE, stringsAsFactors = FALSE,check.names=FALSE)

colnames(sophia) <- c("chr","start","end","Gene_name")
sophia$chr<-gsub('chr','',sophia$chr)

#Gene annotation of database
DUP_gene <-merge(sophia,DUP_db,by.x=c("chr","start","end"),by.y=c("chr","start","end"))
DUP_gene <-DUP_gene[order(DUP_gene$ID, DUP_gene$start),]

DEL_gene <-merge(sophia,DEL_db,by.x=c("chr","start","end"),by.y=c("chr","start","end"))
DEL_gene <-DEL_gene[order(DEL_gene$ID, DEL_gene$start),]

# GO ANNOTATIONS ---------------------------------------------------------------
#annotate GO terms with biomart
#SYMBOL ENTREZID GO EVIDENCE ONTOLOGY
DUP_GO <- bitr(DUP_gene$Gene_name, fromType="SYMBOL", toType=(c("ENTREZID","GO")), 
               OrgDb="org.Hs.eg.db",drop=TRUE)
#delete genes with no asociated terms
DUP_GO <- na.omit(DUP_GO) 

DEL_GO <- bitr(DEL_gene$Gene_name, fromType="SYMBOL", toType=(c("ENTREZID","GO")), 
               OrgDb="org.Hs.eg.db",drop=TRUE)
DEL_GO <- na.omit(DEL_GO) #table of annotated genes

#RVIGO FUNCTION ----------------------------------------------------------------
run_rvgo<- function(enrich_tab, ont){
  #get disease list
  disease_ls <- unique(enrich_tab$disease)
  
  #create GO parents data frame
  #will store the output
  
  df_reduced <- data.frame(matrix(ncol=ncol(enrich_tab),nrow=0))
  colnames(df_reduced) <- colnames(enrich_tab)
  
  #parse diseases
  for (ds in disease_ls) {
    df_ds <- enrich_tab[enrich_tab$disease==ds,]
    
    print(df_ds$GO)
    print(length(df_ds$GO))
    print(ds)
    
    if (length(df_ds$GO)>1){
      
      simMatrix <- calculateSimMatrix(df_ds$GO,
                                      orgdb="org.Hs.eg.db",
                                      ont= ont,
                                      method="Rel")
      
      #rrvgo:the higher score, the better
      #Therefore we must transform p_values with -log(p-value)
      scores <- setNames(-log10(df_ds$p_adjust), df_ds$GO)
      
      reducedTerms <- reduceSimMatrix(simMatrix,
                                      scores,
                                      threshold=0.7,
                                      orgdb="org.Hs.eg.db")
      
      #remove all GOs that are not parents terms
      df_parents_ds <- df_ds[df_ds$GO %in% reducedTerms$parent,]
      
      #add the table to the results df
      df_reduced <- rbind(df_reduced,df_parents_ds)
    }
    
    else{
      df_reduced <- rbind(df_reduced,df_ds)
    }
    
  }
  return(df_reduced)
}


# ENRICHMENT FUNCTION ----------------------------------------------------------
run_GO <- function(gene_df,GO_df,ont){
  
  ##find GO terms with 5-500 genes -----------------------------------------------
  GO_count<-as.data.frame(table(unlist(GO_df[GO_df$ONTOLOGY==ont,]$GO)))
  colnames(GO_count)<-c("GO","number_genes")
  
  GO_filt<- GO_count[GO_count$number_genes>=5 & GO_count$number_genes<=500,]$GO
  
  ## association between GO ids and descriptions ---------------------------------
  TERM_ls<-as.list(GOTERM)
  
  #empty data frame with results of fisher's test
  fisher_df<-data.frame(GO=character(),GOTERM=character(), disease=character(), 
                        gene_count=integer(), gene_total=character(), 
                        gene_ratio=double(), p_value=double())
  
  #diseases with less than 50 cases 
  few_cases_ds <- c("Miopatias","Metabolicas","Inflamatoria","Esterilidad",
                    "Endocrinologica", "Dermatologicas","Prenatal","Cancer","Varios")
  
  ## parse diseases with >50 cases ----------------------------------------------
  
  for (i_ds in seq(11,(ncol(gene_df)-2), by=6)) {
    disease_col=names(gene_df)[i_ds] #get column name
    disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
    
    if (!( disease %in% few_cases_ds)){
      #get AC for disease and pseudocontrols
      aux_ds <- gene_df[, c(4,i_ds)]
      colnames(aux_ds)<-c("Gene_name","AC")
      aux_PC <- gene_df[, c(4,i_ds+3)]
      colnames(aux_PC)<-c("Gene_name","AC")
      
      #aggregate AC by gene
      count_ds <- aggregate(aux_ds$AC, by=list(gene=aux_ds$Gene_name), FUN=sum)
      colnames(count_ds)<-c("Gene_name","AC")
      
      count_PC <- aggregate(aux_PC$AC, by=list(gene=aux_PC$Gene_name), FUN=sum)
      colnames(count_PC)<-c("Gene_name","AC")
      
      #filter tables by AC>2
      count_filt_ds <- count_ds[count_ds$AC>2,]
      count_filt_PC <- count_PC[count_PC$AC>2,]
      
      ## functional annotation of genes -----------------------------------------
      GO_ds <- unique(GO_df[(GO_df$SYMBOL %in% count_filt_ds$Gene_name)
                            & GO_df$ONTOLOGY==ont,]$GO)
      GO_PC <- unique(GO_df[(GO_df$SYMBOL %in% count_filt_PC$Gene_name)
                            & GO_df$ONTOLOGY==ont,]$GO)
      
      ## parse list of GOs ------------------------------------------------------
      for (id in GO_ds){
        
        #run test for GOs with 5-500 genes 
        
        if (id %in% GO_filt){
          
          #get description for GO id
          term_desc<-Term(TERM_ls[[id]])
          
          #list of annotated genes associated with this GO
          gene_GO_list <- unique(GO_df[GO_df$GO==id,]$SYMBOL)
          #list of annotated genes not associated with this GO
          gene_NOGO_list <- unique(GO_df[GO_df$ONTOLOGY==ont & !(GO_df$SYMBOL %in% gene_GO_list) ,]$SYMBOL)
          
          #create contingency matrix
          count_GO_ds <- sum(count_filt_ds[count_filt_ds$Gene_name %in% gene_GO_list,]$AC)
          count_NOGO_ds <- sum(count_filt_ds[count_filt_ds$Gene_name %in% gene_NOGO_list,]$AC)
          count_GO_PC <- sum(count_filt_PC[count_filt_PC$Gene_name %in% gene_GO_list,]$AC)
          count_NOGO_PC <- sum(count_filt_PC[count_filt_PC$Gene_name %in% gene_NOGO_list,]$AC)
          
          cont_mx<-matrix(c(count_GO_ds, count_NOGO_ds, count_GO_PC, count_NOGO_PC),
                          nrow = 2,
                          dimnames = list(c("GO", "NO GO"),
                                          c("Disease", "No disease")))
          
          ## gene count and gene ratio ---------------------------------------------
          gene_count <- length(count_filt_ds[count_filt_ds$Gene_name %in% gene_GO_list,]$Gene_name)
          gene_total <- length(count_filt_ds$Gene_name)
          gene_ratio <- gene_count/gene_total
          
          ## Fisher test: get p-value ----------------------------------------------
          p_value<-fisher.test(cont_mx, alternative="greater")$p.value
          fisher_df[nrow(fisher_df) + 1,] = c(id, term_desc, disease, gene_count, gene_total,gene_ratio, p_value)
          
        }
      }
    }
  }
  ## Benjamini-Hochberg correction of p-value ------------------------------------
  p_adj<-p.adjust(fisher_df$p_value, method="fdr") #adjusted p-value
  fisher_df$p_adjust=p_adj
  
  ## p adjust < 0.05 -------------------------------------------------------------
  fisher_sig <- fisher_df[fisher_df$p_adjust <= 0.05,]
  
  ## Reduce redundat terms with rrvigo -------------------------------------------
  fisher_sig_red <- run_rvgo(fisher_sig,ont)
  
  return(fisher_sig_red)
  
}

#RUN ENRICHMENT FUNCTION AND WRITE OUTPUT --------------------------------------

##run function ----------------------------------------------------------------
DUP_CC_sig_red <- run_GO(DUP_gene,DUP_GO,"CC")
DUP_BP_sig_red <- run_GO(DUP_gene,DUP_GO,"BP")
DUP_MF_sig_red <- run_GO(DUP_gene,DUP_GO,"MF")

DEL_CC_sig_red <- run_GO(DEL_gene,DEL_GO,"CC")
DEL_BP_sig_red <- run_GO(DEL_gene,DEL_GO,"BP")
DEL_MF_sig_red <- run_GO(DEL_gene,DEL_GO,"MF")

# DOTPLOT REPRESENTATION -------------------------------------------------------
prepare_data <- function(enrich_tab, ont){
  #filter by p_adjust
  enrich_filtered<-enrich_tab[enrich_tab$p_adjust<=0.05,]
  #order by p_adjusted
  enrich_filtered <- enrich_filtered[order(enrich_filtered$p_adjust),]
  #show 15 occurrences with highest p-value adjusted
  enrich_filtered <- enrich_filtered[1:15,]
  #delete NA rows (for diseases with less than 15 enriched terms)
  enrich_filtered <- na.omit(enrich_filtered)
  #add ontology column
  enrich_filtered$ontology <- ont
  
  return(enrich_filtered)
}

dotplot_rep <- function(enrich_CC,enrich_BP,enrich_MF,type){
  diseases <- unique(enrich_BP$disease) #list of diseases
  #parse diseases
  for (ds in diseases) {
    enrich_CC_ds<-prepare_data(enrich_CC[enrich_CC$disease==ds,],"CC")
    enrich_BP_ds<-prepare_data(enrich_BP[enrich_BP$disease==ds,],"BP")
    enrich_MF_ds<-prepare_data(enrich_MF[enrich_MF$disease==ds,],"MF")
    
    enrich_filtered <- rbind(enrich_CC_ds,enrich_BP_ds,enrich_MF_ds)
    
    enrich_filtered$gene_count <-as.numeric(as.character(enrich_filtered$gene_count))
    enrich_filtered$gene_ratio <-as.numeric(as.character(enrich_filtered$gene_ratio))
    
    #dotplot
    dp<-ggplot(enrich_filtered, aes(x=gene_ratio, y=reorder(GOTERM, gene_ratio))) +
      geom_point(aes(colour=p_adjust, size=gene_count)) +
      scale_colour_gradient(low = "red", high = "blue") +
      labs(
        x = "Gene ratio",
        y = "",
        colour = "p_adjust",
        size = "Gene count",
        title = ds
      )
    dp_all <- plot(dp + facet_grid(rows = vars(ontology),scales = "free_y"))
  }
}

dotplot_rep(DUP_CC_sig_red, DUP_BP_sig_red, DUP_MF_sig_red,"DUP")
dotplot_rep(DEL_CC_sig_red, DEL_BP_sig_red, DEL_MF_sig_red,"DEL")

#JACCARD INDEX BY ONTOLOGY ----------------------------------------------------
jaccard_index<- function(df){
  ds_list<-unique(df$disease)
  
  #empty matrix
  JI<-matrix(, nrow = length(ds_list), ncol = length(ds_list))
  rownames(JI) <- ds_list
  colnames(JI) <- ds_list
  
  #calculate jaccard index
  for (ds1 in ds_list){
    for (ds2 in ds_list){
      #intersection of d1 and ds2 terms
      intersect_ds<-intersect(df[df$disease==ds1, 1], 
                              df[df$disease==ds2, 1])
      #union of ds1 and ds2 genes
      union_ds<-union(df[df$disease==ds1, 1],  
                      df[df$disease==ds2, 1])
      
      JI[ds1,ds2] <-length(intersect_ds)/length(union_ds)
      JI[ds2,ds1] <- JI[ds1,ds2]
    }
  }
  return(JI)
}

get_lower_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}


plot_JI<- function(df,type,ont){
  JI <-jaccard_index(df)
  
  # Get lower triangle of the correlation matrix
  lower_tri_JI <- get_lower_tri(JI)
  melted_JI <- melt(lower_tri_JI, na.rm=TRUE)
  
  #title function
  title <- paste(type,"_",ont)
  
  p <- ggplot(data = melted_JI, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    ggtitle(title)+
    geom_text(aes(Var1, Var2, label = round(value, digits=2)), color = "black", size = 2) +
    scale_fill_gradient2(low = "white", high = "cornflowerblue",
                         space = "Lab",
                         name="Jaccard Index") +
    theme_minimal() +
    scale_y_discrete(limits=rev)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    coord_fixed() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank())
}

plot_JI(DUP_CC_sig_red,"DUP","CC")
plot_JI(DUP_BP_sig_red,"DUP","BP")
plot_JI(DUP_MF_sig_red,"DUP","MF")

plot_JI(DEL_CC_sig_red,"DEL","CC")
plot_JI(DEL_BP_sig_red,"DEL","BP")
plot_JI(DEL_MF_sig_red,"DEL","MF")

#JACCARD INDEX TOTAL------------------------------------------------------------

jaccard_index<- function(df){
  ds_list<-unique(df$disease)
  
  #empty matrix
  JI<-matrix(, nrow = length(ds_list), ncol = length(ds_list))
  rownames(JI) <- ds_list
  colnames(JI) <- ds_list
  
  #calculate jaccard index
  for (ds1 in ds_list){
    for (ds2 in ds_list){
      #intersection of d1 and ds2 terms
      intersect_ds<-intersect(df[df$disease==ds1, 1], 
                              df[df$disease==ds2, 1])
      #union of ds1 and ds2 genes
      union_ds<-union(df[df$disease==ds1, 1],  
                      df[df$disease==ds2, 1])
      
      JI[ds1,ds2] <-length(intersect_ds)/length(union_ds)
      JI[ds2,ds1] <- JI[ds1,ds2]
    }
  }
  return(JI)
}

get_lower_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}


plot_JI_GO<- function(df,type){
  JI <-jaccard_index(df)
  
  # Get lower triangle of the correlation matrix
  lower_tri_JI <- get_lower_tri(JI)
  melted_JI <- melt(lower_tri_JI, na.rm=TRUE)
  
  #title function
  title <- paste(type)
  
  p <- ggplot(data = melted_JI, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    ggtitle(title)+
    geom_text(aes(Var1, Var2, label = round(value, digits=2)), color = "black", size = 2) +
    scale_fill_gradient2(low = "white", high = "cornflowerblue",
                         space = "Lab",
                         name="Jaccard Index") +
    theme_minimal() +
    scale_y_discrete(limits=rev)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    coord_fixed() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank())
  
}

plot_JI_GO(rbind(DUP_CC_sig_red,DUP_BP_sig_red, DUP_MF_sig_red),"DUP")
plot_JI_GO(rbind(DEL_CC_sig_red,DEL_BP_sig_red, DEL_MF_sig_red),"DEL")
