## IMPORT LIBRARIES -----------------------------------------------------------
library("data.table")
library("ggplot2")
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


#Annotation file

interpro_db <- fread(
  file = 'mart_export.txt', 
  sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names=FALSE,
  select = c("Gene name", "Interpro ID", "Interpro Short Description", "Interpro Description"))

#remove duplicate rows
interpro_db<-interpro_db[!duplicated(interpro_db), ]

#Database with Interpro IDs and entry type

interpro_type <-read.table(
  file = '~/tblab/ana/tests/InterPro/entry.list', 
  sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#get only InterPro Domains
interpro_db <- merge(interpro_db, interpro_type[interpro_type$ENTRY_TYPE =="Domain",],
                     by.x="Interpro ID", by.y="ENTRY_AC")

# FUNCTIONAL ANNOTATION---------------------------------------------------------
DUP_INTERPRO <- merge(DUP_gene, interpro_db, by.x="Gene_name",by.y="Gene name")
DEL_INTERPRO <- merge(DEL_gene, interpro_db, by.x="Gene_name",by.y="Gene name")

# ENRICHMENT FUNCTION ----------------------------------------------------------
run_INTERPRO <- function(INTERPRO_df){
  
  ##find GO terms with 5-500 genes -----------------------------------------------
  gene_INTERPRO<-INTERPRO_df[,c("Gene_name","Interpro ID")]
  INTERPRO_count<-as.data.frame(table(unlist(gene_INTERPRO$`Interpro ID`)))
  colnames(INTERPRO_count)<-c("Interpro ID","number_genes")
  
  INTERPRO_filt<- INTERPRO_count[INTERPRO_count$number_genes>=5 & INTERPRO_count$number_genes<=500,]$`Interpro ID`
  
  
  #empty data frame with results of fisher's test
  fisher_df<-data.frame(ID=character(), description=character(), disease=character(), 
                        gene_count=integer(), gene_total=character(), 
                        gene_ratio=double(), p_value=double())
  
  #diseases with less than 50 cases 
  few_cases_ds <- c("Miopatias","Metabolicas","Inflamatoria","Esterilidad",
                    "Endocrinologica", "Dermatologicas","Prenatal","Cancer","Varios")
  
  ## parse diseases with >50 cases ----------------------------------------------
  
  for (i_ds in seq(11,(ncol(INTERPRO_df)-6), by=6)) {
    disease_col=names(INTERPRO_df)[i_ds] #get column name
    disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
    
    if (!( disease %in% few_cases_ds)){
      #get AC for disease and pseudocontrols
      aux_ds <- INTERPRO_df[, c(1,i_ds)]
      colnames(aux_ds)<-c("Gene_name","AC")
      aux_PC <- INTERPRO_df[, c(1,i_ds+3)]
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
      INTERPRO_ds <- unique(INTERPRO_df[(INTERPRO_df$Gene_name %in% count_filt_ds$Gene_name),]$`Interpro ID`)
      INTERPRO_PC <- unique(INTERPRO_df[(INTERPRO_df$Gene_name %in% count_filt_PC$Gene_name),]$`Interpro ID`)
      
      ## parse list of INTERPRO IDS --------------------------------------------
      for (id in INTERPRO_ds){
        
        #run test for INTERPRO IDs with 5-500 genes 
        
        if (id %in% INTERPRO_filt){
          
          #list of annotated genes associated with this ID
          gene_INTERPRO_list <- unique(INTERPRO_df[INTERPRO_df$`Interpro ID`==id,]$Gene_name)
          #list of annotated genes not associated with this ID
          gene_NOINTERPRO_list <- unique(INTERPRO_df[!(INTERPRO_df$Gene_name %in% gene_INTERPRO_list) ,]$Gene_name)
          
          #create contingency matrix
          count_INTERPRO_ds <- sum(count_filt_ds[count_filt_ds$Gene_name %in% gene_INTERPRO_list,]$AC)
          count_NOINTERPRO_ds <- sum(count_filt_ds[count_filt_ds$Gene_name %in% gene_NOINTERPRO_list,]$AC)
          count_INTERPRO_PC <- sum(count_filt_PC[count_filt_PC$Gene_name %in% gene_INTERPRO_list,]$AC)
          count_NOINTERPRO_PC <- sum(count_filt_PC[count_filt_PC$Gene_name %in% gene_NOINTERPRO_list,]$AC)
          
          cont_mx<-matrix(c(count_INTERPRO_ds, count_NOINTERPRO_ds, count_INTERPRO_PC, count_NOINTERPRO_PC),
                          nrow = 2,
                          dimnames = list(c("INTERPRO", "NO INTERPRO"),
                                          c("Disease", "No disease")))
          
          ## gene count and gene ratio ---------------------------------------------
          gene_count <- length(count_filt_ds[count_filt_ds$Gene_name %in% gene_INTERPRO_list,]$Gene_name)
          gene_total <- length(count_filt_ds$Gene_name)
          gene_ratio <- gene_count/gene_total
          
          #get short description
          term_desc <-unique(INTERPRO_df[INTERPRO_df$`Interpro ID`==id,]$`Interpro Short Description`)
          
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
  
  
  return(fisher_sig)
  
}

# RUN ENRICHMENT --------------------------------------------------------------
DUP_Interpro_fisher_sig <- run_INTERPRO(DUP_INTERPRO)
DEL_Interpro_fisher_sig <- run_INTERPRO(DEL_INTERPRO)


## DOTPLOT REPRESENTATION ------------------------------------------------------
#create function
dotplot_rep <- function(enrich_tab, type){
  diseases <- unique(enrich_tab$disease) #list of diseases
  #parse diseases
  for (ds in diseases) {
    #filter by disease
    enrich_filtered<-enrich_tab[enrich_tab$disease==ds,]
    #show 20 ocurrences with highest p-value adjusted
    enrich_filtered <- enrich_filtered[order(enrich_filtered$p_adjust),]
    enrich_filtered <- enrich_filtered[1:20,]
    enrich_filtered <- na.omit(enrich_filtered)
    
    enrich_filtered$gene_count <-as.numeric(as.character(enrich_filtered$gene_count))
    enrich_filtered$gene_ratio <-as.numeric(as.character(enrich_filtered$gene_ratio))
    
    #dotplot
    dp<-ggplot(enrich_filtered, aes(x=gene_ratio, y=reorder(description, gene_ratio))) +
      geom_point(aes(colour=p_adjust, size=gene_count)) +
      scale_colour_gradient(low = "red", high = "blue") +
      ylab("")+
      ggtitle(ds) 
    
    title=paste(type,"_",ds)
    #plot(dp)
  }
}

#run functions
dotplot_rep(DUP_Interpro_fisher_sig,"DUP")
dotplot_rep(DEL_Interpro_fisher_sig,"DEL")

#JACCARD INDEX------------------------------------------------

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


plot_JI_Interpro<- function(df,type){
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

plot_JI_Interpro(DUP_Interpro_fisher_sig,"DUP")
plot_JI_Interpro(DEL_Interpro_fisher_sig,"DEL")



