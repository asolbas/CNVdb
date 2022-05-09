library("data.table")
library("ggplot2")
theme_set(theme_bw())

#IMPORT DATA FILES -------------------------------------------------------------
#Data bases
DUP_db<- read.table(
  file = '~/tblab/ana/database/DB_DUP_PC.tsv',
  sep = '\t', header = TRUE, stringsAsFactors = FALSE,check.names=FALSE)
#add a unique id for each cnv
DUP_db<-cbind(ID=seq.int(nrow(DUP_db)),DUP_db)

DEL_db<- read.table(
  file = '~/tblab/ana/database/DB_DEL_PC.tsv',
  sep = '\t', header = TRUE, stringsAsFactors = FALSE,check.names=FALSE)
#add unique id for each cnv
DEL_db<-cbind(ID=seq.int(nrow(DEL_db)),DEL_db)

# GENE ANNOTATIONS--------------------------------------------------------------

#Sophia annotation file
sophia <- read.table(
  file = '~/bioinfo/fjd/beds/reanalysisBeds/sophia_clinical_exome_ces_annotated.bed',
  sep = '\t', header = FALSE, stringsAsFactors = FALSE,check.names=FALSE)

colnames(sophia) <- c("chr","start","end","Gene_name")
sophia$chr<-gsub('chr','',sophia$chr)

#Gene annotation of database
DUP_gene <-merge(sophia,DUP_db,by.x=c("chr","start","end"),by.y=c("chr","start","end"))
DUP_gene <-DUP_gene[order(DUP_gene$ID, DUP_gene$start),]

DEL_gene <-merge(sophia,DEL_db,by.x=c("chr","start","end"),by.y=c("chr","start","end"))
DEL_gene <-DEL_gene[order(DEL_gene$ID, DEL_gene$start),]

# KEGG ANNOTATIONS -------------------------------------------------------------
kegg_categories<-read.table(
  file = '~/tblab/ana/tests/kegg/kegg_categories3.tsv', 
  sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#delete disease category rows
kegg_categories<-kegg_categories[kegg_categories$Category != 6,]

#Kegg-gene annotation
kegg_db <- fread(
  file = '~/tblab/ana/tests/dbNSFP4.3_gene.complete.txt', 
  sep = '\t', header = TRUE, stringsAsFactors = FALSE, 
  select = c("Gene_name", "Pathway(KEGG)_id", "Pathway(KEGG)_full"))

#delete genes from kegg db that are not annotated
kegg_db <- kegg_db[!(kegg_db$`Pathway(KEGG)_id`=="."),]
colnames(kegg_db) <- c("Gene_name", "kegg_id", "kegg_term")

#separate genes with >1 term into multiple rows
kegg_db <- kegg_db[, list(Gene_name = Gene_name, 
                          kegg_id = unlist(strsplit(kegg_id, ";")), 
                          kegg_term = unlist(strsplit(kegg_term, ";"))),by=1:nrow(kegg_db)]

#add categories to kegg_db
kegg_cat_db <- merge(kegg_db, kegg_categories, by.x="kegg_id", by.y="KeggID")
kegg_cat_db <- kegg_cat_db[,c("Gene_name","kegg_id","kegg_term","Category.Subcategory","CategoryName"),
                           with=FALSE]

#annotate db
DUP_KEGG <- merge(DUP_gene,kegg_cat_db,by.x="Gene_name", by.y="Gene_name")
DEL_KEGG <- merge(DEL_gene,kegg_cat_db,by.x="Gene_name", by.y="Gene_name")

#ENRICHMENT FUNCTION -----------------------------------------------------------
run_KEGG <- function(KEGG_df){
  
  #empty data frame with results of fisher's test
  fisher_df<-data.frame(category_name=character(), disease=character(), 
                        gene_count=integer(), gene_total=character(), 
                        gene_ratio=double(), p_value=double())
  
  #diseases with less than 50 cases 
  few_cases_ds <- c("Miopatias","Metabolicas","Inflamatoria","Esterilidad",
                    "Endocrinologica", "Dermatologicas","Prenatal","Cancer","Varios")
  
  ## parse diseases with >50 cases ----------------------------------------------
  
  for (i_ds in seq(11,(ncol(KEGG_df)-5), by=6)) {
    disease_col=names(KEGG_df)[i_ds] #get column name
    disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
    
    print(disease)
    
    if (!( disease %in% few_cases_ds)){
      #get AC for disease and pseudocontrols
      aux_ds <- KEGG_df[, c(1,i_ds)]
      colnames(aux_ds)<-c("Gene_name","AC")
      aux_PC <- KEGG_df[, c(1,i_ds+3)]
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
      KEGG_ds <- unique(KEGG_df[(KEGG_df$Gene_name %in% count_filt_ds$Gene_name),]$CategoryName)
      KEGG_PC <- unique(KEGG_df[(KEGG_df$Gene_name %in% count_filt_PC$Gene_name),]$CategoryName)
      
      ## parse list of GOs ------------------------------------------------------
      for (id in KEGG_ds){
        
        #list of annotated genes associated with this KEGG category
        gene_KEGG_list <- unique(KEGG_df[KEGG_df$CategoryName==id,]$Gene_name)
        #list of annotated genes not associated with this KEGG category
        gene_NOKEGG_list <- unique(KEGG_df[!(KEGG_df$CategoryName %in% gene_KEGG_list) ,]$Gene_name)
        
        #create contingency matrix
        count_KEGG_ds <- sum(count_filt_ds[count_filt_ds$Gene_name %in% gene_KEGG_list,]$AC)
        count_NOKEGG_ds <- sum(count_filt_ds[count_filt_ds$Gene_name %in% gene_NOKEGG_list,]$AC)
        count_KEGG_PC <- sum(count_filt_PC[count_filt_PC$Gene_name %in% gene_KEGG_list,]$AC)
        count_NOKEGG_PC <- sum(count_filt_PC[count_filt_PC$Gene_name %in% gene_NOKEGG_list,]$AC)
        
        cont_mx<-matrix(c(count_KEGG_ds, count_NOKEGG_ds, count_KEGG_PC, count_NOKEGG_PC),
                        nrow = 2,
                        dimnames = list(c("KEGG", "NO KEGG"),
                                        c("Disease", "No disease")))
        
        ## gene count and gene ratio ---------------------------------------------
        gene_count <- length(count_filt_ds[count_filt_ds$Gene_name %in% gene_KEGG_list,]$Gene_name)
        gene_total <- length(count_filt_ds$Gene_name)
        gene_ratio <- gene_count/gene_total
        
        ## Fisher test: get p-value ----------------------------------------------
        p_value<-fisher.test(cont_mx, alternative="greater")$p.value
        fisher_df[nrow(fisher_df) + 1,] = c(id, disease, gene_count, gene_total,gene_ratio, p_value)
          
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

#RUN ENRICHMENT TEST AND WRITE OUTPUT-------------------------------------------
DUP_Kegg_fisher_sig <- run_KEGG(DUP_KEGG)
DEL_Kegg_fisher_sig <- run_KEGG(DEL_KEGG)

#DOTPLOT------------------------------------------------------------------------

#create function
dotplot_rep <- function(enrich_filtered, title){
  
  #show 20 ocurrences with highest p-value adjusted
  enrich_filtered <- enrich_filtered[order(enrich_filtered$p_adjust),]
  #enrich_filtered <- enrich_filtered[1:20,]
  #enrich_filtered <- na.omit(enrich_filtered)
  
  enrich_filtered$gene_count <-as.numeric(as.character(enrich_filtered$gene_count))
  enrich_filtered$gene_ratio <-as.numeric(as.character(enrich_filtered$gene_ratio))
  
  #dotplot
  dp<-ggplot(enrich_filtered, aes(x=disease, y=reorder(category_name, gene_count))) +
    geom_point(aes(colour=p_adjust, size=gene_ratio)) +
    scale_colour_gradient(low = "red", high = "blue") +
    ylab("Kegg categories") + xlab("Diseases") +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

}

#run functions
dotplot_rep(DUP_Kegg_fisher_sig,"Duplications")
dotplot_rep(DEL_Kegg_fisher_sig, "Deletions")

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


plot_JI_Kegg<- function(df,type){
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

plot_JI_Kegg(DUP_Kegg_fisher_sig,"DUP")
plot_JI_Kegg(DEL_Kegg_fisher_sig,"DEL")


