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

#HPO ANNOTATIONS ---------------------------------------------------------------
#HPO-gene annotation
hpo_db <- fread(
  file = '~/tblab/ana/tests/dbNSFP4.3_gene.complete.txt', 
  sep = '\t', header = TRUE, stringsAsFactors = FALSE, 
  select = c("Gene_name", "HPO_id", "HPO_name"))

#delete genes from hpo db that are not annotated
hpo_db <- hpo_db[!(hpo_db$HPO_id == "."),]
colnames(hpo_db) <- c("Gene_name", "hpo_id", "hpo_name")

#separate genes with >1 term into multiple rows
hpo_db <- hpo_db[, list(Gene_name = Gene_name, 
                        hpo_id = unlist(strsplit(hpo_id, ";")), 
                        hpo_name = unlist(strsplit(hpo_name, ";"))),by=1:nrow(hpo_db)]

#annotate db
DUP_HPO <- merge(DUP_gene,hpo_db,by.x="Gene_name", by.y="Gene_name")
DEL_HPO <- merge(DEL_gene,hpo_db,by.x="Gene_name", by.y="Gene_name")

#RUN ENRICHMENT ----------------------------------------------------------------
run_HPO <- function(HPO_df){
  
  ##find HPO terms with 5-500 genes -----------------------------------------------
  gene_HPO<-HPO_df[,c("Gene_name","hpo_id")]
  HPO_count<-as.data.frame(table(unlist(gene_HPO$hpo_id)))
  colnames(HPO_count)<-c("hpo_id","number_genes")
  
  HPO_filt<- unique(HPO_count[HPO_count$number_genes>=10 & HPO_count$number_genes<=100,]$hpo_id)
  
  
  #empty data frame with results of fisher's test
  fisher_df<-data.frame(ID=character(), name=character(), disease=character(), 
                        gene_count=integer(), gene_total=character(), 
                        gene_ratio=double(), p_value=double())
  
  #diseases with less than 50 cases 
  few_cases_ds <- c("Miopatias","Metabolicas","Inflamatoria","Esterilidad",
                    "Endocrinologica", "Dermatologicas","Prenatal","Cancer","Varios")
  
  ## parse diseases with >50 cases ----------------------------------------------
  
  for (i_ds in seq(11,(ncol(HPO_df)-6), by=6)) {
    disease_col=names(HPO_df)[i_ds] #get column name
    disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
    
    if (!( disease %in% few_cases_ds)){
      print(disease)
      #get AC for disease and pseudocontrols
      aux_ds <- HPO_df[, c(1,i_ds)]
      colnames(aux_ds)<-c("Gene_name","AC")
      aux_PC <- HPO_df[, c(1,i_ds+3)]
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
      HPO_ds <- unique(HPO_df[(HPO_df$Gene_name %in% count_filt_ds$Gene_name),]$hpo_id)
      HPO_PC <- unique(HPO_df[(HPO_df$Gene_name %in% count_filt_PC$Gene_name),]$hpo_id)
      
      print(length(HPO_ds))
      ## parse list of HPO IDS --------------------------------------------
      for (id in HPO_ds){
        
        #run test for HPO IDs with 10-100 genes 
        
        if (id %in% HPO_filt){
          
          #list of annotated genes associated with this ID
          gene_HPO_list <- unique(HPO_df[HPO_df$hpo_id==id,]$Gene_name)
          #list of annotated genes not associated with this ID
          gene_NOHPO_list <- unique(HPO_df[!(HPO_df$Gene_name %in% gene_HPO_list) ,]$Gene_name)
          
          #create contingency matrix
          count_HPO_ds <- sum(count_filt_ds[count_filt_ds$Gene_name %in% gene_HPO_list,]$AC)
          count_NOHPO_ds <- sum(count_filt_ds[count_filt_ds$Gene_name %in% gene_NOHPO_list,]$AC)
          count_HPO_PC <- sum(count_filt_PC[count_filt_PC$Gene_name %in% gene_HPO_list,]$AC)
          count_NOHPO_PC <- sum(count_filt_PC[count_filt_PC$Gene_name %in% gene_NOHPO_list,]$AC)
          
          cont_mx<-matrix(c(count_HPO_ds, count_NOHPO_ds, count_HPO_PC, count_NOHPO_PC),
                          nrow = 2,
                          dimnames = list(c("HPO", "NO HPO"),
                                          c("Disease", "No disease")))
          
          ## gene count and gene ratio ---------------------------------------------
          gene_count <- length(count_filt_ds[count_filt_ds$Gene_name %in% gene_HPO_list,]$Gene_name)
          gene_total <- length(count_filt_ds$Gene_name)
          gene_ratio <- gene_count/gene_total
          
          #get short description
          name_id <-unique(HPO_df[HPO_df$hpo_id==id,]$hpo_name)
          
          ## Fisher test: get p-value ----------------------------------------------
          p_value<-fisher.test(cont_mx, alternative="greater")$p.value
          fisher_df[nrow(fisher_df) + 1,] = c(id, name_id, disease, gene_count, gene_total,gene_ratio, p_value)
          
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

#RUN TEST ----------------------------------------------------------------------
DUP_hpo_fisher_sig <- run_HPO(DUP_HPO)
DEL_hpo_fisher_sig <- run_HPO(DEL_HPO)


## DOTPLOT REPRESENTATION ------------------------------------------------------
#create function
dotplot_rep <- function(enrich_tab,type){
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
    dp<-ggplot(enrich_filtered, aes(x=gene_ratio, y=reorder(name, gene_ratio))) +
      geom_point(aes(colour=p_adjust, size=gene_count)) +
      scale_colour_gradient(low = "red", high = "blue") +
      ggtitle(ds) 
    
    title=paste(type,"_",ds)
    #plot(dp)
  }
}

#run functions
dotplot_rep(DUP_hpo_fisher_sig,"DUP")
dotplot_rep(DEL_hpo_fisher_sig,"DEL")


