library("data.table")
library("ggplot2")
library("UpSetR")
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
#install_github("hms-dbmi/UpSetR")
library("ComplexHeatmap")
library(reshape2)

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

#Sophia annotation file
sophia <- read.table(
  file = '~/bioinfo/fjd/beds/reanalysisBeds/sophia_clinical_exome_ces_annotated.bed',
  sep = '\t', header = FALSE, stringsAsFactors = FALSE,check.names=FALSE)

colnames(sophia) <- c("chr","start","end","Gene_name")
sophia$chr<-gsub('chr','',sophia$chr)

#Gene annotation of database
DUP_gene<-merge(sophia,DUP_db,by.x=c("chr","start","end"),by.y=c("chr","start","end"))
DUP_gene <-DUP_gene[order(DUP_gene$ID, DUP_gene$start),]

DEL_gene<-merge(sophia,DEL_db,by.x=c("chr","start","end"),by.y=c("chr","start","end"))
DEL_gene <-DEL_gene[order(DEL_gene$ID, DEL_gene$start),]


#GET DISEASES NAMES-------------------------------------------------------------

few_cases_ds <- c("Miopatias","Metabolicas","Inflamatoria","Esterilidad",
                  "Endocrinologica", "Dermatologicas","Prenatal","Cancer","Varios")

ds_list <- vector(length=12)
ds_col_list_DUP <- vector(length=12)
ds_col_list_DEL <- vector(length=12)
#parse diseases
i=1
for (i_ds in seq(11,(ncol(DUP_gene)-2), by=6)) {
  disease_col_DUP=names(DUP_gene)[i_ds] #get column name
  disease_col_DEL=names(DEL_gene)[i_ds] #get column name
  disease<-unlist(strsplit(disease_col_DUP, split='_', fixed=TRUE))[1] #disease name
  
  if (!( disease %in% few_cases_ds)){
    ds_list[i]=disease
    ds_col_list_DUP[i]=disease_col_DUP
    ds_col_list_DEL[i]=disease_col_DEL
    i=i+1
  }
}

#COUNT TOTAL ANNOTATED GENES BY DISEASE----------------------------------------

count_total <- function(genes_df,few_cases_ds){
  
  count_df<-data.frame(disease=character(),count=character(), type=character())
  
  
  for (i in seq(11,(ncol(genes_df)-2), by=6)) {
    disease_col=names(genes_df)[i] #get column name
    disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
    type<-unique(genes_df$SV_type)
    if (!( disease %in% few_cases_ds)){
      genes_df_filt<-genes_df[genes_df[i]>0,]
      count <- length(unique(genes_df_filt$Gene_name))
      count_df[nrow(count_df) + 1,] = c(disease, count,type)
    }
  }
  return(count_df)
}

DUP_count <- count_total(DUP_gene,few_cases_ds)
DEL_count <- count_total(DEL_gene,few_cases_ds)

count_df <- rbind(DUP_count,DEL_count)

#plot gene count

count_df$count <- as.numeric(as.character(count_df$count))

p<- ggplot(data=count_df, aes(x=disease, y=count, fill=type)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=count), vjust=-1, color="black",
            position = position_dodge(0.9), size=3)+
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=c("cornflowerblue","#FF6633"))  +
  ylab("gene count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

plot(p)

# UPSET PLOT GENES BY DISEASE --------------------------------------------------
upset<-function(gene_df,few_cases_ds){
  #Prepare input
  genes_disease<-list() #stores the list of genes of each disease
  
  for (i in seq(11,(ncol(gene_df)-2), by=6)) {
    disease_col=names(gene_df)[i] #get column name
    disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
    #type<-unique(genes_df$SV_type)
    if (!( disease %in% few_cases_ds)){
      genes_filt<-unique(gene_df[gene_df[i]>0,]$Gene_name) #AC>0
      genes_disease[[disease]]<-genes_filt
    }
  }
  #Make combination matrix
  comb_mat<-make_comb_mat(genes_disease)
  return(comb_mat)
}

DUP_mat<-upset(DUP_gene,few_cases_ds)
DEL_mat<-upset(DEL_gene,few_cases_ds)

#Make UpSet plot
DUP_mat<-DUP_mat[comb_degree(DUP_mat)<4]
p <- UpSet(DUP_mat)
plot(p)

DEL_mat<-DEL_mat[comb_degree(DEL_mat)<4]
p<-UpSet(DEL_mat)
plot(p)


#HEATMAP GENES BY DISEASE ------------------------------------------------------
heatmap_genes<-function(gene_df,few_cases_ds){
  #get diseases names
  ds_list <- vector(length=11)
  ds_col_list <- vector(length=11)
  
  #parse diseases
  i=1
  for (i_ds in seq(11,(ncol(gene_df)-2), by=6)) {
    disease_col=names(gene_df)[i_ds] #get column name
    disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
    
    if (!( disease %in% few_cases_ds)){
      ds_list[i]=disease
      ds_col_list[i]=disease_col
      i=i+1
    }
  }
  
  #build heatmap
  heatMap<-matrix(, nrow = length(ds_list), ncol = length(ds_list))
  rownames(heatMap) <- ds_list
  colnames(heatMap) <- ds_list
  
  
  for (i in 1:length(ds_list)){
    for (j in 1:i){
      ds1_col<-which(colnames(gene_df) == ds_col_list[i]) #index column AC disease 1
      ds2_col<-which(colnames(gene_df) == ds_col_list[j]) #index column AC disease 2
      intersect_ds<-intersect(gene_df[gene_df[ds1_col]>=1, 4], #intersection of d1 and ds2 genes
                              gene_df[gene_df[ds2_col]>=1, 4])
      
      heatMap[i,j] <-length(intersect_ds)
      heatMap[j,i] <- heatMap[i,j]
    }
  }
  
  # Get lower triangle of the correlation matrix
  get_upper_tri<-function(cormat){
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
  }
  
  lower_tri <- get_upper_tri(heatMap)
  melted <- melt(lower_tri, na.rm=TRUE)
  
  type=unique(gene_df$SV_type)
  
  #plot heatmap
  ggplot(data = melted, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    ggtitle(type)+
    geom_text(aes(Var1, Var2, label = round(value, digits=2)), color = "black", size = 2.5) +
    scale_fill_gradient2(low = "white", high = "red",
                         space = "Lab",
                         name="gene count") +
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
#plot heatmap
heatmap_genes(DUP_gene,few_cases_ds)
heatmap_genes(DEL_gene,few_cases_ds)

#JACCARD INDEX GENES------------------------------------------------------------
jaccard_index<-function(gene_df,few_cases_ds){
  
  #get diseases names
  ds_list <- vector(length=11)
  ds_col_list <- vector(length=11)
  
  #parse diseases
  i=1
  for (i_ds in seq(11,(ncol(gene_df)-2), by=6)) {
    disease_col=names(gene_df)[i_ds] #get column name
    disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
    
    if (!( disease %in% few_cases_ds)){
      ds_list[i]=disease
      ds_col_list[i]=disease_col
      i=i+1
    }
  }
  
  JI<-matrix(, nrow = length(ds_list), ncol = length(ds_list))
  rownames(JI) <- ds_list
  colnames(JI) <- ds_list
  
  
  for (i in 1:length(ds_list)){
    for (j in 1:i){
      ds1_col<-which(colnames(gene_df) == ds_col_list[i]) #index column AC disease 1
      ds2_col<-which(colnames(gene_df) == ds_col_list[j]) #index column AC disease 2
      intersect_ds<-intersect(gene_df[gene_df[ds1_col]>=1, 4], #intersection of d1 and ds2 genes
                              gene_df[gene_df[ds2_col]>=1, 4])
      union_ds<-union(gene_df[gene_df[ds1_col]>=1, 4],  #union of ds1 and ds2 genes
                      gene_df[gene_df[ds2_col]>=1, 4])
      
      JI[i,j] <-length(intersect_ds)/length(union_ds)
      JI[j,i] <- JI[i,j]
    }
  }
  
  # Get lower triangle of the correlation matrix
  get_upper_tri<-function(cormat){
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
  }
  
  lower_tri_JI <- get_upper_tri(JI)
  melted_JI <- melt(lower_tri_JI, na.rm=TRUE)
  
  type<-unique(gene_df$SV_type)
  
  ggplot(data = melted_JI, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    ggtitle(type)+
    geom_text(aes(Var1, Var2, label = round(value, digits=2)), color = "black", size = 1.75) +
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
#plot jaccard index
jaccard_index(DUP_gene,few_cases_ds)

jaccard_index(DEL_gene,few_cases_ds)

# COUNT NUMBER GENES PER CNV NUMBER (AC) ---------------------------------------
AC_by_gene <- data.frame(gene=character(), AC=character(), type=character(),disease=character())

for (i_ds in seq(11,(ncol(DUP_gene)-2), by=6)) {
  disease_col_DUP=names(DUP_gene)[i_ds] #get column name
  disease_col_DEL=names(DEL_gene)[i_ds] #get column name
  disease<-unlist(strsplit(disease_col_DUP, split='_', fixed=TRUE))[1] #disease name
  
  if (!( disease %in% few_cases_ds)){
    #filter rows with AC>0
    DUP_gene_filt<-DUP_gene[DUP_gene[i_ds]>0,]
    DEL_gene_filt<-DEL_gene[DEL_gene[i_ds]>0,]
    
    #get list with unique genes
    DUP_genes_ls <- unique(DUP_gene_filt$Gene_name)
    DEL_genes_ls <- unique(DEL_gene_filt$Gene_name)
    
    
    for (gene in DUP_genes_ls){
      
      aux <- DUP_gene_filt[DUP_gene_filt$Gene_name==gene,]
      AC <- sum(aux[i_ds])
      AC_by_gene[nrow(AC_by_gene) + 1,] = c(gene, AC, type="DUP",disease)
    }
    
    for (gene in DEL_genes_ls){
      
      aux <- DEL_gene_filt[DEL_gene_filt$Gene_name==gene,]
      AC <- sum(aux[i_ds])
      AC_by_gene[nrow(AC_by_gene) + 1,] = c(gene, AC, type="DEL",disease)
    }
    
    AC_by_gene$AC<-as.numeric(as.character(AC_by_gene$AC))
    
  }
}

ggplot(AC_by_gene,aes(x=AC,,color=disease, linetype=type)) +
  geom_density(aes(y=..count..)) +
  ylab("gene count")+
  scale_x_continuous(breaks = seq(0, 30, 2), limits=c(0,30)) 

# MORE AFFECTED GENES BY PATHOLOGY  --------------------------------------------

#order table by AC in descending order
AC_by_gene_ord <-AC_by_gene[order(-AC_by_gene$AC),]



for (ds in ds_list){
  AC_by_gene_ord_ds <- AC_by_gene_ord[AC_by_gene_ord$disease==ds,]
  
  aux <- rbind(AC_by_gene_ord_ds[AC_by_gene_ord_ds$type=="DEL",][1:15,],
               AC_by_gene_ord_ds[AC_by_gene_ord_ds$type=="DUP",][1:15,])
  
  p <- ggplot(aux, aes(x=reorder(gene,AC), y=AC, fill=type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("cornflowerblue","#FF6633")) +
    xlab("Gen") +
    ggtitle(paste("Genes más afectados en: ",ds)) + 
    theme(legend.title = element_blank()) 

  p <- p + coord_flip()
  
  dp_all <- plot(p + facet_grid(rows = vars(type),scales = "free_y"))
  dp_all
}


