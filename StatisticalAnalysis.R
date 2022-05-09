library("ggvenn")
library("stringr")
library(ggplot2)

#INPUT DATA---------------------------------------------------------------------

#df <- read.table(file = '~/tblab/raquel/CES_relanzados_CNVs/MAF_CNV_FJD_14022020/results/MAF_CNV_Database.tsv', sep = '\t', header = TRUE)
#df <-  read.table(file = '~/bioinfo/NOBACKUP/raquel/MAF_CNV_FJD/results/MAF_CNV_Database_regions.tsv', sep = '\t', header = TRUE)
#df <- read.table(file = '~/bioinfo/fjd/MAF_CNV_FJD/results/MAF_CNV_Database.tsv', sep = '\t', header = TRUE)
df <- read.table(file = '~/bioinfo/fjd/MAF_CNV_FJD/results/MAF_CNV_Database_PC.tsv', sep = '\t', header = TRUE)

#delete X and Y chromosoms
df<- df[!(df$X0=="X" | df$X0=="Y"),]

#add an id column
#each cnv will have an id

df<-cbind(ID=seq.int(nrow(df)),length=df$X2-df$X1,df)

DUP <- df[df$AC_DUP>=1,]
DEL <- df[df$AC_DEL>=1,]
#head(DEL)

#VENN DIAGRAM ------------------------------------------------------------------

#venn diagram: se han considerado solo variantes únicas
x<-list("Duplications"=DUP$ID,"Deletions"=DEL$ID)
ggvenn(x,fill_alpha = 0.4,text_size = 4.5,stroke_size = 0.5, fill_color = c("#FF6600", "#3366FF"))
ggsave("~/tblab/ana/database/Figures/VennDiag.png", bg = "white")

#parsear las enfermedades
i=14
while(i<ncol(DUP)){
  disease<-unlist(strsplit(colnames(DUP[i]), split='_', fixed=TRUE))[1] #disease name
  x<-list("DUP"=DUP[DUP[,i]>=1,]$ID,"DEL"=DEL[DEL[,i+3]>=1,]$ID)
  ggvenn(x,fill_alpha = 0.4,text_size = 4.5,stroke_size = 0.5, fill_color = c("#FF6600", "#3366FF"))
  ggsave(paste("~/tblab/ana/database/Figures/",disease,".png"), bg = "white")
  i=i+14
}

#para considerar el AC (todas las variantes) -->da lo mismo que antes
# DUP_all=rep(DUP$ID,times=DUP$AC_DUP)
# DEL_all=rep(DEL$ID,times=DEL$AC_DEL)
# x<-list("DUP"=DUP_all,"DEL"=DEL_all)
# ggvenn(x,fill_alpha = 0.4,text_size = 4.5,stroke_size = 0.5)

#HISTOGRAM ---------------------------------------------------------------------

#histograma (considera el AC)
AC_df <- data.frame(
  type=c("DUP","DEL"),
  count=c(sum(df$AC_DUP),sum(df$AC_DEL))
)

p<-ggplot(AC_df, aes(x=type, y=count, fill=type)) +
  geom_col(width=0.4) +
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(values=c("cornflowerblue","#FF6633")) +
  geom_text(aes(label=count), vjust=-1, color="black",
            position = position_dodge(0.9), size=2.75)
plot(p)
ggsave("~/tblab/ana/database/Figures/count_total.png")




#Plot histogram by diseases
AC_disease <- data.frame(type=character(), disease=character(),
                         AC=double(), stringsAsFactors=FALSE)

#diseases with less than 100 cases are not considered
few_cases_ds <- c("Miopatias","Metabolicas","Inflamatoria","Esterilidad",
                  "Endocrinologica", "Dermatologicas","Prenatal","Cancer", "Varios")


AC_DUPtot<-0
AC_DELtot<-0

for (i in seq(14,ncol(df), by=14)) {
  disease_col=names(df)[i] #disease col name
  disease<-unlist(strsplit(disease_col, split='_', fixed=TRUE))[1] #disease name
  if (!( disease %in% few_cases_ds)){
    #parse duplications
    AC_disease[nrow(AC_disease) + 1,] <- c("DUP", disease, sum(df[,i]))
    AC_DUPtot<-AC_DUPtot+sum(df[,i])
    #parse deletions
    AC_disease[nrow(AC_disease) + 1,] <- c("DEL", disease, sum(df[,i+3]))
    AC_DELtot<-AC_DELtot+sum(df[,i+3])
  }
}
#Convert AC to numeric value
AC_disease$AC <- as.numeric(as.character(AC_disease$AC))

ggplot(data=AC_disease, aes(x=disease, y=AC, fill=type)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=AC), vjust=-1, color="black",
            position = position_dodge(0.9), size=3)+
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=c("cornflowerblue","#FF6633")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("~/tblab/ana/database/Figures/count_by_disease.png", bg="white")

#print histogram AC total (filtered by diseases with more than 100 cases)

AC_df_filt <- data.frame(
  type=c("DUP","DEL"),
  count=c(AC_DUPtot,AC_DELtot)
)

p<-ggplot(AC_df_filt, aes(x=type, y=count, fill=type)) +
  geom_col(width=0.4) +
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(values=c("cornflowerblue","#FF6633")) +
  geom_text(aes(label=count), vjust=-1, color="black",
            position = position_dodge(0.9), size=3)
plot(p)
ggsave("~/tblab/ana/database/Figures/count_total_filtered_by_disease.png",bg="white")


# REGIONS COUNT BY AC ---------------------------------------------------------
DUP$AC<-DUP$AC_DUP
DUP$type="DUP"
DEL$AC=DEL$AC_DEL
DEL$type="DEL"

regions_by_AC<-rbind(DUP,DEL)

regions_by_AC$AC<-as.numeric(as.character(regions_by_AC$AC))

ggplot(regions_by_AC,aes(x=AC,color=type, linetype=type)) +
  geom_density(aes(y=..count..)) +
  scale_x_continuous(breaks = seq(0, 30, 2), limits=c(1,30))+
  ylab("region count")

ggplot(regions_by_AC,aes(x=AC,color=type, fill=type)) +
  geom_histogram(aes(y=..count..),binwidth=1, alpha=0.2) +
  scale_x_continuous(breaks = seq(0, 30, 2), limits=c(0,30))+
  ylab("region count")

ggsave("~/tblab/ana/database/Figures/regions_by_AC.png",bg="white",width=10,height = 6)



#REGIONS COUNT BY AC BY DISEASE------------------------------------------------

regions_by_AC <- data.frame( AC=character(), number_regions=character(), type=character(),disease=character())

for (i_ds in seq(14,(ncol(df)-2), by=14)) {
  disease_col_DUP=names(df)[i_ds] #get column name
  disease_col_DEL=names(df)[i_ds] #get column name
  disease<-unlist(strsplit(disease_col_DUP, split='_', fixed=TRUE))[1] #disease name
  
  if (!( disease %in% few_cases_ds)){
    
    #Count number of regions by AC number
    DUP_regions<-as.data.frame(table(unlist(df[disease_col_DUP])))
    colnames(DUP_regions) <- c("AC","number_regions")
    DUP_regions$type<-"DUP"
    DUP_regions$disease<-disease
    
    DEL_regions<-as.data.frame(table(unlist(df[disease_col_DEL])))
    colnames(DEL_regions) <- c("AC","number_regions")
    DEL_regions$type<-"DEL"
    DEL_regions$disease<-disease
    
    regions_by_AC<-rbind(regions_by_AC,DUP_regions,DEL_regions)
  }
  
}

regions_by_AC$AC<-as.numeric(as.character(regions_by_AC$AC))

ggplot(regions_by_AC,aes(x=AC,y=number_regions,color=disease, linetype=type)) +
  geom_line() +
  ylab("Number of regions")+
  scale_x_continuous(breaks = seq(0, 30, 2), limits=c(1,30)) +
  ylim(0,3000)

ggsave(paste("~/tblab/ana/database/Figures/regions_by_AC_disease.png"), bg="white",width=10,height=6)






#BOXPLOT-----------------------------------------------------------------------
#import file with confidence variants
#check.names=FALSE to avoid putting an X in front of colnames that are numeric
#variants <- read.table(file = '~/bioinfo/fjd/MAF_CNV_FJD/results/conf_DB_variantsvx.tsv', sep = '\t', header = TRUE,check.names=FALSE)
#variants <- read.table(file = '~/bioinfo/fjd/MAF_CNV_FJD/results/conf_DB_variants_last.tsv', sep = '\t', header = TRUE,check.names=FALSE)
#variants <- read.table(file = '~/bioinfo/fjd/MAF_CNV_FJD/results/conf_DB_variantsv2.tsv', sep = '\t', header = TRUE,check.names=FALSE)
#variants <- read.table(file = '~/bioinfo/fjd/MAF_CNV_FJD/results/conf_DB_variants_filtered_outliers.tsv', sep = '\t', header = TRUE,check.names=FALSE)
variants <- read.table(file = '~/bioinfo/fjd/MAF_CNV_FJD/results/conf_DB_variants_filtered_07_05.tsv', sep = '\t', header = TRUE,check.names=FALSE)

#import file with metadata
#metadata <- read.table(file = '~/bioinfo/fjd/MAF_FJD_v3.0/metadata/date_2021_12_29/mymetadatapathology_uniq_2021_12_29.txt', sep = '\t', header = TRUE)
metadata <- read.table(file = '~/bioinfo/fjd/MAF_FJD_v3.0/metadata/date_2022_03_28/mymetadatapathology_uniq_2022_03_28.txt', sep = '\t', header = TRUE)

new_names<-vector()
for (name in colnames(variants)){
  new_names<-c(new_names,unlist(strsplit(name, split='-CES', fixed=TRUE))[1])
}
colnames(variants)<-new_names

data <- data.frame(SAMPLE=character(), disease=character(), type=character(), AC=character())

#parse variants file
for (i in seq(4,(ncol(variants)), by=1)){
  sample <- colnames(variants[i])
  #sample <- unlist(strsplit(colnames(variants[i]), split='-CES', fixed=TRUE))[1]
  disease <- unique(metadata[metadata$SAMPLE==sample,]$Categoria)
  #AC_DUP <- length(which(variants[i]==1)) 
  #AC_DEL <- length(which(variants[i]==-1))
  AC_DUP <- length(which(variants[i]==1)) + 2*length(which(variants[i]==2))
  AC_DEL <- length(which(variants[i]==-1)) + 2*length(which(variants[i]==-2))
  
  if (! (disease %in% few_cases_ds) ){
    #if ( (length(disease) != 0) & (! (disease %in% few_cases_ds) )){
    data[nrow(data) + 1,] = c(sample, disease, "DUP", AC_DUP)
    data[nrow(data) + 1,] = c(sample, disease, "DEL", AC_DEL)
  }
}

data$AC <- as.numeric(as.character(data$AC))
#data<- na.omit(data)

#outlier.shape=NA to remove outliers
ggplot(data=data, aes(x=disease, y=AC, fill=type)) +
  geom_boxplot() +
  #geom_boxplot(outlier.shape = NA) +
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(values=c("cornflowerblue","#FF6633"))  +
  #ylim(0,30) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


ggsave("~/tblab/ana/database/Figures/boxplot.png", bg="white")

# VIOLIN PLOT -----------------------------------------------------------------
ggplot(data=data, aes(x=disease, y=AC, fill=type)) +
  geom_violin() +
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(values=c("cornflowerblue","#FF6633"))  +
  #ylim(0,270) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("~/tblab/ana/database/Figures/violinplot.png")

# SWARM PLOT ----------------------------------------------------------------
library(ggbeeswarm)
ggplot(data=data, aes(x=disease, y=AC, col=type)) +
  #geom_beeswarm(dodge.width = 0.5) +
  geom_quasirandom(size=1) +
  #scale_color_brewer(palette="Paired")  +
  scale_color_manual(values=c("cornflowerblue","#FF6633"))  +
  #ylim(0,270) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("~/tblab/ana/database/Figures/swarmplot.png")


#STATISTICAL ANALYSIS ---------------------------------------------------------

##Kruskal Wallis test by diseases----------------------------------------------------

#prepare data for DUP and DEL
data_DUP <- data[data$type == "DUP",]
data_DEL <- data[data$type == "DEL",]

#convert diseases into factors
data_DUP$disease <- as.factor(data_DUP$disease)
data_DEL$disease <- as.factor(data_DEL$disease)

#levels(data_DUP$disease)

#Perform test 
kruskal.test(AC~disease, data=data_DUP)
kruskal.test(AC~disease, data=data_DEL)

#Multiple pairwise-comparison between groups 
pairwise.wilcox.test(data_DUP$AC,data_DUP$disease,p.adjust.method = "BH")
pairwise.wilcox.test(data_DEL$AC,data_DEL$disease,p.adjust.method = "BH")

# pairwise.t.test(data_DUP$AC,data_DUP$disease,p.adjust.method = "BH")
# pairwise.t.test(data_DEL$AC,data_DEL$disease,p.adjust.method = "BH")


## Test duplications and deletions ---------------------------------------------
#Boxplot
data_DUP$type <- as.factor(data_DUP$type)
data_DEL$type <- as.factor(data_DEL$type)

ggplot(data=data, aes(x=type, y=AC, fill=type)) +
  geom_boxplot() +
  #geom_boxplot(outlier.shape = NA) +
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(values=c("cornflowerblue","#FF6633"))  +
  #ylim(0,30) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave("~/tblab/ana/database/Figures/boxplot_type.png", bg="white")

ggplot(data=data, aes(x=type, y=AC, fill=type)) +
  geom_violin() +
  #geom_boxplot(outlier.shape = NA) +
  #scale_fill_brewer(palette="Paired") +
  scale_fill_manual(values=c("cornflowerblue","#FF6633"))  +
  #ylim(0,30) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

#run paired t test
t.test(data_DEL$AC, data_DUP$AC, alternative = "two.sided",
       mu = 0, paired = TRUE, conf.level = 0.95) 

## Test dup and del by disease ------------------------------------------------
library(tidyverse)
library(rstatix)
library(ggpubr)
library(data.table)

data$type<-as.factor(data$type)
data$disease<-as.factor(data$disease)


# stat.test <- data %>%
#   group_by(disease) %>%
#   t_test(AC ~ type, alternative = "two.sided",
#          mu = 0, paired = TRUE, conf.level = 0.95) %>%
#   adjust_pvalue(method = "BH") %>%
#   add_significance()
# stat.test

stat.test <- data %>%
  group_by(disease) %>%
  wilcox_test(AC ~ type, alternative = "two.sided",
         paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

##normality test--------------------------------------------------------------
shap.test_DUP <- data_DUP %>%
  group_by(disease) %>%
  shapiro_test(AC)
shap.test_DUP

shap.test_DEL <- data_DEL %>%
  group_by(disease) %>%
  shapiro_test(AC)
shap.test_DEL

shap_type.test <- data %>%
  group_by(type) %>%
  shapiro_test(AC)
shap_type.test

shap.test_ds <- data %>%
  group_by(disease) %>%
  shapiro_test(AC)
shap.test_ds



##Summary satistics ------------------------------------------------------------
summary(data)
summary(data_DUP)
summary(data_DEL)
#count number of samples by disease
table(data_DUP$disease)


