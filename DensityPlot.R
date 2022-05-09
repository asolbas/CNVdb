library(ggplot2)

#BASE DE DATOS NUEVA -----------------------------------------------------------
df <- read.table(file = '~/bioinfo/fjd/MAF_CNV_FJD/results/MAF_CNV_Database_PC.tsv', 
                 sep = '\t', header = TRUE)

#delete X and Y chromosoms
df<- df[!(df$X0=="X" | df$X0=="Y"),]

df<-cbind(length=df$X2-df$X1,df)

p <- ggplot(data=df,aes(x=length)) +
  geom_density(alpha=0.5, color="cornflowerblue", fill="cornflowerblue")
ggsave("~/tblab/ana/database/Figures/density.png",bg="white")
plot(p)

#Log scale
p <- ggplot(data=df,aes(length)) +
  geom_density(alpha=0.5, color="cornflowerblue", fill="cornflowerblue") +
  scale_x_log10() +
  labs(x="size (pb)")

ggsave("~/tblab/ana/database/Figures/density_log.png",bg="white")
plot(p)

df<-df[df$length<1000,]

p <- ggplot(data=df,aes(x=length)) +
  geom_density(alpha=0.5, color="cornflowerblue", fill="cornflowerblue")
ggsave("~/tblab/ana/database/Figures/density_zoom.png",bg="white")
plot(p)



#BASE DE DATOS ANTIGUA ---------------------------------------------------------

df <- read.table(file = '~/tblab/raquel/CES_relanzados_CNVs/MAF_CNV_FJD_14022020/results/MAF_CNV_Database.tsv', sep = '\t', header = TRUE)

df<-cbind(length=df$end-df$start,df)

DUP <- df[df$AC_DUP>=1,]
DEL <- df[df$AC_DEL>=1,]

#he quitado las variantes con un solo nt y más de 5000 nt
DUP <-cbind.data.frame(type="DUP",length=DUP[DUP$length>1 & DUP$length<5000 ,]$length)
DEL <-cbind.data.frame(type="DEL",length=DEL[DEL$length>1 & DEL$length<5000,]$length)
cnvs<-rbind.data.frame(DUP,DEL)

data_cnvs=as.data.frame(cnvs)

ggplot(data=data_cnvs,aes(x=length, group=type,color=type,fill=type)) +
  geom_density(alpha=0.5)

#Probando a quitar las comunes
df <- read.table(file = '~/tblab/raquel/CES_relanzados_CNVs/MAF_CNV_FJD_14022020/results/MAF_CNV_Database.tsv', sep = '\t', header = TRUE)

df<-cbind(length=df$end-df$start,df)

DUP <- df[df$AC_DUP>=1 & df$AC_DEL==0,]
DEL <- df[df$AC_DEL>=1 & df$AC_DUP==0,]

#he quitado las variantes con un solo nt y más de 5000 nt
DUP <-cbind.data.frame(type="DUP",length=DUP[DUP$length>1 & DUP$length<5000 ,]$length)
DEL <-cbind.data.frame(type="DEL",length=DEL[DEL$length>1 & DEL$length<5000,]$length)
cnvs<-rbind.data.frame(DUP,DEL)

data_cnvs=as.data.frame(cnvs)

ggplot(data=data_cnvs,aes(x=length, group=type,color=type,fill=type)) +
  geom_density(alpha=0.5)


den_DUP <- density(DUP[,2]) # returns the density data
plot(den_DUP) # plots the results

hist_DUP <- hist(DUP$length)
plot(hist_DUP)

#quitando los valores mayores de 5000
hist_DUP <- hist(DUP[DUP$length<5000,]$length,breaks=50)
plot(hist_DUP)

#quitando los valores mayores de 5000
den_DUP <- density(DUP[DUP$length<5000,]$length,breaks=50)
plot(den_DUP)

#Considerando el AC

df <- read.table(file = '~/tblab/raquel/CES_relanzados_CNVs/MAF_CNV_FJD_14022020/results/MAF_CNV_Database.tsv', sep = '\t', header = TRUE)
df<-cbind(length=df$end-df$start,df)
DUP <- df[df$AC_DUP>=1 & df$length>1 & df$length<5000,]
DEL <- df[df$AC_DEL>=1 & df$length>1 & df$length<5000,]

DUP <-cbind.data.frame(type="DUP",length=rep(DUP$length,DUP$AC_DUP)) #e.g. si el AC es 3 lo repite 3 veces
DEL <-cbind.data.frame(type="DEL",length=rep(DEL$length,DEL$AC_DEL))

cnvs<-rbind.data.frame(DUP,DEL)
data_cnvs=as.data.frame(cnvs)

ggplot(data=data_cnvs,aes(x=length, group=type,color=type,fill=type)) +
  geom_density(alpha=0.5)

#SIN SEPARAR DUPLICACIONES DE DELECIONES (sin considerar AC)
df <- read.table(file = '~/tblab/raquel/CES_relanzados_CNVs/MAF_CNV_FJD_14022020/results/MAF_CNV_Database.tsv', sep = '\t', header = TRUE)
df<-cbind(length=df$end-df$start,df)
df<-df[df$length>1 & df$length<5000,]

ggplot(data=df,aes(x=length)) +
  geom_density(alpha=0.5, color="cornflowerblue", fill="cornflowerblue")

#BASE DE DATOS ACTUALIZADA
df <-  read.table(file = '~/bioinfo/NOBACKUP/raquel/MAF_CNV_FJD/results/MAF_CNV_Database_regions.tsv', sep = '\t', header = TRUE)

df<-cbind(length=df$end-df$start,df)

DUP <- df[df$AC_DUP>=1,]
DEL <- df[df$AC_DEL>=1,]


DUP <-cbind.data.frame(type="DUP",length=DUP$length)
DEL <-cbind.data.frame(type="DEL",length=DEL$length)
cnvs<-rbind.data.frame(DUP,DEL)

data_cnvs=as.data.frame(cnvs)

#considerando dup y del por separado 
p_sep <- ggplot(data=data_cnvs,aes(x=length, group=type,color=type,fill=type)) +
  geom_density(alpha=0.5)

plot(p_sep)

#eliminando las que tienen más de 5000 nt 
df <-  read.table(file = '~/bioinfo/NOBACKUP/raquel/MAF_CNV_FJD/results/MAF_CNV_Database_regions.tsv', sep = '\t', header = TRUE)

df<-cbind(length=df$end-df$start,df)
DUP <- df[df$length<5000,]
DEL <- df[df$length<5000,]

DUP <-cbind.data.frame(type="DUP",length=DUP$length) 
DEL <-cbind.data.frame(type="DEL",length=DEL$length)

cnvs<-rbind.data.frame(DUP,DEL)
data_cnvs=as.data.frame(cnvs)

p_5000 <- ggplot(data=data_cnvs,aes(x=length, group=type,color=type,fill=type)) +
  geom_density(alpha=0.5)
plot(p_5000)

#Juntando DUP y DEL

df <-  read.table(file = '~/bioinfo/NOBACKUP/raquel/MAF_CNV_FJD/results/MAF_CNV_Database_regions.tsv', sep = '\t', header = TRUE)

df<-cbind(length=df$end-df$start,df)
df<-df[df$length>1 & df$length<5000,]

p <- ggplot(data=df,aes(x=length)) +
  geom_density(alpha=0.5, color="cornflowerblue", fill="cornflowerblue")

plot(p)




