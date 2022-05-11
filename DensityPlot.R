library(ggplot2)

df <- read.table(file = 'MAF_CNV_Database_PC.tsv', 
                 sep = '\t', header = TRUE)

#delete X and Y chromosoms
df<- df[!(df$X0=="X" | df$X0=="Y"),]

df<-cbind(length=df$X2-df$X1,df)

p <- ggplot(data=df,aes(x=length)) +
  geom_density(alpha=0.5, color="cornflowerblue", fill="cornflowerblue")
plot(p)

#Log scale
p <- ggplot(data=df,aes(length)) +
  geom_density(alpha=0.5, color="cornflowerblue", fill="cornflowerblue") +
  scale_x_log10() +
  labs(x="size (pb)")
plot(p)

df<-df[df$length<1000,]

p <- ggplot(data=df,aes(x=length)) +
  geom_density(alpha=0.5, color="cornflowerblue", fill="cornflowerblue")
plot(p)
