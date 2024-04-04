library(ggplot2)
library(plyr)
library(tidyverse)

subtypes <- c('complex','loss','merge','split','strength','zoom')
celltype <- c('K562','DIPG007','DIPGXIII')
svtype <- c('deletion','duplication','dup55','dup33')

K562_del <- c(3,8,0,3,8,1); K562_dup<-c(0,8,4,2,9,1); K562_dup55 <-c(0,1,0,2,6,2); K562_dup33<-c(14,10,4,1,40,3)
DIPG007_del <- c(1,2,0,5,6,6); DIPG007_dup <- c(0,1,0,12,7,13) ; DIPG007_dup55 <- c(1,2,1,6,7,4) ; DIPG007_dup33 <-c(0,0,1,4,4,3)
DIPGXIII_del <- c(0,11,1,0,1,1); DIPGXIII_dup <- c(0,1,2,0,0,2) ; DIPGXIII_dup55 <- c(0,9,1,3,1,0); DIPGXIII_dup33 <- c(0,3,0,0,0,0)

Celltype <- rep(celltype,each = 24)
Svtype <- rep(svtype,each = 6)
Subtypes <- rep(subtypes,times = 4)
Value = c(K562_del,K562_dup,K562_dup55,K562_dup33,DIPG007_del,DIPG007_dup,DIPG007_dup55,DIPG007_dup33,
          DIPGXIII_del,DIPGXIII_dup,DIPGXIII_dup55,DIPGXIII_dup33)

pvalue  <- c()
for(i in 1:length(Value))
{
  if(i <= 6){ pvalue = append(pvalue,Value[i]/sum(K562_del)) }
  else if(i <=12){ pvalue = append(pvalue,Value[i]/sum(K562_dup)) }
  else if(i<=18){ pvalue = append(pvalue,Value[i]/sum(K562_dup55)) }
  else if (i <= 24) { pvalue = append(pvalue,Value[i]/sum(K562_dup33)) }

  else if (i <=30) { pvalue = append(pvalue,Value[i]/sum(DIPG007_del)) }
  else if(i <=36){ pvalue = append(pvalue,Value[i]/sum(DIPG007_dup)) }
  else if(i<=42){ pvalue = append(pvalue,Value[i]/sum(DIPG007_dup55)) }
  else if (i <= 48) { pvalue = append(pvalue,Value[i]/sum(DIPG007_dup33)) }

  else if (i <=54) { pvalue = append(pvalue,Value[i]/sum(DIPGXIII_del)) }
  else if(i <=60){ pvalue = append(pvalue,Value[i]/sum(DIPGXIII_dup)) }
  else if(i<=66){ pvalue = append(pvalue,Value[i]/sum(DIPGXIII_dup55)) }
  else if (i <= 72) { pvalue = append(pvalue,Value[i]/sum(DIPGXIII_dup33)) }
}
pvalue = round(pvalue*100,2)
df <- data.frame(Celltype,Svtype,Subtypes,Value,pvalue)

# Change the order --- use this order
df$Subtypes <- factor(df$Subtypes,level = c('strength','loss','split','zoom','merge','complex'))
df$Celltype <- factor(df$Celltype,levels = c('K562','DIPG007','DIPGXIII'))
df$Svtype <- factor(df$Svtype,levels = c('deletion','duplication','dup55','dup33'))
ggplot(data = df,aes(x = Svtype,y = pvalue,fill = Subtypes)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Celltype) +
  theme(panel.background = element_rect(fill = 'white'))+
  labs(x = 'Svtypes',y = 'percentage(%)') +
  scale_fill_manual(values = c('#619CFF','#B79F00','#00BFC4','#C77CFF','#00BA38','#F8766D'))

"""
# Change the order
df$Subtypes <- factor(df$Subtypes,level = c('strength','loss','merge','complex','split','zoom'))
df$Celltype <- factor(df$Celltype,levels = c('K562','DIPG007','DIPGXIII'))
df$Svtype <- factor(df$Svtype,levels = c('deletion','duplication','dup55','dup33'))
ggplot(data = df,aes(x = Svtype,y = pvalue,fill = Subtypes)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Celltype) +
  theme(panel.background = element_rect(fill = 'white'))+
  labs(x = 'Svtypes',y = 'percentage(%)') +
  scale_fill_manual(values = c('#619CFF','#B79F00','#00BA38','#F8766D','#00BFC4','#C77CFF'))


# X轴 为celltypes - not used
df$Subtypes <- factor(df$Subtypes,level = c('strength','loss','merge','complex','split','zoom'))
ggplot(data = df,aes(x = Celltype,y = pvalue,fill = Subtypes)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Svtype) +
  theme(panel.background = element_rect(fill = 'white'))+
  labs(x = 'Celltypes',y = 'percentage(%)') +
  scale_fill_manual(values = c('#619CFF','#B79F00','#00BA38','#F8766D','#00BFC4','#C77CFF'))
"""
