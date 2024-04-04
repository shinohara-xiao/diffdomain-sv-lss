library(dplyr)

ADDBH <- function(filename)
{
  path = '/lanec2_home/zhangx/DiffCompare/result/'
  tads <- read.table(paste(path,filename,'.txt',sep = ''),header = FALSE,skip=25,fill = TRUE,sep = '\t')
  tads = rename(tads,chr = V1,start = V2,end = V3,region = V4,lambdan = V5,pvalue = V6,nd = V7)
  tads$BH <-p.adjust(tads$pvalue,method  = 'BH')
  write.table(tads,paste(path,filename,'_BH.txt',sep = ''),sep = '\t')
}

ADDBH('NHA_vs_DIPGXIII_10000')
ADDBH('GM12878_vs_K562_10000_f1b5')
ADDBH('NHA_vs_DIPG007_10000')

ADDBH('GM12878_vs_IMR90_25000')
ADDBH('HMEC_vs_NHEK_25000')
ADDBH('GM12878_vs_K562_25000')
ADDBH('NHA_vs_DIPG007_25000')
ADDBH('NHA_vs_DIPGXIII_25000')

ADDBH('GM12878_vs_K562_25000_to1')
ADDBH('NHA_vs_DIPG007_25000_to1')
ADDBH('NHA_vs_DIPGXIII_25000_to1')
