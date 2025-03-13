setwd("~/Desktop/normal/revisit_freq_filtering_May/batch2/jan_batch2_ssm/")
# dir <- list.dirs()[-c(1:3,18:24)]
# dir <- gsub("\\.\\/","",dir)
#vaf_cutoff = 7
library(ggplot2)
library(data.table)
library(plyr)

files <- list.files(pattern =  "_*.ssm")
files <- files[-c(1)]
gene <- c()
loci <- c()

for (i in 1:length(files)) {
  

temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)

locus <- temp_file$locus
names(locus) <- temp_file$name


gene <- c(gene,unique(temp_file$name))
loci <- c(loci,locus)

  
}

a <- data.frame(loci =loci,gene=names(loci))
a_gene <- data.frame(table(gene))

#b <- data.frame(table(gene))
b <- ddply(a,.(loci,gene),nrow)

#70%
c <- b[b$V1>14,]
write.csv(c,file="removed_loci_b2.csv",row.names = F)
c <- as.character(c$loci)




for (i in 1:length(files)) {
  fn <- gsub(".ssm","_freq",files[i])
  temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)
temp_file <- temp_file[!temp_file$locus %in% c,]
temp_file <- temp_file[!temp_file$name == "POLE",]
temp_file <- temp_file[!temp_file$name == "ATR",]
  temp_file$id <- paste("s",seq(0,nrow(temp_file)-1,1),sep = "")
  write.table(temp_file,file = paste(fn,".ssm",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
  
  
}
