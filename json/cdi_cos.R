library(jsonlite)
setwd("~/Desktop/pt/pairtree/final_withoutcnv_4/json/")
files <- list.files(pattern = "*.json")


sink("sdi.txt")
sink("cdi.txt")
sink("cdi_sd.txt")
sink("sdi_sd.txt")
sink("can.txt")
sink('cos.txt')
library(coop)
library(ggplot2)



for (i in 1:length(files)) {
  
json_file=files[i]
a <- fromJSON(json_file,flatten = T)
#grep("phi",a,perl = T)
phi <- round(as.matrix(a[[1]]),2)
phi[phi==0] <- 0.01
sample_name <- gsub("\\..*","",json_file)
names <- as.character(a[["samples"]])
sdi <- as.numeric(a[["sdi"]])
sdi_sd <- sd(sdi)
cdi <- as.numeric(a[["cdi"]])
cdi_sd <- sd(cdi)
# 
# distance <- function(tmp1,tmp2){
# tmp1 <- as.numeric(tmp1)
# tmp2 <- as.numeric(tmp2)
#   for (p in 1:length(tmp1)) {
#     r0 <- abs(tmp1[p]-tmp2[p])/max(tmp1[p],tmp2[p])
#     if(p==1){
#       r=r0
#     }else{
#     r <- r+r0
#   }
#   
#   }
#   return(r)
# }
if(length(sdi)>2){
cos_id <- round(dist(t(phi),method = "canberra"),6)
cos_d <- sd(cos_id)
cos_n <- round(1-cosine(phi),5)[lower.tri(round(1-cosine(phi),5))]
cos_2 <- sd(cos_n)
}else{
  cos_d <- NA
  cos_2 <- NA
}
df1 <- data.frame(samples=names,cd_index=cdi)
df2 <- data.frame(samples=names,sd_index=sdi)

pdf(paste(sample_name,"_plots.pdf",sep = "")) 
a = ggplot(df1,aes(samples,cd_index)) + geom_col() + ggtitle(paste("CDI ",sample_name),subtitle = paste("sd=",cdi_sd))
b = ggplot(df2,aes(samples,sd_index)) + geom_col() + ggtitle(paste("SDI ",sample_name),subtitle = paste("sd=",sdi_sd))
print(a)
print(b)
dev.off()

write(paste(names,sdi,sep = " "),file = "sdi.txt",append = T)

write(paste(names,cdi,sep = " "),file = "cdi.txt",append = T)
write(paste(sample_name,cdi_sd),file = "cdi_sd.txt",append = T)
write(paste(sample_name,sdi_sd),file = "sdi_sd.txt",append = T)
write(paste(sample_name,cos_d),file = "can.txt",append = T)
write(paste(sample_name,cos_2),file = "cos.txt",append = T)

}

sink()
sink()
sink()
sink()
sink()
sink()

rm(list = ls()) 



library(tidyr)
cdi_sd <- read.delim("cdi_sd.txt",sep = " ",header = F)
sdi_sd <- read.delim("sdi_sd.txt",sep = " ",header = F)
cos <- read.delim("can.txt",sep = " ",header = F)
coc <- read.delim("cos.txt",sep = " ",header = F)
df <- cbind(cdi_sd,sdi_sd[,2],cos[,2],coc[,2])
colnames(df) <- c("patient","cdisd","sdisd","can","cos")
df <- na.omit(df)

a <- gather(df,"patient")
colnames(a)[1] <- "label"
a$patient <- rep(df[,1],4)

t <- paste("cdi vs sdi:", round(cor(df$cdisd,df$sdisd),2)," | cdi vs can:",round(cor(df$can,df$cdisd),2), " | sdi vs can",round(cor(df$can,df$sdisd),2))

ggplot(a,aes(x=patient,y=value,color=label))+geom_line(aes(group = label),)  + theme_minimal()+ggtitle("Comparison of standard devistion in cdi, sdi and canberra distance", subtitle = t )

                                                                  