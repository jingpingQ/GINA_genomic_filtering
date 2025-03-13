setwd("~/GINA_genomic_data")
dir <- list.dirs()[-1]
dir <- gsub("\\.\\/","",dir)
library(ggplot2)
library(data.table)


##################################################################################
##  This script is wrote for genomic data (OCAv3) filtering, the data structure is
##  foldered by each of the 34 patients. In each patient's folder, we have one normal
##  sample named with "full" and several other tumor samples. The tumor samples are 
##  labeled with Block number Bxx and foci number Inv/DCIS xx. 
##  This script will loop each of the folder, record the SNPs in normal file, and 
##  substract the SNPS from other tumor samples. For quality control, we set a VAF 
##  filter at 7%, coverage filter at 100 reads, and p-value at 0.05. After filtering, 
##  the output is formatted for PairTree algorithm. See examples in ssm folder.
##  This is the first round of data filtering. Further filtering was done for each 
##  batch using donor affected frequency on each gene.
##  CNV and SNV are recorded separately.

##################################################################################



for (d in 1:length(dir)) {
  
  ll=dir[d]
  
  setwd(paste("~/Desktop/normal/octfull/",ll,sep = ""))
  files <- list.files(pattern =  "*.tsv")
  snps <- list()
  vaf <- list()
  cov <- list()
  #loop number of files in each patient folder.
  for (i in 1:length(files)) {
    #check the file name, if it has "full" in file name, it is germline.
    name <- gsub("\\..*","",files[i])
    #print in console
    cat("sample: ")
    cat(paste(gsub("\\..*","",files[i]),"\n"))

    
    #read in file
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)

    
    if(length(grep("full",name))>0){
      
      #select columns to work with
      temp_file <- temp_file[,c(1,2,5,6,8,9,11,13)]
      temp_file$parse <- temp_file$locus
      #check if in dbSNP ref database
      temp_file[temp_file$dbsnp != "",]$dbsnp<-"YES"

        #file manipulation. some of the entry has two VAFs for polymorphisims, add them all.
        if(length(grep(",",temp_file$VAF))>0){
        for (k in grep(",",temp_file$VAF)) {
          k_temp <- sum(as.numeric(gsub(".*\\=","",unlist(strsplit(temp_file$VAF[k],split = ",")))))
          temp_file$VAF[k] <- k_temp
        }
        
        temp_file$VAF <- as.numeric(temp_file$VAF)
      }
      
      #filter: p value, coverage, variant allele freq.
      temp_file$VAF <- as.numeric(temp_file$VAF)
      temp_file <- temp_file[temp_file$pvalue<0.05,]
      temp_file <- temp_file[temp_file$coverage>100,]
      temp_file <- temp_file[temp_file$VAF>7,]
      temp_file$VAF <- temp_file$VAF/100 #modify VAF
      temp_germ <- temp_file
      
      
      #record the germline.
      germ_snp <- temp_file$parse
      names(germ_snp) <- temp_file$gene
      snps[[i]] <- germ_snp
      names(snps)[i] <- "germ"
      vaf[[i]] <- temp_file$VAF
      names(vaf)[i] <- "germ"
      cov[[i]] <- temp_file$coverage
      names(cov)[i] <- "germ"
      
      next
    }
    
    
    #print for check up.
    cat("_unfiltered total: ")
    cat(nrow(temp_file),"\n")

    
    #remove CNV, CNV will be recorded seperately.
    temp_file$locus <- paste(temp_file$chr,temp_file$start)
    temp_file <- temp_file[temp_file$type != "CNV",]
    temp_file$parse <- temp_file$locus
    temp_file$VAF <- temp_file$VAF/100
    
    
    
    temp_snp <- temp_file$parse[!temp_file$parse %in% unname(unlist(snps["germ"]))]
    names(temp_snp) <- temp_file$gene[!temp_file$parse %in% unname(unlist(snps["germ"]))]
    temp_snp <- c(temp_snp,add_diff)
    snps[[i]] <- temp_snp
    names(snps)[i]<-name
    
    
    #record VAF
    vaf[[i]] <- temp_file$VAF
    names(vaf)[i] <- name
    cov[[i]] <- temp_file$coverage
    names(cov)[i] <- name
    
    
    
  }
  
  #parse germline
  germ <- snps[["germ"]]
  
  snps["germ"] <- NULL
  germ_vaf <- vaf[["germ"]]
  vaf["germ"] <- NULL
  
  
  pdf(file=paste(ll,"pdf",sep = "."))
  
  #plot overall info
  for (v in 1:length(vaf)) {
    if(v==1){
      vaf_d <- data.frame(VAF=vaf[[v]],label=names(vaf)[v])
    }else{
      vaf_temp <- data.frame(VAF=vaf[[v]],label=names(vaf)[v])
      vaf_d<-rbind(vaf_d,vaf_temp)
    }
  }
  
  print(ggplot(vaf_d, aes(x=VAF, color=label)) +
          geom_density()+ggtitle(paste("VAF density plot for patient",ll)))
  
  print(ggplot(vaf_d, aes(x=VAF,color=label,fill=label)) +
          geom_histogram(binwidth = 0.02, position="dodge",alpha=0.5)+ggtitle(paste("VAF hist plot for patient",ll)))
  
  for (c in 1:length(cov)) {
    if(c==1){
      cov_d <- data.frame(COV=cov[[c]],label=names(cov)[c])
    }else{
      cov_temp <- data.frame(COV=cov[[c]],label=names(cov)[c])
      cov_d<-rbind(cov_d,cov_temp)
    }
  }
  
  print(ggplot(cov_d, aes(x=COV, color=label)) +
          geom_density()+ggtitle(paste("coverage density plot for patient",gsub("B.*","",name))))
  
  
  dev.off()
  
  germ_vaf <- vaf[["germ"]]
  vaf["germ"] <- NULL
  

  
  temp_allsnp <- data.frame(loci=gsub("\\#.*","",unlist(snps)),gene=gsub(".*\\.","",names(unlist(snps))),stringsAsFactors = F)
  temp_allsnp <- temp_allsnp[!duplicated(temp_allsnp$loci),]
  #remove suspecious locus filtered by dbSNP and donor affected frequency.
  temp_allsnp <- temp_allsnp[!temp_allsnp$loci %in% c("chr4 55141051","chr5 149433596","chr3 37067097","chr2 48032881","chr12 4383158","chr3 37067097","chr16 89831279","chr3 142266775","chr1 11187893",
                                                      "chr1 40363054",
                                                      "chr11 108183167",
                                                      "chr11 125525195",
                                                      "chr12 133233705",
                                                      "chr12 25368462",
                                                      "chr13 32913055",
                                                     "chr13 32915005",
                                                      "chr13 32929387",
                                                      "chr15 89838236",
                                                      "chr20 36030939",
                                                      "chr4 1807894",
                                                      "chr4 55141051",
                                                      "chr6 152201875",
                                                      "chr7 6036980",
                                                      "chr9 21968199",
                                                      "chr13 28636084",
                                                      "chr7 6026775",
                                                      "chr13 49033747",
                                                      "chr13 32936646"),]
  temp_allsnp <- temp_allsnp[!temp_allsnp$gene %in% c("POLE","ATM"),]

  
  
  
  allsnp <- as.character(temp_allsnp$loci)
  names(allsnp) <- temp_allsnp$gene

  
  
  l <- length(allsnp)-1
  final <- data.frame(id=paste("s",seq(0,l,1),sep = ""),
                      locus = allsnp,
                      name=names(allsnp),
                      var_reads="",
                      total_reads="",
                      var_read_prob=paste(rep("0.5",length(files)-1),collapse = ","),
                      vaf="",
                      pval="")
  final <- final[order(final$locus),]
  names <- c()
  no <- c()
  count <- 0
  
  for (i in 1:length(files)) {

    name <- gsub("\\..*","",files[i])
    
    #read in file
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"),stringsAsFactors=F)

    
    if(length(grep("full",name))>0){
      temp_file <- temp_file[,c(1,2,5,6,8,9,11,13)]
      temp_file$parse <- temp_file$locus
      temp_file[temp_file$dbsnp != "",]$dbsnp<-"YES"
      if(length(unique(temp_file$dbsnp == ""))==2){temp_file[temp_file$dbsnp == "",]$dbsnp<-"NO"}
      cat("sample: ")
      cat(paste(gsub("\\..*","",files[i]),"\n"))
      cat("\n")
      cat("\n")
      #temp_file <- temp_file[temp_file$coverage>100,]
      #temp_file <- temp_file[temp_file$VAF>vaf_cutoff,]
      if(length(grep(",",temp_file$VAF))>0){
        for (k in grep(",",temp_file$VAF)) {
          k_temp <- sum(as.numeric(gsub(".*\\=","",unlist(strsplit(temp_file$VAF[k],split = ",")))))
          temp_file$VAF[k] <- k_temp
        }
        
        temp_file$VAF <- as.numeric(temp_file$VAF)
      }
      
      
      temp_file$VAF <- as.numeric(temp_file$VAF)
      temp_file <- temp_file[temp_file$pvalue<0.05,]
      temp_file <- temp_file[temp_file$coverage>100,]
      temp_file <- temp_file[temp_file$VAF>7,]
      temp_file$VAF <- temp_file$VAF/100
      temp_germ <- temp_file
      temp_germ <- temp_germ[,c(1,2,3,5,6,7)]
      next
    }
    
    
    print(i)
    name <- gsub("\\..*","",files[i])
    names <- c(names,name)
    #read in file
    temp_file <- as.data.frame(fread(file=files[i],header = T,sep = "\t"))
    #temp_file <- temp_file[temp_file$VAF>vaf_cutoff,]
    temp_file$locus <- paste(temp_file$chr,temp_file$start)
    n <- length(temp_file$locus[temp_file$locus %in% allsnp])
    temp_file$`dbSNP151_common_20180423=1` <- gsub("\\_.*","",temp_file$`dbSNP151_common_20180423=1`)
    temp_file$`GNOMAD_EXOME_r2.1.1=1` <- gsub("\\_.*","",temp_file$`GNOMAD_EXOME_r2.1.1=1`)
    
    #check the somatic in db
    extras <- temp_file[temp_file$`dbSNP151_common_20180423=1`=="NO" & temp_file$`GNOMAD_EXOME_r2.1.1=1`=="NO",]
    extras <- extras[extras$pvalue<0.05,]
    extras <- extras[extras$coverage>100,]
    extras <- extras[,c("locus","gene","type","coverage","VAF","ISCN" )]
    extras$VAF <- extras$VAF/100
 

    
    if(length(temp_file$locus[temp_file$locus %in% allsnp]) != 0){
      
      parse_file <- temp_file[temp_file$locus %in% allsnp,]
      parse_file$altread <- round(parse_file$VAF*parse_file$coverage*0.01,0)
      parse_file <- parse_file[,c("locus","gene","coverage","altread","VAF","ISCN","pvalue" )]
      parse_file$VAF <- parse_file$VAF/100

    }else{
      no <- c(no,i)
      next
    }
    
    
    
    # assign(name,parse_file)
    if(length(unique(allsnp %in% temp_file$locus)) == 2){
      temp <- data.frame(locus=allsnp[!allsnp %in% parse_file$locus],
                         gene= names(allsnp[!allsnp %in% parse_file$locus]),
                         coverage=0,
                         altread=1,
                         VAF=0,
                         ISCN=2,
                         pvalue=0.04

      )
      all <- rbind(parse_file,temp)
      all <- all[order(all$locus),]
      #assign(name,all)
    }else{
      all <- parse_file
      all <- all[order(all$locus),]
    }

    all <- all[!duplicated(all$locus),]
    
    count <- count +1
    
    if(count==1){
      #print(final$var_reads)
      final$var_reads <- all$altread
      final$total_reads <- all$coverage
      final$vaf <- all$VAF
      final$pval <- all$pvalue
    }else{
      
      
      final$var_reads <- paste(final$var_reads,all$altread,sep = ",")
      final$total_reads <- paste(final$total_reads,all$coverage,sep = ",")
      final$vaf <- paste(final$vaf,all$VAF,sep = ",")
      final$pval <- paste(final$pval,all$pvalue,sep = ",")
    }

    
  }
  #formatting
  final$var_read_prob <- paste(rep("0.5",length(files)-length(no)-1),collapse = ",")
  cat('filter by normal:',nrow(final),"\n")
  
  #filtering by vaf
  vaf_del<-c()
  for (p in 1:nrow(final)) {
    if(length(grep("\\,",final[1,4]))>0){test_vaf <- unique(unlist(strsplit(final$vaf[p],split = ","))<0.07)
    }else{test_vaf <- final$vaf[p] <0.07}
      
    

    if(isTRUE(test_vaf)){
      vaf_del <- c(vaf_del,as.character(final$locus[p]))
      
    }
    
  }
  
  
  final <- final[!final$locus %in% vaf_del,]
  cat('filter by VAF:',nrow(final),"\n")
  
  
  
  #double check for quality assurance
  #filtering by coverage
  cov_del<-c()
  for (p in 1:nrow(final)) {
    if(length(grep("\\,",final[1,4]))>0){test_cov <- unique(as.numeric(unlist(strsplit(final$total_reads[p],split = ",")))<100)
    }else{test_cov <-final$total_reads[p]<100}
    if(isTRUE(test_cov)){
      cov_del <- c(cov_del,as.character(final$locus[p]))
      
    }
    
  }
  final <- final[!final$locus %in% cov_del,]
  cat('filter by coverage normal:',nrow(final),"\n")
  
  #filtering by p val
  if(length(grep("\\,",final[1,4]))>0){
  cov_del<-c()
  for (p in 1:nrow(final)) {
    test_cov <- unique(as.numeric(unlist(strsplit(final$pval[p],split = ",")))>0.05)
    if(length(test_cov) == 2){
      cov_del <- c(cov_del,as.character(final$locus[p]))
      
    }
    
  }
  final <- final[!final$locus %in% cov_del,]
  }else{final <- final[final$pval<0.05,]}
  cat('filter by p-value:',nrow(final),"\n")
  final <- final[,-c(7,8)]
  
  
  #add coverage 0
  if(length(grep("\\,",final[1,4]))>0){
  for (p in 1:nrow(final)) {
    covs <- as.numeric(unlist(strsplit(final$total_reads[p],split = ",")))
    if(length(unique(covs==0))==2){
      covs[covs==0] <- round(mean(covs[covs!=0]),0)
      final$total_reads[p] <-  paste(covs,collapse = ",") 
    }
  }
  }
  
  final$id <- paste("s",seq(0,nrow(final)-1,1),sep = "")
  write.table(final,file = paste(ll,".ssm",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
  sink(paste(ll,".name.json",sep = ""))
  if(is.null(no)==FALSE){
    new_files <- files[-no]
  }else{
    new_files <- files
  }
  
  new_files <- new_files[-grep("full",new_files)]
  
  cat('{"samples": [')
  for (i in 1:length(new_files)) {
    cat('"')
    cat(paste(gsub(".tsv","",new_files[i])))
    if(i==length(new_files)){
      cat('"')
      break
    }else{
      cat('",')
    }
  }
  #updating for mutation tree
  cat('], "clusters": [',paste("[",'"',final$id,'"',"]",sep = "",collapse = ","),'], "garbage": []}')
  sink() 
  
  
  
  
  
  
  
  
  
  
  
}







