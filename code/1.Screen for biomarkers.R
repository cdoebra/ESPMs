#-------------------------------suvival analysis---------------------------------------------------
cancer=c("HNSC","LIHC","LUAD","LUSC")
set.seed(2023)
#install.packages("tidyverse")
#install.packages("survival")
#install.packages("dplyr")
#install.packages("stringr")
library(tidyverse)
library(survival)
library(dplyr)
library(stringr)
for (i in 1:length(cancer)){
  exp<-read.table(paste("./RNA_seq/",cancer[i],"_RNAseq_TPM_T_classmean.txt",sep = ""),sep = "\t",header = T,check.names = F,row.names = 1)
  exp<-as.data.frame(t(exp))
  t<-read.table(paste("./RNA_seq/",cancer[i],"_suvival_classmean.csv",sep = ""),sep = ",",header = T,check.names = F,row.names = 1)
  eSet<-exp
  cg <- array(dim = nrow(eSet))
  gt<-as.data.frame(array(,dim=c(nrow(eSet),5)))
  for (j in 1:nrow(eSet)) {
    gene_exprs <- eSet[j,]
    gene_exprs <- as.data.frame(t(gene_exprs))
    gene_exprs[,1] <- as.numeric(gene_exprs[,1])
    x<-gene_exprs[,1]
    ###quartile grouping
    gene_exprs$g <- 
      ifelse(x>quantile(x,.75),"high",ifelse(x>quantile(x,.5),"Q2",ifelse(x>quantile(x,.25),"Q3","low")))  
    gene_surv <- t
    gene_surv <- merge(gene_surv,gene_exprs,by.x="case_submitter_id",by.y = "row.names")
    colnames(gene_surv)[1]<-c("Patient")
    gene_surv<-gene_surv %>%filter(g %in% c("high","low")) 
    ###cox proportional-hazards model
    if ("high" %in% gene_surv$g && "low" %in% gene_surv$g){
      try({
        fit<-coxph(Surv(gene_surv$suvival_month,gene_surv$vital_status=="Dead")~gene_surv$g,data = gene_surv)
        gt[j,1:4]<-as.data.frame(summary(fit)$conf.int)[1,]###HR##
        gt[j,5]<-t(as.data.frame(summary(fit)$logtest))[1,3]###pvalue##
      },silent = T)}}
  rownames(gt)<-rownames(eSet)
  colnames(gt)<-c("exp(coef)","exp(-coef)","lower95%","upper95%","p-value")
  write.table(gt,paste("./Survival genes/TCGA_",cancer[i],"_5_suvivalgene_coxph_quartile.txt",sep = ""),sep='\t',quote = FALSE,row.names = TRUE)
}

#Screening for significant survival-related genes
#------------------------------p0.05------------------------------------------------------------------
cancer1=c("LUSC")
for (i in 1:length(cancer1)){
  cox<-read.table(paste("./Survival genes/TCGA_",cancer1[i],"_5_suvivalgene_coxph_quartile.txt",sep = ""),sep="\t",check.names = F)
  cox_diff_genes<- cox%>%filter(cox$`p-value`<0.05)
  write.table(cox_diff_genes,paste("./Survival genes/TCGA_",cancer1[i],"_5_suvivalgene_p0.05_quartile.txt",sep = ""),sep='\t',quote = FALSE)
}
#------------------------------p0.025------------------------------------------------------------------
cancer2=c("HNSC","LUAD")
for (i in 1:length(cancer2)){
  cox<-read.table(paste("./Survival genes/TCGA_",cancer2[i],"_5_suvivalgene_coxph_quartile.txt",sep = ""),sep="\t",check.names = F)
  cox_diff_genes<- cox%>%filter(cox$`p-value`<0.025)
  write.table(cox_diff_genes,paste("./Survival genes/TCGA_",cancer2[i],"_5_suvivalgene_p0.05_quartile.txt",sep = ""),sep='\t',quote = FALSE)
}
#------------------------------p0.001------------------------------------------------------------------
cancer4=c("LIHC")
for (i in 1:length(cancer4)){
  cox<-read.table(paste("./Survival genes/TCGA_",cancer4[i],"_5_suvivalgene_coxph_quartile.txt",sep = ""),sep="\t",check.names = F)
  cox_diff_genes<- cox%>%filter(cox$`p-value`<0.001)
  write.table(cox_diff_genes,paste("./Survival genes/TCGA_",cancer4[i],"_5_suvivalgene_p0.05_quartile.txt",sep = ""),sep='\t',quote = FALSE)
  print(nrow(cox_diff_genes))
}


#------------------------Screening for potential oncology biomarkers with evolutionary features-------------------------------
#loading evolution information
oh<-read.table("./Evolution information/Ohnologs_all_gene_ensg_gene-name.txt",sep = "\t",header = T,check.names = F)
stage<-read.table("./Evolution information/main_HUMAN-gene-name.txt",sep = "\t",header = F,check.names = F)
stage<-stage[stage$V2 %in% c("Eukaryota","Opisthokonta","Eumetazoa"),]
evo_gene<-stage[stage$V3%in%oh$`Associated Gene Name`,]
gene_ratio<-as.data.frame(array(,dim=c(length(cancer),2)))
colnames(gene_ratio)<-c("Evo","Nonevo")
rownames(gene_ratio)<-cancer

for (i in 1:length(cancer)){
  cox_diff_genes<-read.table(paste("./Survival genes/TCGA_",cancer[i],"_5_suvivalgene_p0.05_quartile.txt",sep = ""),sep = '\t',check.names = F)
  filter<-evo_gene[evo_gene$V3%in% rownames(cox_diff_genes),]
  exp<-read.table(paste("./RNA_seq/",cancer[i],"_RNAseq_TPM_T_classmean.txt",sep = ""),sep = "\t",header = T,check.names = F,row.names = 1)
  exp<-exp[,colnames(exp)%in%rownames(cox_diff_genes)]
  t<-read.table(paste("./RNA_seq/",cancer[i],"_suvival_classmean.csv",sep = ""),sep = ",",header = T,check.names = F)
  row.names(t)<-t$case_submitter_id
  t<-t[rownames(t)%in%rownames(exp),]
  evo_exp<-exp[,colnames(exp)%in%filter$V3]
  evo_exp <- evo_exp %>%mutate(class=NA)
  for (j in 1:nrow(evo_exp)){
    evo_exp[j,"class"]<-t[rownames(evo_exp)[j],"class"]
  }
  write.table(evo_exp,paste("./RNA_seq/TCGA_",cancer[i],"_exp_classmean_p0.05_evo_quartile.txt",sep = ""),sep="\t",quote = FALSE)
  
  nonevo_exp<-exp[,!(colnames(exp)%in%filter$V3)]
  nonevo_exp <- nonevo_exp %>%mutate(class=NA)
  for (j in 1:nrow(nonevo_exp)){
    nonevo_exp[j,"class"]<-t[rownames(nonevo_exp)[j],"class"]
  }
  write.table(nonevo_exp,paste("./RNA_seq/TCGA_",cancer[i],"_exp_classmean_p0.05_nonevo_quartile.txt",sep = ""),sep="\t",quote = FALSE)
  gene_ratio[cancer[i],1]<-ncol(evo_exp)-1
  gene_ratio[cancer[i],2]<-ncol(nonevo_exp)-1
}