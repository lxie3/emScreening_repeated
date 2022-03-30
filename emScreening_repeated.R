#joint method ttscreening on binary repeated outcomes for IOW data
#after performing step 1 and step 2 of ttscreening, return the ttscreening passing rate for identified mediators
#this code is intended for DNAm.Res adjusted for 6 cell types and SEX_18(gender) and mother maternal asthma ASTHMY1_0 and Gestational Age (GESTAGE)
#outcome variables are ASCON 1 18months 2 4 10 18 (6 times of observations)
library("MASS")
library(dplyr)
library(lme4)
library(lmerTest)
library(geepack)
library(reshape2)
library(data.table)
library("haven")
#setwd("C:/Users/lxie3/OneDrive - The University of Memphis/Data_Mediation/OneDrive_1_10-13-2020")
#load("DNAm Guthrie 12112019.RData")
load("/home/lxie3/DNAm Guthrie 12112019.RData")

###Functions for the joint method binary repeated outcomes
geeglm.indirect.jointest.cat=function(datanew){
  alpha.p<-summary(lm(datanew[,"M"]~datanew[,"X"]))$coef[2,c("Pr(>|t|)")]
  datanew2=as.data.frame(datanew)
  long <- melt(setDT(datanew2), id.vars = c("ID","X","M"),
               measure.vars=c("Y1","Y2","Y3","Y4", "Y5","Y6"),
               variable.name = "Y",variable.factor=TRUE)
  beta.p =coef(summary(geeglm(value~ M + X +Y,id=ID,
                              family = "binomial",corstr = "AR-1",
                              data =  long)))[2, c("Pr(>|W|)")]
  pvalue=max(alpha.p,beta.p)
  return(pvalue)
}


##Read in phenotype data
#PhenoType=as.matrix(read_sas("C:\\Users\\lxie3\\OneDrive\\mediation\\Data\\OneDrive_1_3-17-2020\\iow_v36a.sas7bdat"))
PhenoType=as.matrix(read_sas("/home/lxie3/iow_v36a.sas7bdat"))
dim(PhenoType)
#[1] 1536 2243
Pool=paste("X_",as.numeric(PhenoType[,"STUDYid"]),sep="")
rownames(PhenoType)=Pool

#screening example:  exposure variable MSMK_03 and outcome variable 1 2 4 10 18 asthma status(5 repeated outcomes)


#prepare the exposure and outcome variables 
PhenoType1=PhenoType[,c("MSMK_03","ASCON_1", "AST_2","ASCON_2", "ASCON_4","M_INVESTIGATORDIAGNOSEDASTHMA_10", "M_INVESTIGATORDIAGNOSEDASTHMA_18",   "SEX_18","ASTHMY1_0","GESTAGE")] #dim 1536 6
PhenoType1_complete <- PhenoType1[complete.cases(PhenoType1), ]
dim(PhenoType1_complete)
#[1] 924 10


##Read in the DNAm data
dim(DNAm)
#[1] 551710    797
DNAm_ID=colnames(DNAm)
sample_ID=intersect(DNAm_ID,rownames(PhenoType1_complete))
length(sample_ID)
#591
rownames(DNAm)=DNAm[,1]
DNAm_sample=DNAm[,sample_ID]
#dim 551710 591
PhenoType.data=PhenoType1[sample_ID,c("MSMK_03","ASCON_1", "AST_2","ASCON_2", "ASCON_4","M_INVESTIGATORDIAGNOSEDASTHMA_10","M_INVESTIGATORDIAGNOSEDASTHMA_18")]
#dim 591 7
table(rownames(PhenoType.data)==colnames(DNAm_sample))
#TRUE 
# 591
PhenoType.adjust=PhenoType1[sample_ID,c("SEX_18","ASTHMY1_0","GESTAGE")]
PhenoType.adjust=as.matrix(PhenoType.adjust)
colnames(PhenoType.adjust)=c("SEX_18","ASTHMY1_0","GESTAGE")
#dim 591 3
#want to use DNAm_sample residual value after regressing cell proportion data and gender and mother asthma

#read in cell proportion data
#celltype=read.csv(file="C:\\Users\\lxie3\\OneDrive\\mediation\\Data\\OneDrive_1_3-17-2020\\F1_guthrie850k_cell_24march2018.csv")
celltype=read.csv(file="/home/lxie3/F1_guthrie850k_cell_24march2018.csv")
dim(celltype)
#[1] 796   8
head(celltype)
rownames(celltype)=paste("X_",as.numeric(celltype[,"studyid"]),sep="")
celltype_sample=celltype[sample_ID,]
dim(celltype_sample)
#[1] 591   8

#need to combine 591*8 cell types and 591*2 gender and asthma into a data frame to compute residual
###Obtain the residuals for the DNAm data
##need HPC
table(rownames(celltype_sample)==colnames(DNAm_sample))
celltype_sample1=as.data.frame(cbind(as.matrix(celltype_sample),PhenoType.adjust))


#TRUE 
#591
library("limma")
res.celltype=function(edata,x){
  log2edata=log2(edata/(1-edata))
x1=as.numeric(x[,1])
x2=as.numeric(x[,2])            
x3=as.numeric(x[,3])         
x4=as.numeric(x[,4])          
x5=as.numeric(x[,5])         
x6=as.numeric(x[,6])
x7=as.numeric(x[,7])
x8=as.numeric(x[,8])
x9=as.numeric(x[,9])
mod1<- model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9)
fit1<-lmFit(log2edata,mod1,method="ls")
fite1<-eBayes(fit1)
res<-residuals(fite1,log2edata)
return(res)
}
DNAm.res=res.celltype(DNAm_sample,celltype_sample1[, c("CD4T","NK","Bcell","Mono","Neu","Eos","SEX_18","ASTHMY1_0","GESTAGE")])
dim(DNAm.res)
#[1] 551710  591
DNAm.res[1:5,1:5]
################################################################################
table(colnames(DNAm.res)==rownames(PhenoType.data))

DNAm.res<-as.data.frame(DNAm.res)


##Start the joint screening method, adjusted by ttscreening 
##two datasets, PhenoType.data, DNAm.res


s=100
p=dim(DNAm.res)[1]#number of potential mediators
n=dim(DNAm.res)[2]#sample size
class(PhenoType.data)="numeric"
joint.out.all=rep(NA, p)


### step 1 :screening each cg
for (i in 1:p){
  ID=1:n
  tempdata <- cbind(ID,X=PhenoType.data[,"MSMK_03"],Y=PhenoType.data[,c("ASCON_1","AST_2","ASCON_2" ,"ASCON_4","M_INVESTIGATORDIAGNOSEDASTHMA_10","M_INVESTIGATORDIAGNOSEDASTHMA_18")], M=t(DNAm.res[i,]))
  colnames(tempdata)=c("ID","X","Y1","Y2","Y3","Y4","Y5","Y6","M")
  joint.out.all[i]=geeglm.indirect.jointest.cat(tempdata)
}


#step 2: use 0.1 as cutoff 
Final_joint=matrix(NA,nrow=p,ncol=2)
colnames(Final_joint)<-c("pass","tt_rate")
rownames(Final_joint)<-rownames(DNAm.res)

for (i in 1:p){
  if (joint.out.all[i] >=0.1) {
    Final_joint[i]=0
  }else{
    ID=1:n
    temp<- cbind(ID,X=PhenoType.data[,"MSMK_03"],Y=PhenoType.data[,c("ASCON_1","AST_2","ASCON_2" ,"ASCON_4","M_INVESTIGATORDIAGNOSEDASTHMA_10","M_INVESTIGATORDIAGNOSEDASTHMA_18")], M=t(DNAm.res[i,]))
    colnames(temp)=c("ID","X","Y1","Y2","Y3","Y4","Y5","Y6","M")
  
    sig_path=rep(NA,s)
    for (k in 1:s) {
      set.seed(k)
      
      sample <- sample.int(n = n, size = floor(.67*n), replace = F)## Every loop generate IDs for train and test datasets
      
      ### we need to try each 
      train <- temp[sample, ]
      test  <- temp[-sample, ]
      Path.out.train=geeglm.indirect.jointest.cat(train)
      Path.out.test=geeglm.indirect.jointest.cat(test)
      
      sig_path[k] <- ifelse (Path.out.train<0.05 &Path.out.test<0.05,1,0)
    }
    
    #Final_joint[i,1]=ifelse (sum(sig_path)>(s/2),1,0)
    #Final_joint[i,1]=ifelse (sum(sig_path)>qs,1,0) # this q is a percentage of passing rate,default value of q = 0.5, 0<q<=1,
    #please write q and s both parameters in user package
    Final_joint[i,1]=ifelse (sum(sig_path)>1,1,0) 
    Final_joint[i,2]=sum(sig_path)/s
  }
  
}


#Final_joint[,1]==1 will give identified mediators from the joint screening method
#will run lavaan package for each of identified mediator

selected<-subset(Final_joint,Final_joint[,1]==1)
No.selected<-dim(selected)[1] #number of identified mediators after ttscreening

selected_ID=rownames(selected)

DNAm_selected<-DNAm.res[selected_ID,] 

selected
 














