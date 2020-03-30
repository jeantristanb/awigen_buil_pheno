library('knitr')
library('kableExtra')
library(VennDiagram)
source("/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/pheno/script/fct.r")
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
Var="meancimtleft"
DirOut="/spaces/jeantristan/GWAS/Romuald/GWAS/buildpheno/"
}else{
Var<-args[1]
DirOut<-args[2]
}
dir.create(DirOut)
## PCA file
### b
## 
FileI="/spaces/jeantristan/GWAS/Romuald/GWAS/buildpheno/Data_cIMT_GWAS_03APR2019.csv"
DirOutF=paste(DirOut,"/report_pdf",sep="")
if(!dir.exists(DirOutF))dir.create(DirOutF)
fileouput=paste(DirOut,"/DataRes_cIMT_GWAS_03APR2019.csv",sep="")

DataI<-read.csv(FileI)
DataI$site_qc<-c("Agincourt", "Dikgale", "Nairobi","Nanoro","Navrongo", "SOWETO")[DataI$site]
DataI$sex_qc<-DataI$sex
DataI$sex<-NA
DataI$sex[DataI$sex_qc=="Women"]<-0
DataI$sex[DataI$sex_qc=="Men"]<-1

names(DataI)[1:2]<-c("FID","IID")
listepca=list()
listenpca=list()

DirPca="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/Pca/prune_50_10_0.1/"
ListeGroupePca=c("All","Females", "Males", "South", "West", "East")

for(HeadPca in ListeGroupePca) {
listepca[[HeadPca]]=paste(DirPca,HeadPca,".eigenvec",sep="")
#listenpca=list(All=8, East=2, West=4, South=3, Women=5, Men=5)
if(length(grep("All", HeadPca))>0)listenpca[[HeadPca]]=8
else if(length(grep("East", HeadPca))>0)listenpca[[HeadPca]]=2
else if(length(grep("West", HeadPca))>0)listenpca[[HeadPca]]=4
else if(length(grep("South", HeadPca))>0)listenpca[[HeadPca]]=3
else if(length(grep("Women", HeadPca))>0)listenpca[[HeadPca]]=5
else if(length(grep("Females", HeadPca))>0)listenpca[[HeadPca]]=5
else if(length(grep("Men", HeadPca))>0)listenpca[[HeadPca]]=5
else if(length(grep("Males", HeadPca))>0)listenpca[[HeadPca]]=5
else {
cat("not found ", HeadPca, "for pca number")
q()
}
if(!file.exists(listepca[[HeadPca]])){
cat("Don't found ", listepca[[HeadPca]])
q()
}
}

DataLim<-read.table("corrvaluelimckd.in", header=T)
infovariable=""

#listenpca=list(All=6, East=4, West=4, South=4)
FamAll="/dataE/AWIGenGWAS/plink/awigen/qc/awigen-qc.fam"
DataPcaAll<-read.table(FamAll)

#Var<-"friedewald_qc"
DeletedSweet=F
VarSite="site_qc"
VarSex<-"sex"
VarId="FID"
VarSweetCtr<-"cohort_id_c";VarSweet<-"SWEET"
Soweto<-"SOWETO"
listM<-list()
ValSexF=0
M111<-NA
M222<-NA
NbSdLim<-NA
# to delete
listcovF<-c("age","sex", "bmi_qc")
Listlistcov<-list(agesexpca=c("age","sex"))
if(Var %in% c("waist_hip_r_c_qc", "visceral_fat_qc")){
Listlistcov[["agesexbmipca"]]=c("age","sex", "bmi_qc")
}

#fcntr<-RankNormTr
#namtr<-"Rank-based inverse normal transformations"
#namtrshort<-"rin"

fcntr<-function(x)x
namtr<-"No transformation"
namtrshort<-"notr"

allcovusedused<-c(listcovF, "hiv_final_status_c")
#Data<-DataI[DataI[,VarId] %in% DataPcaAll[,1] ,]
Data<-DataI
#https://www.ncbi.nlm.nih.gov/pubmed/22898763


if(Var %in% DataLim[,1]){
M111<-DataLim[DataLim[,1]==Var,2]
M222<-DataLim[DataLim[,1]==Var,3]
DeletedSweet=DataLim[DataLim[,1]==Var,4]
if(DataLim[DataLim[,1]==Var,4]==F)listcovF<-c("age","sex")
NbSdLim=DataLim[DataLim[,1]==Var,6]
}
Females<-0
## transform hiv
Data$hiv_final_status_c[Data$hiv_final_status_c==2 & Data$site_qc %in% c("Nanoro","Navrongo")]<-0
Data$age2<-Data$age**2
Data$hiv_final_status_c[Data$hiv_final_status_c==2]<-NA
Data$diabetes_status_c_qc[Data$diabetes_status_c_qc<0]<-NA
Data$htn_jnc7_qc[Data$htn_jnc7_qc<0]<-NA



fileouputlog<-paste(fileouput,".log",sep="")
if(!file.exists(fileouput)){
write.table(Data[,c(names(Data)[c(1,2)], unique(c(VarSite,VarSex, unlist(Listlistcov))))], file=fileouput, quote=F, sep="\t",row.names=F, col.names=T)
writeLines("Var\tVarI\tCov\tNbSd\tRegion\tTrans\tNpca",con=fileouputlog)
}
#FileRnw=paste(DirOut,"/report/",Var, ".Rnw",sep="")
FileRnw="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/pheno/script/TemplateAnalyse.Rnw"
knit(FileRnw, output=paste(DirOutF,"/",Var,".tex",sep=""))
