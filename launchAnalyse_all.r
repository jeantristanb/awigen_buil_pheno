#!/usr/bin/Rscript
library('knitr')
library("optparse")
library('kableExtra')
library(VennDiagram)

dirscript=dirname(commandArgs()[4])
dirscript=strsplit(dirscript,split="=")[[1]][2]
source(paste(dirscript,"script/fct.r",sep='/'))
FileRnw=paste(dirscript,"script/TemplateAnalyse.Rnw",sep='/')

option_list = list(
  make_option("--pheno_file", type="character", default=NULL, 
              help="initial file of phenotype", metavar="character"),
  make_option("--sep", type="character", default=NULL, 
              help="separator for pheno file", metavar="character"),
  make_option("--val_sex", type="character", default="1", 
              help="valeur for females in sex", metavar="character"),
  make_option("--head_id", type="character", default='FID,IID',
              help="head fid", metavar="character"),
  make_option("--head_sex", type="character", 
              help="head sex", metavar="character"),
  make_option("--head_site", type="character", 
              help="head sex", metavar="character"),
  make_option("--params", type="character", default=NULL, 
              help="separator for pheno file", metavar="character"),
  make_option("--tr_var", type="character", default=NULL, 
              help="file for tranformation", metavar="character"),
  make_option("--fam_file", type="character", default=NULL, 
              help="file id to keep", metavar="character"),
  make_option("--dir_out", type="character", default="./", 
              help="dir out file [default= %default]", metavar="character"),
  make_option("--output", type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);
listargcheck=c('dir_out', 'sep', 'pheno_file', 'head_id', 'params', 'fam_file')
for(arg in listargcheck){
if(is.null(args[[arg]])){
cat('not found ',arg, '\n exit\n')
q('no',2)
}
}
#if(length(args)==0){
#Var="meancimtleft"
#DirOut="/spaces/jeantristan/GWAS/Romuald/GWAS/buildpheno/"
#}elseargs{
#Var<-args[1]
#DirOut<-args[2]
#}
dir_out=args[['dir_out']]
file_phenoI=args[['pheno_file']]
sep=GetSep(args[['sep']])
if(is.na(sep)){
cat('\n','not found sep', args[['sep']],'\n')
q('no',2)
}
if(!is.null(args[['params']])){
data_param=read.table(args[['params']] ,header=T, stringsAsFactors=F)
}else{
Var=args[['var']]
q('no', 2)
}

dir.create(dir_out)
## PCA file
### b
## 
DirOutF=paste(dir_out,"/report_pdf",sep="")
if(!dir.exists(DirOutF))dir.create(DirOutF)

fileouput=paste(getwd(),'/',dir_out,"/", args[['output']],sep="")
DataI<-read.csv(file_phenoI, sep=sep)

## must be 2
fidiid=strsplit(args[['head_id']] ,split=',')[[1]]
if(length(fidiid)!=2 | any(!(fidiid %in% names(DataI)))){
cat('\n','problem --head_id  not found', args[['head_id']], '\n')
q('no',2)
}
names(DataI)[fidiid[1]==names(DataI)]<-'FID'
names(DataI)[fidiid[2]==names(DataI)]<-'IID'
if(!is.null(args[['tr_var']])){
if(!file.exists(args[['tr_var']])){
cat('\nFile tr_var:', args[['tr_var']], 'not found \nexit\n')
q('no',2)
}
cat('\ntransformation of variables\n')
tr_var<-read.table(args[['tr_var']],header=T, stringsAsFactors=F,sep=' ')
DataI<-transform_allvar(DataI, tr_var)
}
ListVarId<-c('FID','IID')

DataFam<-read.table(args[['fam_file']])
FamFile=args[['fam_file']]
if(is.null(args[['head_sex']])){
BaliseSex=BalSex=F
VarSex='SexNull'
DataI[,VarSex]='all'
}else{
BaliseSex=T
}
if(is.null(args[['head_site']])){
VarSite='SiteTmp'
DataI[,VarSite]='All'
}else{
VarSite=args[['head_site']]
}
if(!(VarSite %in% names(DataI))){
cat('\n not found header for site', VarSite, '\n exit\n')
q('no', 2)
}
BaliseCheckSoweto=F
ValSexF=args[['val_sex']]

ListlistcovAll<-unique(unlist(strsplit(unlist(strsplit(unlist(strsplit(data_param[!is.na(data_param[,'cov']),'cov'],split=",")),split=',')),split="[:*]")))
cat('begin change variable')
fileouputlog<-paste(fileouput,".log",sep="")
if(!file.exists(fileouput)){
write.table(DataI[,c(ListVarId, unique(c(VarSite,VarSex, unlist(ListlistcovAll))))], file=fileouput, quote=F, sep="\t",row.names=F, col.names=T)
writeLines("Var\tVarI\tCov\tNbSd\tRegion\tTrans\tNpca",con=fileouputlog)
}

for(CmtVar in 1:nrow(data_param)){
## pca
tmp<-getinfopca(data_param[CmtVar,])
listepca=tmp[['pca']]
listenpca=tmp[['npca']]
## 
#UpperLimit=
#fcntr<-function(x)x
#namtr<-"No transformation"
#namtrshort<-"notr"
if(is.na(data_param[CmtVar,'transf']))transf='notr'
else transf=data_param[CmtVar,'transf']
fcntr<-get(transf)
test<-fcntr(1)
namtrshort=transf
namtr<-GetInfoTrans(transf)
Var=as.character(data_param[CmtVar,'variable'])
Data<-DataI
limlowerfor<-data_param[CmtVar,"limlower"]
limlowerrep<-data_param[CmtVar,"limlower_repl"]
limupperfor<-data_param[CmtVar,"limupper"]
limupperrep<-data_param[CmtVar,"limupper_repl"]
NbSdLim=data_param[CmtVar,'sdlim']
Listlistcov<-getinfocov(data_param[CmtVar,])
print(Listlistcov)
print(listepca)
NaValue=data_param[CmtVar,'navalue']
if(!is.na(data_param[CmtVar,'formula']))Data[,Var]=formulastr(Data, data_param[CmtVar,'formula'])
knit(FileRnw, output=paste(DirOutF,"/",Var,".tex",sep=""))
}
#DataI$sex_qc<-DataI$sex
#DataI$sex<-NA
#DataI$sex[DataI$sex_qc=="Women"]<-0
#DataI$sex[DataI$sex_qc=="Men"]<-1



#listepca=list()
#listenpca=list()

#DirPca="/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/Pca/prune_50_10_0.1/"
#ListeGroupePca=c("All","Females", "Males", "South", "West", "East")
#for(HeadPca in ListeGroupePca) {
#listepca[[HeadPca]]=paste(DirPca,HeadPca,".eigenvec",sep="")
##listenpca=list(All=8, East=2, West=4, South=3, Women=5, Men=5)
#if(length(grep("All", HeadPca))>0)listenpca[[HeadPca]]=8
#else if(length(grep("East", HeadPca))>0)listenpca[[HeadPca]]=2
#else if(length(grep("West", HeadPca))>0)listenpca[[HeadPca]]=4
#else if(length(grep("South", HeadPca))>0)listenpca[[HeadPca]]=3
#else if(length(grep("Women", HeadPca))>0)listenpca[[HeadPca]]=5
#else if(length(grep("Females", HeadPca))>0)listenpca[[HeadPca]]=5
#else if(length(grep("Men", HeadPca))>0)listenpca[[HeadPca]]=5
#else if(length(grep("Males", HeadPca))>0)listenpca[[HeadPca]]=5
#else {
#cat("not found ", HeadPca, "for pca number")
#q()
#}
#if(!file.exists(listepca[[HeadPca]])){
#cat("Don't found ", listepca[[HeadPca]])
#q()
#}
#}

#DataLim<-read.table("corrvaluelimckd.in", header=T)
#infovariable=""

#listenpca=list(All=6, East=4, West=4, South=4)

#Var<-"friedewald_qc"
#DeletedSweet=F
#VarSite="site_qc"
#VarSex<-"sex"
#VarId="FID"
#VarSweetCtr<-"cohort_id_c";VarSweet<-"SWEET"
#Soweto<-"SOWETO"
#listM<-list()
#ValSexF=0
#M111<-NA
#M222<-NA
#NbSdLim<-NA
# to delete
#listcovF<-c("age","sex", "bmi_qc")
#Listlistcov<-list(agesexpca=c("age","sex"))
#if(Var %in% c("waist_hip_r_c_qc", "visceral_fat_qc")){
#Listlistcov[["agesexbmipca"]]=c("age","sex", "bmi_qc")
#}

#fcntr<-RankNormTr
#namtr<-"Rank-based inverse normal transformations"
#namtrshort<-"rin"

#fcntr<-function(x)x
#namtr<-"No transformation"
#namtrshort<-"notr"

#allcovusedused<-c(listcovF, "hiv_final_status_c")
#Data<-DataI[DataI[,VarId] %in% DataPcaAll[,1] ,]
#Data<-DataI
#https://www.ncbi.nlm.nih.gov/pubmed/22898763


#if(Var %in% DataLim[,1]){
#M111<-DataLim[DataLim[,1]==Var,2]
#M222<-DataLim[DataLim[,1]==Var,3]
#DeletedSweet=DataLim[DataLim[,1]==Var,4]
#if(DataLim[DataLim[,1]==Var,4]==F)listcovF<-c("age","sex")
#NbSdLim=DataLim[DataLim[,1]==Var,6]
#}
#Females<-0
### transform hiv
#Data$hiv_final_status_c[Data$hiv_final_status_c==2 & Data$site_qc %in% c("Nanoro","Navrongo")]<-0
#Data$age2<-Data$age**2
#Data$hiv_final_status_c[Data$hiv_final_status_c==2]<-NA
#Data$diabetes_status_c_qc[Data$diabetes_status_c_qc<0]<-NA
#Data$htn_jnc7_qc[Data$htn_jnc7_qc<0]<-NA
#
#
#
#if(!file.exists(fileouput)){
#write.table(Data[,c(names(Data)[c(1,2)], unique(c(VarSite,VarSex, unlist(Listlistcov))))], file=fileouput, quote=F, sep="\t",row.names=F, col.names=T)
#writeLines("Var\tVarI\tCov\tNbSd\tRegion\tTrans\tNpca",con=fileouputlog)
#}
