library(scales)
GetSig<-function(lx){
sapply(lx, function(x){
if(is.na(x))return("")
if(x>0.1)return("")
if(x>0.05)return(".")
if(x>0.01)return("*")
if(x>0.001)return("**")
return("***")
})
}

as.characterspe<-function(x,round){
y<-rep(NA, length(x))
bal=!is.na(x) & x>10**-round
y[bal]<-as.character(format(round(x[bal],round), scientific=F))
y[!bal]<-as.character(format(x[!bal], scientific=T,digits=round))
y
}

Dund<-function(x,rep=" ")gsub("_",rep,x)
plotBySite<-function(Data,Var,site, xlab="value", ylab="density", main=""){
        Data<-na.omit(Data[,c(Var,site)])
        ListSit<-unique(Data[,site])
        denslist<-lapply(1:length(ListSit), function(x)density(Data[Data[,site]==ListSit[x],Var]))
        xlim<-range(sapply(1:length(denslist),function(x)return(range(denslist[[x]]$x))))
        ylim<-range(sapply(1:length(denslist),function(x)return(range(denslist[[x]]$y))))
        plot(xlim, ylim,type="n", xlab=xlab,ylab=ylab, main=main)
        Cmt<-1
        for(x in denslist){
                        lines(x,col=Cmt)
                Cmt<-Cmt+1
        }
        legend("topright", legend=ListSit, lty=1, col=1:Cmt,bty='n')
}


plotBySiteSex<-function(Data,Var,site, sex, InfoSite, InfoSex,xlab="value", ylab="density", main=""){
        Data<-na.omit(Data[,c(Var,site,sex)])
        ListSit<-unique(Data[,c(site,sex)])
        ListSit<-merge(merge(ListSit,InfoSite),InfoSex)
        denslist<-lapply(1:nrow(ListSit), function(x)density(Data[Data[,site]==ListSit[x,2] & Data[,sex]==ListSit[x,1],Var]))
        xlim<-range(sapply(1:length(denslist),function(x)return(range(denslist[[x]]$x))))
        ylim<-range(sapply(1:length(denslist),function(x)return(range(denslist[[x]]$y))))
        plot(xlim, ylim,type="n",xlab=xlab,ylab=ylab, main=main)
        Cmt<-1
        for(x in denslist){
                        lines(x,col=ListSit[Cmt,"ColSite"], lty=ListSit[Cmt,"ColSex"])
                Cmt<-Cmt+1
        }
      InfoSexSpe<-unique(ListSit[,c(sex,"ColSex")])
      InfoSiteSpe<-unique(ListSit[,c(site,"ColSite")])
      legend("topleft", legend=c(as.character(InfoSexSpe[,VarSex]), as.character(InfoSiteSpe[,VarSite])), col=c(rep(1,length(InfoSexSpe[,VarSex])), InfoSiteSpe[,"ColSite"]), lty=c(InfoSexSpe[,"ColSex"], rep(1,length(InfoSiteSpe[,"ColSite"]))),bty="n")
}


resumestat<-function(x){
nmiss<-length(x[is.na(x)])
x<-x[!is.na(x)]
res<-c(nmiss,length(x),mean(x), median(x), sd(x),min(x),max(x))
names(res)<-c("nmiss","n","mean","med","sd","min","max")
return(res)
}
GetDistMiss<-function(x){
res<-c(GetPerc(x==-555), GetPerc(x==-999),GetPerc(x==-111) ,GetPerc(x==-222),GetPerc(x>-111 & x<0), GetPerc(x==0), GetPerc(is.na(x)), GetPerc(x>0))
names(res)<-c("-555","-999","-111","-222","<0","==0","NA",">0")
return(res)
}
GetPerc<-function(x)length(which(x))/length(x)*100

PlotPca<-function(Data, Site, listcovF, listepca, listenpca){
DataPca<-read.table(listepca[[Site]])
npca<-listenpca[[Site]]
names(DataPca)<-c("FID","IID",paste("Axis_",1:(ncol(DataPca)-2),sep=""))

DataAndPca<-merge(Data,DataPca,by.x=VarId,by.y="FID")
DataAndPca[,VarSite]<-as.character(DataAndPca[,VarSite])
rownames(DataAndPca)<-DataAndPca[,VarId]
resglm<-glm(as.formula(paste(Var,"~",paste(c(listcovF, paste("Axis_",1:npca,sep="")),collapse="+"),sep="")), data=DataAndPca)
sumresglm<-summary(resglm)
print(kable(sumresglm$coefficients, row.names=T, booktabs = T, caption=paste("summary for glm with pca"),format = "latex")  %>% kable_styling(latex_options = c("striped", "HOLD_position")))
DataAndPca<-merge(DataAndPca,data.frame(Pos=names(resglm$residuals),res_agesexbmipca=resglm$residuals-min(resglm$residuals)+1),by.x=VarId, by.y=1)
boxplot(as.formula(paste("res_agesexbmipca~",VarSex,"+",VarSite)), data=DataAndPca, las=3)

#sizepch=(DataAndPca[,Var]-min(DataAndPca[,Var],na.rm=T))/(max(DataAndPca[,Var])-min(DataAndPca[,Var]))*3+1
DataAndPca<-DataAndPca[!is.na(DataAndPca[,Var]),]
sizepch=(DataAndPca[,Var]-mean(DataAndPca[,Var], na.rm=T))/sd(DataAndPca[,Var],na.rm=T)
sizepch=sizepch-min(sizepch)
sizepch<-order(DataAndPca[,Var])/nrow(DataAndPca)*2.5
Col<-unique(DataAndPca[,c("site","site_qc")])
plot(DataAndPca$Axis_1, DataAndPca$Axis_2,col=alpha(DataAndPca$site,0.4),pch=DataAndPca$sex+15,cex=sizepch,xlab="axis1", ylab="axis2")
legend("topright",legend=Col[,"site_qc"], col=Col[,"site"], pch=1, bty='n')
plot(DataAndPca$Axis_3, DataAndPca$Axis_4,col=alpha(DataAndPca$site,0.4),pch=DataAndPca$sex+15,cex=sizepch,xlab="axis4", ylab="axis4")
legend("topright",legend=Col[,"site_qc"], col=Col[,"site"], pch=1, bty='n')

return(DataAndPca)
}

PlotPca2<-function(Data, Region, listepca, listenpca, Var,VarSite, VarSex,VarId,InfoSite, InfoSex){
DataPca<-read.table(listepca[[Region]])
DataEig<-read.table(gsub("eigenvec","eigenval",listepca[[Region]]))
npca<-listenpca[[Region]]
names(DataPca)<-c("FID","IID",paste("Axis_",1:(ncol(DataPca)-2),sep=""))
DataAndPca<-merge(Data,DataPca,by.x=VarId,by.y="FID")
DataAndPca[,VarSite]<-as.character(DataAndPca[,VarSite])
rownames(DataAndPca)<-DataAndPca[,VarId]
DataAndPca<-DataAndPca[!is.na(DataAndPca[,Var]),]
sizepch=(DataAndPca[,Var]-mean(DataAndPca[,Var], na.rm=T))/sd(DataAndPca[,Var],na.rm=T)
sizepch=sizepch-min(sizepch)
sizepch<-order(DataAndPca[,Var])/nrow(DataAndPca)*2.5
DataAndPca<-merge(merge(DataAndPca,InfoSite,by=VarSite,all.x=T), InfoSex, by=VarSex,all.x=T)
par(mfrow=c(2,2))
plot(DataEig[,1], xlab="eigen num", ylab="variance explained",main="")
plot(DataAndPca$Axis_1, DataAndPca$Axis_2,col=alpha(DataAndPca$ColSite,0.2),pch=DataAndPca$ColSex+15,cex=sizepch,xlab="axis1", ylab="axis2")
legend("topright",legend=InfoSite[,VarSite], col=InfoSite[,"ColSite"], pch=1, bty='n')
plot(DataAndPca$Axis_3, DataAndPca$Axis_4,col=alpha(DataAndPca$ColSite,0.2),pch=DataAndPca$ColSex+15,cex=sizepch,xlab="axis3", ylab="axis4")
plot(DataAndPca$Axis_5, DataAndPca$Axis_6,col=alpha(DataAndPca$ColSite,0.2),pch=DataAndPca$ColSex+15,cex=sizepch,xlab="axis5", ylab="axis6")
return(DataAndPca)
}



RankNormTr<-function(a)qnorm((rank(a)-0.5)/(length(a)-2*0.5+1))

GetLinearisation<-function(Site, listWithPca){
valueI<-function(x)x
squareroot<-function(x)x**2
cuberoot<-function(x){
y<-x**3
y[x<0]<- -y[x<0]
y
}
xval<-listWithPca[[Site]][,"res_agesexbmipca"]
meanval<-mean(xval, na.rm=T)
valsddata<-sd(xval,na.rm=T)*3:6
plot(density(xval, na.rm=T), main="density residuals with pca")
abline(v=meanval+valsddata,col=2:(length(valsddata)+1))
abline(v=meanval-valsddata,col=2:(length(valsddata)+1))
#
listvalsd<-c()
for(sdval in valsddata){
listvalsd<-c(listvalsd, length(xval[xval>meanval+sdval|xval< meanval-sdval]))
}
listvalsd<-matrix(listvalsd, nrow=1)
colnames(listvalsd)<-paste(3:6,"sd")
print(kable(listvalsd, caption=c("values higher than x * sd"), row.names=F, booktabs = T, format = "latex")  %>% kable_styling(latex_options = c("striped", "HOLD_position")))

par(mfrow=c(2,2))
listefunction<-c("valueI", "log", "squareroot", "cuberoot")
Cmt<-1
for(trans in listefunction){
xvaltr<-get(trans)(xval)
plot(density(xvaltr, na.rm=T), main=trans)
kstrs<-ks.test(xvaltr, "pnorm", mean(xvaltr,na.rm=T), sd(xvaltr,na.rm=T))
ksres<-c(kstrs$statistic,kstrs$p.value)
if(Cmt==1)ksresall<-ksres
else ksresall<-rbind(ksresall,ksres)
Cmt<-Cmt+1
}
colnames(ksresall)<-c("D","p.value")
rownames(ksresall)<-listefunction
kable(cbind(Transformation=listefunction,as.data.frame(ksresall)), row.names=T, booktabs = T, caption=paste("ks test for various transformation on data"),format = "latex")  %>% kable_styling(latex_options = c("striped", "HOLD_position"))
}

