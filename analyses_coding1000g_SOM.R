#======================LOAD_INITIAL_FILES====================================================
{
#TAG 'LOAD_INITIAL_FILES'
setwd("/mnt/scratch/fabrizio/LDLD")
data<-read.table("above95/coding1000g/anal/sorted.res")
pairid<-sapply(1:dim(data)[1],function(i) paste(mydata[i,1],mydata[i,2],mydata[i,3],mydata[i,4],sep="."))
data<-cbind(data,pairid)
names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX","pairid")
#save(data,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/rawres.RData")
load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/rawres.RData")
dataFIN<-subset(data,data$pop==7)
ncomp<-read.table("above95/coding1000g/ncomparisons.txt",header=FALSE)
ncomp<-sum(ncomp$V1)
combined_fdr<-p.adjust(p=dchisq(dataFIN$X,df=2*dataFIN$popX), method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
combined_pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)
length(combined_fdr)
dim(dataFIN)
dim(data)
fullsign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
save(fullsign,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/fullsign.RData")
save(combined_pval,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/combined_pval.RData")
#in early code called sign
dim(fullsign)
highlysign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.00001,]
dim(highlysign)

signAbed<-as.data.frame(cbind(as.character(paste0('chr',fullsign$chrA)),fullsign$posA-1,fullsign$posA),stringsAsFactors =FALSE)
signBbed<-as.data.frame(cbind(as.character(paste0('chr',fullsign$chrB)),fullsign$posB-1,fullsign$posB),stringsAsFactors =FALSE)
names(signAbed)<-c("chr","start","end")
names(signBbed)<-c("chr","start","end")
signAbed$start<-as.numeric(signAbed$start)
signAbed$end<-as.numeric(signAbed$end)
signBbed$start<-as.numeric(signBbed$start)
signBbed$end<-as.numeric(signBbed$end)
#create logos
write.table(signAbed,file="above95/coding1000g/anal/snpsA.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file="above95/coding1000g/anal/snpsB.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);

#explorative plots to check consistency of p-values
#custom function to make color transparent
t_col <- function(color, percent = 50, name = NULL) { #color = color name, percent = % transparency,name = an optional name for the color
      rgb.val <- col2rgb(color) # Get RGB values for named color
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], # Make new color using input color as base and alpha set by transparency
               max = 255, alpha = (100-percent)*255/100,names = name)
  # Save the color
  invisible(t.col)
  }

#=============================SIGNIFICANT PAIRS============================================  
#Different ways of obtaining list of significant pairs
#no empirical distribution
combined_T2pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)
combined_T2fdr<-p.adjust(p=combined_T2pval, method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
mydata<-data[data$nA>=0.05*data$Nse,]
mydata<-mydata[mydata$nA<=0.95*mydata$Nse,]
mydata<-mydata[mydata$nB>=0.05*mydata$Nse,]
mydata<-mydata[mydata$nB<=0.95*mydata$Nse,]
range(combined_T2pval) #0.000000e+00 2.496462e-06
head(mydata)

combined_prho2pval<-aggregate(mydata$prho2,by=list(mydata$pairid),FUN=function(x) 1-pchisq(-2*log(prod(x)),df=2*length(x)))
combined_pKulipval<-aggregate(mydata$pKuli,by=list(mydata$pairid),FUN=function(x) 1-pchisq(-2*log(prod(x)),df=2*length(x)))
combined_prho2fdr<-p.adjust(p=combined_prho2pval[,2], method = "fdr", n = ncomp) 
combined_pKulifdr<-p.adjust(p=combined_pKulipval[,2], method = "fdr", n = ncomp) 

mypvalues<-cbind(combined_prho2pval,combined_pKulipval[,2],combined_prho2fdr,combined_pKulifdr)

save(fullsign,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/fullsign.RData")
names(mypvalues)<-c("pairid","combined_prho2","combined_pKuli","combined_prho2_fdr","combined_pKuli_fdr")
mydata<-merge(mydata,mypvalues,by="pairid")
save(mydata,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/mypairs.RData")

  
#plot with pvalues for all populations from snps that are significant in combined p-value based on T2 (without permutations).===========V
#these plots show that various pvalues method correlate pretty well, except for T2 that overestimated, but used for quick calculations.
mydata<-fullsign[fullsign$nA>=0.05*fullsign$Nse,]
mydata<-mydata[mydata$nA<=0.95*mydata$Nse,]
mydata<-mydata[mydata$nB>=0.05*mydata$Nse,]
mydata<-mydata[mydata$nB<=0.95*mydata$Nse,]
plot.new()
mythr<-0.00001
mythr<-0.1
library(ggplot2)
dev.new()
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/LDLDpvalues_distr_signfdrfromT2.pdf")
mydata2<-data.frame(values=c(mydata$pKuli,mydata$prho2,mydata$T2),groups=c(rep("pKuli",length(mydata$pKuli)),rep("prho2",length(mydata$prho2)),rep("T2",length(mydata$prho2))))
ggplot(mydata2, aes(x=values, fill=groups)) +
    geom_histogram(binwidth=.05, position="dodge")+theme_bw()
dev.off()
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/LDLDpvalues_distr_zoom_signfdrfromT2.pdf")
mythr<-0.0001
mydata2<-data.frame(values=c(mydata$pKuli[mydata$pKuli<mythr],mydata$prho2[mydata$prho2<mythr],mydata$T2[mydata$T2<mythr]),groups=c(rep("pKuli",length(mydata$pKuli[mydata$pKuli<mythr])),rep("prho2",length(mydata$prho2[mydata$prho2<mythr])),rep("T2",length(mydata$T2[mydata$T2<mythr]))))
ggplot(mydata2, aes(x=values, fill=groups)) +
    geom_histogram(binwidth=mythr/30, position="dodge")+theme_bw()
dev.off()
#scatterplot of pvalues
dev.new()
library(cowplot)
mysize<-3
myplot<-ggplot(mydata,aes(x=mydata$T2,y=mydata$pKuli))+stat_binhex()+scale_fill_gradient(limits=c(1,10^8), low="white", high="cadetblue4",trans="log")+theme_bw()+theme(legend.position=c(0.75,0.25),legend.title=element_text(size=5),legend.text=element_text(size=5),legend.key.size=unit(0.2,"cm"))
mytest<-cor.test(mydata$T2,mydata$pKuli)
myplot1<-myplot+labs(x="T2",y="pKuli")+annotate("text",x=0.7,y=0.06,label=paste0("p-value=",mytest$p.value),cex=mysize,col="black")+annotate("text",x=0.7,y=0.01,label=paste0("rho=",round(mytest$estimate,digits=4)),cex=mysize,col="black")
myplot<-ggplot(mydata,aes(x=mydata$T2,y=mydata$prho2))+stat_binhex()+scale_fill_gradient(limits=c(1,10^8), low="white", high="cadetblue4",trans="log")+theme_bw()+theme(legend.position="none")
mytest<-cor.test(mydata$T2,mydata$rho2)
myplot2<-myplot+labs(x="T2",y="prho2")+annotate("text",x=0.7,y=0.06,label=paste0("p-value=",mytest$p.value),cex=mysize,col="black")+annotate("text",x=0.7,y=0.01,label=paste0("rho=",round(mytest$estimate,digits=4)),cex=mysize,col="black")
myplot<-ggplot(mydata,aes(x=mydata$T2,y=mydata$D1))+stat_binhex()+ylim(-0.0005,0.0005)+scale_fill_gradient(limits=c(1,10^8), low="white", high="cadetblue4",trans="log")+theme_bw()+theme(legend.position="none")
mytest<-cor.test(mydata$T2[abs(mydata$D1)<1],abs(mydata$D1)[abs(mydata$D1)<1],na.rm=TRUE)
myplot3<-myplot+labs(x="T2",y="D1")+annotate("text",x=0.7,y=-0.00035,label=paste0("p-value=",mytest$p.value),cex=mysize,col="black")+annotate("text",x=0.7,y=-0.00042,label=paste0("rho=",round(mytest$estimate,digits=4)),cex=mysize,col="black")
toprow<-plot_grid(myplot1,myplot2,myplot3,align='h',labels=c('A','B','C'),ncol=3)
myplot<-ggplot(mydata,aes(x=mydata$rho2,y=log(mydata$T2+10^(-16))))+stat_binhex()+scale_fill_gradient(limits=c(1,10^8), low="white", high="cadetblue4",trans="log")+theme_bw()+theme(legend.position="none")
mytest<-cor.test(mydata$rho2,mydata$T2)
myplot1<-myplot+labs(x="rho2",y="log(T2+10^-6)")+annotate("text",x=0.7,y=-1,label=paste0("p-value=",mytest$p.value),cex=mysize,col="black")+annotate("text",x=0.7,y=-3,label=paste0("rho=",round(mytest$estimate,digits=4)),cex=mysize,col="black")
myplot<-ggplot(mydata,aes(x=mydata$rho2,y=log(mydata$pKuli+10^(-16))))+stat_binhex()+scale_fill_gradient(limits=c(1,10^8), low="white", high="cadetblue4",trans="log")+theme_bw()+theme(legend.position="none")
mytest<-cor.test(mydata$rho2,mydata$pKuli)
myplot2<-myplot+labs(x="rho2",y="log(pKuli+10^-6)")+annotate("text",x=0.7,y=-1,label=paste0("p-value=",mytest$p.value),cex=mysize,col="black")+annotate("text",x=0.7,y=-3,label=paste0("rho=",round(mytest$estimate,digits=4)),cex=mysize,col="black")
myplot<-ggplot(mydata,aes(x=mydata$rho2,y=log(mydata$prho2+10^(-16))))+stat_binhex()+scale_fill_gradient(limits=c(1,10^8), low="white", high="cadetblue4",trans="log")+theme_bw()+theme(legend.position="none")
mytest<-cor.test(mydata$rho2,mydata$rho2)
myplot3<-myplot+labs(x="rho2",y="log(prho2+10^-6)")+annotate("text",x=0.7,y=-1,label=paste0("p-value=",mytest$p.value),cex=mysize,col="black")+annotate("text",x=0.7,y=-3,label=paste0("rho=",round(mytest$estimate,digits=4)),cex=mysize,col="black")
bottomrow<-plot_grid(myplot1,myplot2,myplot3,align='h',labels=c('D','E','F'),ncol=3)
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/pvalues_scatter_signfdrfromT2.pdf")
plot_grid(toprow,bottomrow,ncol=1,align='h')
dev.off()
#I checked well p-values and behave well. What I have to notice is that the values for which I have p-value 1 in pkuli and low in rho2 and T2 are those for which
#I really have strange things happening, for example in chr1-chr10, 16890671	47000146, for which in all pops I have all heterozygotes.


#==============================================================================================================^
#===========================================OVERLAP BETWEEN POPS===============================v
#cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/anal/sorted.res | awk '{if (NF==23) {print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/anal/sorted2.res
library(data.table)
data2<-fread("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/anal/sorted2.res") #2G
names(data2)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
#data2<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/anal/sorted2.res")
#cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/*minilog | grep ncomparisons > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/minilog2.tab
#cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/new2/*minilog | grep ncomparisonspair > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/new2/minilog2.tab
ncomppair<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new2/minilog2.tab",header=FALSE)
setwd("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/")
ncomp<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/minilog2.tab",header=FALSE)
ncomppair<-sapply(2:dim(ncomppair)[2], function(x) sum(ncomppair[,x]))
ncomppair<-t(matrix(ncomppair,ncol=12))

mydata<-data2[data2$nA>=0.05*data2$Nse,]
mydata<-mydata[mydata$nA<=0.95*mydata$Nse,]
mydata<-mydata[mydata$nB>=0.05*mydata$Nse,]
mydata<-mydata[mydata$nB<=0.95*mydata$Nse,]

#I could ask: given that significant in one population, what's the chance of being sign in others
#save(mydata,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/mypairs.RData")
#load(file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/mypairs.RData")
rm(data2)
data2_pop_T2fdr<-list()
for (i in 1:12)
{
data2_pop_T2fdr[[i]]<-p.adjust(p=mydata$T2[mydata$pop==i], method = "fdr", n = sum(ncomp[,i+1])) 
}
data2_pop_T2fdr_sign<-list()
#data.table saves your life, so much better, something like 1000 times faster!
for (i in 1:12)
{
  tempdata<-mydata[mydata$pop==i][data2_pop_T2fdr[[i]]<0.05]
  data2_pop_T2fdr_sign[[i]]<-tempdata[, paste(chrA,chrB,posA,posB,sep=".")]
}
#not very informative, I guess p-value would be more
mymat<-matrix(0,nrow=12,ncol=12)
for (i in 1:12)
{
  for (j in 1:12)
  {
  mymat[i,j]<-sum(as.numeric(!is.na(match(data2_pop_T2fdr_sign[[i]],data2_pop_T2fdr_sign[[j]]))))/length(data2_pop_T2fdr_sign[[i]])
  }
}
#how to calculate p-value?
#if I have same sets of pairs N where to start from, then I have:
#ntotA and ntotB, nA and nB, prob that nA and nB overlap more than that. 
#fisher's exact test as described in rpackages.ianhowson.com/bioc/GeneOverlap/man/GeneOverlap.html
#I can do it in two ways,taking only pairs with same freq or all the ones I considered.
#method A) considering all possible comparisons
mymat<-matrix(0,nrow=12,ncol=12)
mymat2<-matrix(0,nrow=12,ncol=12)
for (i in 1:12)
{
  for (j in 1:12)
  { 
  print(c(i,j))
  mynooverlap<-sum(ncomp[,13+1]) - length(union(data2_pop_T2fdr_sign[[j]],data2_pop_T2fdr_sign[[i]]))
  myoverlap<-length(intersect(data2_pop_T2fdr_sign[[j]],data2_pop_T2fdr_sign[[i]]))
  mynooverlapj<-length(setdiff(data2_pop_T2fdr_sign[[j]],data2_pop_T2fdr_sign[[i]]))
  mynooverlapi<-length(setdiff(data2_pop_T2fdr_sign[[i]],data2_pop_T2fdr_sign[[j]]))
  if (myoverlap!=0 && mynooverlapi==0 && mynooverlapj==0) {mymat[i,j]<-0; mymat2[i,j]<-Inf;} else 
  if (myoverlap!=0) {
    mytest<-fisher.test(matrix(c(mynooverlap,mynooverlapj,mynooverlapi,myoverlap),nrow=2))
    mymat[i,j]<-mytest$p.value
    mymat2[i,j]<-mytest$estimate} else {mymat[i,j]<-1; mymat2[i,j]<-1;}
  }
}
poverlap.methodA.p.value<-mymat
poverlap.methodA.odds<-mymat2

#method B) considering only possible pairs for that pop comparisons

mymat<-matrix(0,nrow=12,ncol=12)
mymat2<-matrix(0,nrow=12,ncol=12)
for (i in 1:12)
{
  tempi<-mydata[mydata$pop==i][data2_pop_T2fdr[[i]]<0.05]
  tempi[,pairid:=data2_pop_T2fdr_sign[[i]]]
#  data2_pop_T2fdr_sign[[i]]<-tempdata[, paste(chrA,chrB,posA,posB,sep=".")]
  tempi<-tempi[0.05 <= tempi$nA/tempi$Nse & 0.95 >= tempi$nA/tempi$Nse,]
  tempi<-tempi[0.05 <= tempi$nA/tempi$Nse & 0.95 >= tempi$nA/tempi$Nse,]  
  for (j in 1:12)
  { 
  print(c(i,j))
  tempj<-mydata[mydata$pop==j][data2_pop_T2fdr[[j]]<0.05]
  tempj[,pairid:=data2_pop_T2fdr_sign[[j]]]
#  data2_pop_T2fdr_sign[[i]]<-tempdata[, paste(chrA,chrB,posA,posB,sep=".")]
  tempj<-tempj[0.05 <= tempj$nA/tempj$Nse & 0.95 >= tempj$nA/tempj$Nse,]
  tempj<-tempj[0.05 <= tempj$nA/tempj$Nse & 0.95 >= tempj$nA/tempj$Nse,]  
  mynooverlap<-sum(ncomppair[max(i,j),min(i,j)]) - length(union(tempi$pairid,tempj$pairid))
  myoverlap<-length(intersect(tempj$pairid,tempi$pairid))
  mynooverlapj<-length(setdiff(tempj$pairid,tempi$pairid))
  mynooverlapi<-length(setdiff(tempi$pairid,tempj$pairid))
  if (myoverlap!=0 && mynooverlapi==0 && mynooverlapj==0) {mymat[i,j]<-0; mymat2[i,j]<-Inf;} else 
  if (myoverlap!=0) {
    mytest<-fisher.test(matrix(c(mynooverlap,mynooverlapj,mynooverlapi,myoverlap),nrow=2))
    mymat[i,j]<-mytest$p.value
    mymat2[i,j]<-mytest$estimate} else {mymat[i,j]<-1; mymat2[i,j]<-1;}
  }
}

poverlap.methodB.p.value<-mymat
poverlap.methodB.odds<-mymat2
poverlap<-list(poverlap.methodA.p.value,poverlap.methodA.odds,poverlap.methodB.p.value,poverlap.methodB.odds)
#save(poverlap,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/poverlap.RData")
load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/poverlap.RData")
poverlap.methodB.odds<-poverlap[[4]]
poverlap.methodB.odds[poverlap.methodB.odds==Inf]<--100
poverlap.methodB.odds[poverlap.methodB.odds==-100]<-max(poverlap.methodB.odds)+max(poverlap.methodB.odds)/10
poverlap.methodB.odds<-log(poverlap.methodB.odds,10)
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
rownames(poverlap.methodB.odds)<-popnames
colnames(poverlap.methodB.odds)<-popnames

#show odds because p-values are all <10^-16
#poverlap.methodB.odds<-log(poverlap.methodB.odds)
col_palette<-colorRampPalette(c("papayawhip","yellow","red"))(n=8)
col_breaks<-seq(0,max(poverlap.methodB.odds),(max(poverlap.methodB.odds))/8)
length(col_breaks)
length(col_palette)
library(gplots)
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/heatmap.methodB.odds.pdf")
heatmap.2(poverlap.methodB.odds,trace="none",dendrogram="row",breaks=col_breaks,col=col_palette)
dev.off()

#==============================================================================================================^




























#------------------------------------explorative plots all pops together----------------------------v
colpops<-rep("cadetblue",nspops[1])#TSI
colpops<-c(colpops,rep("cadetblue3",nspops[2]))#IBS
colpops<-c(colpops,rep("chartreuse3",nspops[3]))#PUR
colpops<-c(colpops,rep("khaki1",nspops[4]))#GWD
colpops<-c(colpops,rep("brown1",nspops[5]))#CHB
colpops<-c(colpops,rep("brown",nspops[6]))#JPT
colpops<-c(colpops,rep("firebrick2",nspops[7]))#CHS
colpops<-c(colpops,rep("deepskyblue3",nspops[8]))#FIN
colpops<-c(colpops,rep("darkolivegreen3",nspops[9]))#ACB
colpops<-c(colpops,rep("gold1",nspops[10]))#YRI
colpops<-c(colpops,rep("coral2",nspops[11]))#KHV
colpops<-c(colpops,rep("lightpink",nspops[12]))#STU
isq<-1
seqcenter<-c()
for (isq in 1:12)
{
seqcenter<-c(seqcenter,as.character(infosamples[[isq]]$Main.project.LC.Centers))
}
seqcenter_col<-rep(0,length(seqcenter))
seqcenter_col[seqcenter=='BI']<-"firebrick"
seqcenter_col[seqcenter=='BI,MPIMG']<-"firebrick1"
seqcenter_col[seqcenter=='MPIMG']<-"coral2"
seqcenter_col[seqcenter=='BGI']<-"cadetblue"
seqcenter_col[seqcenter=='BCM,BGI']<-"cadetblue2"
seqcenter_col[seqcenter=='BCM']<-"aquamarine2"
seqcenter_col[seqcenter=='ILLUMINA']<-"chartreuse3"
seqcenter_col[seqcenter=='SC']<-"brown1"
seqcenter_col[seqcenter=='WUGSC']<-"darkgoldenrod1"
seqcenter_col
#pops are separated well when all, but not for significant (that seems to reflect sequencing center problems
pdf("nABallpops_abovefreqthr.pdf")
plot(mylog2[order(mylog2)],col=colpops[order(mylog2)],pch=19,ylab="nAB all",xlab="sample")
dev.off()
pdf("nABallpops_abovefreqthrsign.pdf")
plot(mylog3[order(mylog3)],col=colpops[order(mylog3)],pch=19,ylab="nAB significant",xlab="sample")
dev.off()
pdf("nABallpops_abovefreqthr_SC.pdf")
plot(mylog2[order(mylog2)],col=seqcenter_col[order(mylog2)],pch=19,ylab="nAB all",xlab="sample")
dev.off()
pdf("nABallpops_abovefreqthrsign_SC.pdf")
plot(mylog3[order(mylog3)],col=seqcenter_col[order(mylog3)],pch=19,ylab="nAB significant",xlab="sample")
dev.off()
#---------------------------------------------------------------------------------------------------^
#TAG 'excess of linkage'
nperm<-10
reschr_count<-rep(0,nperm-1)
combined_fdr_l<-list()
for (myperm in 2:nperm)
{
system(paste0("awk '{if (NF==23){print}}' /mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/anal/sorted.res > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/anal/sorted2.res"))
system(paste0("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/*minilog | grep ncomparisons | awk '{print $14}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/ncomparisons.txt"))
}
for (myperm in 2:nperm)
    {
    data<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/anal/sorted2.res"))
    print(dim(data))
    names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
    ncomp<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/ncomparisons.txt"),header=FALSE)
    ncomp<-sum(ncomp$V1)
    dataFIN<-subset(data,data$pop==7)
    combined_fdr<-p.adjust(p=dchisq(dataFIN$X,df=2*dataFIN$popX), method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. #This is why probably at the beginning I was using a smaller cutoff.
    fullsign_reschr<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
    combined_fdr_l[[myperm]]<-combined_fdr
    reschr_count[myperm]<-dim(fullsign_reschr)[1]
    }
save("combined_fdr_l",file="~/Dropbox/LDLD/coding1000g_combined_fdr_l.RData")
load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/fullsign.RData")

pdf("~/Dropbox/LDLD/coding1000g/distr_fdr_datavsempirical.pdf")
histnull<-hist(sqrt(unlist(combined_fdr_l)),xlim=c(0,0.3),breaks=30);
hist(sqrt(fullsign$combined_fdr[fullsign$combined_fdr<0.05]),xlim=c(0,0.3),breaks=histnull$breaks,col=rgb(0.1,0.1,0.1,0.8),main="sqrt(fdr)",xlab="sqrt(fdr)")
histnull$counts<-histnull$counts/(length(combined_fdr_l)-1)
plot(histnull,add=TRUE,col=rgb(0.8,0.8,0.8,0.6),breaks=30)
dev.off()

nperm<-10
reschr_count<-rep(0,nperm-1)
pval_l<-list()
for (myperm in 2:nperm)
    {
    data<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/anal/sorted2.res"))
    print(dim(data))
    names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
    ncomp<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/ncomparisons.txt"),header=FALSE)
    ncomp<-sum(ncomp$V1)
    dataFIN<-subset(data,data$pop==7)
    pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)
    pval<-pval[pval<0.01]
    pval_l[[myperm]]<-pval
    }

save(pval_l,file="~/Dropbox/LDLD/coding1000g/coding1000g_pval_l.RData")
pdf("~/Dropbox/LDLD/coding1000g/distr_combinedpval_datavsempirical.pdf")
mythr<-0.00001
histnull<-hist(sqrt(unlist(pval_l)[unlist(pval_l)<mythr]),xlim=c(0,0.001),breaks=30);
hist(sqrt(combined_pval[combined_pval<mythr]),xlim=c(0,0.001),breaks=histnull$breaks,col=rgb(0.1,0.1,0.1,0.8),main="",xlab="sqrt(combined p-value)")
histnull$counts<-histnull$counts/(length(pval_l)-1)
plot(histnull,col=rgb(0.8,0.8,0.8,0.6),add=TRUE)
dev.off()

fdrres<-empiricall_fdr(combined_pval,unlist(pval_l),length(pval_l)-1)
save(fdrres,file="~/Dropbox/LDLD/coding1000g/coding1000g_fdr.RData")

length(fdrres[fdrres<0.05])
pdf("~/Dropbox/LDLD/coding1000g//fdrplot_datavsempirical.pdf")
plot(fdrres,type="l",xlab="index",ylab="fdr")
dev.off()

length(histnull$breaks)
# TAG 'histograms to illustrate separation method in analysesLDLD_header'------------------------------------------------vvvv

i<-1
if (FALSE)
{
  i<-11
  nAB_l<-unlist(logpop[[i]][1])
  mydata<-nAB_l[order(nAB_l)]
  mod2<-normalmixEM(mydata,lambda=1/2,mu=breakshist(2,mydata),sigma=rep(10,2))
  mod3<-normalmixEM(mydata,lambda=1/3,mu=breakshist(3,mydata),sigma=rep(10,3))
  rescalingfactor<-3*10^4
  pdf("hist_nAB_models.pdf")
  par(mfrow=c(2,3))
  hist(mydata,breaks=20,main="2 k",xlab="nAB")
  x <- seq(min(mydata), max(mydata), length=100)
  curve(dnorm(x, mean=mean(mydata), sd=sd(mydata))*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  hist(mydata,breaks=20,main="4 k",xlab="nAB")
  x <- seq(min(mydata), max(mydata), length=100)
  curve(dnorm(x, mean=mod2$mu[1], sd=mod2$sigma[1])*rescalingfactor, 
	    col="darkblue", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod2$mu[2], sd=mod2$sigma[2])*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  hist(mydata,breaks=20,main="6 k",xlab="nAB")
  x <- seq(min(mydata), max(mydata), length=100)
  curve(dnorm(x, mean=mod3$mu[1], sd=mod3$sigma[1])*rescalingfactor, 
	    col="darkblue", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod3$mu[2], sd=mod3$sigma[2])*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod3$mu[3], sd=mod3$sigma[3])*rescalingfactor, 
	    col="darkgreen", lwd=2, add=TRUE, yaxt="n")
  dev.off()
  i<-9
  nAB_l<-unlist(logpop[[i]][1])
  mydata<-nAB_l[order(nAB_l)]
  mod2<-normalmixEM(mydata,lambda=1/2,mu=breakshist(2,mydata),sigma=rep(10,2))
  mod3<-normalmixEM(mydata,lambda=1/3,mu=breakshist(3,mydata),sigma=rep(10,3))
  mod4<-normalmixEM(mydata,lambda=1/3,mu=breakshist(4,mydata),sigma=rep(10,4))
  mod5<-normalmixEM(mydata,lambda=1/3,mu=breakshist(5,mydata),sigma=rep(10,5))
  mod6<-normalmixEM(mydata,lambda=1/3,mu=breakshist(6,mydata),sigma=rep(10,6))
  rescalingfactor<-10^4
  dev.new()
  pdf("hist_nAB_models.pdf")
  par(mfrow=c(2,3))
  hist(mydata,breaks=20,main="2 k",xlab="nAB")
  x <- seq(min(mydata), max(mydata), length=100)
  curve(dnorm(x, mean=mean(mydata), sd=sd(mydata))*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  hist(mydata,breaks=20,main="4 k",xlab="nAB")
  x <- seq(min(mydata), max(mydata), length=100)
  curve(dnorm(x, mean=mod2$mu[1], sd=mod2$sigma[1])*rescalingfactor, 
	    col="darkblue", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod2$mu[2], sd=mod2$sigma[2])*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  hist(mydata,breaks=20,main="6 k",xlab="nAB")
  curve(dnorm(x, mean=mod3$mu[1], sd=mod3$sigma[1])*rescalingfactor, 
	    col="darkblue", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod3$mu[2], sd=mod3$sigma[2])*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod3$mu[3], sd=mod3$sigma[3])*rescalingfactor, 
	    col="darkgreen", lwd=2, add=TRUE, yaxt="n")
  hist(mydata,breaks=20,main="8 k",xlab="nAB")
  curve(dnorm(x, mean=mod4$mu[1], sd=mod4$sigma[1])*rescalingfactor, 
	    col="darkblue", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod4$mu[2], sd=mod4$sigma[2])*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod4$mu[3], sd=mod4$sigma[3])*rescalingfactor, 
	    col="darkgreen", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod4$mu[4], sd=mod4$sigma[4])*rescalingfactor, 
	    col="gray", lwd=2, add=TRUE, yaxt="n")
  hist(mydata,breaks=20,main="10 k",xlab="nAB")
  curve(dnorm(x, mean=mod5$mu[1], sd=mod5$sigma[1])*rescalingfactor, 
	    col="darkblue", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod5$mu[2], sd=mod5$sigma[2])*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod5$mu[3], sd=mod5$sigma[3])*rescalingfactor, 
	    col="darkgreen", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod5$mu[4], sd=mod5$sigma[4])*rescalingfactor, 
	    col="gray", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod5$mu[5], sd=mod5$sigma[5])*rescalingfactor, 
	    col="orange", lwd=2, add=TRUE, yaxt="n")
  hist(mydata,breaks=20,main="12 k",xlab="nAB")
  curve(dnorm(x, mean=mod6$mu[1], sd=mod6$sigma[1])*rescalingfactor, 
	    col="darkblue", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod6$mu[2], sd=mod6$sigma[2])*rescalingfactor, 
	    col="darkred", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod6$mu[3], sd=mod6$sigma[3])*rescalingfactor, 
	    col="darkgreen", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod6$mu[4], sd=mod6$sigma[4])*rescalingfactor, 
	    col="gray", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod6$mu[5], sd=mod6$sigma[5])*rescalingfactor, 
	    col="orange", lwd=2, add=TRUE, yaxt="n")
  curve(dnorm(x, mean=mod6$mu[6], sd=mod6$sigma[6])*rescalingfactor, 
	    col="pink", lwd=2, add=TRUE, yaxt="n")
  dev.off()
}
#--------------------------------------------------^
#===========================CREATE LOGOS================
#TAG 'CREATE LOGOS'
#cat above95/coding1000g/anal/snpsA.bed above95/coding1000g/anal/snpsB.bed | sort -Vu -k1,1 -k2,2 > above95/coding1000g/anal/snps.bed
data<-read.table("above95/coding1000g/anal/snps.bed",stringsAsFactors =FALSE,header=FALSE)
names(data)<-c("chr","start","end")
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end) 
dim(data) #47193
insnpsbed<-data[data$chr=="chr19",]
insnpsbed[insnpsbed$end==11687195,]
createbedlogo(data,"above95/coding1000g/logosign")

#cat above95/coding1000g/chr*.tab | grep -v '#' | sort -Vu -k1,1 -k2,2 | awk -v OFS='\t' '{print "chr"$1,$2-1,$2}' > above95/coding1000g/chrtab.bed #background from files .tab #NB: a bit more than bed2 in $TARGET. check why!
#bedtools intersect -b above95/coding1000g/anal/snps.bed -a above95/coding1000g/chrtab.bed -v > above95/coding1000g/notsnps.bed
data<-read.table("above95/coding1000g/notsnps.bed",stringsAsFactors =FALSE,header=FALSE)
names(data)<-c("chr","start","end")
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
createbedlogo(data,"above95/coding1000g/logonotsign")

data1<-cbind(data,0,0,0,0)
write.table(myres2,file="above95/coding1000g/logonotsign.temp",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);

#=========================CHECK IMBALANCE in .BAM FILES by using LOGO files==================
#paste <(bedtools intersect -a <( cat above95/coding1000g/anal/snps.bed | sed 's/chr//g' ) -b <( cat above95/coding1000g/chr*.tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2) -loj | awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9,$10}' | sort -Vu -k1,1 -k2,2) above95/coding1000g/logosign | sort -Vu -k1,1 -k2,2 > above95/coding1000g/logosign.bed
#bamlogos for significant snps
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=5000,startfrom=1,filetemp="temp1")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp1")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_1to5k.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=25000,startfrom=20001,filetemp="temp20k")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp20k")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call","tot_samples")
write.table(myres2,file="above95/coding1000g/anal/logo2_20kto25k.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=2000,startfrom=1)
myres2<-as.data.frame(myres,row.names=FALSE)
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_1to2000.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=4000,startfrom=2001)
myres2<-as.data.frame(myres,row.names=FALSE)
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_2001to4000.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=15000,startfrom=10001,filetemp="temp2")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp2")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_10kto15k.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=20000,startfrom=15001,filetemp="temp15k")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp15k")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_15kto20kend.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=40000,startfrom=30001,filetemp="temp30k")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp30k")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call","tot_samples")
write.table(myres2,file="above95/coding1000g/anal/logo2_30kto40kend.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=TRUE,startfrom=40001,filetemp="temp40k")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp40k")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call","tot_samples")
write.table(myres2,file="above95/coding1000g/anal/logo2_40ktoend.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=30000,startfrom=25001,filetemp="temp25k")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp25k")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call","tot_samples")
write.table(myres2,file="above95/coding1000g/anal/logo2_25kto30k.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
data1<-read.table("above95/coding1000g/anal/logo2_1to5k.bed",stringsAsFactors =FALSE,header=TRUE)
data2<-read.table("above95/coding1000g/anal/logo2_5kto10k.bed",stringsAsFactors =FALSE,header=TRUE)
data3<-read.table("above95/coding1000g/anal/logo2_10kto15k.bed",stringsAsFactors =FALSE,header=TRUE)
data4<-read.table("above95/coding1000g/anal/logo2_15kto20kend.bed",stringsAsFactors =FALSE,header=TRUE)
data5<-read.table("above95/coding1000g/anal/logo2_20kto25k.bed",stringsAsFactors =FALSE,header=TRUE)
data6<-read.table("above95/coding1000g/anal/logo2_25kto30k.bed",stringsAsFactors =FALSE,header=TRUE)
data7<-read.table("above95/coding1000g/anal/logo2_30kto40kend.bed",stringsAsFactors =FALSE,header=TRUE)
data8<-read.table("above95/coding1000g/anal/logo2_40ktoend.bed",stringsAsFactors =FALSE,header=TRUE)
data<-rbind(data1[,1:11],data2[,1:11],data3[,1:11],data4[,1:11],data5[,1:11],data6[,1:11],data7[,1:11],data8[,1:11])
write.table(data,file="above95/coding1000g/anal/logo2_all.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE)

#bamlogos for not significant snps
#paste <(bedtools intersect -a <( cat above95/coding1000g/notsnps.bed | sed 's/chr//g' ) -b <( cat above95/coding1000g/chr*.tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2) -loj | awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9,$10}' | sort -Vu -k1,1 -k2,2) above95/coding1000g/logonotsign | sort -Vu -k1,1 -k2,2 > above95/coding1000g/logonotsign.bed
myres<-bamlogo("above95/coding1000g/logonotsign.bed",tilltheend=FALSE,till=5000,startfrom=1,filetemp="temp1")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp1")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_1to5k_notsign.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logonotsign.bed",tilltheend=FALSE,till=10000,startfrom=5001,filetemp="temp5k")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp5k")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_5kto10k_notsign.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logonotsign.bed",tilltheend=FALSE,till=20000,startfrom=10001,filetemp="temp10knot")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp10knot")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_10kto20k_notsign.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
myres<-bamlogo("above95/coding1000g/logonotsign.bed",tilltheend=TRUE,startfrom=20001,filetemp="temp20knot")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp20knot")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_20ktoend_notsign.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);

#========================================================ANALYSES IMBALANCE: not sign

data1<-read.table("above95/coding1000g/anal/logo2_1to5k_notsign.bed",stringsAsFactors =FALSE,header=TRUE)
data2<-read.table("above95/coding1000g/anal/logo2_5kto10k_notsign.bed",stringsAsFactors =FALSE,header=TRUE)
data3<-read.table("above95/coding1000g/anal/logo2_10kto20k_notsign.bed",stringsAsFactors =FALSE,header=TRUE)
data4<-read.table("above95/coding1000g/anal/logo2_20ktoend_notsign.bed",stringsAsFactors =FALSE,header=TRUE)
names(data2)<-names(data1)
names(data3)<-names(data1)
names(data4)<-names(data1)
data<-rbind(data1[,1:11],data2[,1:11],data3[,1:11],data4[,1:11])
write.table(data,file="above95/coding1000g/anal/logo2_allnotsign.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE)
logonotsign<-read.table("above95/coding1000g/anal/logo2_allnotsign.bed",stringsAsFactors =FALSE,header=TRUE)
logonotsign[is.na(logonotsign)] <- 0
logonotsign<-cbind(logonotsign,t(sapply(1:dim(logonotsign)[1],function(x) prHvsE(logonotsign,x))))
dim(logonotsign)
res<-logo_below_thr(logonotsign,0.00001) #0.03417355
res<-logo_below_thr(logonotsign,0.05) #0.1843388



#========================================================ANALYSES IMBALANCE: sign
logosign<-read.table("above95/coding1000g/anal/logo2_all.bed",stringsAsFactors =FALSE,header=TRUE)
logosign[is.na(logosign)] <- 0
logosign<-cbind(logosign,t(sapply(1:dim(logosign)[1],function(x) prHvsE(logosign,x))))
res<-logo_below_thr(logosign,0.00001) #0.05235133
res<-logo_below_thr(logosign,0.05) #0.1960372
logosign_potentialH<-res[[2]]
dataLDLD<-sign
head(res[[1]])[,c(1,3:11,13)]
options(xtable.floating = FALSE)
xtable(head(res[[1]])[,c(1,3:11,13)],digits=8)

system.time(res10k<-ave_LDLD_BED(sign,logosign,myfrom=1,myto=10000))
system.time(res20k<-ave_LDLD_BED(sign,logosign,myfrom=10001,myto=20000))
system.time(res30k<-ave_LDLD_BED(sign,logosign,myfrom=20001,myto=30000))
system.time(res40k<-ave_LDLD_BED(sign,logosign,myfrom=30001,myto=40000))
system.time(resend<-ave_LDLD_BED(sign,logosign,myfrom=40001))


logosign_links<-rbind(res10k,res20k,res30k,res40k,resend)
logosign<-cbind(logosign,logosign_links)
ln<-length(names(logosign))
names(logosign)[ln-2]<-"nlinks"
names(logosign)[ln-1]<-"ave_fdr"
names(logosign)[ln]<-"min_fdr"
#if very low min fdr then high ave fdr..very strange! maybe because more linked snps in total and then ave increases
write.table(logosign,file="above95/coding1000g/anal/logosign_withfdr.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
logosign<-read.table("above95/coding1000g/anal/logosign_withfdr.bed",header=TRUE);
cor(logosign$ave_fdr,logosign$min_fdr,method="spearman") #0.8 #cool, correlate very strongly
cor(logosign$ave_fdr,logosign$nlinks,method="spearman") #-0.09677386
cor(logosign$min_fdr,logosign$nlinks,method="spearman") #-0.4916015 #also, more significant more links

cor(logosign_potentialH$pvalue_no05_alt,logosign_potentialH$nlinks,method="spearman") #more links, the more unlikely to be a real heterozygous #-0.02416674
cor(logosign_potentialH$pvalue_no05_alt,logosign_potentialH$min_fdr,method="spearman") #0.03156202
cor.test(logosign_potentialH$pvalue_no05_alt,logosign_potentialH$min_fdr,method="spearman") #p-value = 0.000000002473

temp<-logosign_potentialH[sample(1:dim(logosign_potentialH)[1],3000),]
plot(temp$pvalue_no05_alt,temp$min_fdr)

#======================================================REMOVE CENTRAL CLUSTER
#TAG 'REMOVE CENTRAL CLUSTER'
myres<-retrieve_partners_collapse(as.character(1.881918),dataLDLD,generatecollapse=TRUE)
newnodes<-retrieve_partners(1,881918,dataLDLD)
chr<-1
pos<-881918
everythingislinked<-run_through_links(dataLDLD$chrA[1],dataLDLD$posA[1],dataLDLD) #nothing is left like this!
everythingislinked_highlysign<-run_through_links(highlysign$chrA[1],highlysign$posA[1],highlysign) #nothing is left like this!
length(everythingislinked[[1]])
dim(everythingislinked[[2]]) #6359
everythingislinked_highlysign
length(everythingislinked_highlysign[[1]]) #807
dim(everythingislinked_highlysign[[2]])[1] #774
everythingislinked[[2]]
dim(dataLDLDcollapsed)
length(dataLDLDlinked)
save.image("LDLD040216.RData")
load("LDLD040216.RData")
everythingislinked[[2]]

rebuild_dataLDLD<-function(removed_snps,dataLDLD,fromremovedsnps=TRUE)
{
  if (fromremovedsnps)
    {
    dataLDLDA<-as.character(paste0(dataLDLD[,1],".",dataLDLD[,3]))
    res<-dataLDLD[(is.na(match(dataLDLDA,removed_snps))),]
    }
  else
    {
    dataLDLDA<-as.character(paste0(dataLDLD[,1],".",dataLDLD[,3]))
    res<-dataLDLD[(!is.na(match(dataLDLDA,removed_snps))),]
    }
return(res)
}
rebuild_dataLDLD<-function(removed_snps,dataLDLD)
{
  dataLDLDA<-as.character(paste0(dataLDLD[,1],".",dataLDLD[,3]))
  #dataLDLDB<-as.character(paste0(dataLDLD$chrB,".",dataLDLD$posB))
  #(!is.na(match(dataLDLDA,removed_snps)))+(!is.na(match(dataLDLDB,removed_snps)))
return(dataLDLD[(is.na(match(dataLDLDA,removed_snps))),])
}


everythingislinked_LDLDleft<-rebuild_dataLDLD(everythingislinked[[1]],dataLDLD)
everythingislinked_logoleft<-rebuild_dataLDLD(everythingislinked[[1]],logosign)
dim(everythingislinked_LDLDleft)[1] #6364
dim(everythingislinked_logoleft)[1] #9143
head(everythingislinked_LDLDleft)

cor(everythingislinked_logoleft$pvalue_no05_alt,everythingislinked_logoleft$nlinks,method="spearman") #-0.02628085
cor(everythingislinked_logoleft$pvalue_no05_alt,everythingislinked_logoleft$min_fdr,method="spearman") #0.01029192 vs 0.03156202
cor.test(everythingislinked_logoleft$pvalue_no05_alt,everythingislinked_logoleft$min_fdr,method="spearman") #0.3251
dim(everythingislinked_logoleft)
res<-logo_below_thr(everythingislinked_logoleft,0.00001) #0.03417355
res<-logo_below_thr(everythingislinked_logoleft,0.05) #0.1843388

unique(everythingislinked[[1]])

#======================================================trios from this
#TAG 'trios from this'
#generate BED files of snps removed.
if (FALSE)
  {
  Cbed<-rebuild_dataLDLD(everythingislinked[[1]][order(everythingislinked[[1]])],logosign,fromremovedsnps=FALSE)
  write.table(Cbed[,1:3],file="above95/coding1000g/anal/Cbed.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
  system("mkdir above95/coding1000g/noC")
  i<-1
  system(paste0("cat above95/coding1000g/anal/Cbed.bed | sed 's/^/chr/' > above95/coding1000g/anal/Cbed.bed2"))
  for (i in 1:22)
    {
    system(paste0("head -1 above95/coding1000g/chr1.tab > above95/coding1000g/noC/chr",i,".tab"))
    system(paste0("cat above95/coding1000g/chr",i,".tab | awk '{if (NR>1){print}}' | sed 's/^/chr/' > above95/coding1000g/noC/chr",i,".tab2"))
    system(paste0("bedtools intersect -a above95/coding1000g/noC/chr",i,".tab2 -b above95/coding1000g/anal/Cbed.bed2 -v | sed 's/chr//' >> above95/coding1000g/noC/chr",i,".tab"))
    system(paste0("cat above95/coding1000g/noC/chr",i,".tab | gzip -9 > above95/coding1000g/noC/chr",i,".coding.vcf.gz"))
    }
  system("nohup ./ICLD11multi.sh /mnt/scratch/fabrizio/LDLD/above95/coding1000g/noC 2 coding 1 12 1 23 nspops.txt")
  }
system("cat above95/coding1000g/noC/*minilog | grep ncomparisons | awk '{print $14}' > above95/coding1000g/noC/ncomparisons.txt")
ncomp<-read.table("above95/coding1000g/noC/ncomparisons.txt",header=FALSE)
ncomp<-sum(ncomp$V1) #528075582
#temp<-everythingislinked_LDLDleft[everythingislinked_LDLDleft$popX!=999,]
temp<-dchisq(everythingislinked_LDLDleft$X,df=2*everythingislinked_LDLDleft$popX)
temp2<-everythingislinked_LDLDleft$T2[everythingislinked_LDLDleft$popX==999]
temp[everythingislinked_LDLDleft$popX==999]<-temp2
combined_fdr<-p.adjust(temp, method = "fdr", n = ncomp)
everythingislinked_LDLDleft$combined_fdr<-combined_fdr
FINsign_noC<-everythingislinked_LDLDleft[combined_fdr<0.05,]
write.table(FINsign_noC,file="above95/coding1000g/anal/sign_noC.res",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE)
write.table(cbind(paste0("chr",FINsign_noC$chrA),FINsign_noC$posA-1,FINsign_noC$posA),file="above95/coding1000g/anal/sign_noCA.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
write.table(cbind(paste0("chr",FINsign_noC$chrB),FINsign_noC$posB-1,FINsign_noC$posB),file="above95/coding1000g/anal/sign_noCB.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
system("sort -Vu -k1,1 -k2,2 above95/coding1000g/anal/sign_noCA.bed > sign_noCA.bed")
system("sort -Vu -k1,1 -k2,2 above95/coding1000g/anal/sign_noCB.bed > sign_noCB.bed")
system("./postproc_bedint_allchr.sh sign_noCA.bed sign_noCA.bed.syn synonymous")
system("./postproc_bedint_allchr.sh sign_noCA.bed sign_noCA.bed.nosyn missense")
system("./postproc_bedint_allchr.sh sign_noCB.bed sign_noCB.bed.syn synonymous")
system("./postproc_bedint_allchr.sh sign_noCB.bed sign_noCB.bed.nosyn missense")
system("bedtools intersect -a above95/coding1000g/anal/sign_noCA.bed -b <( cat sign_noCA.bed.syn sign_noCA.bed.nosyn ) -loj > sign_noCA.ann")
system("bedtools intersect -a above95/coding1000g/anal/sign_noCB.bed -b <( cat sign_noCB.bed.syn sign_noCB.bed.nosyn ) -loj > sign_noCB.ann")
system("awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9,$10,$12,$13,$14,$16}' sign_noCA.ann > sign_noCA.bed")
system("awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9,$10,$12,$13,$14,$16}' sign_noCB.ann > sign_noCB.bed")



dim(FINsign_noC)
FINsign_noC
min(temp2)
sum(temp<min(temp2))
FINsign_noC[,1:4]
FINsign_noC<-FINsign_noC[FINsign_noC$pfisher<0.01,]
dim(FINsign_noC)
FINsign_noC
write.table(FINsign_noC,file="above95/coding1000g/anal/FINsign_noC.res",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE)
mysnps<-rbind(cbind(FINsign_noC$chrA,FINsign_noC$posA),cbind(FINsign_noC$chrB,FINsign_noC$posB))
mysnps<-as.character(unique(paste0(mysnps[,1],".",mysnps[,2])))
mylogo<-logosign[!is.na(match(paste0(logosign[,1],".",logosign[,3]),mysnps)),]
res<-logo_below_thr(mylogo,0.05) #0.3076923
write.table(mylogo,file="above95/coding1000g/anal/FINsign_noC.logo",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE)
res[[1]]
length(mysnps) #21

signAbed<-as.data.frame(cbind(as.character(paste0('chr',FINsign_noC$chrA)),FINsign_noC$posA-1,FINsign_noC$posA),stringsAsFactors =FALSE)
signBbed<-as.data.frame(cbind(as.character(paste0('chr',FINsign_noC$chrB)),FINsign_noC$posB-1,FINsign_noC$posB),stringsAsFactors =FALSE)
names(signAbed)<-c("chr","start","end")
names(signBbed)<-c("chr","start","end")
signAbed$start<-as.numeric(signAbed$start)
signAbed$end<-as.numeric(signAbed$end)
signBbed$start<-as.numeric(signBbed$start)
signBbed$end<-as.numeric(signBbed$end)
setwd("/mnt/scratch/fabrizio/LDLD")
write.table(signAbed,file="above95/coding1000g/anal/snpsA_FIN.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file="above95/coding1000g/anal/snpsB_FIN.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
#sort -V -k1,1 -k2,2 above95/coding1000g/anal/snpsA_FIN.bed > temp
#cat temp > above95/coding1000g/anal/snpsA_FIN.bed
#sort -V -k1,1 -k2,2 above95/coding1000g/anal/snpsB_FIN.bed > temp
#cat temp > above95/coding1000g/anal/snpsB_FIN.bed
bedtools intersect -a above95/coding1000g/anal/snpsB_FIN.bed -b hgncmerged.bed -loj > above95/coding1000g/anal/FINB.bed
bedtools intersect -a above95/coding1000g/anal/snpsA_FIN.bed -b hgncmerged.bed -loj > above95/coding1000g/anal/FINA.bed
FINA<-read.table("FINA.bed",header=FALSE,sep='\t');
FINB<-read.table("FINB.bed",header=FALSE);
FINAB<-cbind(FINA[,7:9],FINB[,7:9])
options(xtable.floating = FALSE)
xtable(FINAB)

head(FINAB)

#======================================================================NETWORK-PLOTS
#TAG NETWORK-PLOTS
library(igraph)

#highlysign_subsampled
highlysign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.000000000000000000000000000000000000000000000000000000000000000000001,]
dim(highlysign)
net<-generate_nework_LDLD(highlysign)
myl_rt <- layout.reingold.tilford(net)
net<-generate_nework_LDLD(highlysign)
myl <- layout.graphopt(net)
par(mfrow=c(1,2))
#plot(net,vertex.label=NA,vertex.size=2,layout=myl,vertex.color=cols)
pdf("connettivity_highlysign_coding.pdf")
plot(net,vertex.label=NA,vertex.size=2,layout=myl)
plot(net,vertex.label=NA,vertex.size=2,layout=myl_rt )
dev.off()

#highlysign_subsampled
highlysign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
highlysign<-highlysign[sample(nrow(highlysign), 1000), ]
dim(highlysign)
net<-generate_nework_LDLD(highlysign)
myl_rt <- layout.reingold.tilford(net)
myl <- layout.graphopt(net)
par(mfrow=c(1,2))
#plot(net,vertex.label=NA,vertex.size=2,layout=myl,vertex.color=cols)
pdf("connettivity_fdr05subsampled1k_coding.pdf")
plot(net,vertex.label=NA,vertex.size=2,layout=myl)
plot(net,vertex.label=NA,vertex.size=2,layout=myl_rt )
dev.off()

#left after removal of big cluster
net<-generate_nework_LDLD(everythingislinked_LDLDleft)
everythingislinked_LDLDleft_subs<-everythingislinked_LDLDleft[sample(nrow(everythingislinked_LDLDleft), 1000), ]
net<-generate_nework_LDLD(everythingislinked_LDLDleft_subs)
myl <- layout.graphopt(net)
myl_rt <- layout.reingold.tilford(net)
pdf("connettivity_nocentersubsampled1k_coding.pdf")
plot(net,vertex.label=NA,vertex.size=2,layout=myl)
plot(net,vertex.label=NA,vertex.size=2,layout=myl_rt )
dev.off()

#left sign after removal of big cluster
FINsign_noC<-everythingislinked_LDLDleft[combined_fdr<0.05,]
dataLDLD<-FINsign_noC
dim(FINsign_noC)
net<-generate_nework_LDLD(dataLDLD)
myl <- layout.graphopt(net)
pal <- rainbow(22, alpha=.5)
cols<-pal[floor(id)]
par(mfrow=c(1,1))
pdf("connettivity_afterfiltering_coding.pdf")
plot(net,vertex.label=NA,vertex.size=4,layout=myl,vertex.color=cols)
dev.off()

myl <- layout.graphopt(net)
#pal <- rainbow(22, alpha=.5)
#cols<-pal[floor(id)]
par(mfrow=c(1,1))
#plot(net,vertex.label=NA,vertex.size=2,layout=myl,vertex.color=cols)
pdf("connettivity_fdr01subsampled1k_coding.pdf")
plot(net,vertex.label=NA,vertex.size=2,layout=myl)
dev.off()

#----TAG:info statistics in 1000genomes--v
if (FALSE)
{ 
  infosamplesFIN<-infosamples[[1]]
  infosamples$In.Final.Phase.Variant.Calling
  infosamplesFIN$ET.Pilot.Platforms
  infosamplesFIN$HC.Pilot.Platforms
  infosamplesFIN$Has.Exome.LOF.Genotypes
  infosamplesFIN$Main.Project..E.Platform
  infosamplesFIN$X..Targets.Covered.to.20x.or.greater
  #infosamplesFIN$EBV.Coverage #Epstein-Barr Virus coverage.
  infosamplesFIN$LC.Non.Duplicated.Aligned.Coverage
  info1000g$In.Final.Phase.Variant.Calling
  info1000g$ET.Pilot.Platforms
  str(info1000g)
  info1000g$ET.Pilot.Centers
  info1000g$EBV.Coverage
  info1000g$EBV.Coverage
  info1000g$Has.Omni.Genotypes
  info1000g$Has.Axiom.Genotypes
  info1000g$Has.Affy.6.0.Genotypes
  infosamplesFIN$Total.Exome.Sequence/1000000
}
#----------------------------------^


#================TAG_MISSENSEBIS===============================

#generate example bamtable from missensebis
if (FALSE) 
{
data<-read.table("above95/missensebis/anal/sorted.res")
names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
dataFIN<-subset(data,data$pop==7)
combined_fdr<-p.adjust(p=dchisq(dataFIN$X,df=2*dataFIN$popX), method = "fdr", n = 689277987)
highlysign<-cbind(dataFIN,combined_fdr)[combined_fdr<10^(-110),][1:20,]
highlysignAbed<-as.data.frame(cbind(as.character(paste0('chr',highlysign$chrA)),highlysign$posA-1,highlysign$posA),stringsAsFactors =FALSE)
highlysignBbed<-as.data.frame(cbind(as.character(paste0('chr',highlysign$chrB)),highlysign$posB-1,highlysign$posB),stringsAsFactors =FALSE)
names(highlysignAbed)<-c("chr","start","end")
names(highlysignBbed)<-c("chr","start","end")
highlysignAbed$start<-as.numeric(highlysignAbed$start)
highlysignAbed$end<-as.numeric(highlysignAbed$end)
highlysignBbed$start<-as.numeric(highlysignBbed$start)
highlysignBbed$end<-as.numeric(highlysignBbed$end)
#create logos
setwd("/mnt/scratch/fabrizio/LDLD")
write.table(highlysignAbed,file="above95/missensebis/anal/hsnpsA.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(highlysignBbed,file="above95/missensebis/anal/hsnpsB.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
system("cat above95/missensebis/anal/hsnpsA.bed above95/missensebis/anal/hsnpsB.bed | sort -Vu -k1,1 -k2,2 > above95/missensebis/anal/hsnps.bed")
data<-read.table("above95/missensebis/anal/hsnps.bed",stringsAsFactors =FALSE,header=FALSE)
names(data)<-c("chr","start","end")
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
createbedlogo(data,"above95/missensebis/logohighlysign")
system("paste <(bedtools intersect -a <( cat above95/missensebis/anal/hsnps.bed | sed 's/chr//g' ) -b <( cat above95/missensebis/chr*.tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2) -loj | awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9,$10}' | sort -Vu -k1,1 -k2,2) above95/missensebis/logohighlysign | sort -Vu -k1,1 -k2,2 > above95/missensebis/logohighlysign.bed")
myres<-bamlogo("above95/missensebis/logohighlysign.bed",tilltheend=FALSE,till=20)
myres2<-as.data.frame(myres,row.names=FALSE)
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/missensebis/anal/logo2_20.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE)
head(myres2)
library(xtable)
options(xtable.floating = FALSE)
xtable(myres2[1:20, ])
}
#================TAG_MISSENSE===============================
#generate example bamtable from missense
if (FALSE) 
{
data<-read.table("above95/missense/anal/sorted.res")
head(data)
names(data)<-c("chrA","chrB","posA","posB","nA","nB","N","nAB","nAA","nBB","D","rAB","Drandom","rABrandom","ppermD","pfisher","pKuli09","pchi","LLDA","LLDB")
combined_fdr<-p.adjust(p=data$pKuli09, method = "fdr", n = 689277987)
highlysign<-cbind(data,combined_fdr)[combined_fdr<0.5,][1:20,]
highlysignAbed<-as.data.frame(cbind(as.character(paste0('chr',highlysign$chrA)),highlysign$posA-1,highlysign$posA),stringsAsFactors =FALSE)
highlysignBbed<-as.data.frame(cbind(as.character(paste0('chr',highlysign$chrB)),highlysign$posB-1,highlysign$posB),stringsAsFactors =FALSE)
names(highlysignAbed)<-c("chr","start","end")
names(highlysignBbed)<-c("chr","start","end")
highlysignAbed$start<-as.numeric(highlysignAbed$start)
highlysignAbed$end<-as.numeric(highlysignAbed$end)
highlysignBbed$start<-as.numeric(highlysignBbed$start)
highlysignBbed$end<-as.numeric(highlysignBbed$end)
#create logos
setwd("/mnt/scratch/fabrizio/LDLD")
write.table(highlysignAbed,file="above95/missense/anal/hsnpsA.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(highlysignBbed,file="above95/missense/anal/hsnpsB.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
system("cat above95/missense/anal/hsnpsA.bed above95/missense/anal/hsnpsB.bed | sort -Vu -k1,1 -k2,2 > above95/missense/anal/hsnps.bed")
data<-read.table("above95/missense/anal/hsnps.bed",stringsAsFactors =FALSE,header=FALSE)
names(data)<-c("chr","start","end")
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
createbedlogo(data,"above95/missense/logohighlysign")
system("paste <(bedtools intersect -a <( cat above95/missense/anal/hsnps.bed | sed 's/chr//g' ) -b <( cat above95/missense/chr*.tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2) -loj | awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9,$10}' | sort -Vu -k1,1 -k2,2) above95/missense/logohighlysign | sort -Vu -k1,1 -k2,2 > above95/missense/logohighlysign.bed")
myres<-bamlogo("above95/missense/logohighlysign.bed",tilltheend=FALSE,till=15)
myres2<-as.data.frame(myres,row.names=FALSE)
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/missense/anal/logo2_20.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE)
head(myres2)
library(xtable)
options(xtable.floating = FALSE)
xtable(myres2[1:15, ])
}

highlysign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.01,]
dim(highlysign)
highlysign<-highlysign[sample(nrow(highlysign), 2000), ]
library(igraph)
node1<-as.numeric(paste0(highlysign$chrA,".",highlysign$posA))
node2<-as.numeric(paste0(highlysign$chrB,".",highlysign$posB))
id<-unique(c(node1,node2))
length(id)
floor(id)
nodes<-data.frame(id)
links<-data.frame(node1,node2,highlysign$combined_fdr)
dim(links)#149
#nodes <- read.csv("~/Dropbox/dropbox_tablet/manuals/R/GeneNet_works/Data/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
#links <- read.csv("~/Dropbox/dropbox_tablet/manuals/R/GeneNet_works/Data/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
#dim(nodes)
#dim(links)
#head(links)
#net <- graph.data.frame(links, nodes, directed=T)
names(links)<-c("from","to","weight")
net <- graph.data.frame(links, nodes, directed=F)
#myl <- layout.kamada.kawai(net)
#myl <- layout.random(net)
#myl <- layout.spring(net, mass=1)
#myl <- layout.fruchterman.reingold(net, repulserad=vcount(net)^3,area=vcount(net)^2.4)
#myl <- fruchterman.reingold.grid(net)
#myl <- layout.grid(net)
myl <- layout.graphopt(net)
#pal <- rainbow(22, alpha=.5)
#cols<-pal[floor(id)]
par(mfrow=c(1,1))
#plot(net,vertex.label=NA,vertex.size=2,layout=myl,vertex.color=cols)
pdf("connettivity_fdr01subsampled1k_coding.pdf")
plot(net,vertex.label=NA,vertex.size=2,layout=myl)
dev.off()
