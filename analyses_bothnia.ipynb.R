setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
options(scipen=999)
setwd("/mnt/scratch/fabrizio/LDLD")

# Prepare files

#./postproc_split2.sh botnia/FIN2
#./postproc_split2.sh botnia/FIN2/reschr1
#./postproc_split2.sh botnia/FIN2/reschr2
#./postproc_split2.sh botnia/FIN2/reschr3
#./postproc_split2.sh botnia/FIN2/reschr4
#./postproc_split2.sh botnia/FIN2/reschr5
#system(paste0("cat /mnt/scratch/fabrizio/LDLD/botnia/FIN2/*minilog | grep ncomparisons | awk '{print $2}' > /mnt/scratch/fabrizio/LDLD/botnia/FIN2/ncomparisons.txt"))

#---import files--------------------------------------------------------------------------------------v
samplesbotnia<-system("head -1 botnia/FIN2/chr1.tab | awk '
{
    for (i = 10; i <= NF; i++)
        print $i
}' ",intern=TRUE)
nspops<-as.numeric(system("cat /mnt/scratch/fabrizio/LDLD/nspopsbotnia.txt",intern=TRUE))
ncomp<-read.table("/mnt/scratch/fabrizio/LDLD/botnia/FIN2/ncomparisons.txt",header=FALSE)
ncomp<-sum(ncomp$V1)

#for factorial bigger than 90 (present in botnia) ICLD6multi fails because of some strange limits of mpg. 
#I have to use log factorial to make it manageable.
databot<-read.table("botnia/FIN2/anal/sorted.res")
names(databot)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","T2","Xtot","X","pop","popX")

databotcombined_pval<-dchisq(databot$X,df=2*databot$popX) 
databotcombined_fdr<-p.adjust(p=databotcombined_pval, method = "fdr", n = ncomp) #NB: i did not correct ncomp, so they might be a bit more
databot_sign<-databot[databotcombined_fdr<0.05,]
databot<-databot[order(databotcombined_fdr),]

#sum(as.numeric(databotcombined_fdr[databotcombined_fdr<0.05])) #0

nperm<-3
reschr_count<-rep(0,nperm)
combined_fdr_l<-list()
combined_pval_l<-list()
for (myperm in 1:nperm)
    {
    data<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/botnia/FIN2/reschr",myperm,"/anal/sorted.res"))
    names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","T2","Xtot","X","pop","popX")
    system(paste0("cat /mnt/scratch/fabrizio/LDLD/botnia/FIN2/reschr",myperm,"/*minilog | grep ncomparisons | awk '{print $2}' > /mnt/scratch/fabrizio/LDLD/botnia/FIN2/reschr",myperm,"/ncomparisons.txt"))
    ncomp<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/botnia/FIN2/reschr",myperm,"/ncomparisons.txt"),header=FALSE)
    ncomp<-sum(ncomp$V1)
    combined_pval<-dchisq(data$X,df=2*data$popX)  #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
    combined_fdr<-p.adjust(p=combined_pval, method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
    fullsign_reschr<-cbind(data,combined_fdr)[combined_fdr<0.05,]
    combined_fdr_l[[myperm]]<-combined_fdr
    combined_pval_l[[myperm]]<-combined_pval
    reschr_count[myperm]<-dim(fullsign_reschr)[1]
    }

    
res<-empiricall_fdr(databotcombined_fdr,c(combined_fdr_l[[1]],combined_fdr_l[[2]],combined_fdr_l[[3]]),3)
respval<-empiricall_fdr(databotcombined_pval,c(combined_pval_l[[1]],combined_pval_l[[2]],combined_pval_l[[3]]),3)

myhist<-hist(log(databotcombined_pval[databotcombined_pval<0.05]+10^(-20),10),breaks=30,col="gold",main="T2 vs T2 reshuffled",xlim=c(-5,-3),xlab="log(pvalue(T2)+10^-20)")
myhist$counts<-myhist$counts*10

plot(myhist,add=TRUE)


myhist<-hist(log(unlist(combined_pval_l)[unlist(combined_pval_l)<=0.001]+10^(-20),10),breaks=30,col="cadetblue",add=TRUE,xlim=c(-5,-3))
myhist$counts<-myhist$counts/3
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/removal1/empfdrvsT2_hist.pdf")
hist(log(databotcombined_pval[databotcombined_pval<=0.001]+10^(-20),10),breaks=30,col="gold",main="T2 vs T2 reshuffled",xlim=c(-5,-3),xlab="log(pvalue(T2)+10^-20)")
plot(myhist,add=TRUE,col="cadetblue")
dev.off()

hist(log(combined_pval_l[[1]][combined_pval_l[[1]]<0.05]+10^(-20),10),breaks=30,col="cadetblue",add=TRUE,xlim=c(-5,-3))

databotcombined_pval

sum(as.numeric(respval<0.05)) #0



res<-empiricall_fdr(databotcombined_fdr,c(combined_fdr_l[[1]],combined_fdr_l[[2]],combined_fdr_l[[3]],combined_fdr_l[[4]],combined_fdr_l[[5]]),5)
respval<-empiricall_fdr(databotcombined_pval,c(combined_pval_l[[1]],combined_pval_l[[2]],combined_pval_l[[3]],combined_pval_l[[4]],combined_pval_l[[5]]),5)

res<-res[order(res)]

length(res)
length(res[res<0.05])
plot(res,type="l",xlab="index",ylab="fdr")
abline(0.05,0,col="firebrick",lty=2)

databot_sign<-databot[1:length(res[res<0.05]),]
dim(databot_sign)

FIN2A<-system("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINA.bed | sed 's/chr//' | awk -v OFS='.' '{print $1,$3}'", intern=TRUE)
FIN2B<-system("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINB.bed | sed 's/chr//' | awk -v OFS='.' '{print $1,$3}'", intern=TRUE)
mylinks<-sapply(1:length(FIN2A),function(x) paste0(FIN2A[x],".",FIN2B[x]))
#FIN2B<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINB.bed")

head(cbind(FIN2A,FIN2B))
dim(cbind(FIN2A,FIN2B))

length(unique(c(FIN2A,FIN2B)))
nsnpsFIN2<-as.numeric(system("cat botnia/FIN2/*tab |wc -l",intern=TRUE))-231
nsnpsFIN2

dataLDLD<-databot_sign
dataLDLDA<-as.character(paste0(dataLDLD[,1],".",dataLDLD[,3],".",dataLDLD[,2],".",dataLDLD[,4]))
dataLDLD_overlap<-dataLDLD[!is.na(match(dataLDLDA,mylinks)),]

mylinks[order(mylinks)][1:10]

dbinom(prob=(nsnpsFIN2/2)/ncomp,size=length(dataLDLDA),x=0)
dbinom(prob=(nsnpsFIN2/2)/ncomp,size=length(dataLDLDA),x=1)
dbinom(prob=(nsnpsFIN2/2)/ncomp,size=length(dataLDLDA),x=2)
dbinom(prob=(nsnpsFIN2/2)/ncomp,size=length(dataLDLDA),x=3)
dbinom(prob=(nsnpsFIN2/2)/ncomp,size=length(dataLDLDA),x=4)

overlap is not significant

# TRIOS

There is always very little overlap with the scan, even when using trios information in order to avoid substructure.

see ~/Dropbox/LDLD/botnia/epistatictrasmission.R

setwd("~/Dropbox/LDLD/botnia")
load("restrios_botnia.RData")

pos1<-restrios[,1][order(restrios[,7])[1:2000]]
pos2<-restrios[,2][order(restrios[,7])[1:2000]]
top5klinkstrios<-sapply(1:2000,function(x) paste0(pos1[x],".",pos2[x]))
top5klinkstrios[1:10]

dataLDLD<-top5klinkstrios
dataLDLD_overlap<-dataLDLD[!is.na(match(dataLDLD,mylinks))]
dataLDLD_overlap

par(mfrow=c(2,2))
hist(as.numeric(restrios[,3]),main="all")
hist(as.numeric(restrios[,4]),main="sum")
hist(as.numeric(restrios[,5]),main="mother vs child")
hist(as.numeric(restrios[,6]),main="father vs child")

par(mfrow=c(1,1))
hist(as.numeric(restrios[,7]),main="father vs mother")

par(mfrow=c(2,2))
hist(as.numeric(restrios[,3])[as.numeric(restrios[,3])<0.1],main="all")
hist(as.numeric(restrios[,4])[as.numeric(restrios[,4])<0.1],main="sum")
hist(as.numeric(restrios[,5])[as.numeric(restrios[,5])<0.1],main="mother vs child")
hist(as.numeric(restrios[,6])[as.numeric(restrios[,6])<0.1],main="father vs child")

par(mfrow=c(1,1))
hist(as.numeric(restrios[,7])[as.numeric(restrios[,7])<0.1],main="father vs mother")