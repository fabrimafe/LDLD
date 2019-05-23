#==========reanalyses of data starting on 15-3-2017================
#======================LDLD shell calculations==================
{
#see LDLD_1000g_preparation_files.sh
}
#====================R initial variables collections=======================
{
if (TRUE) {
mynamefreq<-5
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
computed_exact_pvalues<-TRUE
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
npopulations<-12
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE
nperm=1
myrepls<-"merged" #myrepls<-3:5
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/"
freqfiltering<-mynamefreq/100
myrepl<-myrepls[1]
myfolder<-paste0(myfolder_root,myrepl,"/min",mynamefreq)
myfolderdropbox<-paste0(myfolderdropbox_root,myrepl,"/min",mynamefreq)
}

if (TRUE) {
mynamefreq<-1
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
computed_exact_pvalues<-TRUE
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
npopulations<-12
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE
nperm=1
myrepls<-"merged" #myrepls<-3:5
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/"
freqfiltering<-mynamefreq/100
myrepl<-myrepls[1]
myfolder<-paste0(myfolder_root,myrepl,"/min",mynamefreq)
myfolderdropbox<-paste0(myfolderdropbox_root,myrepl,"/min",mynamefreq)
}


if (TRUE) {
mynamefreq<-5
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
computed_exact_pvalues<-TRUE
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
npopulations<-12
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE
nperm=1
myrepls<-"merged" #myrepls<-3:5
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/introns.1000g/"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/introns.1000g/"
freqfiltering<-mynamefreq/100
myrepl<-myrepls[1]
myfolder<-paste0(myfolder_root,myrepl,"/min",mynamefreq)
myfolderdropbox<-paste0(myfolderdropbox_root,myrepl,"/min",mynamefreq)
}


if (TRUE) {
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/introns.1000g/repl"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/introns.1000g/repl"
myrepls<-1
npopulations<-12
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
PREPROCESSING<-FALSE
computed_exact_pvalues<-TRUE
using_sign_for_nAB<-TRUE
nperm=1
mynamefreq<-5
freqfiltering<-mynamefreq/100
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
}

if (TRUE) {
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/repl"
myrepls<-1:3
npopulations<-12
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
computed_exact_pvalues<-TRUE
PREPROCESSING<-TRUE
using_sign_for_nAB<-TRUE
nperm=1
mynamefreq<-1
freqfiltering<-mynamefreq/100
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
}


if (TRUE) {
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/repl"
myrepls<-1:10
npopulations<-12
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
computed_exact_pvalues<-TRUE
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE
nperm=1
mynamefreq<-5
freqfiltering<-mynamefreq/100
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
}

if (TRUE) {
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl"
myrepls<-2
npopulations<-12
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
computed_exact_pvalues<-TRUE
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE
nperm=1
mynamefreq<-1
freqfiltering<-mynamefreq/100
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
}

if (TRUE) {
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/repl"
myrepls<-1:3
npopulations<-12
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
computed_exact_pvalues<-TRUE
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE
nperm=1
mynamefreq<-1
freqfiltering<-mynamefreq/100
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
}


setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
source("~/Dropbox/general_utils/general_functions.R")
options(scipen=999)
setwd("/mnt/scratch/fabrizio/LDLD")
library(data.table)
library(lme4)
for ( myrepl in myrepls){
myfolder<-paste0(myfolder_root,myrepl,"/min",mynamefreq)
myfolderdropbox<-paste0(myfolderdropbox_root,myrepl,"/min",mynamefreq)
}

#===================SIMULATE FALSE POSITIVES==================================================
#mytab<-fread(paste0("cat ",myfolder,"/reschr1/chr",i,".tab | sed 's/|//g' "))


}
#=============================================================================================
#=============STEP0: PREPROCESSING and IMPORTING==============================================
#=============================================================================================
{

for ( myrepl in myrepls){
myfolder<-paste0(myfolder_root,myrepl,"/min",mynamefreq)
myfolderdropbox<-paste0(myfolderdropbox_root,myrepl,"/min",mynamefreq)

#======================LOAD_INITIAL_FILES====================================================
system(paste0("mkdir -p ",myfolderdropbox))
#TAG 'LOAD_INITIAL_FILES'
if (TRUE)
{
if (PREPROCESSING) {
    print("aggregating res files")
    system(paste0("mkdir -p ",myfolder,"/anal/"))
    system(paste0("~/Dropbox/LDLD/scripts/postproc_split.sh ",myfolder))
    print("aggregating log files for each pops")
    print("aggregating all log files")
    organize_log_files(myfolder,npopulations) #generates nAB files i.e. ./logs/all.popX.log.gz
    system(paste0("cat ",myfolder,"/*minilog | grep ncomparisonspair  > ",myfolder,"/ncomparisons_pairwise.txt"))
    system(paste0("cat ",myfolder,"/*minilog | grep ncomparisons | grep -v pair | awk -v npops=", npopulations ,"  '{for (i=2;i<=(npops+2);i++){printf \"%d\\t\", $i }; printf \"\\n\"}' > ",myfolder,"/ncomparisons.txt"))
    for (myperm in 1:nperm)
        {
        print(paste("permutation: ",myperm))
        system(paste0("cat ",myfolder,"/reschr",myperm,"/*minilog | grep ncomparisons | grep -v pair | awk '{print $14}' > ",myfolder,"/reschr",myperm,"/ncomparisons.txt"))
        print("aggregating res files for permutations")
        system(paste0("mkdir -p ",myfolder,"/reschr",myperm,"/logs"))
        system(paste0("mkdir -p ",myfolder,"/reschr",myperm,"/anal"))
        system(paste0("~/Dropbox/LDLD/scripts/postproc_split.sh ",myfolder,"/reschr",myperm))
        print("aggregating log files for permutations")
        organize_log_files(paste0(myfolder,"/reschr",myperm),12)
        }
    print("generating mutation files")
    setwd("~/Dropbox/LDLD/scripts/")
    system("gcc ~/Dropbox/LDLD/scripts/LDLD_filter5.c -lgmp -lm -fno-stack-protector -o ~/Dropbox/LDLD/scripts/filterbyfreq.out")
    system(paste0("rm ",myfolder,"/logs/all.freqlog"))
    for (i in 1:22)
        {
        system(paste0("nohup ~/Dropbox/LDLD/scripts/filterbyfreq.out ",myfolder,"/chr",i,".tab ",myfolder,"/chr",i,".freqlog ",myfolder,"/chr",i,".mutlog ", freqfiltering," 1220"))
        system(paste0("nohup ~/Dropbox/LDLD/scripts/filterbyfreq_multi_sharing.out ",myfolder,"/chr",i,".tab nspops.txt ",freqfiltering," ",npopulations," 1 > ",myfolder,"/chr",i,".myfreq.sharing"))
        system(paste0("cat ",myfolder,"/chr",i,".freqlog | awk '{printf \"%d\\t\",",i,"; print $0}' >> ",myfolder,"/logs/all.freqlog"))
        #for (j in 1:nperm)
        #   {
        #   system(paste0("nohup ~/Dropbox/LDLD/scripts/filterbyfreq.out ",myfolder,"/reschr",j,"/chr",i,".tab ",myfolder,"/reschr",j,"/chr",i,".freqlog ",myfolder,"/reschr",j,"/chr",i,".mutlog ", freqfiltering," 1220"))
        #   }
        }
    system(paste0("gzip -f ",myfolder,"/logs/all.freqlog"))
    setwd("/mnt/scratch/fabrizio/LDLD")
    system(paste0("cat ",myfolder,"/chr*.myfreq.sharing | grep -v [an] | awk '{for (i=1;i<=12;i++){if (i<=NF){printf \"%d\\t\",$i} else {printf  \"-1\\t\"};};printf \"\\n\"}' > ",myfolder,"/myfreq.sharing"))
    system(paste0("mkdir ",myfolderdropbox))
}


#combining p-values theory
if (FALSE)
{
#Stouffer's method
#Fisher's exact test is also good (actually it is naturally) as 1-sided p-value test.
#However Kulinskaya correction makes it 2 sided.
#Notice I can always correct fdr independently for negative and positive using empirical null. However combined p-value is both positive and negative.
#One way would be to do fisher's exact test without Kuli, so simply at values of .: http://mathworld.wolfram.com/FishersExactTest.html
#Now, comparison with null distribution to get a false discovery rate can be done also with any statistics. However in that case there would be the underlying assumption that significant.
#notice that fisher's exact method and stouffer's give basically the same results. So I keep method as it is and work later in R.

library(metap)
sumz(c(0.05,0.05)) #Stouffer's
sumlog(c(0.05,0.05)) #Fisher's #more conservative
z_value<-function(pvalue) {
#qnorm(1-pvalue/2)
qnorm(1-pvalue)
}
mypvalues<-c(0.05,0.05)
sum(sapply(mypvalues,function(x) z_value(x)))
qnorm(pnorm(0.5)) #0
1-pnorm(sum(sapply(mypvalues,function(x) z_value(x)))/sqrt(2))

}
mydata<-fread(paste0("zcat ",myfolder,"/anal/sorted.res.gz"))
names23<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","pFisherpos","pFisherneg","pKuli","T2","Xtot","X","pop","popX")
#names23<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
names20<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","T2","Xtot","X","pop","popX")
if (computed_exact_pvalues) { mynames<-names23 } else { mynames<-names20 }
setnames(mydata,mynames)
#save(mydata,file=paste0(myfolder,"/mydata.RData"))
ncomp<-read.table(paste0(myfolder,"/ncomparisons.txt"),header=FALSE)
ncomp<-sum(as.numeric(ncomp$V1))

#T2 combined p_values as calculated by LDLD
{
dataFIN<-subset(mydata,mydata$pop==index_focal_pop)
combined_pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)
combined_fdr<-p.adjust(p=combined_pval, method = "fdr", n = ncomp)
length(combined_fdr)
dim(dataFIN)
dim(mydata)
fullsign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
save(dataFIN,file=paste0(myfolder,"/dataFIN.RData"))
save(fullsign,file=paste0(myfolder,"/fullsign.RData"))
save(combined_pval,file=paste0(myfolder,"/combined_pval.RData"))
highlysign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.00001,]
save(combined_fdr,file=paste0(myfolder,"/combined_fdr.RData"))
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
write.table(signAbed,file=paste0(myfolder,"/anal/snpsA.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file=paste0(myfolder,"/anal/snpsB.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
}
#directional p_values
{
mydata$T2[mydata$T2==0]<-10^(-21)#min(mydata$pFisherpos[mydata$pFisherpos!=0])
mydata$T2[mydata$T2>1]<-1#min(mydata$pFisherneg[mydata$pFisherneg!=0])
mydata$pKuli[mydata$pKuli==0]<-10^(-21)#min(mydata$pFisherpos[mydata$pFisherpos!=0])
mydata$pFisherpos[mydata$pFisherpos==0]<-10^(-21)#min(mydata$pFisherpos[mydata$pFisherpos!=0])
mydata$pFisherneg[mydata$pFisherneg==0]<-10^(-21)#min(mydata$pFisherneg[mydata$pFisherneg!=0])
mydata$pKuli[mydata$pKuli>1]<-1#min(mydata$pFisherpos[mydata$pFisherpos!=0])
mydata$pFisherpos[mydata$pFisherpos>1]<-1#min(mydata$pFisherpos[mydata$pFisherpos!=0])
mydata$pFisherneg[mydata$pFisherneg>1]<-1#min(mydata$pFisherneg[mydata$pFisherneg!=0])
temp<-mydata[mydata$nA/mydata$Nse>=0.01 & mydata$nB/mydata$Nse>=0.01 & mydata$nA/mydata$Nse<=0.99 & mydata$nB/mydata$Nse<=0.99,]
save(temp,file=paste0(myfolder,"/temp.RData"))
rm(mydata)
#combined pvalues for positive and negative
require(metap)
fishersmethod<-function(z) if ( length(z)>=2) {return(sumlog(z)$p)} else return(z)
pairid<-temp[, paste(chrA,chrB,posA,posB,sep=".")]
combined_pval_T2<-aggregate(temp$T2,by=list(pairid),FUN=function(x) fishersmethod(x))
combined_fdr_T2<-p.adjust(p=combined_pval_T2$x, method = "fdr", n = ncomp)
if (computed_exact_pvalues) { 
combined_pval_pos<-aggregate(temp$pFisherpos,by=list(pairid),FUN=function(x) fishersmethod(x))
combined_pval_neg<-aggregate(temp$pFisherneg,by=list(pairid),FUN=function(x) fishersmethod(x))
combined_pval<-aggregate(temp$pKuli,by=list(pairid),FUN=function(x) fishersmethod(x))
combined_fdr_pos<-p.adjust(p=combined_pval_pos$x, method = "fdr", n = ncomp)
combined_fdr_neg<-p.adjust(p=combined_pval_neg$x, method = "fdr", n = ncomp)
combined_fdr<-p.adjust(p=combined_pval$x, method = "fdr", n = ncomp)
combined_pvals<-list(combined_pval,combined_pval_pos,combined_pval_neg,combined_pval_T2,combined_fdr,combined_fdr_pos,combined_fdr_neg,combined_fdr_T2)
names(combined_pvals)<-c("combined_pval","combined_pval_pos","combined_pval_neg","combined_pval_T2","combined_fdr","combined_fdr_pos","combined_fdr_neg","combined_fdr_T2")
} else {
combined_pvals<-list(combined_pval_T2,combined_fdr_T2)
names(combined_pvals)<-c("combined_pval","combined_fdr")
}
save(combined_pvals,file=paste0(myfolder,"/combined_pvals_theoretical.RData"))


}
}
#================================================================================================^
#===============IDENTIFYING SIGNIFICANT PAIRS BASED ON THEORETICAL P-VALUES======================
#================================================================================================
if(FALSE)
{
#Different ways of obtaining list of significant pairs
#T2
combined_T2pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)
combined_T2fdr<-p.adjust(p=combined_T2pval, method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
mydata<-mydata[mydata$nA>=0.05*mydata$Nse,]
mydata<-mydata[mydata$nA<=0.95*mydata$Nse,]
mydata<-mydata[mydata$nB>=0.05*mydata$Nse,]
mydata<-mydata[mydata$nB<=0.95*mydata$Nse,]
range(combined_T2pval) #0.000000e+00 2.496462e-06


#Below only do it if also exact p-values are computed
if (computed_exact_pvalues)
{
pairid<-mydata[, paste(chrA,chrB,posA,posB,sep=".")]
mydata[,pairid:=pairid]
dim(mydata)
length(pairid)
combined_prho2pval<-aggregate(mydata$prho2,by=list(pairid),FUN=function(x) 1-pchisq(-2*log(prod(x)),df=2*length(x)))
combined_pKulipval<-aggregate(mydata$pKuli,by=list(pairid),FUN=function(x) 1-pchisq(-2*log(prod(x)),df=2*length(x)))
combined_prho2fdr<-p.adjust(p=combined_prho2pval[,2], method = "fdr", n = ncomp) 
combined_pKulifdr<-p.adjust(p=combined_pKulipval[,2], method = "fdr", n = ncomp) 

mypvalues<-cbind(combined_prho2pval,combined_pKulipval[,2],combined_prho2fdr,combined_pKulifdr)
sum(as.numeric(combined_prho2fdr<0.05))
sum(as.numeric(combined_pKulifdr<0.05))
sum(as.numeric(combined_T2fdr<0.05))

#save(fullsign,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/fullsign.RData")
names(mypvalues)<-c("pairid","combined_prho2","combined_pKuli","combined_prho2_fdr","combined_pKuli_fdr")
mydata<-merge(mydata,mypvalues,by="pairid")
save(mydata,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/mypairs.RData")

load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/mypairs.RData")

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
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/LDLDpvalues_distr_signfdrfromT2_new.pdf")
mydata2<-data.frame(values=c(mydata$pKuli,mydata$prho2,mydata$T2),groups=c(rep("pKuli",length(mydata$pKuli)),rep("prho2",length(mydata$prho2)),rep("T2",length(mydata$prho2))))
ggplot(mydata2, aes(x=values, fill=groups)) +
    geom_histogram(binwidth=.05, position="dodge")+theme_bw()
dev.off()
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/LDLDpvalues_distr_zoom_signfdrfromT2_new.pdf")
mythr<-0.0001
mydata2<-data.frame(values=c(mydata$pKuli[mydata$pKuli<mythr],mydata$prho2[mydata$prho2<mythr],mydata$T2[mydata$T2<mythr]),groups=c(rep("pKuli",length(mydata$pKuli[mydata$pKuli<mythr])),rep("prho2",length(mydata$prho2[mydata$prho2<mythr])),rep("T2",length(mydata$T2[mydata$T2<mythr]))))
ggplot(mydata2, aes(x=values, fill=groups)) +
    geom_histogram(binwidth=mythr/30, position="dodge")+theme_bw()+xlim(0,0.0001)
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


#==============================================================================================================v
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
ncomppair

mydata<-data2[data2$nA>=0.05*data2$Nse,]
mydata<-mydata[mydata$nA<=0.95*mydata$Nse,]
mydata<-mydata[mydata$nB>=0.05*mydata$Nse,]
mydata<-mydata[mydata$nB<=0.95*mydata$Nse,]

#I could ask: given that significant in one population, what's the chance of being sign in others
#save(mydata,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/mypairs.RData")
#load(mydata,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/mypairs.RData")

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
save("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/mymat")

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
  mynooverlap<-sum(ncomp[,13+1]) - length(union(tempi$pairid,tempj$pairid))
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


data2_pop_T2fdr_sign[[1]]
save("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/mymat")
}
}
#================================================================================================^
#===============IDENTIFYING SIGNIFICANT PAIRS BASED ON EMPIRICAL DISTRIBUTION====================
#================================================================================================
if (TRUE)
{
#compute fdr
{
#load(paste0(myfolder,"/combined_pvals.RData"))
#load(paste0(myfolder,"/fullsign.RData"))
#combined_pval_sign<-combined_pvals
#combined_fdr_sign<-fullsign$combined_fdr
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
nperm<-1
reschr_count<-rep(0,nperm)
combined_pval_l<-list();dataFIN_l<-list();pop_pval_l<-list();pop_pvalpos_l<-list();pop_pvalneg_l<-list()
for (myperm in 1:nperm)
    {
    print(myperm)
    print("import")
    mydata<-fread(paste0("zcat ",myfolder,"/reschr",myperm,"/anal/sorted.res.gz"))
#    data<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/reschr",myperm,"/anal/sorted.res"))
    names23<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","pFisherpos","pFisherneg","pKuli","T2","Xtot","X","pop","popX")
    names20<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","T2","Xtot","X","pop","popX")
    if (computed_exact_pvalues) { mynames<-names23 } else { mynames<-names20 }
    #names(data)<-mynames #if data.frame
    setnames(mydata,mynames)
    #if (dim(mydata)[2]==20){names(mydata)<-names20} else if (dim(mydata)[2]==23){names(mydata)<-names23}
    ncomp<-read.table(paste0(myfolder,"/reschr",myperm,"/ncomparisons.txt"),header=FALSE)
    ncomp<-sum(as.numeric(ncomp$V1))
    dataFIN<-subset(mydata,mydata$pop==index_focal_pop)
    combined_pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)    
    combined_fdr<-p.adjust(p=combined_pval, method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
    fullsign_reschr<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
    combined_pval_l[[myperm]]<-combined_pval
    dataFIN_l[[myperm]]<-dataFIN
    reschr_count[myperm]<-dim(fullsign_reschr)[1]
        #directional p_values
        {
        print("calculate pvalues")
        mydata$T2[mydata$T2==0]<-10^(-21)#min(mydata$pFisherpos[mydata$pFisherpos!=0])
        mydata$T2[mydata$T2>1]<-1#min(mydata$pFisherneg[mydata$pFisherneg!=0])
        mydata$pKuli[mydata$pKuli==0]<-10^(-21)#min(mydata$pFisherpos[mydata$pFisherpos!=0])
        mydata$pFisherpos[mydata$pFisherpos==0]<-10^(-21)#min(mydata$pFisherpos[mydata$pFisherpos!=0])
        mydata$pFisherneg[mydata$pFisherneg==0]<-10^(-21)#min(mydata$pFisherneg[mydata$pFisherneg!=0])
        mydata$pKuli[mydata$pKuli>1]<-1#min(mydata$pFisherpos[mydata$pFisherpos!=0])
        mydata$pFisherpos[mydata$pFisherpos>1]<-1#min(mydata$pFisherpos[mydata$pFisherpos!=0])
        mydata$pFisherneg[mydata$pFisherneg>1]<-1#min(mydata$pFisherneg[mydata$pFisherneg!=0])
        temp<-mydata[mydata$nA/mydata$Nse>=0.01 & mydata$nB/mydata$Nse>=0.01 & mydata$nA/mydata$Nse<=0.99 & mydata$nB/mydata$Nse<=0.99,]
        rm(mydata)
        #combined pvalues for positive and negative
        fishersmethod<-function(z) if ( length(z)>=2) {return(sumlog(z)$p)} else return(z)
        pairid<-temp[, paste(chrA,chrB,posA,posB,sep=".")]
        combined_pval_T2<-aggregate(temp$T2,by=list(pairid),FUN=function(x) fishersmethod(x))
        combined_fdr_T2<-p.adjust(p=combined_pval_T2$x, method = "fdr", n = ncomp)
        if (computed_exact_pvalues) { 
        combined_pval_pos<-aggregate(temp$pFisherpos,by=list(pairid),FUN=function(x) fishersmethod(x))
        combined_pval_neg<-aggregate(temp$pFisherneg,by=list(pairid),FUN=function(x) fishersmethod(x))
        combined_pval<-aggregate(temp$pKuli,by=list(pairid),FUN=function(x) fishersmethod(x))
        combined_fdr_pos<-p.adjust(p=combined_pval_pos$x, method = "fdr", n = ncomp)
        combined_fdr_neg<-p.adjust(p=combined_pval_neg$x, method = "fdr", n = ncomp)
        combined_fdr<-p.adjust(p=combined_pval$x, method = "fdr", n = ncomp)
        combined_pvals<-list(combined_pval,combined_pval_pos,combined_pval_neg,combined_pval_T2,combined_fdr,combined_fdr_pos,combined_fdr_neg,combined_fdr_T2)
        pop_pval_l[[myperm]]<-lapply(0:(npopulations-1),function(x) temp[temp$pop==x]$pKuli)
        pop_pvalpos_l[[myperm]]<-lapply(0:(npopulations-1),function(x) temp[temp$pop==x]$pFisherpos)
        pop_pvalneg_l[[myperm]]<-lapply(0:(npopulations-1),function(x) temp[temp$pop==x]$pFisherneg)
        names(combined_pvals)<-c("combined_pval","combined_pval_pos","combined_pval_neg","combined_pval_T2","combined_fdr","combined_fdr_pos","combined_fdr_neg","combined_fdr_T2")
        } else {
        combined_pvals<-list(combined_pval_T2,combined_fdr_T2)
        names(combined_pvals)<-c("combined_pval","combined_fdr")
        }
        combined_pval_l[[myperm]]<-combined_pvals
        }
    }
#    save(combined_pval_l,file=paste0(myfolder,"/combined_pval_l.RData"))
#    save(dataFIN_l,file=paste0(myfolder,"/dataFIN_l.RData"))
    save(combined_pval_l,file=paste0(myfolder,"/combined_pval_l.RData"))
    save(pop_pval_l,file=paste0(myfolder,"/pop_pval_l.RData"))
    load(paste0(myfolder,"/combined_pvals_theoretical.RData"))
    empfdr_T2<-empiricall_fdr(combined_pvals[["combined_pval_T2"]]$x,unlist(sapply(combined_pval_l, function(x) x[["combined_pval_T2"]]$x)),nperm)
    combined_pvals[["combined_fdr_T2"]]<-empfdr_T2
    if (computed_exact_pvalues) {     
    empfdr_pos<-empiricall_fdr(combined_pvals[["combined_pval_pos"]]$x,unlist(sapply(combined_pval_l, function(x) x[["combined_pval_pos"]]$x)),nperm)
    empfdr_neg<-empiricall_fdr(combined_pvals[["combined_pval_neg"]]$x,unlist(sapply(combined_pval_l, function(x) x[["combined_pval_neg"]]$x)),nperm)
    empfdr<-empiricall_fdr(combined_pvals[["combined_pval"]]$x,unlist(sapply(combined_pval_l, function(x) x[["combined_pval"]]$x)),nperm)
    #pop_empfdr<-empiricall_fdr(combined_pvals[["combined_pval_neg"]]$x,unlist(sapply(combined_pval_l, function(x) x[["combined_pval_neg"]]$x)),nperm)
    load(paste0(myfolder,"/temp.RData"))
    pop_fdr_l<-list();pop_fdrpos_l<-list();pop_fdrneg_l<-list();
    for (ipop in 1:npopulations) { pop_fdr_l[[ipop]]<-empiricall_fdr(temp$pKuli[temp$pop==(ipop-1)],c(unlist(sapply(1:nperm,function(x) pop_pval_l[[x]][[ipop]]))),nperm) }
    for (ipop in 1:npopulations) { pop_fdrpos_l[[ipop]]<-empiricall_fdr(temp$pFisherpos[temp$pop==(ipop-1)],c(unlist(sapply(1:nperm,function(x) pop_pvalpos_l[[x]][[ipop]]))),nperm) }
    for (ipop in 1:npopulations) { pop_fdrneg_l[[ipop]]<-empiricall_fdr(temp$pFisherneg[temp$pop==(ipop-1)],c(unlist(sapply(1:nperm,function(x) pop_pvalneg_l[[x]][[ipop]]))),nperm) }
    combined_pvals[["combined_fdr"]]<-empfdr
    combined_pvals[["combined_fdr_pos"]]<-empfdr_pos
    combined_pvals[["combined_fdr_neg"]]<-empfdr_neg
    combined_pvals[["pop_fdr"]]<-pop_fdr_l
    combined_pvals[["poppos_fdr"]]<-pop_fdrpos_l
    combined_pvals[["popneg_fdr"]]<-pop_fdrneg_l
}
    save(combined_pvals,file=paste0(myfolder,"/combined_pvals.RData"))
    load(paste0(myfolder,"/combined_pvals.RData"))

if (FALSE)
{
    sum(pop_pval_l[[1]][[ipop]]<0.05)
    mybars<-rbind(sapply(1:npopulations,function(ipop) sum(temp$pFisherpos[temp$pop==(ipop-1)]<10^(-6))),    sapply(1:npopulations,function(ipop) sum(pop_pvalpos_l[[1]][[ipop]]<10^(-6))))
    pdf(paste0(myfolderdropbox,"/barplot_sign_vs_null_106.pdf"))
    barplot(mybars,names.arg=popnames,beside=TRUE,cex.names=0.8,col=c("cadetblue","gold"),ylab="n.links",space=c(0,0.01))
    dev.off()
    pdf(paste0(myfolderdropbox,"/barplot_sign_vs_null_logfold_106.pdf"))
    barplot(mybars[1,]/mybars[2,],names.arg=popnames,cex.names=0.8,col=c("cadetblue"),ylab="fold enrichment")
    dev.off()
    mybars<-sapply(1:npopulations,function(ipop) sum(combined_pvals$poppos_fdr[[ipop]]<0.05))
    pdf(paste0(myfolderdropbox,"/barplot_sign_links.pdf"))
    x<-barplot(mybars,names.arg=popnames,col=c("cadetblue"),ylab="n.links",xaxt="n",cex.names=1.5,cex.axis=1.5)
    text(cex=1.5, x=x-.25, y=-max(mybars)/15, popnames, xpd=TRUE, srt=45)
    dev.off()
}
    
if (FALSE){
    length(combined_pvals$combined_pval_pos$Group.1[combined_pvals$combined_fdr_pos<0.05])
    length(combined_pvals$combined_pval_pos$Group.1)
    
    
    for (ipop in 1:npopulations) { pop_fdr_l[[ipop]]<-empiricall_fdr(temp$pKuli[temp$pop==(ipop-1)],c(unlist(sapply(1:nperm,function(x) pop_pval_l[[x]][[ipop]]))),nperm) }
    for (ipop in 1:npopulations) { pop_fdrpos_l[[ipop]]<-empiricall_fdr(temp$pFisherpos[temp$pop==(ipop-1)],c(unlist(sapply(1:nperm,function(x) pop_pvalpos_l[[x]][[ipop]]))),nperm) }
    for (ipop in 1:npopulations) { pop_fdrneg_l[[ipop]]<-empiricall_fdr(temp$pFisherneg[temp$pop==(ipop-1)],c(unlist(sapply(1:nperm,function(x) pop_pvalneg_l[[x]][[ipop]]))),nperm) }
    }
    

}
}
}

#check individuals repl
if (FALSE)
{
#I notice diff in coding repl1 and repl2, where they should be the same: more sign in repl2, with more discoveries.
#from figure barplot_logabsolutenAB.pdf it seems I used to do better. Also I get so much more significant (this is why probably I separate less, more false positive. Why? Starting file is smaller. Are reshuffled ok?
#do individual pairs get the same pvalue?
#look at in repl4/min5/chr1.tabchr13.tab.res line3 : 915490  20641422        rs138027487     rs170758. Yes, linkage values are the same. #both positions are present in repl1 vcf files but not shown in res!!
#ok, let's start from repl1 since smaller. First link in repl1/min5/chr1.tabchr13.tab.res #linkage values are the same
#this at least means that it is not fucked up. Good.
#repl2 data is overall larger, although starting files are not. This seems unlikely coming from threshold pvalue (although it might). 
#Frequency threshold used in the past could have been files with min1, or higher filtering for total freq. If that s the case, I would have more rare alleles/alleles in just one pop, that could be more easily significant (also by chance). 
#This would explain everything, although in theory it should be controlled by reshuffling. Then there is the reshuffling. Maybe I should have more replicates for coding. To exclude that reshuffling let's calculate theoretical fdr.
#what can I look at now?
#links that are present in one and not the other, what features do they have?
#cool, the problem it is solved: the difference is truly in LDLD: now I also output cases in which last cols in res is 1.
#what does that mean? That I also print cases in which significant in only 1 pop but not overall. I feel these cases should be excluded, or not?
#actually yes. However, I used to avoid printinf also those that significant in only 1 pop, but that is the only one for which above freq threshold. These should be kept. So ok as I do now.With the only thing that among them I probably have
#more false positive.

myrepls=1:10
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/repl"

load(paste0(myfolder,"/combined_pvals.RData"))
sum(combined_pvals$combined_fdr_pos<0.05)
sum(combined_pvals$combined_pval_pos<0.00005)
sum(combined_pvals$combined_pval_T2<0.00005)
load(paste0(myfolder,"/temp.RData"))
dim(temp) #more significant in 2, but temp smaller
mydata<-fread(paste0("zcat ",myfolder,"/anal/sorted.res.gz")) #data are definitely bigger in repl1
all.allele.freqsA<-dim(matrix(mydata$V7,ncol=12,byrow=T)/matrix(mydata$V9,ncol=12,byrow=T))

#
myrepls=1:5
myfolder_root<-"/mnt/scratch/fabrizio/LDLD/above95/introns.1000g/repl"
myfolderdropbox_root<-"~/Dropbox/LDLD/ipynb/figs/introns.1000g/repl"
for ( myrepl in myrepls){
myfolder<-paste0(myfolder_root,myrepl,"/min",mynamefreq)
myfolderdropbox<-paste0(myfolderdropbox_root,myrepl,"/min",mynamefreq)
load(paste0(myfolder,"/combined_pvals.RData"))
print(c(sum(combined_pvals$combined_fdr_pos<0.05),sum(combined_pvals$combined_fdr_neg<0.05)))
}

}

}

#================================================================================================^
#======================PLOT SIGNIFICANT PAIRS BASED ON EMPIRICAL DISTRIBUTION====================
#================================================================================================
{
#plots of fdr
{
#        load(paste0(myfolder,"/combined_pvals_theoretical.RData"))
#        combined_pvals_theoretical<-combined_pvals
        load(paste0(myfolder,"/combined_pvals.RData"))
        load(paste0(myfolder,"/combined_pval_l.RData"))
        #save(empfdr,file=paste0(myfolder,"/empfdr.RData"))
        load(paste0(myfolder,"/combined_fdr.RData"))
        if (FALSE) {
        plot(empfdr,type="l",xlab="index",ylab="fdr")
        pdf(paste0(myfolderdropbox,"/empfdrvsT2.pdf"))
        plot(combined_fdr[order(combined_fdr)],type="l",xlab="index",ylab="fdr",ylim=c(0,0.05),xlim=c(1,150000),col="cadetblue",lwd=2.5)
        lines(empfdr[order(empfdr)],col="gold",lwd=2.5)
        abline(h=0.05,col="black",lty=2,lwd=2.5)
        legend(0,0.03,c("T2 fdr","T2 empirical"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        lty=c(1,1), # gives the legend appropriate symbols (lines)
        lwd=c(2.5,2.5),col=c("cadetblue","gold"))
        dev.off()
        #NB: empirical curves obtained with #new (high res scan) and the lower resolution one are exactly the same, just that now I can explore higher fdr (that is useless anyway, apart for this plot)
        }
if (computed_exact_pvalues){
        mybreaks<-unlist(c(0,sapply(20:5,function(x) 10^(-x))))
        myticks<-sort(unique(c(mybreaks,mybreaks/2)))
        mylabs<-c(0.00000000000000000001,0.00000000000000001000,0.00000000000001000000,0.00000000001000000000,0.00000001000000000000,0.00001000000000000000)
        myxlabs<-list(pos=which(myticks %in% (mylabs/2)),lab=c(expression(paste('<',10^'-20')),expression(paste(10^'-17')),expression(paste(10^'-14')),expression(paste(10^'-11')),expression(paste(10^'-8')),expression(paste(10^'-5'))))
        pdf(paste0(myfolderdropbox,"/empfdr_hist.pdf"))
        myhistdata<-myhist_f(mybreaks,combined_pvals[["combined_pval"]]$x[combined_pvals[["combined_pval"]]$x<10^(-5)])
        myhistdata_reschr<-myhist_f(mybreaks,combined_pval_l[[1]][["combined_pval"]]$x[combined_pval_l[[1]][["combined_pval"]]$x<10^(-5)])
        myhistdata<-c(rbind(myhistdata,myhistdata_reschr))
        mybarplot_f(1:length(myticks),log(myhistdata+1,10),myylim=0,mycols=rep(c("gold","cadetblue"),length(mybreaks)),ticks_are_breaks=TRUE,myxlabel='p.value',myylabel=expression(paste(log["10"],"(count+1)")),add=FALSE,xaxis=myxlabs)
        legend("top",c("observed","reshuffled"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        pch=15,col=c("gold","cadetblue"))
        dev.off()
        pdf(paste0(myfolderdropbox,"/empfdr_pos_hist.pdf"))
        myhistdata<-myhist_f(mybreaks,combined_pvals[["combined_pval_pos"]]$x[combined_pvals[["combined_pval_pos"]]$x<10^(-5)])
        myhistdata_reschr<-myhist_f(mybreaks,combined_pval_l[[1]][["combined_pval_pos"]]$x[combined_pval_l[[1]][["combined_pval_pos"]]$x<10^(-5)])
        myhistdata<-c(rbind(myhistdata,myhistdata_reschr))
        mybarplot_f(1:length(myticks),log(myhistdata+1,10),myylim=0,mycols=rep(c("gold","cadetblue"),length(mybreaks)),ticks_are_breaks=TRUE,myxlabel='p.value',myylabel=expression(paste(log["10"],"(count+1)")),add=FALSE,xaxis=myxlabs)
        legend("top",c("observed","reshuffled"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        pch=15,col=c("gold","cadetblue"))
        dev.off()
        pdf(paste0(myfolderdropbox,"/empfdr_neg_hist.pdf"))
        myhistdata<-myhist_f(mybreaks,combined_pvals[["combined_pval_neg"]]$x[combined_pvals[["combined_pval_neg"]]$x<10^(-5)])
        myhistdata_reschr<-myhist_f(mybreaks,combined_pval_l[[1]][["combined_pval_neg"]]$x[combined_pval_l[[1]][["combined_pval_neg"]]$x<10^(-5)])
        myhistdata<-c(rbind(myhistdata,myhistdata_reschr))
        mybarplot_f(1:length(myticks),log(myhistdata+1,10),myylim=0,mycols=rep(c("gold","cadetblue"),length(mybreaks)),ticks_are_breaks=TRUE,myxlabel='p.value',myylabel=expression(paste(log["10"],"(count+1)")),add=FALSE,xaxis=myxlabs)
        legend("top",c("observed","reshuffled"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        pch=15,col=c("gold","cadetblue"))
        dev.off()
        pdf(paste0(myfolderdropbox,"/empfdrT2_hist.pdf"))
        myhistdata<-myhist_f(mybreaks,combined_pvals[["combined_pval_T2"]]$x[combined_pvals[["combined_pval_T2"]]$x<10^(-5)])
        myhistdata_reschr<-myhist_f(mybreaks,combined_pval_l[[1]][["combined_pval_T2"]]$x[combined_pval_l[[1]][["combined_pval_T2"]]$x<10^(-5)])
        myhistdata<-c(rbind(myhistdata,myhistdata_reschr))
        mybarplot_f(1:length(myticks),log(myhistdata+1,10),myylim=0,mycols=rep(c("gold","cadetblue"),length(mybreaks)),ticks_are_breaks=TRUE,myxlabel='p.value',myylabel=expression(paste(log["10"],"(count+1)")),add=FALSE,xaxis=myxlabs)
        legend("top",c("observed","reshuffled"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        pch=15,col=c("gold","cadetblue"))
        dev.off()
 
        load(paste0(myfolder,"/combined_pvals_theoretical.RData"))
        combined_pvals_theoretical<-combined_pvals
        load(paste0(myfolder,"/combined_pvals.RData"))
        myempfdr<-c(
        sum(combined_pvals[["combined_fdr_T2"]]<0.05),
        sum(combined_pvals[["combined_fdr"]]<0.05),
        sum(combined_pvals[["combined_fdr_pos"]]<0.05),
        sum(combined_pvals[["combined_fdr_neg"]]<0.05))
        myfdr<-c(
        sum(combined_pvals_theoretical[["combined_fdr_T2"]]<0.05),
        sum(combined_pvals_theoretical[["combined_fdr"]]<0.05),
        sum(combined_pvals_theoretical[["combined_fdr_pos"]]<0.05),
        sum(combined_pvals_theoretical[["combined_fdr_neg"]]<0.05))
        mymat<-rbind(myempfdr,myfdr)
        colnames(mymat)<-c("T2","pKuli","D'>0","D'<0")
        pdf(paste0(myfolderdropbox,"/empfdr_vs_theoreticalfdr_hist.pdf"))
        barplot(mymat,beside=TRUE,col=c("floralwhite","azure3"),ylab="count")
        legend("topright",fill=c("floralwhite","azure3"), legend=c("empirical","theoretical"))
        dev.off()
}


}
#write snps files
{
snps_sign<-unique(combined_pvals[["combined_pval"]]$Group.1[combined_pvals[["combined_fdr"]]<0.05])
if (length(snps_sign)>=1){
mysnps<-matrix(unlist(sapply(snps_sign,function(x) strsplit(x,split="[.]"))),ncol=4,byrow=T)
signAbed<-cbind(mysnps[,1],as.numeric(mysnps[,3])-1,mysnps[,3])
signBbed<-cbind(mysnps[,2],as.numeric(mysnps[,4])-1,mysnps[,4])
write.table(signAbed,file=paste0(myfolder,"/anal/snpsA.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file=paste0(myfolder,"/anal/snpsB.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
}
snps_sign<-unique(combined_pvals[["combined_pval_pos"]]$Group.1[combined_pvals[["combined_fdr_pos"]]<0.05])
if (length(snps_sign)>=1){
mysnps_pos<-matrix(unlist(sapply(snps_sign,function(x) strsplit(x,split="[.]"))),ncol=4,byrow=T)
signAbed<-cbind(mysnps_pos[,1],as.numeric(mysnps_pos[,3])-1,mysnps_pos[,3])
signBbed<-cbind(mysnps_pos[,2],as.numeric(mysnps_pos[,4])-1,mysnps_pos[,4])
write.table(signAbed,file=paste0(myfolder,"/anal/snpsA_pos.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file=paste0(myfolder,"/anal/snpsB_pos.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
}
snps_sign<-unique(combined_pvals[["combined_pval_neg"]]$Group.1[combined_pvals[["combined_fdr_neg"]]<0.05])
if (length(snps_sign)>=1){
mysnps_neg<-matrix(unlist(sapply(snps_sign,function(x) strsplit(x,split="[.]"))),ncol=4,byrow=T)
signAbed<-cbind(mysnps_neg[,1],as.numeric(mysnps_neg[,3])-1,mysnps_neg[,3])
signBbed<-cbind(mysnps_neg[,2],as.numeric(mysnps_neg[,4])-1,mysnps_neg[,4])
write.table(signAbed,file=paste0(myfolder,"/anal/snpsA_neg.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file=paste0(myfolder,"/anal/snpsB_neg.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
}
#load(paste0(myfolder,"/dataFIN.RData"))
#dataFIN_sign<-dataFIN[order(combined_pval_sign)<(sum(empfdr<0.05)+0.5),]
#save(dataFIN_sign,file=paste0(myfolder,"/dataFIN_sign.RData"))
#dataFIN_top10000<-dataFIN[order(combined_pval_sign)[1:10000],]
#save(dataFIN_top10000,file=paste0(myfolder,"/dataFIN_top10000.RData"))
}                                   
#individual p-value per population and sharings
{
#I can look at sharing in many different ways:
#I have pop i and pop j; links_i and links_j, snp_i and snp_j (if significant denoted by *). 
#-fractionsharing of snps that are shared between pops snp_i^snp_j / (snp_i) and snp_i^snp_j / (snp_j)
#-fractionsharing of links that are shared between pops link_i^link_j / (link_i) and link_i^link_j / (link_j)
#restricting to only those that are present in both populations equals restricting on snp_i^snp_j or link_i^link_j
# fraction of snps that are significant among shared snps: #(snp_i*^snp_j*)/#(snp_i^snp_j) #~very low
# fraction of links that are significant among shared links: #(link_i*^link_j*)/#(link_i^link_j) #~very low
# prob of a link being significant in i| sign in j: #(link_i*^link_j*)/#link_j*
# prob of a shared shared link being significant in i| sign in j: #(link_i*^link_j*)/#link_j*
load(paste0(myfolder,"/combined_pvals.RData"))
load(paste0(myfolder,"/temp.RData"))
library(gplots)

links_sign_by_pop<-sapply(1:12,function(ipop) sum(combined_pvals[["poppos_fdr"]][[ipop]]<0.05)) #count of significant links per pop
links_sign_by_pop_neg<-sapply(1:12,function(ipop) sum(combined_pvals[["popneg_fdr"]][[ipop]]<0.05)) #count of significant links per pop

pdf(paste0(myfolderdropbox,"/barplot_sign_links.pdf"))
barplot(links_sign_by_pop,names.arg=popnames,cex.names=1,col="cadetblue",ylab="n.links")
dev.off()
pairid_sign_by_pop<-lapply(1:npopulations,function(ipop) data2pairid(temp[temp$pop==(ipop-1)])[combined_pvals[["poppos_fdr"]][[ipop]]<0.05]) #significant links per pop
ncompair<-read.table(paste0(myfolder,"/ncomparisons_pairwise.txt"),header=FALSE)
ncompair<-matrix(apply(ncompair[,-1]/1000000,MARGIN=2,sum),ncol=12,byrow=T) #/1000000 to avoid overflow
nshared_signlinks<-sapply(1:npopulations, function(i) sapply(1:npopulations,function(j) sum(!is.na(match(pairid_sign_by_pop[[i]],pairid_sign_by_pop[[j]])))))
ncompair[ncompair==0]<-1
ncompair[upper.tri(ncompair)]<-t(ncompair)[upper.tri(t(ncompair))] #making it square
(nshared_signlinks/ncompair)*1000000 #rows indicate prob for a link of being significant in both population 
nshared_signlinks/links_sign_by_pop #rows indicate fraction of links that are shared given a pop in order of links_sign_by_pop
pop_significant<-links_sign_by_pop!=0
if (sum(pop_significant)>0){
popnames_temp<-popnames[pop_significant]
ncompair<-ncompair[pop_significant,pop_significant]
nshared_signlinks<-nshared_signlinks[pop_significant,pop_significant] #observed number of shared links
links_sign_by_pop<-links_sign_by_pop[pop_significant]

expected_shared_signlinks<-sapply(1:sum(pop_significant), function(x) sapply(1:sum(pop_significant), function(z) (links_sign_by_pop[x]/ncompair[x,z])*1000000*links_sign_by_pop[z])) #expected number of shared links
excess_sharing_links<-nshared_signlinks/expected_shared_signlinks
excess_sharing_links[upper.tri(excess_sharing_links)]<-1
excess_sharing_links[excess_sharing_links==0]<-10^(-9)
for (z in 1:sum(pop_significant)) { excess_sharing_links[z,z]<-1}
excess_sharing_links<-log(excess_sharing_links)
rownames(excess_sharing_links)<-popnames_temp
colnames(excess_sharing_links)<-popnames_temp
excess_sharing_links[upper.tri(excess_sharing_links)]<-t(excess_sharing_links)[upper.tri(t(excess_sharing_links))] #making it square

myodds_sharing<-(nshared_signlinks/links_sign_by_pop)/((nshared_signlinks/ncompair)*1000000) #odds of being shared significant given that present in both (measure of how much sharing is overlapping) 

col_palette<-colorRampPalette(c("papayawhip","yellow","red"))(n=8)
col_breaks<-seq(min(excess_sharing_links,na.rm=TRUE),max(excess_sharing_links,na.rm=TRUE),(max(excess_sharing_links,na.rm=TRUE)-min(excess_sharing_links,na.rm=TRUE))/8)
length(col_breaks)
length(col_palette)
pdf(paste0(myfolderdropbox,"/heatmap.sharing_signlinks.pdf")) #normalized by shared alleles
heatmap.2(excess_sharing_links,trace="none",dendrogram="row",breaks=col_breaks,col=col_palette,main="linkage sharing")
dev.off()

fraction_shared_links<-nshared_signlinks/links_sign_by_pop
col_palette<-colorRampPalette(c("papayawhip","yellow","red"))(n=8)
col_breaks<-seq(0,max(fraction_shared_links),(max(fraction_shared_links))/8)
rownames(fraction_shared_links)<-popnames_temp
colnames(fraction_shared_links)<-popnames_temp
length(col_breaks)
length(col_palette)
myheat<-heatmap.2(fraction_shared_links,trace="none",Rowv=TRUE,Colv=TRUE,dendrogram="row",breaks=col_breaks,col=col_palette,main="fraction shared links",symm=TRUE)
pdf(paste0(myfolderdropbox,"/heatmap.fraction_sharing_signlinks.pdf")) #not normalized by shared alleles
heatmap.2(fraction_shared_links,trace="none",Rowv=myheat$rowInd,Colv=myheat$rowInd,dendrogram="row",breaks=col_breaks,col=col_palette,main="fraction shared links",symm=TRUE)
dev.off()
}

#check better, it seems really strange that less sharing: how come less sharing of snps when I have more
if (FALSE){
pop_snps_sign<-lapply(1:12,function(ipop) unique(c(data2snpA(temp[temp$pop==(ipop-1)][combined_pvals[["poppos_fdr"]][[ipop]]<0.05]),data2snpB(temp[temp$pop==(ipop-1)][combined_pvals[["poppos_fdr"]][[ipop]]<0.05])))) #snps in significant links
nshared_signsnps<-sapply(1:12, function(i) sapply(1:12,function(j) sum(!is.na(match(pop_snps_sign[[i]],pop_snps_sign[[j]])))))
nshared_signsnps/sapply(pop_snps_sign,function(x) length(x))
shared_snps<-read.table(paste0(myfolder,"/myfreq.sharing"))
matrix_sharing_snps<-sapply(0:(npopulations-1), function(i) sapply(0:(npopulations-1),function(j) {
sharing_vec<-apply(shared_snps==i | shared_snps==j,MARGIN=1,sum);
return(sum(sharing_vec==2)/sum(sharing_vec==2 | sharing_vec==1))
}))
for (z in 1:npopulations) { matrix_sharing_snps[z,z]<-1}
population_specific_snps<-(nshared_signsnps/sapply(pop_snps_sign,function(x) length(x)))/matrix_sharing_snps
rownames(population_specific_snps)<-popnames
colnames(population_specific_snps)<-popnames
col_breaks<-seq(0,2,0.2)#seq(0,max(population_specific_snps),(max(population_specific_snps))/8)
col_palette<-colorRampPalette(c("red","yellow","deepskyblue"))(n=(length(col_breaks)-1))
length(col_breaks)
length(col_palette)
pdf(paste0(myfolderdropbox,"/heatmap.sharedsnps_enrichment.pdf"))
heatmap.2(population_specific_snps,trace="none",dendrogram="row",breaks=col_breaks,col=col_palette,main="shared snps enrichment")
dev.off()
}
}
}
#================================================================================================^
#===============PER SAMPLE ANALYSES: IDENTIFYING BIASES==========================================
#================================================================================================
#---import files--------------------------------------------------------------------------------------v
#import infosamples
if (TRUE)
{
info1000g<-read.table("/mnt/scratch/fabrizio/LDLD/20130606_sample_info.txt",header=TRUE,sep='\t')
samples1000g<-system("zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | 
head -300 | grep '#CHROM' | head -1 | awk '{for (i=10; i<=NF; i++) print $i}'",intern=TRUE)
samplesabove95<-system("cat ~/workspace/1000genomes/above95.unrelated.samples",intern=TRUE)
nspops<-as.numeric(system("cat /mnt/scratch/fabrizio/LDLD/nspops.txt",intern=TRUE))
samples_order1000g_above95<-intersect(samplesabove95,samples1000g)
infosamples<-lapply(1:12,function(x) c())
for (i in 1:12)
{
  mysamples<-samples_imypop(samples_order1000g_above95,nspops,i)
  infosamples[[i]]<-info1000g[match(mysamples,info1000g$Sample),]
}
popsamplesi<-sapply(1:12,function(x) samples_imypopi(nspops,x))
}
#import other files (nA) to test mutations
if (FALSE)
{
#system("cat above95/coding1000g/chr*.tabchr*.minilog | grep nAB_above_freq_threshold | grep -v and_significant > above95/coding1000g/tot.nABabovethr")
#system("cat above95/coding1000g/chr*.tabchr*.minilog | grep nAB_above_freq_threshold | grep and_significant > above95/coding1000g/tot.nABabovethrsign")
#mylog2<-read.table("above95/coding1000g/tot.nABabovethr")
#mylog2<-sapply(2:dim(mylog2)[2],function(x) sum(mylog2[,x]))
#mylog3<-read.table("above95/coding1000g/tot.nABabovethrsign")
#mylog3<-sapply(2:dim(mylog3)[2],function(x) sum(mylog3[,x]))
#-------------------- create whole genome freqlog file (genotype per sample) ---
#rm ${myfolder}/all.freqlog
#for i in `seq 1 22`; do 
#cat ${myfolder}/chr$i.freqlog | awk -v FS='\t' -v OFS='\t' -v CHR=${i} '{if (NR>1){print CHR,$0}}' >> ${myfolder}/all.freqlog
#done
#gzip ${myfolder}/all.freqlog
#--------------------- creates polymorphism per individual (mutlog) file (newer)   ------------------v
#myfolder=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5
#gcc ~/Dropbox/LDLD/scripts/LDLD_filter5.c -lgmp -lm -fno-stack-protector -o filterbyfreq.out
#for j in "/reschr1/" "/reschr2/"; do #"/"
#for i in `seq 1 22`; do 
#./filterbyfreq.out ${myfolder}${j}chr$i.tab ${myfolder}${j}chr$i.freqlog ${myfolder}${j}chr$i.tab.mutlog 0.05 1220 & 
#done;done

#import nA
for (i in 1:22)
{
if (i==1){mutperssample<-as.numeric(read.table(paste0(myfolder,"/chr",i,".mutlog")))}
else {res2<-as.numeric(read.table(paste0(myfolder,"/chr",i,".mutlog")));mutperssample<-rbind(mutperssample,res2)}
}
l_mutsamples<-lapply(1:12,function(x) c())
for (i in 1:12)
{
  l_mutsamples[[i]]<-apply(mutperssample[,samples_imypopi(nspops,i)],2,sum)
}
save(l_mutsamples,file=paste0(myfolder,"/l_mutsamples.RData"))
#import nA for reshuffled
if (FALSE){
l_mutsamples_reschr<-list()
for (myperm in 1:nperm){for (i in 1:22)
{
if (i==1){mutperssample<-as.numeric(read.table(paste0(myfolder,"/reschr",myperm,"/chr",i,".mutlog")))}
else {res2<-as.numeric(read.table(paste0(myfolder,"/reschr",myperm,"/chr",i,".mutlog")));mutperssample<-rbind(mutperssample,res2)}
}
l_mutsamples_reschr[[myperm]]<-lapply(1:12,function(x) c())
for (i in 1:12)
{
l_mutsamples_reschr[[myperm]][[i]]<-apply(mutperssample[,samples_imypopi(nspops,i)],2,sum)
}}
}
} 
#----------------------------------------------------------------------------------------------------^
#=======================
#done by taking only significant links 
#=======================
#compute nAB
{
#using_sign_for_nAB<-FALSE
#if (!using_sign_for_nAB){load(paste0(myfolder,"/dataFIN_top10000.RData"));dataFIN_sign<-dataFIN_top10000}
if (exists("mydata")){rm(mydata)}
if (exists("temp")){rm(temp)}
load(paste0(myfolder,"/combined_pvals.RData"))
load(paste0(myfolder,"/combined_pval_l.RData"))
#calculate nAB with all minor alleles
{
mylogpop<-list();mylogpop_pos<-list();mylogpop_neg<-list();mylogpop_pos20<-list();mylogpop_pos01<-list();logpop_nAB<-list();logpop_nAB_pos<-list();logpop_nAB_neg<-list();logpop_nAB_pos20<-list();logpop_nAB_pos01<-list();
for ( ipop in 0:11)    
{
    print(ipop)
    mydata_temp<-fread(paste0("zcat ", myfolder,"/logs/all.pop",ipop,".log.gz"))
    #for some fields the output is not 0-9, but either triallelic or something, so I should remove those
    #real problem comes however when repeats. I updated files to fix, but check. otherwise: #for i in `seq 0 11`;do zcat all.pop${i}.log.gz | sort -Vu | gzip -f > temp; mv temp all.pop${i}.log.gz;done
    mydata_temp<-mydata_temp[apply(mydata_temp[,5:(ncol(mydata_temp)-3),with=FALSE],MARGIN=1,function(x) all(x<10)),]
    #mysign<-apply(dataFIN[empfdr<0.05,c(1,3,2,4),with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse="."))
       if (sum(combined_pvals[["combined_fdr"]]<0.05)>0){
           mylogpop[[ipop+1]]<-intersect_links(mydata_temp,combined_pvals[["combined_pval"]]$Group.1[combined_pvals[["combined_fdr"]]<0.05],B.is.pairid=TRUE,myorder=c(1,2,3,4))
           logpop_nAB[[ipop+1]]<-apply(mylogpop[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
           } else {mylogpop[[ipop+1]]<-0;logpop_nAB[[ipop+1]]<-0}
      if (sum(combined_pvals[["combined_fdr_neg"]]<0.05)>0){
          mylogpop_neg[[ipop+1]]<-intersect_links(mydata_temp,combined_pvals[["combined_pval_neg"]]$Group.1[combined_pvals[["combined_fdr_neg"]]<0.05],B.is.pairid=TRUE,myorder=c(1,2,3,4))
          logpop_nAB_neg[[ipop+1]]<-apply(mylogpop_neg[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
      } else {mylogpop_neg[[ipop+1]]<-0;logpop_nAB_neg[[ipop+1]]<-0}
       if (sum(combined_pvals[["combined_fdr_pos"]]<0.05)>0){
           mylogpop_pos[[ipop+1]]<-intersect_links(mydata_temp,combined_pvals[["combined_pval_pos"]]$Group.1[combined_pvals[["combined_fdr_pos"]]<0.05],B.is.pairid=TRUE,myorder=c(1,2,3,4))
           logpop_nAB_pos[[ipop+1]]<-apply(mylogpop_pos[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
           } else {mylogpop_pos[[ipop+1]]<-0;logpop_nAB_pos[[ipop+1]]<-0}
      if (sum(combined_pvals[["combined_fdr_pos"]]<0.2)>0){
          mylogpop_pos20[[ipop+1]]<-intersect_links(mydata_temp,combined_pvals[["combined_pval_pos"]]$Group.1[combined_pvals[["combined_fdr_pos"]]<0.2],B.is.pairid=TRUE,myorder=c(1,2,3,4))
          logpop_nAB_pos20[[ipop+1]]<-apply(mylogpop_pos20[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
          } else {mylogpop_pos20[[ipop+1]]<-0;logpop_nAB_pos20[[ipop+1]]<-0}
     if (sum(combined_pvals[["combined_fdr_pos"]]<0.0001)>0){
         mylogpop_pos01[[ipop+1]]<-intersect_links(mydata_temp,combined_pvals[["combined_pval_pos"]]$Group.1[combined_pvals[["combined_fdr_pos"]]<0.0001],B.is.pairid=TRUE,myorder=c(1,2,3,4))
         logpop_nAB_pos01[[ipop+1]]<-apply(mylogpop_pos01[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
         } else {mylogpop_pos01[[ipop+1]]<-0;logpop_nAB_pos01[[ipop+1]]<-0}
}

#certain samples have crazy amount of double hets -> check in ipop=8 who is V78
#> sum(mylogpop_pos20[[ipop+1]]$V78==5)
#[1] 13298
#> sum(mylogpop_pos20[[ipop+1]]$V79==5)
#[1] 71

dataLDLD<-mylogpop_pos20[[ipop+1]]
dataLDLDA<-apply(dataLDLD[,1:4,with=FALSE],MARGIN=1,function(x) paste0(x,collapse='.'))
dataLDLDB<-apply(mydata_temp[,1:4,with=FALSE],MARGIN=1,function(x) paste0(x,collapse='.'))
length(dataLDLDA)
length(unique(dataLDLDA))
length(dataLDLDB)
length(unique(dataLDLDB))


save(mylogpop,file=paste0(myfolder,"/mylogpop.RData"))
logpop_nAB<-sapply(logpop_nAB,function(x) x[c(-length(x),-length(x)+1)]);save(logpop_nAB,file=paste0(myfolder,"/logpop_nAB.RData"))
save(mylogpop_pos,file=paste0(myfolder,"/mylogpop_pos.RData"))
logpop_nAB_pos<-sapply(logpop_nAB_pos,function(x) x[c(-length(x),-length(x)+1)]);save(logpop_nAB_pos,file=paste0(myfolder,"/logpop_nAB_pos.RData"))
save(mylogpop_neg,file=paste0(myfolder,"/mylogpop_neg.RData"))
logpop_nAB_neg<-sapply(logpop_nAB_neg,function(x) x[c(-length(x),-length(x)+1)]);save(logpop_nAB_neg,file=paste0(myfolder,"/logpop_nAB_neg.RData"))
save(mylogpop_pos20,file=paste0(myfolder,"/mylogpop_pos20.RData"))
logpop_nAB_pos20<-sapply(logpop_nAB_pos20,function(x) x[c(-length(x),-length(x)+1)]);save(logpop_nAB_pos20,file=paste0(myfolder,"/logpop_nAB_pos20.RData"))
save(mylogpop_pos01,file=paste0(myfolder,"/mylogpop_pos01.RData"))
logpop_nAB_pos01<-sapply(logpop_nAB_pos01,function(x) x[c(-length(x),-length(x)+1)]);save(logpop_nAB_pos01,file=paste0(myfolder,"/logpop_nAB_pos01.RData"))
}

#null models for nAB
{
niterations<-1    
n_sign_links<-sum(combined_pvals[["combined_fdr"]]<0.05)
if (n_sign_links>0){
for ( myperm in 1:nperm){ #total nAB kept fixed
mylogpop_reschr<-list();logpop_reschr_nAB<-list()
for ( ipop in 0:(npopulations-1))
    {
    fraction_nlinks<-1
    mydata_temp<-fread(paste0("zcat ",myfolder,"/reschr",myperm,"/logs/all.pop",ipop,".log.gz"))
    print(ipop)
    for (i in 1:niterations)
        {
        sign_links_reschr<-combined_pval_l[[myperm]][["combined_pval"]]$Group.1[order(combined_pval_l[[myperm]][["combined_pval"]]$x)<=round(n_sign_links/fraction_nlinks)]
        mylogpop_reschr[[ipop+1]]<-intersect_links(mydata_temp,sign_links_reschr,B.is.pairid=TRUE,myorder=c(1,2,3,4))
        #combined_pval_l[[myperm]][["combined_pval"]]$Group.1[combined_pval_l[[myperm]][["combined_fdr"]]<0.05],B.is.pairid=TRUE,myorder=c(1,2,3,4))
        #intersect_links(mydata_temp,sign_links_reschr)
        logpop_reschr_nAB[[ipop+1]]<-apply(mylogpop_reschr[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
        fraction_nlinks<-fraction_nlinks*(sum(logpop_reschr_nAB[[ipop+1]])/sum(logpop_nAB[[ipop+1]]))
        }
    }
save(logpop_reschr_nAB,file=paste0(myfolder,"/logpop_reschr",myperm,"_nAB.RData"))
save(mylogpop_reschr,file=paste0(myfolder,"/mylogpop_reschr",myperm,".RData"))
}}
n_sign_links<-sum(combined_pvals[["combined_fdr_neg"]]<0.05)
if (n_sign_links>0){
for ( myperm in 1:nperm){ #total nAB kept fixed
mylogpop_reschr<-list();logpop_reschr_nAB<-list()
for ( ipop in 0:(npopulations-1))
    {
    fraction_nlinks<-1
    mydata_temp<-fread(paste0("zcat ",myfolder,"/reschr",myperm,"/logs/all.pop",ipop,".log.gz"))
    print(ipop)
    for (i in 1:niterations)
        {
        sign_links_reschr<-combined_pval_l[[myperm]][["combined_pval_neg"]]$Group.1[order(combined_pval_l[[myperm]][["combined_pval_neg"]]$x)<=round(n_sign_links/fraction_nlinks)]
        mylogpop_reschr[[ipop+1]]<-intersect_links(mydata_temp,sign_links_reschr,B.is.pairid=TRUE,myorder=c(1,2,3,4))
        #combined_pval_l[[myperm]][["combined_pval"]]$Group.1[combined_pval_l[[myperm]][["combined_fdr"]]<0.05],B.is.pairid=TRUE,myorder=c(1,2,3,4))
        #intersect_links(mydata_temp,sign_links_reschr)
        logpop_reschr_nAB[[ipop+1]]<-apply(mylogpop_reschr[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
        fraction_nlinks<-fraction_nlinks*(sum(logpop_reschr_nAB[[ipop+1]])/sum(logpop_nAB_neg[[ipop+1]]))
        }
    }
    mylogpop_neg_reschr<-mylogpop_reschr
    logpop_nAB_neg_reschr<-logpop_reschr_nAB
save(logpop_nAB_neg_reschr,file=paste0(myfolder,"/logpop_nAB_neg_reschr",myperm,"_nAB.RData"))
save(mylogpop_neg_reschr,file=paste0(myfolder,"/mylogpop_neg_reschr",myperm,".RData"))
}}
n_sign_links<-sum(combined_pvals[["combined_fdr_pos"]]<0.05)
if (n_sign_links>0){
for ( myperm in 1:nperm){ #total nAB kept fixed
mylogpop_reschr<-list();logpop_reschr_nAB<-list()
for ( ipop in 0:(npopulations-1))
    {
    fraction_nlinks<-1
    mydata_temp<-fread(paste0("zcat ",myfolder,"/reschr",myperm,"/logs/all.pop",ipop,".log.gz"))
    print(ipop)
    for (i in 1:niterations)
        {
        sign_links_reschr<-combined_pval_l[[myperm]][["combined_pval"]]$Group.1[order(combined_pval_l[[myperm]][["combined_pval"]]$x)<=round(n_sign_links/fraction_nlinks)]
        mylogpop_reschr[[ipop+1]]<-intersect_links(mydata_temp,sign_links_reschr,B.is.pairid=TRUE,myorder=c(1,2,3,4))
        #combined_pval_l[[myperm]][["combined_pval"]]$Group.1[combined_pval_l[[myperm]][["combined_fdr"]]<0.05],B.is.pairid=TRUE,myorder=c(1,2,3,4))
        #intersect_links(mydata_temp,sign_links_reschr)
        logpop_reschr_nAB[[ipop+1]]<-apply(mylogpop_reschr[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
        fraction_nlinks<-fraction_nlinks*(sum(logpop_reschr_nAB[[ipop+1]])/sum(logpop_nAB_pos[[ipop+1]]))
        }
    }
    mylogpop_pos_reschr<-mylogpop_reschr
    logpop_nAB_pos_reschr<-logpop_reschr_nAB
save(logpop_nAB_pos_reschr,file=paste0(myfolder,"/logpop_nAB_pos_reschr",myperm,"_nAB.RData"))
save(mylogpop_pos_reschr,file=paste0(myfolder,"/mylogpop_pos_reschr",myperm,".RData"))
}}
n_sign_links<-sum(combined_pvals[["combined_fdr_pos"]]<0.2)
if (n_sign_links>0){
for ( myperm in 1:nperm){ #total nAB kept fixed
mylogpop_reschr<-list();logpop_reschr_nAB<-list()
for ( ipop in 0:(npopulations-1))
    {
    fraction_nlinks<-1
    mydata_temp<-fread(paste0("zcat ",myfolder,"/reschr",myperm,"/logs/all.pop",ipop,".log.gz"))
    print(ipop)
    for (i in 1:niterations)
        {
        sign_links_reschr<-combined_pval_l[[myperm]][["combined_pval"]]$Group.1[order(combined_pval_l[[myperm]][["combined_pval"]]$x)<=round(n_sign_links/fraction_nlinks)]
        mylogpop_reschr[[ipop+1]]<-intersect_links(mydata_temp,sign_links_reschr,B.is.pairid=TRUE,myorder=c(1,2,3,4))
        #combined_pval_l[[myperm]][["combined_pval"]]$Group.1[combined_pval_l[[myperm]][["combined_fdr"]]<0.05],B.is.pairid=TRUE,myorder=c(1,2,3,4))
        #intersect_links(mydata_temp,sign_links_reschr)
        logpop_reschr_nAB[[ipop+1]]<-apply(mylogpop_reschr[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
        fraction_nlinks<-fraction_nlinks*(sum(logpop_reschr_nAB[[ipop+1]])/sum(logpop_nAB_pos20[[ipop+1]]))
        }
    }
    mylogpop_pos20_reschr<-mylogpop_reschr
    logpop_nAB_pos20_reschr<-logpop_reschr_nAB
save(logpop_nAB_pos20_reschr,file=paste0(myfolder,"/logpop_nAB_pos20_reschr",myperm,"_nAB.RData"))
save(mylogpop_pos20_reschr,file=paste0(myfolder,"/mylogpop_pos20_reschr",myperm,".RData"))
}}
}


}
#build infoplots
{
#->log tables with only significant links

require(vioplot)
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
load(paste0(myfolder,"/logpop_nAB_pos20.RData"))
load(paste0(myfolder,"/logpop_nAB_neg.RData"))
load(paste0(myfolder,"/logpop_nAB.RData"))
logpop_reschr_nAB_l<-list();logpop_reschr_nAB_neg_l<-list();logpop_reschr_nAB_pos_l<-list();logpop_reschr_nAB_pos20_l<-list();

suppressWarnings(try(for (z in 1:nperm) { load(paste0(myfolder,"/logpop_reschr",z,"_nAB.RData")); logpop_reschr_nAB_l[[z]]<-logpop_reschr_nAB; },silent=TRUE))
suppressWarnings(try(for (z in 1:nperm) { load(paste0(myfolder,"/logpop_nAB_neg_reschr",z,"_nAB.RData")); logpop_reschr_nAB_neg_l[[z]]<-logpop_nAB_neg_reschr; },silent=TRUE))
suppressWarnings(try(for (z in 1:nperm) { load(paste0(myfolder,"/logpop_nAB_pos_reschr",z,"_nAB.RData")); logpop_reschr_nAB_pos_l[[z]]<-logpop_nAB_pos_reschr; },silent=TRUE))
suppressWarnings(try(for (z in 1:nperm) { load(paste0(myfolder,"/logpop_nAB_pos20_reschr",z,"_nAB.RData")); logpop_reschr_nAB_pos20_l[[z]]<-logpop_nAB_pos_reschr; },silent=TRUE))


#short infoplots
{
#violin plots
if (FALSE)
{
myvio<-function(x) vioplot(x)

pdf(paste0(myfolderdropbox,"/infoplot_violin.pdf"))
#par(mfrow=c(1,12),mar=c(0.5,0.5,0.5,0.5))
plot(-1:(-3),xlim=c(0,12),ylim=c(0,1.1))
for (i in 1:12) { vioplot(logpop_nAB[[i]]/max(logpop_nAB[[i]]),add=TRUE,at=i,range=c(0.5,0.6))}
dev.off()


myhist<-hist(logpop_nAB[[i]]/max(logpop_nAB[[i]]))

quantile_intervals_col<-seq(0,1,0.1)
quantile_intervals<-seq(0,1,0.1)
myshades_tot<-c()
for (i in 1:12) { 
mydata_pop<-logpop_nAB[[i]]/max(logpop_nAB[[i]])
mycolor_variable<-infosamples[[i]]$Total.Exome.Sequence
myshades_tot<-c(myshades_tot,sapply(1:(length(quantile_intervals)-1), function(x) mean(mycolor_variable[mydata_pop>quantile(mydata_pop,quantile_intervals)[x] & mydata_pop<quantile(mydata_pop,quantile_intervals)[x+1] ])))
}
#myshades<-sapply(1:(length(quantile_intervals)-1), function(x) mean(mycolor_variable[mydata_pop>quantile(mydata_pop,quantile_intervals)[x] & mydata_pop<quantile(mydata_pop,quantile_intervals)[x+1] ]))
#range_mycolor_variable<-range(unlist(sapply(infosamples,function(x) x$Total.Exome.Sequence)))
range_mycolor_variable<-range(myshades_tot)
mycolintervals<-seq(range_mycolor_variable[1],range_mycolor_variable[2],(range_mycolor_variable[2]-range_mycolor_variable[1])/length(seq(0,1,0.1)))

if (FALSE){
myhist<-hist(rnorm(100),breaks=10)
mydata<-rnorm(50)
myshades<-sapply(quantile(mydata,quantile_intervals), function(x) mean(mydata[mydata>x]))[-length(quantile_intervals)]
mypalette<-colorRampPalette(c('gray','white'))(ncolclasses)
myrangevalues<-range(mydata)
mycolintervals<-seq(min(mydata),max(mydata),(max(mydata)-min(mydata))/(ncolclasses))
mycols<-sapply(myshades,function(x) mypalette[sum(x>=mycolintervals)])
}

vioplot_f<-function(mydata,nbreaks=FALSE,myadd=FALSE,mycol="gray",at=0,widthx=1){ 
if (!is.numeric(nbreaks)){nbreaks<-quantile(mydata,seq(0,1,0.1))} 
myhist<-hist(mydata,breaks=nbreaks,plot=FALSE);
if ( myadd) { points(myhist$density,myhist$mids,type="n") } else { plot(myhist$density,myhist$mids,ylim=c(min(myhist$breaks),max(myhist$breaks)),xlim=c(-max(myhist$density),max(myhist$density)),type="n") }
#lo <- loess(c(0,myhist$density,0)~c(myhist$breaks[1],myhist$mids,myhist$breaks[length(myhist$breaks)]))
#xcoords<-predict(lo,ycoords)
ycoords<-sort(c(myhist$breaks,myhist$mids))
coords<-approx(c(myhist$breaks[1],myhist$mids,myhist$breaks[length(myhist$breaks)]),c(0,myhist$density,0),xout=ycoords)
xcoords<-coords$y; ycoords[ycoords<0]<-0
xcoords<-xcoords*widthx/max(xcoords);
ycoords<-coords$x
if (length(mycol)==1) { polygon(c(xcoords,-rev(xcoords)),c(ycoords,rev(ycoords)),col=mycol) } else
    {
    for ( z in 0:(length(nbreaks)-1)){ #((0:(length(nbreaks)-1)))){ 
    polygon(c(xcoords[(1+2*z):(1+2*z+2)]-0.04,-rev(xcoords[(1+2*z):(1+2*z+2)])+0.04)+at,c(ycoords[(1+2*z):(1+2*z+2)],rev(ycoords[(1+2*z):(1+2*z+2)])),col=mycol[z+1],lty=0)
    }
    polygon(c(xcoords,-rev(xcoords))+at,c(ycoords,rev(ycoords)))
}
}

ncolclasses<-10
mydata_pop<-log(logpop_nAB[[i]]+0.1)/log(max(logpop_nAB[[i]]+0.1))
ord<-order(mydata_pop)[-(1:2)]
infosamples[[i]]$Main.project.LC.Centers[ord]
mycolor_variable<-infosamples[[i]]$Total.Exome.Sequence
mycolintervals<-seq(min(mycolor_variable),max(mycolor_variable),(max(mycolor_variable)-min(mycolor_variable))/ncolclasses)
myshades<-sapply(1:(length(quantile_intervals)-1), function(x) mean(mycolor_variable[mydata_pop>quantile(mydata_pop,quantile_intervals)[x] & mydata_pop<quantile(mydata_pop,quantile_intervals)[x+1] ]))
mycols<-sapply(myshades,function(x) sum(x>=mycolintervals))

}

#infoplot barplots
for (mysign_logpop_nAB in c("","_pos","_neg","_pos20"))
{
mylogpopnAB<-get(paste0("logpop_nAB",mysign_logpop_nAB))
mylogpopnAB_reschr<-get(paste0("logpop_reschr_nAB",mysign_logpop_nAB,"_l"))
if (sum(unlist(mylogpopnAB))>0){
#infoplot barplots normalized
{
add_null_distr<-FALSE
mystep<-0.025
mycolors_sq<-c("gold1","deepskyblue1","coral","darkslategray4","forestgreen","khaki4","tan","lightsalmon","lemonchiffon3","lightblue1","firebrick1","mediumpurple1","mediumaquamarine","red4","seagreen3","azure1","peachpuff","orangered2","limegreen","hotpink3")
mylevels<-unique(unlist(sapply(1:npopulations,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
allpops<-unlist(sapply(1:npopulations,function(x) unlist(mylogpopnAB[[x]])[1:(length(mylogpopnAB[[x]])-2)] ))
allpops<-log(allpops+1)
labs<-c(min(allpops),min(allpops)+(max(allpops)-min(allpops))/2,max(allpops)) #labels in absolute scale
allpops<-allpops-min(allpops); allpops<-allpops/max(allpops);
table(unlist(sapply(1:npopulations,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
myticks<-seq(0,1,1/30) #seq(floor(min(allpops-0.02)/mystep)*mystep,1,mystep)
pdf(paste0(myfolderdropbox,"/barplot_relativenAB",mysign_logpop_nAB,".pdf"))
if (add_null_distr) { 
    par(mfrow=c(4,13),mar=c(0.5,0.1,0.1,0.1));
    layout(matrix(1:(14*4), 4, 14, byrow = TRUE),widths=c(0.5,rep(1,14*4-1)))
    plot.new()
    #layout(matrix(1:24, 2, 12, byrow=TRUE))
    i<-10
    mydata_pop<-log(mylogpopnAB_reschr[[i]]+0.1)/log(max(mylogpopnAB_reschr[[i]]+0.1)); 
    mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
    mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks))
    barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
    #axis(2,at=1:length(myticks),myticks,cex.axis=0.9)
    text(1*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/2,0.1,paste0("null ",popnames[[i]]))
    box("plot", col="black")  
    } else 
    {
    par(mfrow=c(4,npopulations+1),mar=c(0.5,0.1,0.1,0.1))
    layout(matrix(1:((npopulations+1)*4), 4, (npopulations+1), byrow = TRUE),widths=c(0.5,rep(1,(npopulations+1)*4-1)))
    plot.new()
    }
for (i in 1:npopulations){
mydata_pop<-log(mylogpopnAB[[i]]+1)
mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
mydata_pop<-(mydata_pop-min(labs))/(max(labs)-min(labs))
#mydata_pop<-log(mylogpopnAB[[i]]+0.1)/log(max(mylogpopnAB[[i]]+0.1))
#mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
#ord<-order(mydata_pop)[-(1:2)]
mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks))
barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
if (i==1){
axis(2,at=c(0,length(myticks)/2,length(myticks)),round(labs,1),cex.axis=0.9)
}
text(2*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/3,0.1,popnames[[i]])
#rect( grconvertX(0.005, from='ndc'), grconvertY(0.505, from='ndc'),
#     grconvertX(0.495, from='ndc'), grconvertY(0.995, from='ndc'))
box("plot", col="black")  
}
#box("figure", col="green")  
#box('outer')
dev.off()
}
if (FALSE){
#infoplot barplots non normalized
{
add_null_distr<-FALSE
mystep<-0.025
mycolors_sq<-c("gold1","deepskyblue1","coral","darkslategray4","forestgreen","khaki4","tan","lightsalmon","lemonchiffon3","lightblue1","firebrick1","mediumpurple1","mediumaquamarine","red4","seagreen3","azure1","peachpuff","orangered2","limegreen","hotpink3")
mylevels<-unique(unlist(sapply(1:12,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
allpops<-unlist(sapply(1:12,function(x) unlist(mylogpopnAB[[x]][1:(length(mylogpopnAB[[x]])-2)] )))
table(unlist(sapply(1:12,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
mystep<-(max(allpops)-min(allpops))/30
myticks<-seq(floor(min(allpops-0.02)/mystep)*mystep,floor(min(allpops-0.02)/mystep)*mystep+34*mystep ,mystep)
pdf(paste0(myfolderdropbox,"/barplot_absolutenAB",mysign_logpop_nAB,".pdf"))
if (add_null_distr) { 
    par(mfrow=c(4,13),mar=c(0.5,0.1,0.1,0.1));
    layout(matrix(1:(14*4), 4, 14, byrow = TRUE),widths=c(0.5,rep(1,14*4-1)))
    plot.new()
    #layout(matrix(1:24, 2, 12, byrow=TRUE))
    i<-10
    mydata_pop<-mylogpopnAB_reschr[[i]]; 
    mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
    mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks))
    barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
    #axis(2,at=1:length(myticks),myticks,cex.axis=0.9)
    text(1*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/2,0.1,paste0("null ",popnames[[i]]))
    box("plot", col="black")  
    } else 
    {
    par(mfrow=c(4,npopulations+1),mar=c(0.5,0.1,0.1,0.1))
    layout(matrix(1:((npopulations+1)*4), 4, (npopulations+1), byrow = TRUE),widths=c(0.5,rep(1,(npopulations+1)*4-1)))
    plot.new()
    }
for (i in 1:12){
mydata_pop<-mylogpopnAB[[i]]
mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
#ord<-order(mydata_pop)[-(1:2)]
mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks))
barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
text(2*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/3,0.1,popnames[[i]])
#rect( grconvertX(0.005, from='ndc'), grconvertY(0.505, from='ndc'),
#     grconvertX(0.495, from='ndc'), grconvertY(0.995, from='ndc'))
box("plot", col="black")  
}
#box("figure", col="green")  
#box('outer')
dev.off()
myhist_f(1:10,mybreaks=0:10)
}
#infoplot barplots non normalized with log
{
add_null_distr<-FALSE
mystep<-0.025
mycolors_sq<-c("gold1","deepskyblue1","coral","darkslategray4","forestgreen","khaki4","tan","lightsalmon","lemonchiffon3","lightblue1","firebrick1","mediumpurple1","mediumaquamarine","red4","seagreen3","azure1","peachpuff","orangered2","limegreen","hotpink3")
mylevels<-unique(unlist(sapply(1:12,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
allpops<-unlist(sapply(1:12,function(x) unlist(log(mylogpopnAB[[x]])[1:(length(mylogpopnAB[[x]])-2)] )))
if (sum(allpops==-Inf)>0) { 
for ( x in 1:length(mylogpopnAB)){mylogpopnAB[[x]]<-mylogpopnAB[[x]]+1}
}
table(unlist(sapply(1:12,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
allpops<-unlist(sapply(1:12,function(x) unlist(log(mylogpopnAB[[x]])[1:(length(mylogpopnAB[[x]])-2)] )))
mystep<-(max(allpops)-min(allpops))/30
myticks<-seq(floor(min(allpops-0.02)/mystep)*mystep,floor(min(allpops-0.02)/mystep)*mystep+34*mystep ,mystep)
if (sum(allpops==-Inf)>0) { 
pdf(paste0(myfolderdropbox,"/barplot_logplus1absolutenAB",mysign_logpop_nAB,".pdf"))
} else {pdf(paste0(myfolderdropbox,"/barplot_logabsolutenAB",mysign_logpop_nAB,".pdf"))}
if (add_null_distr) { 
    par(mfrow=c(4,13),mar=c(0.5,0.1,0.1,0.1));
    layout(matrix(1:(14*4), 4, 14, byrow = TRUE),widths=c(0.5,rep(1,14*4-1)))
    plot.new()
    #layout(matrix(1:24, 2, 12, byrow=TRUE))
    i<-10
    mydata_pop<-mylogpopnAB_reschr[[i]]; 
    mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
    mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks))
    barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
    #axis(2,at=1:length(myticks),myticks,cex.axis=0.9)
    text(1*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/2,0.1,paste0("null ",popnames[[i]]))
    box("plot", col="black")  
    } else 
    {
    par(mfrow=c(4,npopulations+1),mar=c(0.5,0.1,0.1,0.1))
    layout(matrix(1:((npopulations+1)*4), 4, (npopulations+1), byrow = TRUE),widths=c(0.5,rep(1,(npopulations+1)*4-1)))
    plot.new()
    }
for (i in 1:12){
mydata_pop<-log(mylogpopnAB[[i]])
mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
#ord<-order(mydata_pop)[-(1:2)]
mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks))
barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
text(2*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/3,0.1,popnames[[i]])
#rect( grconvertX(0.005, from='ndc'), grconvertY(0.505, from='ndc'),
#     grconvertX(0.495, from='ndc'), grconvertY(0.995, from='ndc'))
box("plot", col="black")  
}
#box("figure", col="green")  
#box('outer')
dev.off()
myhist_f(1:10,mybreaks=0:10)
}
}
}
}

}


load(paste0(myfolder,"/logpop_nAB_pos.RData"))
load(paste0(myfolder,"/logpop_nAB_neg.RData"))
load(paste0(myfolder,"/logpop_nAB.RData"))
#violin infoplots that look like butterflies
{

vioplot_f<-function(mydata,nbreaks=FALSE,myadd=FALSE,mycol="gray",at=0,widthx=1){ 
if (!is.numeric(nbreaks)){nbreaks<-quantile(mydata,seq(0,1,0.1))} else {nbreaks<-quantile(mydata,nbreaks)}
myhist<-hist(mydata,breaks=nbreaks,plot=FALSE);
if ( myadd) { points(myhist$density,myhist$mids,type="n") } else { plot(myhist$density,myhist$mids,ylim=c(min(myhist$breaks),max(myhist$breaks)),xlim=c(-max(myhist$density),max(myhist$density)),type="n") }
#lo <- loess(c(0,myhist$density,0)~c(myhist$breaks[1],myhist$mids,myhist$breaks[length(myhist$breaks)]))
#xcoords<-predict(lo,ycoords)
ycoords<-sort(c(myhist$breaks,myhist$mids))
coords<-approx(c(myhist$breaks[1],myhist$mids,myhist$breaks[length(myhist$breaks)]),c(0,myhist$density,0),xout=ycoords)
xcoords<-coords$y; ycoords[ycoords<0]<-0
xcoords<-xcoords*widthx/max(xcoords);
ycoords<-coords$x
if (length(mycol)==1) { polygon(c(xcoords,-rev(xcoords)),c(ycoords,rev(ycoords)),col=mycol) } else
    {
    for ( z in 0:(length(nbreaks)-1)){ #((0:(length(nbreaks)-1)))){ 
    polygon(c(xcoords[(1+2*z):(1+2*z+2)]+0.04,-rev(xcoords[(1+2*z):(1+2*z+2)])-0.04)+at,c(ycoords[(1+2*z):(1+2*z+2)],rev(ycoords[(1+2*z):(1+2*z+2)])),col=mycol[z+1],lty=0)
    }
    polygon(c(xcoords+0.04,-rev(xcoords)-0.04)+at,c(ycoords,rev(ycoords)))
}
}


mycolors_sq<-c("lightsalmon","tan","coral","darkslategray4","lemonchiffon3","lightblue1","gold1","firebrick1","mediumpurple1","mediumaquamarine","red4","seagreen3","azure1","deepskyblue1","forestgreen","peachpuff","orangered2","limegreen","hotpink3","khaki4")
quantile_intervals<-seq(0,1,0.05)
ncolclasses<-10
#mypalette<-colorRampPalette(c('black','gray',"white","white","white"))(ncolclasses)
#mypalette<-colorRampPalette(c('red','orange',"yellow","white"))(ncolclasses)
#mypalette<-colorRampPalette(c("gray8","gray78","gray92","gray100","white"))(ncolclasses)
mypalette<-colorRampPalette(c("gray8","gray78","gray100","white"))(ncolclasses)
pdf(paste0(myfolderdropbox,"/infoplot_violin.pdf"))
par(mfrow=c(2,1))
plot(-1:(-3),xlim=c(0,12),ylim=c(0,1.1),ylab="nAB/max(nAB)",xlab="population",xaxt='n',yaxt='n')
axis(side=1, at=0:12, labels=c("null",popnames),cex.axis=0.5)
for (i in 1:12) { 
#mydata_pop<-logpop_nAB[[i]]/max(logpop_nAB[[i]]+0.1)
mydata_pop<-log(logpop_nAB[[i]]+0.1)/log(max(logpop_nAB[[i]]+0.1))
ord<-order(mydata_pop)[-(1:2)]
infosamples[[i]]$Main.project.LC.Centers[ord]
mycolor_variable<-infosamples[[i]]$Total.Exome.Sequence
mycolintervals<-seq(min(mycolor_variable),max(mycolor_variable),(max(mycolor_variable)-min(mycolor_variable))/ncolclasses)
myshades<-sapply(1:(length(quantile_intervals)-1), function(x) mean(mycolor_variable[mydata_pop>quantile(mydata_pop,quantile_intervals)[x] & mydata_pop<quantile(mydata_pop,quantile_intervals)[x+1] ]))
mycols<-sapply(myshades,function(x) sum(x>=mycolintervals))
#print(mycols)
vioplot_f(mydata_pop[ord],mycol=mypalette[mycols],at=i,widthx=0.4,myadd=TRUE,nbreaks=quantile_intervals)
points(rep(i,length(ord)),mydata_pop[ord], col=mycolors_sq[as.numeric(unlist(sapply(1:12, function(z) infosamples[[z]]$Main.project.LC.Centers)))[samples_imypopi(nspops,i)][ord]],pch=19,cex=0.4)
#vioplot_f(mydata_pop,mycol=1:20,at=i,widthx=0.45,myadd=TRUE,nbreaks=quantile_intervals) 
#plot(1:10,col=mypalette[1:ncolclasses],pch=19)
}
i<-1
mydata_pop<-log(logpop_reschr_nAB_l[[1]][[i]]+0.1)/log(max(logpop_reschr_nAB_l[[1]][[i]]+0.1))
ord<-order(mydata_pop)[-(1:2)]
infosamples[[i]]$Main.project.LC.Centers[ord]
mycolor_variable<-infosamples[[i]]$Total.Exome.Sequence
mycolintervals<-seq(min(mycolor_variable),max(mycolor_variable),(max(mycolor_variable)-min(mycolor_variable))/ncolclasses)
myshades<-sapply(1:(length(quantile_intervals)-1), function(x) mean(mycolor_variable[mydata_pop>quantile(mydata_pop,quantile_intervals)[x] & mydata_pop<quantile(mydata_pop,quantile_intervals)[x+1] ]))
mycols<-sapply(myshades,function(x) sum(x>=mycolintervals))
#print(mycols)
vioplot_f(mydata_pop[ord],mycol=mypalette[mycols],at=0,widthx=0.5,myadd=TRUE,nbreaks=quantile_intervals)
points(0,mean(mydata_pop[1]),col="gray",pch=19,cex=0.4)
#points(rep(0,length(ord)),mydata_pop[ord], col=mycolors_sq[as.numeric(unlist(sapply(1:12, function(z) infosamples[[z]]$Main.project.LC.Centers)))[samples_imypopi(nspops,i)][ord]],pch=19,cex=0.4)
dev.off()
table(print(as.numeric(unlist(sapply(1:12, function(z) infosamples[[z]]$Main.project.LC.Centers)))))
mycolors_sq[unique(as.numeric(unlist(sapply(1:12, function(z) infosamples[[z]]$Main.project.LC.Centers))))]
}
#nAB histograms with gaussian clustering overlayed
{
ipop<-5
popnames[[ipop]]
plot_posterior=FALSE
plot_null=FALSE
mycols=c("red","orange","gold","yellow","yellowgreen")
for (ipop in 1:npopulations){
if (plot_posterior) { pdf(paste0(myfolderdropbox,"/barplot_posterior_",popnames[ipop],".pdf")) } else { pdf(paste0(myfolderdropbox,"/barplot_",popnames[ipop],".pdf")) }
Sys.sleep(1)
print(paste0("pop ",ipop))
    mydata_pop<-logpop_nAB[[ipop]]; 
    mydata_pop<-log(mydata_pop[1:(length(mydata_pop)-2)])
if (plot_null)
    {
    par(mfrow=c(2,2),mar=c(1,1,1,1));
    mydata_pop_null<-log(logpop_reschr_nAB[[ipop]]+1); 
    mydata_pop_null<-mydata_pop_null[1:(length(mydata_pop_null)-2)]
    allpops<-c(mydata_pop_null,mydata_pop)
    mystep<-(max(allpops)-min(allpops))/30
    myticks<-seq(min(allpops),max(allpops),mystep)
    x<-seq(myticks[1],myticks[length(myticks)],(myticks[length(myticks)]-myticks[1])/1000)
    mybarplot_f(myticks,myhist_f(myticks,mydata_pop_null),ticks_are_breaks=TRUE,myxlabel='log(nAB)',myylabel='number of individuals')
    myp<-dnorm(x,mean(mydata_pop_null),sd=sd(mydata_pop_null))
    my_rescale_posteriors<-max(mydata_pop_null)/max(myp)
    lines(x,myp, col=mycols)
    } else 
    {
    allpops<-mydata_pop;
    }
mystep<-(max(allpops)-min(allpops))/30
myticks<-seq(min(allpops),max(allpops),mystep)
x<-seq(myticks[1],myticks[length(myticks)],(myticks[length(myticks)]-myticks[1])/1000)
myhistdata<-myhist_f(myticks,mydata_pop)
mybarplot_f(myticks,myhistdata,ticks_are_breaks=TRUE,myxlabel='log(nAB)',myylabel='number of individuals')
mymodel<-FALSE; class(mymodel) <- "try-error";
mykk<-10
while (class(mymodel) == "try-error" )
{
print(mykk)
if( mykk >= 2 ) { 
mymodel<-try ( partition_samples_nAB(mydata_pop,maxk=mykk))
mykk<-mykk - 1; } else {mymodel=list(k_from_AICc=1,model_fromAICc=list(mu=mean(mydata_pop),sigma=sd(mydata_pop),lambda=1))}
if (class(mymodel) != "try-error"){if (length(mymodel$model_fromAICc)==1){class(mymodel) = "try-error"}}
};
myp<-sapply(1:mymodel$k_from_AICc, function(ik) dnorm(x,mean=mymodel$model_fromAICc$mu[ik],sd=mymodel$model_fromAICc$sigma[ik])*mymodel$model_fromAICc$lambda[ik])
if (plot_posterior) {marginal_sum_probs<-apply(myp,MARGIN=1,sum)} else {marginal_sum_probs<-1}
my_rescale_posteriors<-max(myhistdata)/max(myp/marginal_sum_probs)
for (ik in 1:mymodel$k_from_AICc)
{
lines(x,myp[,ik]*my_rescale_posteriors/marginal_sum_probs,col=rgb(ik/mymodel$k_from_AICc,0.6,1-ik/mymodel$k_from_AICc),lwd=3)
}
#barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
dev.off()
}


}
#combining infoplot (barplots non normalized with log) and clustering
{
add_null_distr<-FALSE
mystep<-0.025
mycolors_sq<-c("gold1","deepskyblue1","coral","darkslategray4","forestgreen","khaki4","tan","lightsalmon","lemonchiffon3","lightblue1","firebrick1","mediumpurple1","mediumaquamarine","red4","seagreen3","azure1","peachpuff","orangered2","limegreen","hotpink3")
mylevels<-unique(unlist(sapply(1:12,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
allpops<-unlist(sapply(1:12,function(x) unlist(log(logpop_nAB[[x]])[1:(length(logpop_nAB[[x]])-2)] )))
if (sum(allpops==-Inf)>0) { 
for ( x in 1:length(logpop_nAB)){logpop_nAB[[x]]<-logpop_nAB[[x]]+1}
}
table(unlist(sapply(1:12,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
allpops<-unlist(sapply(1:12,function(x) unlist(log(logpop_nAB[[x]])[1:(length(logpop_nAB[[x]])-2)] )))
mystep<-(max(allpops)-min(allpops))/30
myticks_bars<-seq(floor(min(allpops-0.02)/mystep)*mystep,floor(min(allpops-0.02)/mystep)*mystep+34*mystep ,mystep)
if (sum(allpops==-Inf)>0) { 
pdf(paste0(myfolderdropbox,"/barplot_logplus1absolutenAB_curves.pdf"))
} else {pdf(paste0(myfolderdropbox,"/barplot_logabsolutenAB_curves.pdf"))}
if (add_null_distr) { 
    par(mfrow=c(4,13),mar=c(0.5,0.1,0.1,0.1));
    layout(matrix(1:(14*4), 4, 14, byrow = TRUE),widths=c(0.5,rep(1,14*4-1)))
    plot.new()
    #layout(matrix(1:24, 2, 12, byrow=TRUE))
    i<-10
    mydata_pop<-logpop_reschr_nAB[[i]]; 
    mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
    mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks_bars))
    barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
    #axis(2,at=1:length(myticks_bars),myticks_bars,cex.axis=0.9)
    text(1*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/2,0.1,paste0("null ",popnames[[i]]))
    box("plot", col="black")  
    } else 
    {
    par(mfrow=c(4,npopulations+1),mar=c(0.5,0.1,0.1,0.1))
    layout(matrix(1:((npopulations+1)*4), 4, (npopulations+1), byrow = TRUE),widths=c(0.5,rep(1,(npopulations+1)*4-1)))
    plot.new()
    }
for (i in 1:12){
mydata_pop<-log(logpop_nAB[[i]])
mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
#ord<-order(mydata_pop)[-(1:2)]
mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks_bars))
barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col=mycolors_sq,yaxt='n',xaxt='n')
text(2*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/3,0.1,popnames[[i]])
#rect( grconvertX(0.005, from='ndc'), grconvertY(0.505, from='ndc'),
#     grconvertX(0.495, from='ndc'), grconvertY(0.995, from='ndc'))
myticks_bars<-seq(floor(min(allpops-0.02)/mystep)*mystep,floor(min(allpops-0.02)/mystep)*mystep+34*mystep ,mystep)
#myticks_curves<-seq(min(allpops),max(allpops),mystep)
x<-seq(myticks_bars[1],myticks_bars[length(myticks_bars)],(myticks_bars[length(myticks_bars)]-myticks_bars[1])/1000)
myhistdata<-myhist_f(myticks_bars,mydata_pop)
#mybarplot_f(myticks_bars,myhistdata,ticks_are_breaks=TRUE,myxlabel='log(nAB)',myylabel='number of individuals')
mymodel<-FALSE; class(mymodel) <- "try-error";
mykk<-10
while (class(mymodel) == "try-error" )
{
print(mykk)
if( mykk >= 2 ) { 
mymodel<-try ( partition_samples_nAB(mydata_pop,maxk=mykk,phtr=0.001))
mykk<-mykk - 1; } else {mymodel=list(k_from_AICc=1,model_fromAICc=list(mu=mean(mydata_pop),sigma=sd(mydata_pop),lambda=1))}
if (class(mymodel) != "try-error"){if (length(mymodel$model_fromAICc)==1){class(mymodel) = "try-error"}}
};
myp<-sapply(1:mymodel$k_from_AICc, function(ik) dnorm(x,mean=mymodel$model_fromAICc$mu[ik],sd=mymodel$model_fromAICc$sigma[ik])*mymodel$model_fromAICc$lambda[ik])
if (plot_posterior) {marginal_sum_probs<-apply(myp,MARGIN=1,sum)} else {marginal_sum_probs<-1}
my_rescale_posteriors<-max(t(mydata_pop_bySC))/max(myp/marginal_sum_probs)
for (ik in 1:mymodel$k_from_AICc)
{
#colored lines
#lines(myp[,ik]*my_rescale_posteriors/marginal_sum_probs,34*(x-min(x))/max(x-min(x)),col=rgb(ik/mymodel$k_from_AICc,0.6,1-ik/mymodel$k_from_AICc),lwd=3)
#black and red lines
lines(myp[,ik]*my_rescale_posteriors/marginal_sum_probs,34*(x-min(x))/max(x-min(x)),col=rep(c("black","red"),5)[ik],lwd=1.5,lty=1)
}
box("plot", col="black")  
}
#box("figure", col="green")  
#box('outer')
dev.off()
}

#combining infoplot (barplots non normalized with log but no colors) and clustering 
{
add_null_distr<-FALSE
mystep<-0.025
mycolors_sq<-c("gold1","deepskyblue1","coral","darkslategray4","forestgreen","khaki4","tan","lightsalmon","lemonchiffon3","lightblue1","firebrick1","mediumpurple1","mediumaquamarine","red4","seagreen3","azure1","peachpuff","orangered2","limegreen","hotpink3")
mylevels<-unique(unlist(sapply(1:12,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
allpops<-unlist(sapply(1:12,function(x) unlist(log(logpop_nAB[[x]])[1:(length(logpop_nAB[[x]])-2)] )))
if (sum(allpops==-Inf)>0) { 
for ( x in 1:length(logpop_nAB)){logpop_nAB[[x]]<-logpop_nAB[[x]]+1}
}
table(unlist(sapply(1:12,function(x) as.character(infosamples[[x]]$Main.project.LC.Centers))))
allpops<-unlist(sapply(1:12,function(x) unlist(log(logpop_nAB[[x]])[1:(length(logpop_nAB[[x]])-2)] )))
mystep<-(max(allpops)-min(allpops))/30
myticks_bars<-seq(floor(min(allpops-0.02)/mystep)*mystep,floor(min(allpops-0.02)/mystep)*mystep+34*mystep ,mystep)
if (sum(allpops==-Inf)>0) { 
pdf(paste0(myfolderdropbox,"/barplot_logplus1absolutenAB_noinfo_curves.pdf"))
} else {pdf(paste0(myfolderdropbox,"/barplot_logabsolutenAB_noinfo_curves.pdf"))}
if (add_null_distr) { 
    par(mfrow=c(4,13),mar=c(0.5,0.1,0.1,0.1));
    layout(matrix(1:(14*4), 4, 14, byrow = TRUE),widths=c(0.5,rep(1,14*4-1)))
    plot.new()
    #layout(matrix(1:24, 2, 12, byrow=TRUE))
    i<-10
    mydata_pop<-logpop_reschr_nAB[[i]]; 
    mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
    mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks_bars))
    mydata_pop_bySC<-apply(mydata_pop_bySC,MARGIN=1,sum)
    barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col='gray94',yaxt='n',xaxt='n')
    #axis(2,at=1:length(myticks_bars),myticks_bars,cex.axis=0.9)
    text(1*(max(apply(mydata_pop_bySC,MARGIN=1,sum))-min(apply(mydata_pop_bySC,MARGIN=1,sum)))/2,0.1,paste0("null ",popnames[[i]]))
    box("plot", col="black")  
    } else 
    {
    par(mfrow=c(4,npopulations+1),mar=c(0.5,0.1,0.1,0.1))
    layout(matrix(1:((npopulations+1)*4), 4, (npopulations+1), byrow = TRUE),widths=c(0.5,rep(1,(npopulations+1)*4-1)))
    plot.new()
    }
for (i in 1:12){
mydata_pop<-log(logpop_nAB[[i]])
mydata_pop<-mydata_pop[1:(length(mydata_pop)-2)]
#ord<-order(mydata_pop)[-(1:2)]
mydata_pop_bySC<-sapply(mylevels,function(x) myhist_f(mydata_pop[as.character(infosamples[[i]]$Main.project.LC.Centers)==x],mybreaks=myticks_bars))
mydata_pop_bySC<-apply(mydata_pop_bySC,MARGIN=1,sum)
barplot(t(mydata_pop_bySC),horiz=TRUE,space=0,col='gray94',yaxt='n',xaxt='n')
text(2*(max(mydata_pop_bySC)-min(mydata_pop_bySC))/3,0.1,popnames[[i]])
#rect( grconvertX(0.005, from='ndc'), grconvertY(0.505, from='ndc'),
#     grconvertX(0.495, from='ndc'), grconvertY(0.995, from='ndc'))
myticks_bars<-seq(floor(min(allpops-0.02)/mystep)*mystep,floor(min(allpops-0.02)/mystep)*mystep+34*mystep ,mystep)
#myticks_curves<-seq(min(allpops),max(allpops),mystep)
x<-seq(myticks_bars[1],myticks_bars[length(myticks_bars)],(myticks_bars[length(myticks_bars)]-myticks_bars[1])/1000)
myhistdata<-myhist_f(myticks_bars,mydata_pop)
#mybarplot_f(myticks_bars,myhistdata,ticks_are_breaks=TRUE,myxlabel='log(nAB)',myylabel='number of individuals')
mymodel<-FALSE; class(mymodel) <- "try-error";
mykk<-10
while (class(mymodel) == "try-error" )
{
print(mykk)
if( mykk >= 2 ) { 
mymodel<-try ( partition_samples_nAB(mydata_pop,maxk=mykk,phtr=0.001))
mykk<-mykk - 1; } else {mymodel=list(k_from_AICc=1,model_fromAICc=list(mu=mean(mydata_pop),sigma=sd(mydata_pop),lambda=1))}
if (class(mymodel) != "try-error"){if (length(mymodel$model_fromAICc)==1){class(mymodel) = "try-error"}}
};
myp<-sapply(1:mymodel$k_from_AICc, function(ik) dnorm(x,mean=mymodel$model_fromAICc$mu[ik],sd=mymodel$model_fromAICc$sigma[ik])*mymodel$model_fromAICc$lambda[ik])
if (plot_posterior) {marginal_sum_probs<-apply(myp,MARGIN=1,sum)} else {marginal_sum_probs<-1}
my_rescale_posteriors<-max(t(mydata_pop_bySC))/max(myp/marginal_sum_probs)
for (ik in 1:mymodel$k_from_AICc)
{
#colored lines
#lines(myp[,ik]*my_rescale_posteriors/marginal_sum_probs,34*(x-min(x))/max(x-min(x)),col=rgb(ik/mymodel$k_from_AICc,0.6,1-ik/mymodel$k_from_AICc),lwd=3)
#black and red lines
lines(myp[,ik]*my_rescale_posteriors/marginal_sum_probs,34*(x-min(x))/max(x-min(x)),col=rep(c("black","red"),5)[ik],lwd=1.5,lty=1)
}
box("plot", col="black")  
}
#box("figure", col="green")  
#box('outer')
dev.off()
}


#explorations prior to infoplot
if (FALSE)
{
#> sum(logpop_reschr[[ipop+1]])
#[1] 168814
#> sum(logpop_nAB[[ipop+1]])
#[1] 134130.5
#> ipop
#[1] 11
#> ipop<-0
#> sum(logpop_nAB[[ipop+1]])
#[1] 190377
#> sum(logpop_reschr[[ipop+1]])
#[1] 174556.5


load(paste0(myfolder,"/logpop_nAB_pos.RData"))
load(paste0(myfolder,"/logpop_nAB_neg.RData"))
load(paste0(myfolder,"/logpop_nAB_neg_reschr",myperm,"_nAB.RData"))
load(paste0(myfolder,"/logpop_nAB_pos_reschr",myperm,"_nAB.RData"))
load(paste0(myfolder,"/mylogpop_neg_reschr",myperm,".RData"))
load(paste0(myfolder,"/mylogpop_pos_reschr",myperm,".RData"))

require(vioplot)
load(paste0(myfolder,"/logpop_nAB.RData"))
logpop_reschr_nAB_l<-list()
for (z in 1:nperm) { load(paste0(myfolder,"/logpop_reschr",z,"_nAB.RData")); logpop_reschr_nAB_l[[z]]<-logpop_reschr_nAB; }


pdf("~/Dropbox/LDLD/temp.pdf")
plot(nAB_l[ord],pch=19,xlab="YRI samples",ylab="nAB")
dev.off()
myreschrnAB<-sapply(1:nperm, function(x) {nAB_lr<-logpop_reschr_nAB_l[[x]][[ipop]]; ordr<-order(nAB_lr[1:length(nAB_lr)]);return(nAB_lr[ordr])})
myreschr_mean<-apply(myreschrnAB,MARGIN=1,mean)
myreschr_sd<-apply(myreschrnAB,MARGIN=1,sd)
reshufflings<-sapply(1:length(myreschr_mean), function(z) rnorm(mynrepl,myreschr_mean[z],myreschr_sd[z]))
cdat <- as.list(as.data.frame(reshufflings))
names(cdat)[1] <- "x"  # vioplot() needs the first element to be called 'x'
pdf("~/Dropbox/LDLD/temp.pdf")
do.call(vioplot,c(cdat,list(col="cadetblue",colMed="cadetblue")))
do.call(vioplot,c(cdat,list(col="cadetblue",colMed="cadetblue",add=TRUE)))
dev.off()
nAB_lr<-logpop_reschr_nAB[[ipop]]
ordr<-order(nAB_lr[1:length(nAB_lr)])
plot(1:length(nAB_l),nAB_l[1:length(nAB_l)][ord]/max(nAB_l))
points(1:length(nAB_lr),nAB_lr[1:length(nAB_lr)][ordr]/max(nAB_l),col="red")

length(logpop_reschr_nAB_l[[1]][[1]])

mean(logpop_reschr_nAB_l[[1]][[1]])
sd(logpop_reschr_nAB_l[[1]][[1]])
ipop<-5;nperm<-2
logpop_reschr_nAB_l_shufflefake<-logpop_reschr_nAB_l
for (myperm in 1:nperm)
{
logpop_reschr_nAB_l_shufflefake[[myperm]][[ipop]]<-rnorm(length(logpop_reschr_nAB_l[[myperm]][[ipop]]),mean(logpop_reschr_nAB_l[[myperm]][[ipop]]),sd(logpop_reschr_nAB_l[[myperm]][[ipop]])*2.5)
}

pdf(paste0(myfolderdropbox,"/infoplot.fake.",ipop,".pdf"))
infoplot(ipop,mymaxk=6,mypthr=0.05,mynABdefault=FALSE,mynAB=logpop_nAB,showreal_confidence_interval=TRUE,logpop_reschr_nAB_l=logpop_reschr_nAB_l_shufflefake,reshufflings=0,mylog=FALSE,showmutperssample=FALSE,mypopname=popnames[ipop],finite_mixture_criterion=2) #version with real reshuffling of chromosomes
dev.off()
}

#infoplots with full information
if (FALSE)
{

#all nAB #ipop=1
for (ipop in 1:12) 
{
pdf(paste0(myfolderdropbox,"/infoplot.",ipop,".pdf"))
infoplot(ipop,mymaxk=6,mypthr=0.05,mynABdefault=FALSE,mynAB=logpop_nAB,showreal_confidence_interval=TRUE,logpop_reschr_nAB_l=logpop_reschr_nAB_l,reshufflings=0,mylog=FALSE,showmutperssample=FALSE,mypopname=popnames[ipop],finite_mixture_criterion=2,continuous_fields=c("X..Targets.Covered.to.20x.or.greater","Total.LC.Sequence","LC.Non.Duplicated.Aligned.Coverage")) #version with real reshuffling of chromosomes
dev.off()
}

#panel pos
my.locations<-rep("topleft",6)
my.labels<-c("(a)","(b)","(c)","(d)","(e)","(f)")
mypops_plot<-1:6#1:6
pdf(paste0(myfolderdropbox,"/infoplot_pos_panel_",paste0(mypops_plot,collapse=''),".pdf"))
par(mfrow=c(3,2))
par(cex = 0.5)
#par(mar = c(3, 2, 1, 1), oma = c(5, 5, 2, 2))
par(mar = c(4, 4, 2, 2), oma = c(4, 4, 2, 2))
counter<-1
for (ipop in mypops_plot) #7:12 
{
infoplot(ipop,mymaxk=6,mypthr=0.05,mynABdefault=FALSE,mynAB=logpop_nAB_pos,showreal_confidence_interval=TRUE,logpop_reschr_nAB_l=logpop_reschr_nAB_pos_l,reshufflings=0,mylog=FALSE,showmutperssample=FALSE,mypopname=popnames[ipop],finite_mixture_criterion=2,myylab="            nAB          ") #version with real reshuffling of chromosomes
put.fig.letter(label=my.labels[counter], location=my.locations[counter], font=2);counter<-counter+1
}
dev.off()
#panel neg
my.locations<-rep("topleft",6)
my.labels<-c("(a)","(b)","(c)","(d)","(e)","(f)")
mypops_plot<-7:12#1:6#7:12
paste0(mypops_plot,collapse='')
pdf(paste0(myfolderdropbox,"/infoplot_neg_panel_",paste0(mypops_plot,collapse=''),".pdf"))
par(mfrow=c(3,2))
par(cex = 0.5)
#par(mar = c(3, 2, 1, 1), oma = c(5, 5, 2, 2))
par(mar = c(4, 4, 2, 2), oma = c(4, 4, 2, 2))
counter<-1
for (ipop in mypops_plot) #7:12 
{
infoplot(ipop,mymaxk=6,mypthr=0.05,mynABdefault=FALSE,mynAB=logpop_nAB_neg,showreal_confidence_interval=TRUE,logpop_reschr_nAB_l=logpop_reschr_nAB_neg_l,reshufflings=0,mylog=FALSE,showmutperssample=FALSE,mypopname=popnames[ipop],finite_mixture_criterion=2,myylab="            nAB          ") #version with real reshuffling of chromosomes
put.fig.letter(label=my.labels[counter], location=my.locations[counter], font=2);counter<-counter+1
}
dev.off()


#pos
for (ipop in 1:12) 
{
pdf(paste0(myfolderdropbox,"/infoplot_pos.",ipop,".pdf"))
infoplot(ipop,mymaxk=6,mypthr=0.05,mynABdefault=FALSE,mynAB=logpop_nAB_pos,showreal_confidence_interval=TRUE,logpop_reschr_nAB_l=logpop_reschr_nAB_pos_l,reshufflings=0,mylog=FALSE,showmutperssample=FALSE,mypopname=popnames[ipop],finite_mixture_criterion=2) #version with real reshuffling of chromosomes
dev.off()
}
#neg
for (ipop in 1:12) #neg
{
pdf(paste0(myfolderdropbox,"/infoplot_neg.",ipop,".pdf"))
infoplot(ipop,mymaxk=6,mypthr=0.05,mynABdefault=FALSE,mynAB=logpop_nAB_neg,showreal_confidence_interval=TRUE,logpop_reschr_nAB_l=logpop_reschr_nAB_neg_l,reshufflings=0,mylog=FALSE,showmutperssample=FALSE,mypopname=popnames[ipop],finite_mixture_criterion=2) #version with real reshuffling of chromosomes
dev.off()
}


#nAB plots with nA in gray
for (ipop in 1:12)
{
pdf(paste0(myfolderdropbox,"/infoplot_mut.",ipop,".pdf"))
infoplot(ipop,mymaxk=6,mypthr=0.05,mynABdefault=FALSE,mynAB=logpop_nAB,showreal_confidence_interval=TRUE,logpop_reschr_nAB_l=logpop_reschr_nAB_l,reshufflings=0,mylog=FALSE,showmutperssample=TRUE,mypopname=popnames[ipop],finite_mixture_criterion=2,mutpersample=l_mutsamples[[ipop]]) #version with real reshuffling of chromosomes
dev.off()
}

#nA plots with in nAB gray
for (ipop in 1:12)
{
pdf(paste0(myfolderdropbox,"/infoplot_nA.",ipop,".pdf"))
infoplot(ipop,mymaxk=6,mypthr=0.05,mynABdefault=FALSE,mynAB=l_mutsamples,showreal_confidence_interval=FALSE,logpop_reschr_nAB_l=0,reshufflings=0,mylog=FALSE,showmutperssample=TRUE,mypopname=popnames[ipop],finite_mixture_criterion=2,mutpersample=logpop_nAB[[ipop]],myylab="nA") #version with real reshuffling of chromosomes
dev.off()
}

}


}
#test effects of variables
if (FALSE)
{

unlist_info<-function(myinfo) { res<-as.character(unlist(sapply(1:12,function(i) infosamples[[i]][[myinfo]])));if (length(is.na(res))>1){ res[is.na(res)]<-0;res[res==""]<-0};return(res) }
#random intercepts and slopes (1|<name>)+(0+<cov_1>|<name>)+(0+<cov_2>|<name>) 
library(lme4)
contr=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000),calc.derivs=F)

#only when need to model it
#Main.Project..E.Platform: all illumina

#for EXOME
mydata_model<-data.frame(nAB=scale(log(as.numeric(unlist(logpop_nAB_pos)))),Population=as.factor(unlist_info("Population")),LC.Pilot.Platforms=as.factor(unlist_info("LC.Pilot.Platforms")),LC.Pilot.Centers=as.factor(unlist_info("LC.Pilot.Centers")),
In.High.Coverage.Pilot=as.factor(unlist_info("In.High.Coverage.Pilot")),HC.Pilot.Centers=as.factor(unlist_info("HC.Pilot.Centers")),
In.Exon.Targetted.Pilot=as.factor(unlist_info("In.Exon.Targetted.Pilot")),ET.Pilot.Platforms=as.factor(unlist_info("ET.Pilot.Platforms")),ET.Pilot.Centers=as.factor(unlist_info("ET.Pilot.Centers")),
Has.Sequence.in.Phase1=as.factor(unlist_info("Has.Sequence.in.Phase1")),Main.project.LC.Centers=as.factor(unlist_info("Main.project.LC.Centers")),
Main.project.LC.platform=as.factor(unlist_info("Main.project.LC.platform")),Total.LC.Sequence=scale(as.numeric(unlist_info("Total.LC.Sequence"))),LC.Non.Duplicated.Aligned.Coverage=scale(as.numeric(unlist_info("LC.Non.Duplicated.Aligned.Coverage"))),
Main.Project.E.Centers=as.factor(unlist_info("Main.Project.E.Centers")),Main.Project..E.Platform=as.factor(unlist_info("Main.Project..E.Platform")),Total.Exome.Sequence=scale(as.numeric(unlist_info("Total.Exome.Sequence"))),X..Targets.Covered.to.20x.or.greater=scale(as.numeric(unlist_info("X..Targets.Covered.to.20x.or.greater"))),
Has.Omni.Genotypes=as.factor(unlist_info("Has.Omni.Genotypes")),Has.Axiom.Genotypes=as.factor(unlist_info("Has.Axiom.Genotypes")),Has.Affy.6.0.Genotypes=as.factor(unlist_info("Has.Affy.6.0.Genotypes")))

res_Ec=lmer(nAB ~  Main.Project.E.Centers+(1+Main.Project.E.Centers|Population),data=mydata_model,REML=FALSE,control=contr)
res_allcov_Ec=lmer(nAB ~ Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.Project.E.Centers+ (1+Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.Project.E.Centers|Population),data=mydata_model,control=contr,REML=F)
res_covnoE_Ec=lmer(nAB ~  X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.Project.E.Centers+(1+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.Project.E.Centers|Population),data=mydata_model,control=contr,REML=F)
res_allcovnoL_Ec=lmer(nAB ~  Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+Main.Project.E.Centers+(1+Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+Main.Project.E.Centers|Population),data=mydata_model,control=contr,REML=F)
res_allcovnoT_Ec=lmer(nAB ~  X..Targets.Covered.to.20x.or.greater+Main.Project.E.Centers+(1+X..Targets.Covered.to.20x.or.greater+Main.Project.E.Centers|Population),data=mydata_model,control=contr,REML=F)
res_allcov_Ec_chips=lmer(nAB ~  Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.Project.E.Centers+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes+(1+Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.Project.E.Centers+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes|Population),data=mydata_model,control=contr,REML=F)
res_allcov_chips=lmer(nAB ~  Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes+(1+Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes|Population),data=mydata_model,control=contr,REML=F)
res_EC_chips=lmer(nAB ~  Main.Project.E.Centers+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes+(1+Main.Project.E.Centers+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes|Population),data=mydata_model,control=contr,REML=F)
summary(res_allcov_Ec_chips)$coefficients
anova(res_allcov_Ec_chips, res_allcov_Ec, test="Chisq")
anova(res_allcov_Ec_chips, res_allcov_chips, test="Chisq")
anova(res_allcov_Ec_chips, res_EC_chips, test="Chisq")


#for INTERGENIC - only difference is that Main.project.LC.Centers rather than Main.Project.E.Centers
mydata_model<-data.frame(nAB=scale(log(as.numeric(unlist(logpop_nAB_pos20)))),Population=as.factor(unlist_info("Population")),LC.Pilot.Platforms=as.factor(unlist_info("LC.Pilot.Platforms")),LC.Pilot.Centers=as.factor(unlist_info("LC.Pilot.Centers")),
In.High.Coverage.Pilot=as.factor(unlist_info("In.High.Coverage.Pilot")),HC.Pilot.Centers=as.factor(unlist_info("HC.Pilot.Centers")),
In.Exon.Targetted.Pilot=as.factor(unlist_info("In.Exon.Targetted.Pilot")),ET.Pilot.Platforms=as.factor(unlist_info("ET.Pilot.Platforms")),ET.Pilot.Centers=as.factor(unlist_info("ET.Pilot.Centers")),
Has.Sequence.in.Phase1=as.factor(unlist_info("Has.Sequence.in.Phase1")),Main.project.LC.Centers=as.factor(unlist_info("Main.project.LC.Centers")),
Main.project.LC.platform=as.factor(unlist_info("Main.project.LC.platform")),Total.LC.Sequence=scale(as.numeric(unlist_info("Total.LC.Sequence"))),LC.Non.Duplicated.Aligned.Coverage=scale(as.numeric(unlist_info("LC.Non.Duplicated.Aligned.Coverage"))),
Main.Project.E.Centers=as.factor(unlist_info("Main.Project.E.Centers")),Main.Project..E.Platform=as.factor(unlist_info("Main.Project..E.Platform")),Total.Exome.Sequence=scale(as.numeric(unlist_info("Total.Exome.Sequence"))),X..Targets.Covered.to.20x.or.greater=scale(as.numeric(unlist_info("X..Targets.Covered.to.20x.or.greater"))),
Has.Omni.Genotypes=as.factor(unlist_info("Has.Omni.Genotypes")),Has.Axiom.Genotypes=as.factor(unlist_info("Has.Axiom.Genotypes")),Has.Affy.6.0.Genotypes=as.factor(unlist_info("Has.Affy.6.0.Genotypes")))
res_Ec=lmer(nAB ~  Main.project.LC.Centers+ (1+Main.project.LC.Centers|Population),data=mydata_model,REML=FALSE,control=contr)
res_allcov_Ec=lmer(nAB ~Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.project.LC.Centers+  (1+Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.project.LC.Centers|Population),data=mydata_model,control=contr,REML=F)
res_covnoE_Ec=lmer(nAB ~  X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.project.LC.Centers+(1+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.project.LC.Centers|Population),data=mydata_model,control=contr,REML=F)
res_allcovnoL_Ec=lmer(nAB ~  Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+Main.project.LC.Centers+(1+Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+Main.project.LC.Centers|Population),data=mydata_model,control=contr,REML=F)
res_allcovnoT_Ec=lmer(nAB ~  X..Targets.Covered.to.20x.or.greater+Main.project.LC.Centers+(1+X..Targets.Covered.to.20x.or.greater+Main.project.LC.Centers|Population),data=mydata_model,control=contr,REML=F)
res_allcov_Ec_chips=lmer(nAB ~  Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.project.LC.Centers+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes+(1+Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Main.project.LC.Centers+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes|Population),data=mydata_model,control=contr,REML=F)
res_allcov_chips=lmer(nAB ~  Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes+(1+Total.Exome.Sequence+X..Targets.Covered.to.20x.or.greater+LC.Non.Duplicated.Aligned.Coverage+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes|Population),data=mydata_model,control=contr,REML=F)
res_EC_chips=lmer(nAB ~  Main.project.LC.Centers+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes+(1+Main.project.LC.Centers+Has.Affy.6.0.Genotypes+Has.Axiom.Genotypes+Has.Affy.6.0.Genotypes|Population),data=mydata_model,control=contr,REML=F)

summary(res_allcov_Ec_chips)$coefficients
anova(res_allcov_Ec_chips, res_allcov_Ec, test="Chisq")
anova(res_allcov_Ec_chips, res_allcov_chips, test="Chisq")
anova(res_allcov_Ec_chips, res_EC_chips, test="Chisq")
#outputs save in ~/Dropbox/LDLD/ms/Ajhg/table_model_namesSeqCenters.odt
}

#reshuffle lognAB
{
for (i in 1:10){
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
logpop_nAB_pos<-lapply(1:npopulations, function(x) sample(logpop_nAB_pos[[x]]))
save(logpop_nAB_pos,file=paste0(myfolder,"/logpop_nAB_pos_sim",i,".RData"))
}
}
#verify correlation nAB across iterates
if (FALSE)
{
npopulations<-12
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/introns.1000g/repl3/min1"
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
intro3<-sapply(logpop_nAB_pos, function(i) i[1:(length(i)-2)])
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/introns.1000g/repl2/min1"
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
intro2<-sapply(logpop_nAB_pos, function(i) i[1:(length(i)-2)])
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl3/min1"
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
inter3<-sapply(logpop_nAB_pos, function(i) i[1:(length(i)-2)])
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl2/min1"
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
inter2<-sapply(logpop_nAB_pos, function(i) i[1:(length(i)-2)])
#true overall, but intergenic and intronic separated
cor.test(unname(unlist(inter2)),unname(unlist(inter3))) #p-value = 7.565e-07
cor.test(unname(unlist(intro2)),unname(unlist(intro3))) #p-value < 2.2e-16
cor.test(unname(unlist(inter2)),unname(unlist(intro2))) #p-value = 0.05454
cor.test(unname(unlist(inter3)),unname(unlist(intro3))) #p-value = 0.1706
#check if true even by population
sapply(1:npopulations, function(i) cor.test(inter2[[i]],inter3[[i]])$p.value) #5 out of 12 correlated but not highly, no pattern
sapply(1:npopulations, function(i) cor.test(intro2[[i]],intro3[[i]])$p.value) #all highly correlated except for STU and GWD (probably effect of exome) 
intergenic_nAB<-list();introns_nAB<-list()
for (ipop in 1:npopulations) {intergenic_nAB[[ipop]]<-inter2[[ipop]]+inter3[[ipop]]}
for (ipop in 1:npopulations) {introns_nAB[[ipop]]<-intro2[[ipop]]+intro3[[ipop]]}
#unname(unlist(introns_nAB))
#unname(unlist(intergenic_nAB))
save(intergenic_nAB,file=paste0("/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/intergenic_nAB.RData"))
save(introns_nAB,file=paste0("/mnt/scratch/fabrizio/LDLD/above95/introns.1000g/introns_nAB.RData"))

  mychr<-22
  myorderedsamples_file <-"~/workspace/1000genomes/above95.unrelated.samples"
  mysource_file <-paste0("/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr",mychr,".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz")
  npops_file<-"/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/nspops.txt"
  targetsamples<-system(paste0("cat ",myorderedsamples_file),intern=TRUE) #previously called samplesabove95
  sourcesamples<-system(paste0("zcat ",mysource_file," | head -1500 | grep '#CHROM' | tail -1 | awk '{for (i=10; i<=NF; i++) print $i}'"),intern=TRUE) #previously called samples1000g
#reordered_samples<-paste0(intersect(targetsamples,sourcesamples),collapse=',') #previously called samples_order1000g_above95
l_npops<-as.numeric(system(paste0("cat ",npops_file),intern=TRUE))
reordered_samples<-intersect(targetsamples,sourcesamples)
cmd<-paste0(c(1:9,match(reordered_samples,sourcesamples)+9),collapse=',$')
cmd<-paste0("grep -v '##' | awk -v OFS='\t' '{print $",cmd,"}' ")
cmd2<-paste0("zcat ",mysource_file," |",cmd," | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/nspops.txt 0.01 12 0 | awk '{for (i=1;i<=NF;i++){if ((i<=2) || (i>=9)){printf \"%s\\t\",$i}};printf \"\n\"}' | sed 's/|/\t/g'")

cmd2<-paste0("zcat ",mysource_file," |",cmd," | sed 's/|/\t/g'| less")
cmd<-paste0("zcat ",mysource_file, "| grep -v '##' | awk -v OFS='\t' '{print $",cmd,"}' >> ",mydestination_file)

}



#save.image("~/Dropbox/LDLD/temp.RData")
#test if nAB explains differences in mutational load
if (FALSE)
{
load(paste0(myfolder,"/logpop_reschr",myperm,"_nAB.RData"))
load(paste0(myfolder,"/mylogpop_reschr",myperm,".RData"))
logpop_reschr_nAB_l<-list();for (z in 1:nperm) { load(paste0(myfolder,"/logpop_reschr",z,"_nAB.RData")); logpop_reschr_nAB_l[[z]]<-logpop_reschr_nAB; }

cor_nAvsnAB_reschr<-lapply(1:nperm, function(myperm) lapply(1:12,function(x) cor.test(l_mutsamples_reschr[[myperm]][[x]],logpop_reschr_nAB_l[[myperm]][[x]]))) #never significant
cor_nAvsnAB<-lapply(1:12,function(x) cor.test(l_mutsamples[[x]],logpop_nAB[[x]]))

sapply(1:12,function(i) cor_nAvsnAB[[i]]$p.value)
sapply(1:12,function(i) cor_nAvsnAB_reschr[[myperm]][[i]]$p.value)

#weird, there is more correlation
pdf(paste0(myfolderdropbox,"/nAvsnAV.pdf"))
#l_mutsamples_reschr
par(mfrow=c(1,2))
plot(rbind(c(0,0.9),c(3.1,1.1)),pch=19,type="n",xlab="nA",ylab="nAB",)
for (i in 1:12)
{
points(logpop_nAB[[i]]/mean(logpop_nAB[[i]]),l_mutsamples[[i]]/mean(l_mutsamples[[i]]),pch=20,col=rainbow(12)[i])
reg1 <- lm(unname(l_mutsamples[[i]]/mean(l_mutsamples[[i]]))~unname(logpop_nAB[[i]]/mean(logpop_nAB[[i]])))
mytest<-cor.test(l_mutsamples[[i]]/mean(l_mutsamples[[i]]),logpop_nAB[[i]]/mean(logpop_nAB[[i]]))
abline(reg1,lwd=1,col=rainbow(12)[i],lty=1 )
}
myperm<-1
plot(rbind(c(0.5,0.9),c(1.5,1.1)),pch=19,type="n",xlab="nA",ylab="nAB",)
for (i in 1:12)
{
points(logpop_reschr_nAB_l[[myperm]][[i]]/mean(logpop_reschr_nAB_l[[myperm]][[i]]),l_mutsamples_reschr[[myperm]][[i]]/mean(l_mutsamples_reschr[[myperm]][[i]]),pch=20,col=rainbow(12)[i])
reg1 <- lm(unname(l_mutsamples_reschr[[myperm]][[i]]/mean(l_mutsamples_reschr[[myperm]][[i]]))~unname(logpop_reschr_nAB_l[[myperm]][[i]]/mean(logpop_reschr_nAB_l[[myperm]][[i]])))
abline(reg1,lwd=1,col=rainbow(12)[i],lty=1 )
}
dev.off()
gdriveLDLD=0B3oqoLi0PJOlbUZOVEp2UURpcHc


mytest<-cor.test(mymutsamples_rel,mylogpop_rel)

}
#test if nAB inhomogeneous
if (FALSE)
{


load(paste0(myfolder,"/logpop_reschr",1,"_nAB.RData"))
load(paste0(myfolder,"/logpop_nAB.RData"))
sapply(1:12, function(x) var(logpop_nAB[[x]])/var(logpop_reschr_nAB[[x]]))

res<-sapply(1:12, function(x) var(logpop_nAB[[x]])/var(logpop_reschr_nAB[[x]]))
mean(res[res<900]) #52.1 86.2 161.6719
mean(c(52.1, 86.2, 121.6719))
res1<-sapply(1:12, function(x) var(logpop_nAB[[x]]))
res0<-sapply(1:12, function(x) var(logpop_reschr_nAB[[x]]))
wilcox.test(res1,res0)
#/var(logpop_reschr_nAB[[x]])


load(paste0(myfolder,"/logpop_nAB_nlinks.RData"))
load(paste0(myfolder,"/logpop_reschr",myperm,"_nAB.RData"))
load(paste0(myfolder,"/mylogpop_reschr",myperm,".RData"))
logpop_reschr_nAB_l<-list();for (z in 1:nperm) { load(paste0(myfolder,"/logpop_reschr",z,"_nAB.RData")); logpop_reschr_nAB_l[[z]]<-logpop_reschr_nAB; }

pvalues_nAB<-test_inhomogeneous_nAB(logpop_nAB,logpop_reschr_nAB_l)
#> pvalues_nAB$AICc_relativelik_t.test
# [1] 0.018248843278179695 0.013753703574454275 0.001054902690828902
# [4] 0.000000000007970619 0.025873817767663655 0.005963726755546480
# [7] 0.015156462811308773 0.003462870156458958 0.013026951614287258
#[10] 0.021476524776942826 0.016003271369754338 0.263458266128182350
#> pvalues_nAB$nAB_sd_t.test
# [1] 0.0033317057 0.0033242548 0.0005135633 0.0275048145 0.0003884504
# [6] 0.0001859795 0.0007861456 0.0041047774 0.0045937317 0.0006946674
# [11] 0.0005388422 0.0651578585
}
#correlation with mutational load at rare sites
{
doubletonsfile="~/Dropbox/LDLD/1000genomes/mutload/maxdoubletons.gw.tsv"
myfolder="/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5"
load(paste0(myfolder,"/logpop_nAB_pos20.RData"))
mylognAB<-logpop_nAB_pos20

myfolder="/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1"
load(paste0(myfolder,"/logpop_nAB_pos20.RData"))
mylognAB<-logpop_nAB_pos20

myfolder="/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1"
load(paste0(myfolder,"/logpop_nAB.RData"))
mylognAB<-logpop_nAB

doubletonsfile="~/Dropbox/LDLD/1000genomes/mutload/maxdoubletons.exons.tsv"
myfolder="/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5"
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
mylognAB<-logpop_nAB_pos

mydata<-read.table(doubletonsfile)
mydata<-apply(mydata,MARGIN=2,sum)

init<-1;end<-0
nspops<-sapply(mylognAB,function(x) length(x))
mynA<-list()
for (ipop in 1:length(nspops)){
end<-end+nspops[ipop]
print(c(init,end))
mynA[[ipop]]<-mydata[init:end]
init<-init+nspops[ipop]
}

sapply(1:length(nspops), function(ipop) cor.test(mylognAB[[ipop]],mynA[[ipop]]))

}

#=======================
#done by taking significance threshold 
#=======================
if (FALSE)
{
load("LDLD170216.RData")
load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/logpoptot_fromfilelog.RData")
par(mfrow=c(1,1))
if (FALSE) #nAB explains differences in mutational load
{
length(l_mutsamples[[1]])
mylogpop_rel<-c()
mymutsamples_rel<-c()
for (i in 1:12)
{
mylogpop_rel<-c(mylogpop_rel,logpop[[i]][[1]]/mean(logpop[[i]][[1]]))
mymutsamples_rel<-c(mymutsamples_rel,l_mutsamples[[i]]/mean(l_mutsamples[[i]]))
}
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/nABvsnA.pdf")
mytest<-cor.test(mymutsamples_rel,mylogpop_rel)
plot(mylogpop_rel,mymutsamples_rel,pch=19,col="darkblue",xlab="nAB",ylab="nA")
reg1 <- lm(mymutsamples_rel~mylogpop_rel)
abline(reg1,lwd=2,col="red" )
text(4,1.1,labels=paste0("r=",mytest$estimate))
text(4,1.07,labels=paste0("p.value=",mytest$p.value))
dev.off()

length(l_mutsamples[[1]])
mylogpop_rel<-c()
mymutsamples_rel<-c()

#always increase by definition if I increase exponent, and anyway small difference, so it does not really make sense changing everything to have it squared.
#anyway check again being sure that mutations are global per individual while nAB only from significant, otherwise a bit trivial that for those there are more.
i<-2
cor.test(logpop[[i]][[1]],l_mutsamples[[i]])
cor.test(logpop[[i]][[1]],l_mutsamples[[i]]^2) 

for (i in 1:12)
{
mylogpop_rel<-c(mylogpop_rel,logpop[[i]][[1]]/mean(logpop[[i]][[1]]))
mymutsamples_rel<-c(mymutsamples_rel,l_mutsamples[[i]]/mean(l_mutsamples[[i]]))
}



#Furthermore as well as more variance mutational load in pops with more linkage uneveness!!!
sd(l_mutsamples[[3]])/mean(l_mutsamples[[3]])
sd(l_mutsamples[[5]])/mean(l_mutsamples[[5]])
sd(l_mutsamples[[1]])/mean(l_mutsamples[[1]])
sd(l_mutsamples[[8]])/mean(l_mutsamples[[8]])
sd(l_mutsamples[[4]])/mean(l_mutsamples[[4]])
sd(l_mutsamples[[12]])/mean(l_mutsamples[[12]])
}

#Developing null distribution for nAB clustering (OBSOLETE because now not in blocks anymore)-------------------------------v
if (FALSE)
{
#-Thoughts on block separations:
#Other solution is using pairs of chromosome as independent blocks. This is conservative in finding differences between distributions. 
#However, it should still have the power.
#Also, should I consider the contribution like I do now, just sum. Well, it depends on the way of permuting it. If this way it is fine.
#see it like that however: when many pairs most likely I have enough power to consider even chromosomes as single blocks.
#when lower amount of pairs, separating becomes meaningful.

#elements in contingency tables (logpoptot)
mylog<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/logs/all.pop0.log")) #any pop is the same->to fix

myblocks<-create_blocks(mylog)
permbychr<-permute_chrAB_byblock(mylog,myblocks,100000)
myblocks2<-create_blocks(mylog,space_between_indep_blocks=1000000)
permbychrblocks2<-permute_chrAB_byblock(mylog,myblocks2,100000)
#save.image("permbychrblocks2.RData")
str(permbychrblocks)
sum(permbychrblocks2[,1]) #428177
sum(permbychr[,1]) #428177
var(permbychr[,1]) #137336.6
var(permbychrblocks2[,1]) #9129.772
#as expected less variance in fragmenting more

setwd("/mnt/scratch/fabrizio/LDLD")
load("nulldistr1iteration.RData")
load("permbychrblocks2.RData")
load("permbychr.RData")
load("logpop_fromfilelog.RData")
#always reload last version of functions in analysesLDLD_header.R after loading workspace images
source("~/Dropbox/LDLD/analysesLDLD_header.R") 

#ok, so let's say that in case of this population we want to go forward in removing individuals

#system.time(permbychr<-permute_chrAB(mylog,10000)) #   user   system  elapsed 3885.680    0.410 3886.155
#save.image("permbychr.RData")
head(permbychr)
for (i in 1:10000)
{permbychr[,i]<-permbychr[,i][order(permbychr[,i])]}
for (i in 1:10000)
{permbychrblocks2[,i]<-permbychrblocks2[,i][order(permbychrblocks2[,i])]}

mydata<-permbychr[,1]
sum(permbychr[,10000]) #428177
#system.time(nulldistr1iteration<-generate_nulldistr_mixtures(permbychr,1000)) #time 448.160
#save(nulldistr1iteration,file="nulldistr1iteration.RData")
head(as.numeric(nulldistr1iteration[,1]))
head(as.numeric(nulldistr1iteration[,2]))
head(as.numeric(nulldistr1iteration[,3]))
head(as.numeric(nulldistr1iteration[,4]))

#compare value with distribution
setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
mypop<-partition_samples_nAB(logpop[[1]][[1]])[c(1,4,7,8)]
setwd("/mnt/scratch/fabrizio/LDLD")
sum(as.numeric(as.numeric(nulldistr1iteration[,1])>=mypop[1]))/1000 # 0.053 it seems that by taking only number of groups, both with LR and AIC I don't even get significant, although almost
sum(as.numeric(as.numeric(nulldistr1iteration[,2])>=mypop[2]))/1000 # 0.115
sum(as.numeric(as.numeric(nulldistr1iteration[,3])<=mypop[3]))/1000 # 0 instead by taking p-values or RL I get extremely significant.
sum(as.numeric(as.numeric(nulldistr1iteration[,4])<=mypop[4]))/1000 #0

system.time(testfull<-test_nulldistr_mixtures(12,1000,logpoptot,logpop,myblocks))
#save(testfull,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/testfull.RData")
}
load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/testfull.RData")
#plot tests and infoplots
if (FALSE)
{
sapply(1:12,function(x) testfull[[x]]$pval_sd)
sapply(1:12,function(x) testfull[[x]]$pval_nLR)
sapply(1:12,function(x) testfull[[x]]$pval_nAIC)
sapply(1:12,function(x) testfull[[x]]$pval_pLR)
sapply(1:12,function(x) testfull[[x]]$pval_pAIC)


library(vioplot)
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/empiricaldistrperm_pop1.pdf")
par(mfrow=c(2,2))
permvar<-sapply(1:10000,function(x) var(permbychr[,x]))
hist(log(permvar),xlim=c(log(min(permvar)),log(max(max(permvar),var(logpop[[1]][[1]]))+10)),col="gray")
abline(v=log(var(logpop[[1]][[1]])),col="red",lwd=1.5)
mypop<-partition_samples_nAB(logpop[[1]][[1]])[c(1,4,7,8)]
hist(as.numeric(nulldistr1iteration[,1]),breaks=0:6,main="nblocks LR 1000 permutations",xlab="nblocks",col="gray")
abline(v=as.numeric(mypop[1])-0.5,col="red",lwd=1.5)
hist(as.numeric(nulldistr1iteration[,3]),breaks=50,main="LR 1000 permutations",xlab="pvalue LR",col="gray")
abline(v=as.numeric(mypop[3]),col="red",lwd=1.5)
hist(as.numeric(nulldistr1iteration[,4]),breaks=50,main="RL AICc 1000 permutations",xlab="RL",col="gray")
abline(v=as.numeric(mypop[4]),col="red",lwd=1.5)
dev.off()

par(mfrow=c(1,1))
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/empiricaldistrperm_nAB_pop1.pdf")
vioplot(permbychr[1,],permbychr[10,],permbychr[20,],permbychr[30,],permbychr[40,],permbychr[50,],permbychr[60,],permbychr[70,],permbychr[80,],permbychr[90,],permbychr[100,],permbychr[106,],ylim=c(0,25000),col="firebrick",colMed="firebrick2")
vioplot(permbychrblocks2[1,],permbychrblocks2[10,],permbychrblocks2[20,],permbychrblocks2[30,],permbychrblocks2[40,],permbychrblocks2[50,],permbychrblocks2[60,],permbychrblocks2[70,],permbychrblocks2[80,],permbychrblocks2[90,],permbychrblocks2[100,],permbychrblocks2[106,],add=TRUE,col="cadetblue3",colMed="cadetblue4")
mypop<-logpop[[1]][[1]][order(logpop[[1]][[1]])]
points(c(mypop[1],mypop[10],mypop[20],mypop[30],mypop[40],mypop[50],mypop[60],mypop[70],mypop[80],mypop[90],mypop[100],mypop[106]),pch=19,col="darkblue")
dev.off()

cdat <- as.list(as.data.frame(t(as.matrix(permbychr))))
names(cdat)[1] <- "x"  # vioplot() needs the first element to be called 'x'
do.call(vioplot,c(cdat,list(col="firebrick",colMed="firebrick",add=TRUE)))

load("LDLD170216.RData")
setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
if (FALSE)
{
    par(mfrow=c(1,1))
    for (i in 1:12)
    {
        pdf(paste0("~/Dropbox/LDLD/ipynb/figs/coding1000g/logplot",i,".pdf"))
        infoplot(i,mymaxk=5,reshufflings=testfull[[i]]$raw)
        dev.off()
    }
    pdf("logplot_panel.pdf")
    par(mfrow=c(2,2))
    for (i in 1:4)
        {
        infoplot(i)
        }
    dev.off()
}
}
}
#================================================================================================^
#===============PER SAMPLE PER SNP ANALYSES: IDENTIFYING BIASED SNPs -> removal1 ================
#================================================================================================
#import genotypes (nAall)
{
load(paste0(myfolder,"/l_mutsamples.RData"))
nAmutations<-c(rbind(unlist(l_mutsamples),unlist(l_mutsamples))) #mutational load per individual
load(paste0(myfolder,"/combined_pvals.RData"))
nspops<-as.numeric(system(paste0("cat ",myfolder,"/nspops.txt"),intern=TRUE))
nAgenotypes<-fread(paste0("zcat ",myfolder,"/logs/all.freqlog.gz")) #nAs per snp; I have mutlog per pop; I have nAB
if (is.na(nAgenotypes[1,dim(nAgenotypes)[2],with=FALSE])) {nAgenotypes<-nAgenotypes[,1:(dim(nAgenotypes)[2]-1),with=FALSE]}
nAgenotypes_l<-list()
for (i in 1:12)
{
nAgenotypes_l[[i]]<-nAgenotypes[,samples_imypopi(nspops,i)+2,with=FALSE]
}

nAgenotypes<-subset(nAgenotypes,c(nAgenotypes[,2,with=FALSE]<10^11))
nAgenotypes$V2<-as.integer(nAgenotypes$V2)
snpsnames<-apply(as.data.frame(cbind(nAgenotypes$V1,nAgenotypes$V2)),MARGIN=1, function(x) paste0(x,collapse="."))

#checks: in first row, I trasform genotypes nAgenotypes in single chromosome nAall.
#positions where in nAall where I would have a 1 for first row
#genotypes correspond
#which(nAgenotypes[1,-c(1,2)]==2)
#which(nAgenotypes[1,-c(1,2)]==1)
#sort(c((which(nAgenotypes[1,-c(1,2)]==2))*2-1,(which(nAgenotypes[1,-c(1,2)]==2))*2))
#which(nAgenotypes[1,-c(1,2)]==1)*2-1

nsamples<-ncol(nAgenotypes)-2
nsites<-nrow(nAgenotypes)
#if small substructure, can I take that into account? if sample very different, more linkage. therefore it makes no sense to correct by distance. help comes from multiple populations. And if single pop I might use derived/ancestral rather than minor/major!
nAall2b<-floor(nAgenotypes[,3:ncol(nAgenotypes),with=FALSE]/2) #2->1,1->0
nAall1b<-ceiling(nAgenotypes[,3:ncol(nAgenotypes),with=FALSE]/2) #2->1
nAall<-cbind(nAall1b,nAall2b)
myindices<-c(rbind(1:(ncol(nAall)/2),(1+ncol(nAall)/2):ncol(nAall)))
nAall<-nAall[,myindices,with=FALSE]


#tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr1.vcf.gz 1:69427-69428 |  ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/nspops.txt 0.05 12 0 | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk '{for (i=10;i<=NF;i++){if ($i==1){printf "%s\t",i-9 }};printf "\n" }' #agrees with shell scan.
#tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr1.vcf.gz 1:69427-69428 |  ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/nspops.txt 0.05 12 0 | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk '{for (i=10;i<=NF;i++){if ($i==2){printf "%s\t",i-9 }};printf "\n" }'
#zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/logs/all.freqlog.gz | grep 9427508 | awk '{for (i=3;i<=NF;i++){if ($i==1){printf "%s\t",(i-2) }};printf "\n" }'

#tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr1.vcf.gz 1:9427507-9427508 |  ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/nspops.txt 0.05 12 0 | awk '{for (i=1;i<NF;i++){if ((i<=2) || (i>9)){printf "%s\t",$i}};printf "%s\n",$NF}' | sed 's/|/\t/g' > ~/Dropbox/LDLD/scripts/temp.tab
}
#snp discovery with glm
if (FALSE)
{
#load(paste0(myfolder,"/dataFIN.RData"))
load(paste0(myfolder,"/combined_pval.RData"))
load(paste0(myfolder,"/empfdr.RData"))
load(paste0(myfolder,"/dataFIN_l.RData"))
load(paste0(myfolder,"/combined_pval_l.RData"))
load(paste0(myfolder,"/logpop_nAB_nlinks.RData"))
#logisticpvals_nAvsnAB_l<-corr_nAvsnAB(nAall_l,logpop_nAB,l_mutsamples)
#save(logisticpvals_nAvsnAB_l,file=paste0(myfolder,"/logisticpvals_nAvsnAB_l.RData"))
load(paste0(myfolder,"/logisticpvals_nAvsnAB_l.RData"))
load(paste0(myfolder,"/dataFIN_sign.RData"))

#in all SNPS #553
#GLM
{
logisticpval_nAvsnAB<-apply(sapply(1:12,function(z) logisticpvals_nAvsnAB_l[[z]] ), MARGIN=1, function(y) pchisq( -2*sum(log(y)), df=2*length(y), lower.tail=FALSE))
logisticpfdr_nAvsnAB_all<-p.adjust(p=logisticpval_nAvsnAB, method = "fdr", n = length(logisticpval_nAvsnAB))
#In SNPS in significant pairs #550
fullsign_snps<-unique(c(apply(cbind(dataFIN_sign$chrA,dataFIN_sign$posA),MARGIN=1,function(x) paste0(x,collapse=".")),
apply(cbind(dataFIN_sign$chrB,dataFIN_sign$posB),MARGIN=1,function(x) paste0(x,collapse="."))))
logisticpval_nAvsnAB<-logisticpval_nAvsnAB[!is.na(match(snpsnames,fullsign_snps))]
logisticpfdr_nAvsnAB_sign<-p.adjust(p=logisticpval_nAvsnAB, method = "fdr", n = length(logisticpval_nAvsnAB))
sum(logisticpfdr_nAvsnAB_all<0.05) 
sum(logisticpfdr_nAvsnAB_sign<0.05) 
#I guess this means that I will keep doing it in all SNPs
removed_snps<-snpsnames[logisticpfdr_nAvsnAB_all<0.05]
#create BED file with removed snps
write.table(stringsnp2bed(as.matrix(removed_snps)),file=paste0(myfolder,"/removed1_glm.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
#write.table(stringsnp2bed(as.matrix(removed_snps)),file=paste0(myfolder,"/removed1.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)

}

}
#snp discovery with glmm
{
library(lme4)
for (mysign_logpop_nAB in c("_pos"))#c("_pos01","_pos","_neg","_pos20"))
{
load(paste0(myfolder,"/logpop_nAB",mysign_logpop_nAB,".RData"))
logpopnAB<-get(paste0("logpop_nAB",mysign_logpop_nAB))
#mylogpopnAB_reschr<-get(paste0("logpop_reschr_nAB",mysign_logpop_nAB,"_l"))

#logpopnAB_clean<-sapply(mylogpopnAB, function(x) x[1:(length(x)-2)])
data0<-data.table(sampleID=1:length(unlist(logpopnAB)),population=unlist(sapply(1:length(logpopnAB),function(x) rep(x,length(logpopnAB[[x]])))),nAB=unlist(logpopnAB))
data0<-rbind(data0,data0)[myindices,]
data0$sampleID<-as.factor(data0$sampleID)
data0$population<-as.factor(data0$population)
data0$nAmutations<-nAmutations
setnames(nAall,paste0("V",1:ncol(nAall)))
contr=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000),calc.derivs=F)
#nAall<-nAall[,order(data0$population),with=FALSE]
#data0<-data0[order(data0$population),]
data0$nAB<- scale(data0$nAB)
data0$nAmut<- scale(data0$nAmutations)
    
#onlychr1    
time.myf<-system.time(bad_snps_pvalues<-sapply(1:9594,function(x) {print(x);myfx(x,myfreq=freqfiltering)}))
save(bad_snps_pvalues,file=paste0(myfolder,"/bad_snps_pvalues_chr1",mysign_logpop_nAB,".RData"))
bad_snps_fdr<-p.adjust(p=unlist(bad_snps_pvalues),method='fdr',n=9594)
nsnps_sign<-sum(bad_snps_fdr<0.05)
if (nsnps_sign>0){
    bad_snps_pvalues_1<-sapply(bad_snps_pvalues,function(x) if (length(x)==0) return(1) else return(x))
    bad_snps<-snpsnames[bad_snps_pvalues_1<=sort(bad_snps_pvalues_1)[nsnps_sign]]
    bad_snps.df<-stringsnp2bed(bad_snps)
    write.table(bad_snps.df,file=paste0(myfolder,"/removed1_chr1",mysign_logpop_nAB,".bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
}
#get many snps but bad, via command line the same thing gives me few snps. Are they good at least?

#let's check snp 1
iline<-1
data_temp<-cbind(nA=t(nAall[iline,]),data0)
#setnames(data_temp,c('nA','sampleID','population','nAB'))
#data_temp<-data_temp[unlist(unname(aggregate(data_temp$nA,by=list(data_temp$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
data_temp_snp1_R<-data_temp
snpsnames[iline] #"1.69428"

myfx<-function(iline,myfreq=0.05) {
data_temp0<-cbind(nA=t(nAall[iline,]),data0)
setnames(data_temp0,c('nA','sampleID','population','nAB'))
data_temp<-data_temp0[unlist(unname(aggregate(data_temp0$nA,by=list(data_temp0$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
if ( nrow(data_temp)>0 ){
glmres<-preliminary_check_glm_nAvsnAB(data_temp$nA,data_temp$nAB,data_temp$population)
if (glmres<1 ) { if (length(unique(data_temp$population))>1) {return(c(fx(data_temp),fx(data_temp0),glmres))} else {return(c(fx1(data_temp),fx(data_temp0),glmres))} }}
}

myfx_parse_mut(iline) #0.7506533 #and also clearly glm is much more powerful

system.time(temp<-sapply(c(1:10,771),function(x) myfx_parse_mut(x)));
print(temp)

system.time(temp<-sapply(c(1:10,771),function(x) myfx_parse_mut(x)));
print(temp)

system.time(temp<-sapply(c(1000:1005,771),function(x) myfx_parse_mut(x)));
print(temp)

sapply(800:820,function(x) myfx_parse(x));
system.time(temp<-sapply(c(905,1717,25001,25006,1005,771),function(x) myfx_parse(x)));
print(temp)

system.time(temp<-sapply(c(905,1717,25001,25006,1005,771),function(x) myfx_parse_mut(x)));
print(temp)


sapply(c(1:10,771),function(x) myfx_parse_mut(x))
#          [,1]      [,2]       [,3]      [,4]       [,5]       [,6]       [,7]
#[1,] 0.7506533 0.8961193 0.03593717 0.3370001 0.87920527 0.87920527 0.87920527
#[2,] 0.5296736 0.1097622 0.01210315 0.1678743 0.03137448 0.03137448 0.03137448
#           [,8]      [,9]    [,10]
#[1,] 0.87920527 0.6043772 0.331842
#[2,] 0.03137448 0.5277385 0.293897


MYCHR=1
TARGETVCF=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz
myRData=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/logpop_nAB_pos.RData
mytag=logpop_nAB_pos
myinit=69428
myend=69430
tabix $TARGETVCF ${MYCHR}:${myinit}-${myend} |  ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/nspops.txt 0.05 12 0 | awk '{for (i=1;i<NF;i++){if ((i<=2) || (i>9)){printf "%s\t",$i}};printf "%s\n",$NF}' | sed 's/|/\t/g'  | ~/Dropbox/LDLD/scripts/nAvsnAB_scan_header.R $myRData $mytag >> ${destination_folder}/chr${MYCHR}.${myinit}.${myend}.$mytag
#output
#[1] "1 69428 0.750653326231213 FALSE FALSE"

MYCHR=1
TARGETVCF=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz
myRData=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/logpop_nAB_pos.RData
mytag=logpop_nAB_pos
myinit=69428
myend=69430
tabix $TARGETVCF ${MYCHR}:${myinit}-${myend} |  ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/nspops.txt 0.05 12 0 | awk '{for (i=1;i<NF;i++){if ((i<=2) || (i>9)){printf "%s\t",$i}};printf "%s\n",$NF}' | sed 's/|/\t/g'  | ~/Dropbox/LDLD/scripts/nAvsnAB_scan_header_mut.R $myRData $mytag >> ${destination_folder}/chr${MYCHR}.${myinit}.${myend}.$mytag


    
time.myf<-system.time(bad_snps_pvalues<-sapply(1:nrow(nAall),function(x) {print(x);myfx(x,myfreq=freqfiltering)}))
save(bad_snps_pvalues,file=paste0(myfolder,"/bad_snps_pvalues",mysign_logpop_nAB,".RData"))
load(paste0(myfolder,"/bad_snps_pvalues",mysign_logpop_nAB,".RData"))


bad_snps_fdr<-p.adjust(p=unlist(bad_snps_pvalues),method='fdr',n=nrow(nAall))
nsnps_sign<-sum(bad_snps_fdr<0.05)
if (nsnps_sign>0){
    bad_snps_pvalues_1<-sapply(bad_snps_pvalues,function(x) if (length(x)==0) return(1) else return(x))
    bad_snps<-snpsnames[bad_snps_pvalues_1<=sort(bad_snps_pvalues_1)[nsnps_sign]]
    bad_snps.df<-stringsnp2bed(bad_snps)
    write.table(bad_snps.df,file=paste0(myfolder,"/removed1",mysign_logpop_nAB,".bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
}

#for (i in names(nAall1)) { nAall1[eval(parse(text=i))==2,i]<-1;print(i) }
#for (i in names(nAall2)) { nAall2[eval(parse(text=i))==1,i]<-0;print(i) }
#for (i in names(nAall2)) { nAall2[eval(parse(text=i))==2,i]<-1;print(i) }

nAall<-cbind(nAall1,nAall2)
setnames(nAall,paste0("V",1:ncol(nAall)))
contr=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000),calc.derivs=F)
nAall<-nAall[,order(data0$population),with=FALSE]
data0<-data0[order(data0$population),]
data0$nAB<- scale(data0$nAB)


#unlist(unname(aggregate(data_temp$nA,by=list(data_temp$population),function(x) if (sum(x)/length(x)>0.05 & sum(x)/length(x)<0.95){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2]))

#aggregate(data_temp$nA,by=list(data_temp$population),function(x) sum(x)/length(x))
#aggregate(data_temp$nA,by=list(data_temp$population),function(x) length(x))

list_data0_nAll<-list(data0,nAall)
#save(list_data0_nAll,file=paste0(myfolder,"/list_data0_nAll.RData"))
data0<-list_data0_nAll[[1]]
nAall<-list_data0_nAll[[2]]

myfx<-function(data_temp,myfreq=0.01) {
data_temp<-data_temp[unlist(unname(aggregate(data_temp$nA,by=list(data_temp$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
#print(paste0("pop above thr: ",nrow(data_temp)))
if ( nrow(data_temp)>0 ){
if (preliminary_check_glm_nAvsnAB(data_temp$nA,data_temp$nAB,data_temp$population)<0.01 ) { if (length(unique(data_temp$population))>1) {return(fx(data_temp))} else {return(fx1(data_temp))} }}
}

#myfolder='/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5'
#mysign_logpop_nAB='_pos20'
time.myf<-system.time(bad_snps_pvalues<-sapply(1:nrow(nAall),function(x) {print(x);myfx(x,myfreq=freqfiltering)}))
save(bad_snps_pvalues,file=paste0(myfolder,"/bad_snps_pvalues",mysign_logpop_nAB,".RData"))
mysign_logpop_nAB='_pos'
load(paste0(myfolder,"/bad_snps_pvalues",mysign_logpop_nAB,".RData"))


bad_snps_fdr<-p.adjust(p=unlist(bad_snps_pvalues),method='fdr',n=nrow(nAall))
nsnps_sign<-sum(bad_snps_fdr<0.05)
if (nsnps_sign>0){
    bad_snps_pvalues_1<-sapply(bad_snps_pvalues,function(x) if (length(x)==0) return(1) else return(x))
    bad_snps<-snpsnames[bad_snps_pvalues_1<=sort(bad_snps_pvalues_1)[nsnps_sign]]
    bad_snps.df<-stringsnp2bed(bad_snps)
    write.table(bad_snps.df,file=paste0(myfolder,"/removed1",mysign_logpop_nAB,".bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
}
#check a single variant. For example
{
mybadsnp<-paste(bad_snps.df$chr[1],as.character(bad_snps.df$end)[1],sep='.')
mylen<-length(nAall[which(snpsnames==mybadsnp),])
mygenotypes<-unlist(unname(nAall[which(snpsnames==mybadsnp),1:(mylen/2)]+nAall[which(snpsnames==mybadsnp),(mylen/2+1):mylen]))

#tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr1.vcf.gz 1:9427507-9427508 |  ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/nspops.txt 0.05 12 0 | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk '{for (i=10;i<=NF;i++){if ($i==1){printf "%s\t",i-9 }};printf "\n" }' #agrees with shell scan.
#zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/logs/all.freqlog.gz | grep 9427508 | awk '{for (i=3;i<=NF;i++){if ($i==1){printf "%s\t",(i-2) }};printf "\n" }'

#tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr1.vcf.gz 1:9427507-9427508 |  ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/nspops.txt 0.05 12 0 | awk '{for (i=1;i<NF;i++){if ((i<=2) || (i>9)){printf "%s\t",$i}};printf "%s\n",$NF}' | sed 's/|/\t/g' > ~/Dropbox/LDLD/scripts/temp.tab
mynA<-read.table("~/Dropbox/LDLD/scripts/temp.tab")
mynA<-mynA[,-(1:2)]
unname(which(mynA==1))
nAll2<-mynA[which(1:ncol(mynA) %% 2==0)];nAll1<-mynA[which(1:ncol(mynA) %% 2==1)]
mygenotypes<-unlist(unname(nAll1+nAll2))


}

}

}

#tests
if (FALSE)
{
#tests:
{
#very cool: first of all myfx_glmm_nofilter on all populations (even below freq threshold) and myfx_glmm superstrongly correlate. 
#second, myf is very fast and at the same time retain all highly significant cases.
time.glm<-system.time(res_glm<-sapply(1:200,function(x) {myfx_glm(x)}))
time.glmm<-system.time(res_glmm<-sapply(1:200,function(x) {myfx_glmm(x)}))
time.myf<-system.time(res<-sapply(1:200,function(x) {myfx(x)}))
time.glmm.nof<-system.time(res_glmm.nof<-sapply(1:200,function(x) {myfx_glmm_nofilter(x)}))
cbind(time.glm,time.glmm,time.glmm.nof,time.myf)
#           time.glm time.glmm time.glmm.nof time.myf
#user.self    25.144   180.512       338.580   51.480
cor.test(res_glm,res_glmm) #0.4767534
cor.test(res_glmm.nof,res_glmm) #0.95
cor.test(res_glmm.nof,res_glm) #0.41
#> as.numeric(res_glm<0.01)[res_glm<0.05]
# [1] 0 0 1 1 1 1 0 0 1 1 0 1 1 1 0 0 0 0 1 0 0 1 0 0 0 1 0 1 0 1
#> as.numeric(res_glm<0.01)[res_glm>0.05]
#  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [38] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [75] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[112] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[149] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
}
{
#iline where it stops = 16272
#25648421 3 in pop0 from sorted.res.gz should be 
#but it is 2 from
myfreq<-freqfiltering
iline<-7661
data_temp<-cbind(nA=t(nAall[iline]),data0)
setnames(data_temp,c('nA','sampleID','population','nAB'))
aggregate(data_temp$nA,by=list(data_temp$population),function(x) sum(x))
data_temp<-data_temp[unlist(unname(aggregate(data_temp$nA,by=list(data_temp$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
population<-data_temp$population
mypop<-unique(population)[1]
mynAB<-data_temp$nAB[data_temp$population==mypop]
mynA<-data_temp$nA[data_temp$population==mypop]
model <- glm( mynA ~ mynAB,family=binomial(link='logit'))
model0 <- glm( mynA ~ 1,family=binomial(link='logit'))
pvalue<-anova(model0, model, test = "Chisq")[[5]][2]

#empty

mypos<-unname(sapply(snpsnames, function(x) strsplit(x,split='[.]')[[1]][2]))
mypos[[iline]]
data_temp<-cbind(nA=t(nAall[which(mypos==mypos[[iline]])]),data0)
setnames(data_temp,c('nA','sampleID','population','nAB'))
aggregate(data_temp$nA,by=list(data_temp$population),function(x) sum(x)) #do not coincide. 

sum(nAall[iline,]) #4
sum(nAall[which(mypos==mypos[[iline]]),]) #13
#overall frequencies match, that means that conversion must have worked fine, and positions are the right ones

sum(unlist(unname(nAall[iline,]))[data0$population==1]) #2 this is what observed in R

data0<-data.table(sampleID=1:length(unlist(logpopnAB_clean)),population=unlist(sapply(1:length(logpopnAB_clean),function(x) rep(x,length(logpopnAB_clean[[x]])))),nAB=unlist(logpopnAB_clean))
data0<-rbind(data0,data0)
data0$sampleID<-as.factor(data0$sampleID)
data0$population<-as.factor(data0$population)

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #3 this is what observed in LDLD

data0<-data0[order(data0$population),]

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #2 this is what observed in LDLD

nAall<-nAall[,order(data0$population),with=FALSE]

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #2 this is what observed in LDLD

data0<-data.table(sampleID=1:length(unlist(logpopnAB_clean)),population=unlist(sapply(1:length(logpopnAB_clean),function(x) rep(x,length(logpopnAB_clean[[x]])))),nAB=unlist(logpopnAB_clean))
data0<-rbind(data0,data0)
data0$sampleID<-as.factor(data0$sampleID)
data0$population<-as.factor(data0$population)
nAall<-cbind(nAall1,nAall2)
setnames(nAall,paste0("V",1:ncol(nAall)))
contr=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000),calc.derivs=F)
nAall<-nAall[,order(data0$population),with=FALSE]
data0<-data0[order(data0$population),]
data0$nAB<- scale(data0$nAB)

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #2 this is what observed in LDLD
}
#check nAall
{
#iline where it stops = 2503
#25648421 3 in pop0 from sorted.res.gz should be 
#but it is 2 from 
data_temp<-cbind(nA=t(nAall[2503]),data0)
setnames(data_temp,c('nA','sampleID','population','nAB'))
aggregate(data_temp$nA,by=list(data_temp$population),function(x) sum(x))
data_temp<-data_temp[unlist(unname(aggregate(data_temp$nA,by=list(data_temp$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
#empty

mypos<-unname(sapply(snpsnames, function(x) strsplit(x,split='[.]')[[1]][2]))
data_temp<-cbind(nA=t(nAall[which(mypos==120197784)]),data0)
setnames(data_temp,c('nA','sampleID','population','nAB'))
aggregate(data_temp$nA,by=list(data_temp$population),function(x) sum(x)) #do not coincide. 

sum(nAall[2503,]) #4
sum(nAall[which(mypos==120197784),]) #13
#overall frequencies match, that means that conversion must have worked fine, and positions are the right ones

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #2 this is what observed in R

data0<-data.table(sampleID=1:length(unlist(logpopnAB_clean)),population=unlist(sapply(1:length(logpopnAB_clean),function(x) rep(x,length(logpopnAB_clean[[x]])))),nAB=unlist(logpopnAB_clean))
data0<-rbind(data0,data0)
data0$sampleID<-as.factor(data0$sampleID)
data0$population<-as.factor(data0$population)

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #3 this is what observed in LDLD

data0<-data0[order(data0$population),]

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #2 this is what observed in LDLD

nAall<-nAall[,order(data0$population),with=FALSE]

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #2 this is what observed in LDLD

data0<-data.table(sampleID=1:length(unlist(logpopnAB_clean)),population=unlist(sapply(1:length(logpopnAB_clean),function(x) rep(x,length(logpopnAB_clean[[x]])))),nAB=unlist(logpopnAB_clean))
data0<-rbind(data0,data0)
data0$sampleID<-as.factor(data0$sampleID)
data0$population<-as.factor(data0$population)
nAall<-cbind(nAall1,nAall2)
setnames(nAall,paste0("V",1:ncol(nAall)))
contr=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000),calc.derivs=F)
nAall<-nAall[,order(data0$population),with=FALSE]
data0<-data0[order(data0$population),]
data0$nAB<- scale(data0$nAB)

sum(unlist(unname(nAall[2503,]))[data0$population==1]) #2 this is what observed in LDLD
}
}


mypvals<-10^-runif(1143,5,16)
write.table(file='~/Dropbox/LDLD/figs/pvals.tab',cbind(mypvals,p.adjust(mypvals,method='fdr',n=128000)),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)

#check on bad-snps
{
#============================
#===== TESTS before revisions ========
#============================

MYFILE=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/bad_snps_glmm/bad_snps_logpop_nAB_pos_glmm_all.bed
MYFILE=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/bad_snps_glmm/bad_snps_logpop_nAB_pos_glmm_all.bed
MYFILE=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/bad_snps_glmm/bad_snps_logpop_nAB_pos20_glmm_all.bed
MYFILE=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/bad_snps_glmm/bad_snps_logpop_nAB_pos20_glmm_all.bed
cat ${MYFILE} | uniq > temp; mv temp ${MYFILE}

#GWAS
{
#very few in GWAS -> no table in paper
cd ~/Dropbox/general_utils/
wget 'https://www.ebi.ac.uk/gwas/api/search/downloads/full'
wget 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative'
cat ~/Dropbox/general_utils/full | awk -v FS='\t' -v OFS='\t' '{print "chr"$12,$13-1,$13,$4,$5,$6,$8,$22,$25,$28,$33}' | sort -k1,1 -k2,2n > ~/Dropbox/general_utils/gwas_catalog.bed
#to liftover from hg38 to hg19, then I can overlap

cat gwas_catalog_v1.0-associations_e89_r2017-07-31.tsv | awk -v FS='\t' -v OFS='\t' '{print "chr"$12,$13-1,$13}' | grep -v '[x;-]' | sort -k1,1 -k2,2n > ~/Dropbox/general_utils/gwas_catalog_v1.0-associations_e89_r2017-07-31.bed
bedtools intersect -b /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed -a <( cat ~/Dropbox/general_utils/gwas_catalog_v1.0-associations_e89_r2017-07-31_hg19.bed | sed 's/chr//g' )
#2	27587723	27587724 from coding #in hg38 27364857
#2013-06-14	23382691	Lauc G	2013-01-31	PLoS Genet	www.ncbi.nlm.nih.gov/pubmed/23382691	Loci associated with N-glycosylation of human immunoglobulin G show pleiotropy with autoimmune diseases and haematological cancers.	IgG glycosylation	2,247 European ancestry individuals	NA	2p23.3	2	27364857	NR	EIF2B4			8890			rs1058065-G	rs1058065	0	1058065	synonymous_variant	0	0.97552739311682	4E-6	5.3979400086720375	(IGP56)	0.54	[0.31-0.77] unit increase	Illumina [~ 2500000] (imputed)	N
#22	23299159	23299160 from intergenic
bedtools intersect -a <( sortBed -i /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/removed1_pos20.bed ) -b <( cat ~/Dropbox/general_utils/gwas_catalog_v1.0-associations_e89_r2017-07-31_hg19.bed | sed 's/chr//g' ) | uniq
#in intergenic
bedtools intersect -b /home/fabrizio_mafessoni/Dropbox/LDLD/ms/Ajhg/TableS5.bed.txt -a <( cat ~/Dropbox/general_utils/gwas_catalog_v1.0-associations_e89_r2017-07-31_hg19.bed | sed 's/chr//g' )

}
#quality 1000 genomes
{
cd /mnt/scratch/fabrizio/LDLD/above95/downloads
wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.pilot_mask.autosomes.bed' #in release..supporting
for mychr in `seq 1 22`;do
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/20140520.chr${mychr}.strict_mask.fasta.gz"
done
flatten () {
awk '{if (substr($1,1,1)==">"){convert=0;for (k = 1; k <= 9; k++){if (substr($1,2,1)==k || substr($1,2,1)=="X"){convert=convert+1;counter=0;}}; 
if (convert>0){print $1}} else if (convert>0){for (i = 1; i <= length($1); i++){counter=counter+1;print substr($1,i,1);}}}'
}
for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
for mychr in `seq 1 22`;do
zcat /mnt/scratch/fabrizio/LDLD/above95/downloads/20140520.chr${mychr}.strict_mask.fasta.gz | flatten | awk -v OFS='\t' -v MYCHR=${mychr} 'BEGIN{COUNTER=0}{if (substr($1,1,1)!=">"){COUNTER=COUNTER+1;print MYCHR,COUNTER-1,COUNTER,$1}}' | grep [${letter}] | gzip -f  > /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.${letter}.bed.gz &
bedtools merge -i /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.${letter}.bed.gz | gzip -f > /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz
done;done

for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
for mychr in `seq 1 22`;do
echo $mychr
bedtools intersect -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz -sorted >> /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.${letter}.bed
done;done
for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
for mychr in `seq 1 22`;do
echo $letter $mychr
bedtools intersect -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/chr${mychr}.vcf.gz -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz -sorted |wc -l >> /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/overlap.all.$letter
done;done

for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/overlap.all.$letter | awk -v let=$letter 'BEGIN{co=0}{co=co+$1}END{print let,co}'
done
#H 3781
#P 230977
#L 12059
#N 0
#ZQ 70508
wc -l /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap*
   24 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.H.bed
  745 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.P.bed
     67 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.L.bed
    0 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.N.bed
  306 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.ZQ.bed
 1142 total

 obs<-c(3781,230977,12059,70508)
 candidates<-c(24,745,67,306)
 obs/sum(obs)
candidates/sum(candidates)

#>  obs/sum(obs)
#[1] 0.01191523 0.72788781 0.03800205 0.22219491
#> candidates/sum(candidates)
#[1] 0.02101576 0.65236427 0.05866900 0.26795096
#we see slight enrichment in H
}
#Complete Genomics
{

#generate files
#CGBED=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/#Complete_Public_Genomes_69genomes_all_listvariants.bed.gz
#bzcat /mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/#Complete_Public_Genomes_69genomes_all_listvariants.tsv.bz2 | grep -v chromosome | awk -v OFS='\t' '{print $2,$3,$4}' | gzip -f > temp.bed.gz
#bedtools sort -i temp.bed.gz | gzip -f > $CGBED #maybe I should only do snps. Now it is snps and things falling within deletions. no insertions.
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
zcat $CGBED | awk -v mychrom=chr$MYCHROM '{if (($1==mychrom) && ($3>$2)){print}}' | sed 's/chr//g' | gzip -f > $CGBED_CHR &
done


#====SET UP VARIABLES======
{
#non merged------------------------
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5

MYFOLDERBACKGROUND=$MYFOLDER
MYBADFOLDER=${MYFOLDER}/bad_snps_glmm
MYBADSNPS=bad_snps_logpop_nAB_pos_glmm_filt.bed
MYBADSNPS=bad_snps_logpop_nAB_pos_glmm_all.bed
MYBADSNPS=bad_snps_logpop_nAB_pos20_glmm_all.bed
Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab


#merged----------------------------
MYMIN=1 #MYMIN=5
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min${MYMIN}
MYFOLDERBACKGROUND=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl1/min${MYMIN}
MYBADFOLDER=${MYFOLDER}/bad_snps_glmm
MYBADSNPS=bad_snps_logpop_nAB_pos20_glmm_filt.bed
MYBADSNPS=bad_snps_logpop_nAB_pos20_glmm_all.bed
MYBADSNPS=bad_snps_logpop_nAB_pos20_glmm_all_top.bed
MYDESTINATION_FOLDER=${MYFOLDER}/anal
mkdir -p $MYDESTINATION_FOLDER

Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

sortBed -i ${MYBADFOLDER}/${MYBADSNPS} > temp
mv temp ${MYBADFOLDER}/${MYBADSNPS}


}
#1) create background set (stratify background snps by frequency )
{
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
mkdir -p $MYFOLDER/anal
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
MYFILE=chr${MYCHROM}.vcf.gz
echo $MYFILE
IPS1=$(( $IPS + 1 ))
bcftools view -q ${MYPS[$IPS]}:minor -Q ${MYPS[$IPS1]}:minor  $MYFOLDERBACKGROUND/$MYFILE  | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2}'  >> $MYFOLDER/anal/byfreq_${IPS}.bed
done;
gzip -f $MYFOLDER/anal/byfreq_${IPS}.bed
done
}
#2) create background set ( overlap stratified background snps with CG)
{
#COMPLETE GENOMICS
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
bedtools intersect -a <( zcat $MYFOLDER/anal/byfreq_${IPS}.bed.gz | sort -k1,1 -k2,2n ) -b $CGBED_CHR -sorted | uniq >> $MYFOLDER/anal/byfreqCG_${IPS}.bed
done;
#gzip -f $MYFOLDER/anal/byfreqCG_${IPS}.bed
bedtools sort -i $MYFOLDER/anal/byfreqCG_${IPS}.bed | gzip -f > $MYFOLDER/anal/byfreqCG_${IPS}.bed.gz
rm $MYFOLDER/anal/byfreqCG_${IPS}.bed
done




}
#step 3 and 4 differ if merged
#---non merged---all bad snps are included in background vcf files
{
#3) overlap 
#sort -k1,1n -k2,2n ${MYFOLDER}/${MYBADSNPS} > temp
#cp temp ${MYFOLDER}/${MYBADSNPS}
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
for IPS in `seq 0 11`;do
echo $IPS
bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b $MYFOLDER/anal/byfreqCG_${IPS}.bed.gz | gzip -f > $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz
bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b $MYFOLDER/anal/byfreq_${IPS}.bed.gz | gzip -f > $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz
done

rm $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
for IPS in `seq 0 11`;do
TOT_SNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}.bed.gz | wc -l )
TOT_SNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}.bed.gz | wc -l )
TOT_BADSNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | wc -l )
TOT_BADSNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz | wc -l )
echo ${MYPS[${IPS}]} ${TOT_SNPS_CG[$IPS]} ${TOT_SNPS[$IPS]} ${TOT_BADSNPS_CG[$IPS]} ${TOT_BADSNPS[$IPS]} >> $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
done
Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

}
#---merged
{
#for intergenic I have to first extract my variants from 1000g

for MYCHROM in `seq 1 22`;do
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -R ${MYBADFOLDER}/${MYBADSNPS} -h | gzip -f > ${MYDESTINATION_FOLDER}/chr${MYCHROM}.vcf.gz
done

#3) stratify badsnps by frequency 
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
MYFILE=chr${MYCHROM}.vcf.gz
echo $MYFILE
IPS1=$(( $IPS + 1 ))
bcftools view -q ${MYPS[$IPS]}:minor -Q ${MYPS[$IPS1]}:minor  ${MYDESTINATION_FOLDER}/$MYFILE  | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2}'  >> ${MYFOLDER}/anal/byfreq_${IPS}_${MYBADSNPS}
done;
gzip -f ${MYFOLDER}/anal/byfreq_${IPS}_${MYBADSNPS}
done


#4) overlap background set ( overlap stratified background snps with CG)
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
bedtools intersect -a <( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | sort -k1,1 -k2,2n ) -b $CGBED_CHR -sorted >> $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}
done;
gzip -f $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}
done

rm $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
for IPS in `seq 0 11`;do
TOT_SNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}.bed.gz | wc -l )
TOT_SNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}.bed.gz | wc -l )
TOT_BADSNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | wc -l )
TOT_BADSNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz | wc -l )
echo ${MYPS[${IPS}]} ${TOT_SNPS_CG[$IPS]} ${TOT_SNPS[$IPS]} ${TOT_BADSNPS_CG[$IPS]} ${TOT_BADSNPS[$IPS]} >> $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
done
Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab


}
#plots
{
#cp plots of overlap with CG
cp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/anal/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab_sharingCG.pdf ~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min1
cp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/anal/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab_sharingCG.pdf ~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min5
cp /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/anal/byfreq_summary_bad_snps_logpop_nAB_pos_glmm_all.bed.tab_sharingCG.pdf ~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min1
cp /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal/byfreq_summary_bad_snps_logpop_nAB_pos_glmm_all.bed.tab_sharingCG.pdf ~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min5

#cp 
cp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/bad_snps_glmm/*bed ~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min1
cp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/bad_snps_glmm/*bed ~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min5
cp /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/bad_snps_glmm/*bed ~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min1
cp /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/bad_snps_glmm/*bed ~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min5

#generate tables
cat ~/Dropbox/LDLD/ms/genomebiology/TableS2.bed | grep '#' > ~/Dropbox/LDLD/ms/genomebiology/TableS2.bed
cat ~/Dropbox/LDLD/ms/genomebiology/TableS2.bed | grep '#' > ~/Dropbox/LDLD/ms/genomebiology/TableS3.bed
cat ~/Dropbox/LDLD/ms/genomebiology/TableS2.bed | grep '#' > ~/Dropbox/LDLD/ms/genomebiology/TableS4.bed
cat ~/Dropbox/LDLD/ms/genomebiology/TableS2.bed | grep '#' > ~/Dropbox/LDLD/ms/genomebiology/TableS5.bed
sortBed -i ~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min5/bad_snps_logpop_nAB_pos_glmm_all.bed >> ~/Dropbox/LDLD/ms/genomebiology/TableS2.bed
sortBed -i ~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min1/bad_snps_logpop_nAB_pos_glmm_all.bed >> ~/Dropbox/LDLD/ms/genomebiology/TableS3.bed
sortBed -i ~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min5/bad_snps_logpop_nAB_pos20_glmm_all.bed >> ~/Dropbox/LDLD/ms/genomebiology/TableS4.bed
sortBed -i ~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min1/bad_snps_logpop_nAB_pos20_glmm_all.bed >> ~/Dropbox/LDLD/ms/genomebiology/TableS5.bed


#4) estimate enrichment
#chisquare
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal/byfreq_summary_GoNL_bad_snps_logpop_nAB_pos_glmm_all.bed.tab #p-value = 2.2e-16 #fdr_total 0.1657119
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/anal/byfreq_summary_GoNL_bad_snps_logpop_nAB_pos_glmm_all.bed.tab #p-value = 2.2e-16 #fdr_total 0.1745332
#/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/anal/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab #p-value < 2.2e-16 #fdr_total: 0.670791
#/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/anal/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab #p-value < 6.059e-14  #fdr total: 0.5478318


}
}
#GoNL vs Complete Genomics
{
#====SET UP VARIABLES======
{

zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*.vcf.gz | grep -v '#' | uniq |wc -l #87663
CO=0
for i in `seq 1 22`;do
CO1=$( bedtools intersect -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${i}.vcf.gz -b /mnt/sequencedb/GoNL/gonl.chr${i}.snps_indels.r5.vcf.gz -sorted | uniq | wc -l )
CO=$(( $CO + $CO1 ))
echo $CO
done
#53742
zcat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl2/min5/chr*.vcf.gz | grep -v '#' | uniq |wc -l #87953
zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*.vcf.gz | grep -v '#' | uniq |wc -l #87663
CO=0
for i in `seq 1 22`;do
CO1=$( bedtools intersect -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${i}.vcf.gz -b /mnt/sequencedb/GoNL/gonl.chr${i}.snps_indels.r5.vcf.gz -sorted | uniq | wc -l )
CO=$(( $CO + $CO1 ))
echo $CO
done
CO=0
for MYCHROM in `seq 1 22`;do
CO1=$( bedtools intersect -a /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -b /mnt/sequencedb/GoNL/gonl.chr${MYCHROM}.snps_indels.r5.vcf.gz -sorted | wc -l )
CO=$(( $CO + $CO1 ))
echo $CO
done

#non merged------------------------
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5 

MYFOLDERBACKGROUND=$MYFOLDER
MYBADFOLDER=${MYFOLDER}/bad_snps_glmm
MYBADSNPS=bad_snps_logpop_nAB_pos_glmm_filt.bed
MYBADSNPS=bad_snps_logpop_nAB_pos_glmm_all.bed

for i in `seq 1 22`;do
bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b /mnt/sequencedb/GoNL/gonl.chr${i}.snps_indels.r5.vcf.gz -wo | grep PASS | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' >> ${MYBADFOLDER}/GoNL_${MYBADSNPS}
done
sortBed -i ${MYBADFOLDER}/GoNL_${MYBADSNPS} > temp
mv temp ${MYBADFOLDER}/GoNL_${MYBADSNPS}
MYBADSNPS=GoNL_${MYBADSNPS}
#merged----------------------------


MYMIN=1 #MYMIN=5
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min${MYMIN}
MYFOLDERBACKGROUND=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl1/min${MYMIN}
MYBADFOLDER=${MYFOLDER}/bad_snps_glmm
MYBADSNPS=bad_snps_logpop_nAB_pos20_glmm_filt.bed
MYBADSNPS=bad_snps_logpop_nAB_pos20_glmm_all.bed
MYDESTINATION_FOLDER=${MYFOLDER}/anal
mkdir -p $MYDESTINATION_FOLDER

#Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

sortBed -i ${MYBADFOLDER}/${MYBADSNPS} > temp
mv temp ${MYBADFOLDER}/${MYBADSNPS}

for i in `seq 1 22`;do
bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b /mnt/sequencedb/GoNL/gonl.chr${i}.snps_indels.r5.vcf.gz -wo | grep PASS | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' >> ${MYBADFOLDER}/GoNL_${MYBADSNPS}
done
sortBed -i ${MYBADFOLDER}/GoNL_${MYBADSNPS} > temp
mv temp ${MYBADFOLDER}/GoNL_${MYBADSNPS}
MYBADSNPS=GoNL_${MYBADSNPS}

}
#1) create background set (stratify background snps by frequency )
{
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
mkdir -p $MYFOLDER/anal
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
MYFILE=chr${MYCHROM}.vcf.gz
echo $MYFILE
IPS1=$(( $IPS + 1 ))
bedtools intersect -a $MYFOLDER/anal/byfreq_${IPS}.bed.gz -b /mnt/sequencedb/GoNL/gonl.chr${MYCHROM}.snps_indels.r5.vcf.gz -wo | grep PASS | awk -v OFS='\t' '{print $1,$2,$3}' >> $MYFOLDER/anal/byfreq_${IPS}_GoNL.bed
done;
gzip -f $MYFOLDER/anal/byfreq_${IPS}_GoNL.bed
done
}
#2) create background set ( overlap stratified background snps with CG)
{
#COMPLETE GENOMICS
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
bedtools intersect -a <( zcat $MYFOLDER/anal/byfreq_${IPS}_GoNL.bed.gz | sort -k1,1 -k2,2n ) -b $CGBED_CHR -sorted | uniq >> $MYFOLDER/anal/byfreqCG_${IPS}_GoNL.bed
done;
#gzip -f $MYFOLDER/anal/byfreqCG_${IPS}.bed
bedtools sort -i $MYFOLDER/anal/byfreqCG_${IPS}_GoNL.bed  | uniq | gzip -f > $MYFOLDER/anal/byfreqCG_${IPS}_GoNL.bed.gz
rm $MYFOLDER/anal/byfreqCG_${IPS}.bed
done
}
#step 3 and 4 differ if merged
#---non merged---all bad snps are included in background vcf files
{
#3) overlap 
#sort -k1,1n -k2,2n ${MYFOLDER}/${MYBADSNPS} > temp
#cp temp ${MYFOLDER}/${MYBADSNPS}
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
for IPS in `seq 0 11`;do
echo $IPS
bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b $MYFOLDER/anal/byfreqCG_${IPS}_GoNL.bed.gz | uniq | gzip -f > $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz
bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b $MYFOLDER/anal/byfreq_${IPS}_GoNL.bed.gz | uniq | gzip -f > $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz
done

rm $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
for IPS in `seq 0 11`;do
TOT_SNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_GoNL.bed.gz | uniq | wc -l  )
TOT_SNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_GoNL.bed.gz | uniq | wc -l  )
TOT_BADSNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | uniq | wc -l  )
TOT_BADSNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz | uniq | wc -l  )
echo ${MYPS[${IPS}]} ${TOT_SNPS_CG[$IPS]} ${TOT_SNPS[$IPS]} ${TOT_BADSNPS_CG[$IPS]} ${TOT_BADSNPS[$IPS]} >> $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
done
Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

}
#---merged
{
#for intergenic I have to first extract my variants from 1000g

for MYCHROM in `seq 1 22`;do
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -R ${MYBADFOLDER}/${MYBADSNPS} -h | gzip -f > ${MYDESTINATION_FOLDER}/chr${MYCHROM}_GoNL.vcf.gz
done

#3) stratify badsnps by frequency 
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
MYFILE=chr${MYCHROM}_GoNL.vcf.gz
echo $MYFILE
IPS1=$(( $IPS + 1 ))
bcftools view -q ${MYPS[$IPS]}:minor -Q ${MYPS[$IPS1]}:minor  ${MYDESTINATION_FOLDER}/$MYFILE  | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2}'  >> ${MYFOLDER}/anal/byfreq_${IPS}_${MYBADSNPS}
done;
gzip -f ${MYFOLDER}/anal/byfreq_${IPS}_${MYBADSNPS}
done


#4) overlap background set ( overlap stratified background snps with CG)
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
bedtools intersect -a <( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | sort -k1,1 -k2,2n ) -b $CGBED_CHR -sorted >> $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}
done;
gzip -f $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}
done

rm $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
for IPS in `seq 0 11`;do
TOT_SNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_GoNL.bed.gz | uniq | wc -l  )
TOT_SNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_GoNL.bed.gz | uniq | wc -l  )
TOT_BADSNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | uniq | wc -l  )
TOT_BADSNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz | uniq | wc -l  )
echo ${MYPS[${IPS}]} ${TOT_SNPS_CG[$IPS]} ${TOT_SNPS[$IPS]} ${TOT_BADSNPS_CG[$IPS]} ${TOT_BADSNPS[$IPS]} >> $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
done
Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

#chisquare
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal/byfreq_summary_GoNL_bad_snps_logpop_nAB_pos_glmm_all.bed.tab #p-value = 0.04575 #fdr_total 0.5132398
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/anal/byfreq_summary_GoNL_bad_snps_logpop_nAB_pos_glmm_all.bed.tab #p-value = 0.03506 #fdr_total 0.4596427
#/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/anal/byfreq_summary_GoNL_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab #p-value < 2.2e-16 #fdr_total: 0.8359624
#/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/anal/byfreq_summary_GoNL_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab #p-value < 2.2e-16  #fdr total: 0.7324062


}
#calculate power of better datasets in removing bad snps
{

MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
mkdir -p $MYFOLDER/anal
for IPS in `seq 0 11`;do
echo $IPS
IPS1=$(( $IPS + 1 ))
zcat $MYFOLDER/anal/byfreq_${IPS}.bed.gz | grep -v '#' | uniq | wc -l  > $MYFOLDER/anal/byfreq_${IPS}_tot
done;

for IPS in `seq 0 11`;do
TOT_SNPS[$IPS]=$( cat $MYFOLDER/anal/byfreq_${IPS}_tot )
TOT_SNPS_OVERLAP[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_GoNL.bed.gz | uniq | wc -l  )
TOT_SNPS[$IPS]=$(( ${TOT_SNPS[$IPS]} - ${TOT_SNPS_OVERLAP[$IPS]} ))
TOT_BADSNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | uniq | wc -l  )
BADSNPS_OVERLAP[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_GoNL_${MYBADSNPS}.gz | uniq | wc -l )
TOT_BADSNPS[$IPS]=$(( ${TOT_BADSNPS[$IPS]} - ${BADSNPS_OVERLAP[$IPS]} ))
echo ${MYPS[${IPS}]} ${TOT_SNPS_OVERLAP[$IPS]} ${TOT_SNPS[$IPS]} ${BADSNPS_OVERLAP[$IPS]} ${TOT_BADSNPS[$IPS]} >> temp
done

for IPS in `seq 0 11`;do
TOT_SNPS[$IPS]=$( cat $MYFOLDER/anal/byfreq_${IPS}_tot )
TOT_SNPS_OVERLAP[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_HRC.bed.gz | uniq | wc -l  )
TOT_SNPS[$IPS]=$(( ${TOT_SNPS[$IPS]} - ${TOT_SNPS_OVERLAP[$IPS]} ))
TOT_BADSNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | uniq | wc -l  )
BADSNPS_OVERLAP[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_overlapHRC.bed.gz | uniq | wc -l )
TOT_BADSNPS[$IPS]=$(( ${TOT_BADSNPS[$IPS]} - ${BADSNPS_OVERLAP[$IPS]} ))
echo ${MYPS[${IPS}]} ${TOT_SNPS_OVERLAP[$IPS]} ${TOT_SNPS[$IPS]} ${BADSNPS_OVERLAP[$IPS]} ${TOT_BADSNPS[$IPS]} >> temp
done
Rscript ~/Dropbox/LDLD/scripts/overlap_badsnps.R temp


}
}
#HRC vs Complete Genomics
{

#bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz > overlapHRC_exons_min1.bed #overlapHRC_exons_min5.bed overlapHRC_intergenic_min1.bed overlapHRC_intergenic_min5.bed
cp /mnt/sequencedb/HRC/overlapHRC_exons_min1.bed /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/bad_snps_glmm/overlapHRC.bed
cp /mnt/sequencedb/HRC/overlapHRC_exons_min5.bed /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/bad_snps_glmm/overlapHRC.bed
cp /mnt/sequencedb/HRC/overlapHRC_intergenic_min1.bed /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/bad_snps_glmm/overlapHRC.bed
cp /mnt/sequencedb/HRC/overlapHRC_intergenic_min5.bed /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/bad_snps_glmm/overlapHRC.bed

#====SET UP VARIABLES======
{
#non merged------------------------
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5 

MYFOLDERBACKGROUND=$MYFOLDER
MYBADFOLDER=${MYFOLDER}/bad_snps_glmm
MYBADSNPS=overlapHRC.bed

sortBed -i ${MYBADFOLDER}/${MYBADSNPS} > temp
mv temp ${MYBADFOLDER}/${MYBADSNPS}

#merged----------------------------
MYMIN=1 #MYMIN=5
MYFOLDER=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min${MYMIN}
MYFOLDERBACKGROUND=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl1/min${MYMIN}
MYBADFOLDER=${MYFOLDER}/bad_snps_glmm
MYBADSNPS=overlapHRC.bed
MYDESTINATION_FOLDER=${MYFOLDER}/anal
mkdir -p $MYDESTINATION_FOLDER

#Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

sortBed -i ${MYBADFOLDER}/${MYBADSNPS} > temp
mv temp ${MYBADFOLDER}/${MYBADSNPS}

cd /mnt/sequencedb/HRC/

}
#1) create background set (stratify background snps by frequency )
{
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
mkdir -p $MYFOLDER/anal
for IPS in `seq 0 11`;do
for MYCHROM in `seq 1 22`;do
echo $IPS
IPS1=$(( $IPS + 1 ))
bedtools intersect -a $MYFOLDER/anal/byfreq_${IPS}.bed.gz -b <( tabix /mnt/sequencedb/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz ${MYCHROM} | awk -v OFS='\t' '{print $1,$2-1,$2}' ) | awk -v OFS='\t' '{print $1,$2,$3}' >> $MYFOLDER/anal/byfreq_${IPS}_HRC.bed
done
sortBed -i $MYFOLDER/anal/byfreq_${IPS}_HRC.bed | gzip -f > $MYFOLDER/anal/byfreq_${IPS}_HRC.bed.gz
done
}
#2) create background set ( overlap stratified background snps with CG)
{
#COMPLETE GENOMICS
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
bedtools intersect -a <( zcat $MYFOLDER/anal/byfreq_${IPS}_HRC.bed.gz | sort -k1,1 -k2,2n ) -b $CGBED_CHR -sorted | uniq >> $MYFOLDER/anal/byfreqCG_${IPS}_HRC.bed
done;
#gzip -f $MYFOLDER/anal/byfreqCG_${IPS}.bed
bedtools sort -i $MYFOLDER/anal/byfreqCG_${IPS}_HRC.bed  | uniq | gzip -f > $MYFOLDER/anal/byfreqCG_${IPS}_HRC.bed.gz
#rm $MYFOLDER/anal/byfreqCG_${IPS}.bed
done
}
#step 3 and 4 differ 
#---non merged---all bad snps are included in background vcf files
{
#3) overlap 
#sort -k1,1n -k2,2n ${MYFOLDER}/${MYBADSNPS} > temp
#cp temp ${MYFOLDER}/${MYBADSNPS}
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
for IPS in `seq 0 11`;do
echo $IPS
bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b $MYFOLDER/anal/byfreqCG_${IPS}_HRC.bed.gz | uniq | gzip -f > $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz
bedtools intersect -a ${MYBADFOLDER}/${MYBADSNPS} -b $MYFOLDER/anal/byfreq_${IPS}_HRC.bed.gz | uniq | gzip -f > $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz
done

rm $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
for IPS in `seq 0 11`;do
TOT_SNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_HRC.bed.gz | wc -l )
TOT_SNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_HRC.bed.gz | wc -l )
TOT_BADSNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | wc -l )
TOT_BADSNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz | wc -l )
echo ${MYPS[${IPS}]} ${TOT_SNPS_CG[$IPS]} ${TOT_SNPS[$IPS]} ${TOT_BADSNPS_CG[$IPS]} ${TOT_BADSNPS[$IPS]} >> $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
done
Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

}
#---merged
{
#for intergenic I have to first extract my variants from 1000g

for MYCHROM in `seq 1 22`;do
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -R ${MYBADFOLDER}/${MYBADSNPS} -h | gzip -f > ${MYDESTINATION_FOLDER}/chr${MYCHROM}_HRC.vcf.gz
done

#3) stratify badsnps by frequency 
MYPS[0]=0
MYPS[1]=0.01
MYPS[2]=0.02
MYPS[3]=0.03 
MYPS[4]=0.05 
MYPS[5]=0.075 
MYPS[6]=0.1 
MYPS[7]=0.15 
MYPS[8]=0.2 
MYPS[9]=0.25 
MYPS[10]=0.3 
MYPS[11]=0.4 
MYPS[12]=0.5
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
MYFILE=chr${MYCHROM}_HRC.vcf.gz
echo $MYFILE
IPS1=$(( $IPS + 1 ))
bcftools view -q ${MYPS[$IPS]}:minor -Q ${MYPS[$IPS1]}:minor  ${MYDESTINATION_FOLDER}/$MYFILE  | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2}'  >> ${MYFOLDER}/anal/byfreq_${IPS}_${MYBADSNPS}
done;
gzip -f ${MYFOLDER}/anal/byfreq_${IPS}_${MYBADSNPS}
done


#4) overlap background set ( overlap stratified background snps with CG)
for IPS in `seq 0 11`;do
echo $IPS
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
bedtools intersect -a <( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | sort -k1,1 -k2,2n ) -b $CGBED_CHR -sorted >> $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}
done;
gzip -f $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}
done

rm $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
for IPS in `seq 0 11`;do
TOT_SNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}.bed.gz | wc -l )
TOT_SNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}.bed.gz | wc -l )
TOT_BADSNPS[$IPS]=$( zcat $MYFOLDER/anal/byfreq_${IPS}_${MYBADSNPS}.gz | wc -l )
TOT_BADSNPS_CG[$IPS]=$( zcat $MYFOLDER/anal/byfreqCG_${IPS}_${MYBADSNPS}.gz | wc -l )
echo ${MYPS[${IPS}]} ${TOT_SNPS_CG[$IPS]} ${TOT_SNPS[$IPS]} ${TOT_BADSNPS_CG[$IPS]} ${TOT_BADSNPS[$IPS]} >> $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab
done
Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

#chisquare
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal/byfreq_summary_overlapHRC.bed.tab #p-value = 8.993e-06 #total_fdr 0.3159992
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/anal/byfreq_summary_overlapHRC.bed.tab #p-value = p-value = 0.00263 #total_fdr 0.4131174
#/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/anal/byfreq_summary_overlapHRC.bed.tab #pvalue 2.2e-16 #fdr_total 0.7424441
#/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/anal/byfreq_summary_overlapHRC.bed.tab #p-value = 2.2e-16 #total_fdr 0.7869323



}
#plots
{
}















}
#additional functions
if (FALSE)
{
# I could have trick: check if at least half 2x difference between nAB upper and lower half and calculate for those.

#system.time(res1<-sapply(1:500,function(i) myfx(i))) #800 sec: I have to shorten it.
#save(res1,file=paste0(myfolder,"/res1.RData"))
#load(paste0(myfolder,"/res1.RData"))

#cor.test gives similar results but not faster
logreg_snps2nAB<-function(mygen,mynAB)#calculate logistic regression to see snps that contribute the most to individuals with highest nAB
{   #roger suggests on comparing the models because even if "perfect separation" and p-value with Wald approximation is useless, still the fit is very good and likelihood ratio works perfectly
    if ( sum(mygen)/length(mygen) > 0.05 & sum(mygen)/length(mygen) < 0.95 ) {
    model <- glm( mygen ~ mynAB,family=binomial(link='logit'))
    model0 <- glm( mygen ~ 1,family=binomial(link='logit'))
    return(anova(model0, model, test = "Chisq")[[5]][2])
    } else return(NA)
}    
    
logreg_snps2nAB<-function(mygen,mynAB,mymut)#calculate logistic regression to see snps that contribute the most to individuals with highest nAB
{   #roger suggests on comparing the models because even if "perfect separation" and p-value with Wald approximation is useless, still the fit is very good and likelihood ratio works perfectly
    model <- glm( mygen ~ mynAB + mymut,family=binomial(link='logit'))
    model0 <- glm( mygen ~ mymut,family=binomial(link='logit'))
    return(anova(model0, model, test = "Chisq")[[5]][2])
}


check_nAvsnAB<-function(mygen,mynAB){
    if ( sum(mygen)/length(mygen) > 0.05 & sum(mygen)/length(mygen) < 0.95 ) {
    abs(mean(mygen[order(mynAB)][1:round(length(mynAB)/2)])-mean(mygen))/mean(mygen)
} else return(NA)
}

}
if (FALSE){
if (FALSE)
{
#old version
#load("fullsign.RData")
fullsign_snps<-unique(paste0(c(fullsign[,1],fullsign[,2]),".",c(fullsign[,3],fullsign[,4])));
length(fullsign_snps)
resres_fullsign<-list()
for (ipop in 1:12)
    {
    resres_fullsign[[ipop]]<-subset(resres[[ipop]],!is.na(match(resresnames,fullsign_snps))) #only significant    
    }
#identify biased SNPs 
#(SNPs mostly present in nAB rich individuals)
tempres<-corr_nAvsnAB(resres,logpop,l_mutsamples) #this is done by looking at all SNPs.
tempres_sign<-corr_nAvsnAB(resres_fullsign,logpop,l_mutsamples,snpsnames=fullsign_snps) #this is done by looking only at SNPs that are linked.
#save.image("resres_complete.RData")

#load("resres_complete.RData")

#more or less the same number of biased snps from first filtering in a way or the other. Difference might become
#more important when less total significant SNPs are considered
dim(tempres)
dim(tempres_sign)
dim(tempres[tempres[,3]<=0.05,])
dim(tempres[tempres[,2]<=0.05,])
removed_snps<-tempres[tempres[,3]<=0.05,]
#dim(tempres_sign[tempres_sign[,3]<=0.05,])
#dim(tempres_sign[tempres_sign[,2]<=0.05,])
#removed_snpssign<-tempres_sign[tempres_sign[,3]<=0.05,]

#create BED file with removed snps
write.table(stringsnp2bed(removed_snps),file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed1.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)

}
#remove snps from data for next cycle
system(paste0("mkdir ",myfolder,"/removal1"))
for (mychr in 1:22)
  {
  system(paste0("bedtools intersect -v -b ",myfolder,"/removed1.bed -a ",myfolder,"/chr",mychr,".tab -header | gzip > ",myfolder,"/removal1/chr",mychr,".coding.vcf.gz"))
  }
  rm removed1dbsnps.bed
  touch removed1dbsnps.bed
  for i in `seq 1 22`; do 
  echo $i
  zcat /mnt/sequencedb/ucsc/goldenPath/hg19/snp142.bed.gz | grep chr$i$'\t' > /mnt/scratch/fabrizio/LDLD/temp.bed
  bedtools intersect -a <(sort -Vu -k1,1 -k2,2 ${myfolder}/removed1.bed | sed 's/^/chr/g' | grep chr${i}) -b <(sort -Vu -k1,1 -k2,2 /mnt/scratch/fabrizio/LDLD/temp.bed | grep chr${i} ) -sorted -wb | awk -v OFS='\t' '{print $4,$5,$6,$7}' >>  ${myfolder}/removed1dbsnps.bed
  done
  #bedtools intersect -a <(cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed1.bed | sed 's/^/chr/g') -b /mnt/sequencedb/ucsc/goldenPath/hg19/snp142.bed.gz -sorted | awk -v OFS='\t' '{print $4,$5,$6,$7}' >  /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed.bed
}
#table comparing tests before and after removal
if (FALSE)
{
library(xtable)
load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/testfull.RData")
testfull_mat<-rbind(popnames,
c("p.sd",sapply(1:12,function(x) testfull[[x]]$pval_sd)),
c("p.nLR",sapply(1:12,function(x) testfull[[x]]$pval_nLR)),
c("p.nAIC",sapply(1:12,function(x) testfull[[x]]$pval_nAIC)),
c("p.pLR",sapply(1:12,function(x) testfull[[x]]$pval_pLR)),
c("p.pAIC",sapply(1:12,function(x) testfull[[x]]$pval_pAIC)))

load("removal1_test_bis.RData")
test_removal1_mat<-rbind(popnames,
c("p.sd",sapply(1:12,function(x) test_removal1[[x]]$pval_sd)),
c("p.nLR",sapply(1:12,function(x) test_removal1[[x]]$pval_nLR)),
c("p.nAIC",sapply(1:12,function(x) test_removal1[[x]]$pval_nAIC)),
c("p.pLR",sapply(1:12,function(x) test_removal1[[x]]$pval_pLR)),
c("p.pAIC",sapply(1:12,function(x) test_removal1[[x]]$pval_pAIC)))

xtable(testfull_mat)
xtable(test_removal1_mat)
}
#reads imbalance
{
#TAG 'CREATE LOGOS'
#cat above95/coding1000g/anal/snpsA.bed above95/coding1000g/anal/snpsB.bed | sort -Vu -k1,1 -k2,2 > above95/coding1000g/anal/snps.bed
mydata<-read.table("~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt",stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
dim(data) #47193
createbedlogo(mydata,paste0(myfolder,"/basessign"))
paste <( bedtools intersect -b ~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 ) ) /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign | awk '{if (length($4)==1 && length($5)==1){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign.bed
bamlogo_sign<-bamlogo("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign.bed",tilltheend=TRUE,till=100,startfrom=1,filetemp="temp")
#save(bamlogo_sign,file=paste0(myfolder,"bamlogo_sign.RData"))
load(paste0(myfolder,"bamlogo_sign.RData"))
bamlogo_sign<-as.data.frame(bamlogo_sign)
names(bamlogo_sign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
bamlogo_sign<-bamlogo_sign[bamlogo_sign$A!="NaN",]
bamlogo_sign$A<-as.numeric(as.character(bamlogo_sign$A))
bamlogo_sign$C<-as.numeric(as.character(bamlogo_sign$C))
bamlogo_sign$G<-as.numeric(as.character(bamlogo_sign$G))
bamlogo_sign$T<-as.numeric(as.character(bamlogo_sign$T))
bamlogo_sign$tot_second_call<-as.numeric(as.character(bamlogo_sign$tot_second_call))
bamlogo_sign$ref<-as.character(bamlogo_sign$ref)
bamlogo_sign$alt<-as.character(bamlogo_sign$alt)
imbalance_sign<-cbind(bamlogo_sign,t(sapply(1:dim(bamlogo_sign)[1],function(x) prHvsE(bamlogo_sign,x))))
imbalance_sign$pvalue_no05_alt[imbalance_sign$pvalue_no05_alt<=1]

dim(logonotsign)
res<-logo_below_thr(logonotsign,0.00001) #0.03417355


#background set
bedtools subtract -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' ) -b ~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt | shuf -n 1000 | sort -Vu -k1,1 -k2,2 > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign.bed 
mydata<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign.bed",stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
createbedlogo(mydata,paste0(myfolder,"/basesnotsign"))
#to do
paste <( bedtools intersect -b /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign.bed -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 ) ) /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign | awk '{if (length($4)==1 && length($5)==1){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign.bed
bamlogo_notsign<-bamlogo("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign.bed",tilltheend=TRUE,till=100,startfrom=1,filetemp="tempnot")
#save(bamlogo_notsign,file=paste0(myfolder,"bamlogo_notsign.RData"))
load(paste0(myfolder,"bamlogo_notsign.RData"))
bamlogo_notsign<-as.data.frame(bamlogo_notsign)
names(bamlogo_notsign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
bamlogo_notsign<-bamlogo_notsign[bamlogo_notsign$A!="NaN",]
bamlogo_notsign$A<-as.numeric(as.character(bamlogo_notsign$A))
bamlogo_notsign$C<-as.numeric(as.character(bamlogo_notsign$C))
bamlogo_notsign$G<-as.numeric(as.character(bamlogo_notsign$G))
bamlogo_notsign$T<-as.numeric(as.character(bamlogo_notsign$T))
bamlogo_notsign$tot_second_call<-as.numeric(as.character(bamlogo_notsign$tot_second_call))
bamlogo_notsign$ref<-as.character(bamlogo_notsign$ref)
bamlogo_notsign$alt<-as.character(bamlogo_notsign$alt)
imbalance_notsign<-cbind(bamlogo_notsign,t(sapply(1:dim(bamlogo_notsign)[1],function(x) prHvsE(bamlogo_notsign,x))))

imbalance_notsign<-imbalance_notsign[imbalance_notsign$pvalue_no05_alt<=1,]
imbalance_sign<-imbalance_sign[imbalance_sign$pvalue_no05_alt<=1,]
dat_sign<-log(imbalance_sign$pvalue_no05_alt)
dat_sign[dat_sign<(-250)]<-(-250)
dat_notsign<-log(imbalance_notsign$pvalue_no05_alt)
dat_notsign[dat_notsign<(-250)]<-(-250)
mybreaks<-seq(-250,0,25)
myhistsign<-hist(dat_sign,breaks=mybreaks)
myhistnotsign<-hist(dat_notsign,breaks=mybreaks)
myhistsign$counts<-myhistsign$counts/sum(myhistsign$counts)
myhistnotsign$counts<-myhistnotsign$counts/sum(myhistnotsign$counts)
pdf(paste0(myfolderdropbox,"/hist_allelic_imbalance.pdf"))
data_temp<-c(rbind(myhistsign$count,myhistnotsign$count))
myax=data.frame(pos=seq(1,19,2),lab=mybreaks[1:10])
mybarplot_f(1:length(data_temp),data_temp,myylim=c(0,1),mycols=rep(c("gold","cadetblue"),length(data_temp)),myxlabel='log(p-value)', myylabel='density',add=FALSE,xaxis=myax,cex.axis=0.9)
legend(0,0.6,c("candidates","all"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        pch=(c(15,15)),col=c("gold","cadetblue"))
dev.off()


myhistsign<-hist(imbalance_sign$pvalue_no05_alt,breaks=mybreaks)
myhistnotsign<-hist(imbalance_notsign$pvalue_no05_alt,breaks=mybreaks)


mybreaks<-seq(0,1,0.05)

dat_notsign


myhistsign<-hist(imbalance_sign$pvalue_no05_alt,breaks=mybreaks)
myhistnotsign<-hist(imbalance_notsign$pvalue_no05_alt,breaks=mybreaks)

dev.off()
wilcox.test(imbalance_sign$pvalue_no05_alt,imbalance_notsign$pvalue_no05_alt)

for i in /mnt/454/HGDP/genomes_Bteam/HG19Align/*.bam /mnt/sequencedb/11men/martin/hg19_1000g/BAM/*.bam ; do /r1/people/pruefer/bin/BamTable2 -r ~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt;done

logonotsign[is.na(logonotsign)] <- 0
logonotsign<-cbind(logonotsign,t(sapply(1:dim(logonotsign)[1],function(x) prHvsE(logonotsign,x))))
dim(logonotsign)
res<-logo_below_thr(logonotsign,0.00001) #0.03417355
res<-logo_below_thr(logonotsign,0.05) #0.1843388


}

#============================
#===== TESTS for GBE revisions - 12/04/18==
#============================
{

#total sum of sites in frequencies
{

mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_totfreq
for mysupp in 6 7 4 5;do
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
count_n=0;count_n20=0
for MYCHR in `seq 1 22`;do
echo $MYCHR
source ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g_2awk.sh
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
count_n1=$( zcat $my1000gfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq |wc -l)
count_n=$(( $count_n + $count_n1 ));
done
echo "count in supp " $mysupp " is " $count_n >> $mylogfile
done

#bed to speed up next analyses
for myminfreq in 1 5;do
for MYCHR in `seq 1 22`;do
echo $MYCHR
source ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g_2awk.sh
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
zcat $my1000gfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mask.freq${myminfreq}.chr${MYCHR}.bed.gz
done
done

#check that TableS4.bed is from freq 5% and TableS5.bed from 1%
for mysupp in 6 7;do
for myminfreq in 1 5;do
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz
temp1000gfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/temp.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed
tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min${myminfreq}/chr${MYCHR}.vcf.gz -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5}' | sort -k 2,2n | uniq > $mybedfile
bedtools intersect -a $mybedfile -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed > temp; mv temp $mybedfile
tabix $my1000gfile -R $mybedfile | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8}' > $temp1000gfile
cat $temp1000gfile |wc -l
done;done
#139
#139
#139
#135
#that shows that in TableS5 there are some snps that no in minfreq 5

MYCHR=1
for mysupp in 6 7 4 5;do
for myminfreq in 1 5;do
mybedfile=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mask.freq${myminfreq}.chr${MYCHR}.bed.gz 
bedtools intersect -a $mybedfile -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed |wc -l
done;done
#from this it seems that TableS6 is 1% and TableS7 5%, while TableS4 is 5% and TableS5 is 1%. This means that I wrote well on the manuscript and inverted the files -> I fix by reinverting them so that TableS6 has 5% and viceversa
#415
#414
#417
#417
#102
#102
#105
#102



}

#N/S
{
#I directly use BED files attached in papers
#/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE
#TableS4.bed and TableS5.bed are exons
#TableS6.bed and TableS7.bed are intergenic
cd /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE
wc -l TableS4.bed #700
#cat /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS4.bed | grep -v '#' | awk '{print $2}' | uniq |wc -l #696
#clean up
for i in TableS4.bed TableS5.bed TableS6.bed TableS7.bed;do
#cat $i | grep '#' > temp;
cp revisions/header temp
#cat $i | grep -v '#' | awk 'BEGIN{var==0;var1==0}{if ( $1==1 && var==1){var1=1};if ($1!=1){var=1}; if (var1==1){print} }' >> temp
cat $i | grep -v '#' | sort -k1,1 -k2,2n | awk 'BEGIN{var=0}{if ($2!=var){print};var=$2}' | uniq >> temp
mv temp $i
done
#now no doubles



#notice that in this vcf ( /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz ) there are multiple alternative alleles
#in my analyses I always only analyzed the first (1 in vcf), thus I only extract this
for MYCHR in `seq 1 22`;do 
done

count_n=0;count_s=0;count_nc=0;
for MYCHR in `seq 1 22`;do
echo $MYCHR
tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS4.bed | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5}' | sort -k 2,2n | uniq | sed 's/,.*//g' > ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed; 
count_n1=$( tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz -R ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed | grep missense | wc -l ) 
count_s1=$( tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz -R ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed | grep -v missense | grep synonymous | wc -l ) 
count_nc1=$( tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz -R ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed | grep -v missense | grep -v synonymous | grep coding | wc -l ) 
count_n=$(( $count_n + $count_n1 ))
count_s=$(( $count_s + $count_s1 ))
count_nc=$(( $count_nc + $count_nc1 ))
done
echo $count_n $count_s $count_nc
# numbers do not sum to 696 because more variants at same position in functional_annotation file of 1000 genomes

#I notice very strange phenomenon for which tabix takes some extra variants, and I am sure that the bed file is correct.
#tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS4.bed | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5}' | sort -k 2,2n | uniq > $mybedfile
#so I take them and then parse it again.
#here the criterion is most severe
for mysupp in 4 5;do
if [ $mysupp -eq 4 ];then
    myminfreq=5
else 
    myminfreq=1
fi
count_n=0;count_s=0;count_nc=0;count_tot=0
count_n20=0;count_s20=0;count_nc20=0;count_tot20=0
for MYCHR in `seq 1 22`;do
echo $MYCHR
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz
temp1000gfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/temp.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed
tabix /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min${myminfreq}/chr${MYCHR}.vcf.gz -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5}' | sort -k 2,2n | uniq > $mybedfile
bedtools intersect -a $mybedfile -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed > temp; mv temp $mybedfile
tabix $my1000gfile -R $mybedfile | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8}' > $temp1000gfile
count_n1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep missense | wc -l ) 
count_s1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep synonymous | wc -l ) 
count_nc1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep coding | wc -l )
count_n2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep missense | grep -v synonymous | wc -l ) 
count_s2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep synonymous | wc -l ) 
count_nc2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep coding | wc -l )
bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep -v coding >> revisions/bed_ns/others
count_tot1=$( cat ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed  |wc -l  )
count_n=$(( $count_n + $count_n1 )); count_s=$(( $count_s + $count_s1 )); count_nc=$(( $count_nc + $count_nc1 ));count_tot=$(( $count_tot + $count_tot1 ))
count_n20=$(( $count_n20 + $count_n2 )); count_s20=$(( $count_s20 + $count_s2 )); count_nc20=$(( $count_nc20 + $count_nc2 ));count_tot20=$(( $count_tot20 + $count_tot1 ))
done
echo $count_n $count_s $count_nc $count_tot #467 166 38 696
echo $count_n20 $count_s20 $count_nc20 $count_tot20 #461 166 38 696
done

for mysupp in 4 5;do
if [ $mysupp -eq 4 ];then
    myminfreq=5
else 
    myminfreq=1
fi
count_n=0;count_s=0;count_nc=0;count_tot=0
count_n20=0;count_s20=0;count_nc20=0;count_tot20=0
for MYCHR in `seq 1 22`;do
echo $MYCHR
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz
temp1000gfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/temp.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed
tabix $my1000gfile -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5}' | sort -k 2,2n | uniq > $mybedfile
bedtools intersect -a $mybedfile -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed > temp; mv temp $mybedfile
tabix $my1000gfile -R $mybedfile | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8}' > $temp1000gfile
count_n1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep missense | wc -l ) 
count_s1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep synonymous | wc -l ) 
count_nc1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep coding | wc -l )
count_n2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep missense | grep -v synonymous | wc -l ) 
count_s2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep synonymous | wc -l ) 
count_nc2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep coding | wc -l )
bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep -v coding >> revisions/bed_ns/others
count_tot1=$( cat ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed  |wc -l  )
count_n=$(( $count_n + $count_n1 )); count_s=$(( $count_s + $count_s1 )); count_nc=$(( $count_nc + $count_nc1 ));count_tot=$(( $count_tot + $count_tot1 ))
count_n20=$(( $count_n20 + $count_n2 )); count_s20=$(( $count_s20 + $count_s2 )); count_nc20=$(( $count_nc20 + $count_nc2 ));count_tot20=$(( $count_tot20 + $count_tot1 ))
done
echo $count_n $count_s $count_nc $count_tot 
echo $count_n20 $count_s20 $count_nc20 $count_tot20 
done
#467 166 38 696
#461 166 38 696
#511 171 37 744
#501 171 37 744
}

#mappability
{
for mymap in 50 99;do
for mysupp in `seq 4 7`;do
mymappabilityfile=/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_${mymap}.bed.gz
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed
tabix $mymappabilityfile -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | sort -k 1,1 -k 2,2n | uniq > $mybedfile
bedtools intersect -a $mybedfile -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l
done;done
#531
#572
#3141
#10661
#456
#498
#1673
#6192

#from parser of 1000 genomes for my populations ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g.sh
mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_map
for mymap in 50 99;do
for mysupp in 6 7 4 5;do
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
count_n=0;count_n20=0
for MYCHR in `seq 1 22`;do
echo $MYCHR
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed.gz
source ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g_2awk.sh
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
mymappabilityfile=/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_${mymap}.bed.gz
#tabix $my1000gfile -R $mymappabilityfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
bedtools intersect -a $my1000gfile -b <( tabix $mymappabilityfile $MYCHR ) -sorted | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
count_n2=$(cat $mybedfile | wc -l )
count_n1=$( bedtools intersect -a $mybedfile -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l )
count_n=$(( $count_n + $count_n1 ));
count_n20=$(( $count_n20 + $count_n2 )); 
echo "count is " $count_n2 >> $mylogfile
done
echo "n overlaps with map in TableS" $mysupp "with map" $mymap " are" $count_n $count_n20 >> $mylogfile.$mymap
done
done

}

#accessibility masks 1000 genomes
{
cd /mnt/scratch/fabrizio/LDLD/above95/downloads
wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.pilot_mask.autosomes.bed' #in release..supporting
for mychr in `seq 1 22`;do
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/20140520.chr${mychr}.strict_mask.fasta.gz"
done
for mychr in `seq 1 22`;do
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/PilotMask/20140520.chr${mychr}.pilot_mask.fasta.gz"
done

flatten () {
awk '{if (substr($1,1,1)==">"){convert=0;for (k = 1; k <= 9; k++){if (substr($1,2,1)==k || substr($1,2,1)=="X"){convert=convert+1;counter=0;}}; 
if (convert>0){print $1}} else if (convert>0){for (i = 1; i <= length($1); i++){counter=counter+1;print substr($1,i,1);}}}'
}

awk -v myvar=${letter}  '{if ($4==myvar){print}}' |

for letter in 'H' 'P' 'L' 'N' 'ZQ';do  
for mychr in `seq 1 22`;do
zcat /mnt/scratch/fabrizio/LDLD/above95/downloads/20140520.chr${mychr}.strict_mask.fasta.gz | flatten | awk -v OFS='\t' -v MYCHR=${mychr} 'BEGIN{COUNTER=0}{if (substr($1,1,1)!=">"){COUNTER=COUNTER+1;print MYCHR,COUNTER-1,COUNTER,$1}}' |  grep [${letter}] | gzip -f  > /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.${letter}.bed.gz &
bedtools merge -i /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.${letter}.bed.gz | gzip -f > /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz
done;done

for letter in 'H' 'P' 'L' 'N' 'ZQ';do  
for mychr in `seq 1 22`;do
zcat /mnt/scratch/fabrizio/LDLD/above95/downloads/20140520.chr${mychr}.pilot_mask.fasta.gz | flatten | awk -v OFS='\t' -v MYCHR=${mychr} 'BEGIN{COUNTER=0}{if (substr($1,1,1)!=">"){COUNTER=COUNTER+1;print MYCHR,COUNTER-1,COUNTER,$1}}' | grep [${letter}] |  gzip -f  > /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.pilot.${letter}.bed.gz
bedtools merge -i /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.pilot.${letter}.bed.gz | gzip -f > /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.pilot.mask.${letter}.bed.gz &
done;done


mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_accessibilitymap_pilot
for mysupp in 6 7;do
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
count_n=0;count_tot=0
for mychr in `seq 1 22`;do
echo $mychr
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mask.freq${myminfreq}.chr${mychr}.bed.gz
count_tot1=$( bedtools intersect -a $mybedfile -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.pilot.mask.${letter}.bed.gz | wc -l )
count_n1=$( bedtools intersect -a ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.pilot.mask.${letter}.bed.gz |wc -l )
count_tot=$(( $count_tot + count_tot1 ))
count_n=$(( $count_n + count_n1 ))
done;
echo $mysupp $letter $count_n $count_tot >> $mylogfile
done;done

mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_accessibilitymap
for mysupp in 6 7;do
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
count_n=0;count_tot=0
for mychr in `seq 1 22`;do
echo $mychr
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mask.freq${myminfreq}.chr${mychr}.bed.gz
count_tot1=$( bedtools intersect -a $mybedfile -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz | wc -l )
count_n1=$( bedtools intersect -a ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz |wc -l )
count_tot=$(( $count_tot + count_tot1 ))
count_n=$(( $count_n + count_n1 ))
done;
echo $mysupp $letter $count_n $count_tot >> $mylogfile
done;done


mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_accessibilitymap_pilot
for mysupp in 4 5;do
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
count_n=0;count_tot=0
for mychr in `seq 1 22`;do
echo $mychr
mybedfile=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr${mychr}.vcf.gz
count_tot1=$( bedtools intersect -a $mybedfile -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.pilot.mask.${letter}.bed.gz | wc -l )
count_n1=$( bedtools intersect -a ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.pilot.mask.${letter}.bed.gz |wc -l )
count_tot=$(( $count_tot + count_tot1 ))
count_n=$(( $count_n + count_n1 ))
done;
echo $mysupp $letter $count_n $count_tot >> $mylogfile
done;done


mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_accessibilitymap
for mysupp in 4 5;do
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
count_n=0;count_tot=0
for mychr in `seq 1 22`;do
echo $mychr
mybedfile=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr${mychr}.vcf.gz
count_tot1=$( bedtools intersect -a $mybedfile -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz | wc -l )
count_n1=$( bedtools intersect -a ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz |wc -l )
count_tot=$(( $count_tot + count_tot1 ))
count_n=$(( $count_n + count_n1 ))
done;
echo $mysupp $letter $count_n $count_tot >> $mylogfile
done;done



for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
for mychr in `seq 1 22`;do
echo $letter $mychr
bedtools intersect -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/chr${mychr}.vcf.gz -b /mnt/scratch/fabrizio/LDLD/above95/downloads/chr${mychr}.mask.${letter}.bed.gz -sorted |wc -l >> /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/overlap.all.$letter
done;done

for letter in 'H' 'P' 'L' 'N' 'ZQ';do 
cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/overlap.all.$letter | awk -v let=$letter 'BEGIN{co=0}{co=co+$1}END{print let,co}'
done
#H 3781
#P 230977
#L 12059
#N 0
#ZQ 70508
wc -l /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap*
   24 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.H.bed
  745 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.P.bed
     67 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.L.bed
    0 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.N.bed
  306 /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1/removed1.bed.overlap.ZQ.bed
 1142 total

 obs<-c(3781,230977,12059,70508)
 candidates<-c(24,745,67,306)
 obs/sum(obs)
candidates/sum(candidates)

#>  obs/sum(obs)
#[1] 0.01191523 0.72788781 0.03800205 0.22219491
#> candidates/sum(candidates)
#[1] 0.02101576 0.65236427 0.05866900 0.26795096
#we see slight enrichment in H
}

#base composition
{
#I found files /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5bamlogo_notsign.RData and min5bamlogo_sign.RData but not code that I used to generate them. I should look in my personal pc
setwd("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2")
load("min5bamlogo_sign.RData")
load("min5bamlogo_notsign.RData")
> table(bamlogo_sign[,4])/sum(table(bamlogo_sign[,4]))

        A         C         G         T 
0.2552239 0.2268657 0.2373134 0.2805970 
> table(bamlogo_sign[,5])/sum(table(bamlogo_sign[,5]))

        A         C         G         T 
0.1447761 0.3656716 0.3388060 0.1507463 
> table(bamlogo_notsign[,4])/sum(table(bamlogo_notsign[,4]))

        A         C         G         T 
0.1856540 0.2964135 0.3491561 0.1687764 
> table(bamlogo_notsign[,5])/sum(table(bamlogo_notsign[,5]))

        A         C         G         T 
0.3217300 0.1940928 0.1993671 0.2848101 
> table(bamlogo_notsign[,5])/sum(table(bamlogo_notsign[,5]))+table(bamlogo_notsign[,4])/sum(table(bamlogo_notsign[,4]))

        A         C         G         T 
0.5073840 0.4905063 0.5485232 0.4535865 
> table(bamlogo_sign[,4])/sum(table(bamlogo_sign[,4]))+table(bamlogo_sign[,5])/sum(table(bamlogo_sign[,5]))

        A         C         G         T 
0.4000000 0.5925373 0.5761194 0.4313433 

#TAG 'CREATE LOGOS'
#cat above95/coding1000g/anal/snpsA.bed above95/coding1000g/anal/snpsB.bed | sort -Vu -k1,1 -k2,2 > above95/coding1000g/anal/snps.bed
setwd(myfolder)
mydata<-read.table("~/Dropbox/LDLD/ms/GBE/TableS4.bed",stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
dim(data) #47193
createbedlogo(mydata,paste0(myfolder,"/basessign"))
#paste <( bedtools intersect -b ~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 ) ) /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign | awk '{if (length($4)==1 && length($5)==1){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/TableS4.bed
myfolder=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5
paste <( bedtools intersect -b $mybedfile -a <( cat ${myfolder}/chr*tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 ) ) ${myfolder}/basessign | awk '{if (length($4)==1 && length($5)==1){print}}' > ${myfolder}/basessign.bed
bamlogo_sign<-bamlogo(paste0(myfolder,"/basessign.bed"),tilltheend=TRUE,till=100,startfrom=1,filetemp="temp3")
save(bamlogo_sign,file=paste0(myfolder,"bamlogo_sign.RData"))
load(paste0(myfolder,"bamlogo_sign.RData"))
bamlogo_sign<-as.data.frame(bamlogo_sign)
names(bamlogo_sign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
bamlogo_sign<-bamlogo_sign[bamlogo_sign$A!="NaN",]
bamlogo_sign$A<-as.numeric(as.character(bamlogo_sign$A))
bamlogo_sign$C<-as.numeric(as.character(bamlogo_sign$C))
bamlogo_sign$G<-as.numeric(as.character(bamlogo_sign$G))
bamlogo_sign$T<-as.numeric(as.character(bamlogo_sign$T))
bamlogo_sign$tot_second_call<-as.numeric(as.character(bamlogo_sign$tot_second_call))
bamlogo_sign$ref<-as.character(bamlogo_sign$ref)
bamlogo_sign$alt<-as.character(bamlogo_sign$alt)
imbalance_sign<-cbind(bamlogo_sign,t(sapply(1:dim(bamlogo_sign)[1],function(x) prHvsE(bamlogo_sign,x))))
imbalance_sign$pvalue_no05_alt[imbalance_sign$pvalue_no05_alt<=1]

dim(logonotsign)
res<-logo_below_thr(logonotsign,0.00001) #0.03417355


#background set
bedtools subtract -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' ) -b ~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt | shuf -n 1000 | sort -Vu -k1,1 -k2,2 > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign.bed 
mydata<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign.bed",stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
createbedlogo(mydata,paste0(myfolder,"/basesnotsign"))
#to do
paste <( bedtools intersect -b /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign.bed -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 ) ) /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign | awk '{if (length($4)==1 && length($5)==1){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign.bed
bamlogo_notsign<-bamlogo("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign.bed",tilltheend=TRUE,till=100,startfrom=1,filetemp="tempnot")
#save(bamlogo_notsign,file=paste0(myfolder,"bamlogo_notsign.RData"))
load(paste0(myfolder,"bamlogo_notsign.RData"))
bamlogo_notsign<-as.data.frame(bamlogo_notsign)
names(bamlogo_notsign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
bamlogo_notsign<-bamlogo_notsign[bamlogo_notsign$A!="NaN",]
bamlogo_notsign$A<-as.numeric(as.character(bamlogo_notsign$A))
bamlogo_notsign$C<-as.numeric(as.character(bamlogo_notsign$C))
bamlogo_notsign$G<-as.numeric(as.character(bamlogo_notsign$G))
bamlogo_notsign$T<-as.numeric(as.character(bamlogo_notsign$T))
bamlogo_notsign$tot_second_call<-as.numeric(as.character(bamlogo_notsign$tot_second_call))
bamlogo_notsign$ref<-as.character(bamlogo_notsign$ref)
bamlogo_notsign$alt<-as.character(bamlogo_notsign$alt)
imbalance_notsign<-cbind(bamlogo_notsign,t(sapply(1:dim(bamlogo_notsign)[1],function(x) prHvsE(bamlogo_notsign,x))))

imbalance_notsign<-imbalance_notsign[imbalance_notsign$pvalue_no05_alt<=1,]
imbalance_sign<-imbalance_sign[imbalance_sign$pvalue_no05_alt<=1,]
dat_sign<-log(imbalance_sign$pvalue_no05_alt)
dat_sign[dat_sign<(-250)]<-(-250)
dat_notsign<-log(imbalance_notsign$pvalue_no05_alt)
dat_notsign[dat_notsign<(-250)]<-(-250)
mybreaks<-seq(-250,0,25)
myhistsign<-hist(dat_sign,breaks=mybreaks)
myhistnotsign<-hist(dat_notsign,breaks=mybreaks)
myhistsign$counts<-myhistsign$counts/sum(myhistsign$counts)
myhistnotsign$counts<-myhistnotsign$counts/sum(myhistnotsign$counts)
pdf(paste0(myfolderdropbox,"/hist_allelic_imbalance.pdf"))
data_temp<-c(rbind(myhistsign$count,myhistnotsign$count))
myax=data.frame(pos=seq(1,19,2),lab=mybreaks[1:10])
mybarplot_f(1:length(data_temp),data_temp,myylim=c(0,1),mycols=rep(c("gold","cadetblue"),length(data_temp)),myxlabel='log(p-value)', myylabel='density',add=FALSE,xaxis=myax,cex.axis=0.9)
legend(0,0.6,c("candidates","all"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        pch=(c(15,15)),col=c("gold","cadetblue"))
dev.off()


myhistsign<-hist(imbalance_sign$pvalue_no05_alt,breaks=mybreaks)
myhistnotsign<-hist(imbalance_notsign$pvalue_no05_alt,breaks=mybreaks)


mybreaks<-seq(0,1,0.05)

dat_notsign


myhistsign<-hist(imbalance_sign$pvalue_no05_alt,breaks=mybreaks)
myhistnotsign<-hist(imbalance_notsign$pvalue_no05_alt,breaks=mybreaks)

dev.off()
wilcox.test(imbalance_sign$pvalue_no05_alt,imbalance_notsign$pvalue_no05_alt)

for i in /mnt/454/HGDP/genomes_Bteam/HG19Align/*.bam /mnt/sequencedb/11men/martin/hg19_1000g/BAM/*.bam ; do /r1/people/pruefer/bin/BamTable2 -r ~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt;done

logonotsign[is.na(logonotsign)] <- 0
logonotsign<-cbind(logonotsign,t(sapply(1:dim(logonotsign)[1],function(x) prHvsE(logonotsign,x))))
dim(logonotsign)
res<-logo_below_thr(logonotsign,0.00001) #0.03417355
res<-logo_below_thr(logonotsign,0.05) #0.1843388





}

#hardy-weinberg
{
MYCHR=1
#from parser of 1000 genomes for my populations ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g.sh
mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_map
for mymap in 50 99;do
for mysupp in 6 7 4 5;do
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
count_n=0;count_n20=0
for MYCHR in `seq 1 22`;do
echo $MYCHR
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed.gz
source ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g_2awk.sh
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
mymappabilityfile=/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_${mymap}.bed.gz
#tabix $my1000gfile -R $mymappabilityfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
bedtools intersect -a $my1000gfile -b <( tabix $mymappabilityfile $MYCHR ) -sorted | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
count_n2=$(cat $mybedfile | wc -l )
count_n1=$( bedtools intersect -a $mybedfile -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l )
count_n=$(( $count_n + $count_n1 ));
count_n20=$(( $count_n20 + $count_n2 )); 
echo "count is " $count_n2 >> $mylogfile
done
echo "n overlaps with map in TableS" $mysupp "with map" $mymap " are" $count_n $count_n20 >> $mylogfile
done
done


for mysupp in 4 5;do
myminfreq=1
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
for MYCHR in `seq 1 22`;do
myhwtabfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/output.coding.min${myminfreq}.chr${MYCHR}.tab
vcftools --gzvcf /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr${MYCHR}.vcf.gz --hardy --out $myhwtabfile
bedtools intersect -a <( cat $myhwtabfile.hwe | grep -v CHR | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8}' ) -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | sort -k2,2n | uniq >>  ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_hwe.tableS${mysupp}
bedtools subtract -a <( cat $myhwtabfile.hwe | grep -v CHR | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8}' ) -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed  | sort -k2,2n | uniq >>  ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_hwe.notableS${mysupp}
done;done

mydata_sign<-read.table("~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_hwe.tableS4")
mydata_nosign<-read.table("~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_hwe.notableS4")
names(mydata_sign)<-c("CHR","init","end","OBS(HOM1/HET/HOM2)","E(HOM1/HET/HOM2)","ChiSq_HWE","P_HWE","P_HET_DEFICIT","P_HET_EXCESS")
names(mydata_nosign)<-c("CHR","init","end","OBS(HOM1/HET/HOM2)","E(HOM1/HET/HOM2)","ChiSq_HWE","P_HWE","P_HET_DEFICIT","P_HET_EXCESS")
#pvalues: strong excess of hets
mean(mydata_sign$P_HET_DEFICIT)
mean(mydata_nosign$P_HET_DEFICIT)
mean(mydata_sign$P_HET_EXCESS)
mean(mydata_nosign$P_HET_EXCESS)
> mean(mydata_sign$P_HET_DEFICIT)
[1] 0.8778417
> mean(mydata_nosign$P_HET_DEFICIT)
[1] 0.2537137
> mean(mydata_sign$P_HET_EXCESS)
[1] 0.2178896
> mean(mydata_nosign$P_HET_EXCESS)
[1] 0.9046385

}
#simulations
{
#1% of sites susceptible to errors
#of these, I would take a part of the individuals and I would add mutations. Just hets since I find this huge excess.
#add invented snps rather than modifying errors, at low frequencies. Maybe I should just add them with the same frequency distribution of other snps.
#actually sims should be simple.
#
cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5
zcat chr22.vcf.gz | grep -v '#' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/min5/chr22.vcf 

library(data.table)
myfreq<-0.05
nsamplesout<-100
#I sample from STU with replacement, because population with smallest bias
mynspos<-read.table("~/Dropbox/LDLD/nspops.txt")
mydata<-fread("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/min5/chr22.vcf")

myneworder<-sample((sum(mynspos$V1[1:11])+10):sum(mynspos$V1[1:12]),nsamplesout,replace=TRUE)
mydata_temp<-mydata[,myneworder,with=FALSE]
myf<-function(x){
temp<-unlist(strsplit(unlist(x),split="|"))
sum(temp==1)/(sum(temp==0)+sum(temp==1))
}
myfreqs<-apply(mydata_temp, MARGIN=1,FUN=function(x) myf(x)) 
mydata_temp<-mydata_temp[myfreqs>=myfreq & myfreqs<=(1-myfreq),]
mydatapos_temp<-mydata$V2[myfreqs>=myfreq & myfreqs<=(1-myfreq)]
n_erroneous_snps<-round(nrow(mydata_temp)/100)
id_erroneous_snps<-sample(1:nrow(mydata_temp),n_erroneous_snps)

for ( id_error in id_erroneous_snps){
mydata_temp[id_error,]
}


}


}


}
#==================================================================================
#=============================ANALYSES after removal1==============================
#==================================================================================
#load
if (TRUE)
{
setwd("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/")
data<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/sorted2.res")
names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
system("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/*minilog | grep ncomparisons | awk '{print $14}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/ncomparisons.txt")
ncomp<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/ncomparisons.txt",header= FALSE )
ncomp<-sum(ncomp$V1)
dataFIN<-subset(data,data$pop==7)
dataFIN$popX[dataFIN$popX==999]<-1
combined_fdr<-p.adjust(p=dchisq(dataFIN$X,df=2*dataFIN$popX), method = "fdr", n = ncomp) #not thebest approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
length(combined_fdr)
dim(dataFIN)
dim(data)
fullsign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
#save(fullsign,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/fullsign.RData")
#in early code called sign
rm(data,dataFIN)
dim(fullsign)
#str(fullsign)
#rm(chrA,chrB)
signAbed<-as.data.frame(cbind(as.character(paste0('chr',fullsign$chrA)),fullsign$posA-1,fullsign$posA),stringsAsFactors = FALSE )
signBbed<-as.data.frame(cbind(as.character(paste0('chr',fullsign$chrB)),fullsign$posB-1,fullsign$posB),stringsAsFactors = FALSE )
names(signAbed)<-c("chr","start","end")
names(signBbed)<-c("chr","start","end")
signAbed$start<-as.numeric(signAbed$start)
signAbed$end<-as.numeric(signAbed$end)
signBbed$start<-as.numeric(signBbed$start)
signBbed$end<-as.numeric(signBbed$end)
removal1_combined_fdr<-combined_fdr
}
#remove biased snps
if (TRUE)
{
fullsign_snps<-unique(paste0(c(fullsign[,1],fullsign[,2]),".",c(fullsign[,3],fullsign[,4])));
length(fullsign_snps)
resres_fullsign<-list()
for (ipop in 1:12)
{
resres_fullsign[[ipop]]<-subset(removal1_resresgen[[ipop]],!is.na(match(resresnames,fullsign_s
nps))) #only significant
}

#remove snps that contribute the most to nAB
removal1_removesnps<-corr_nAvsnAB(removal1_resresgen,logpop,l_mutsamples)
removal1_removesnps_sign<-corr_nAvsnAB(resres_fullsign,logpop,l_mutsamples,snpsnames=fullsign_snps)
#save(removal1_removesnps,file="removal1_removesnps.RData")
#save(removal1_removesnps_sign,file="removal1_removesnps_sign.RData")


#load("~/Dropbox/LDLD/removal1_removesnps.RData")
load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/removal1_removesnps_sign.RData")
#load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/removal1_removesnps.RData") #I have lost this file, maybe in crash of scratch! regenerate!

#almost empty! possibly because I try to clean whole data-sets, not restricting to the linked ones, so that I have a very big n for the fdr correction. #let's try restricting it.
#removal1_removesnps[removal1_removesnps[,2]<=0.05,]
removal1_removesnps_sign[removal1_removesnps_sign[,3]<=0.05,]


#fdr for real data
dim(fullsign)
#res<-rebuild_dataLDLD(removal1_removesnps[removal1_removesnps[,3]<0.05,1],fullsign)
res_sign<-rebuild_dataLDLD(removal1_removesnps_sign[removal1_removesnps_sign[,3]<0.05,1],fullsign)
#dim(res)
dim(res_sign)
#combined_fdr<-p.adjust(p=dchisq(res$X,df=2*res$popX), method = "fdr", n = ncomp) #NB: i did not correct ncomp, so they might be a bit more
combined_pval_sign<-dchisq(res_sign$X,df=2*res_sign$popX)
combined_fdr_sign<-p.adjust(p=combined_pval_sign, method = "fdr", n = ncomp)
#signleft<-res[combined_fdr<0.05,]
signleft_sign<-res_sign[combined_fdr_sign<0.05,]
#dim(signleft) #16578
dim(signleft_sign) #16578
}

#fdr and pval for null distribution
if (FALSE)
{
nperm<-3
reschr_count<-rep(0,nperm)
combined_fdr_l<-list()
combined_pval_l<-list()
for (myperm in 1:nperm)
{
data<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/reschr",myperm,"/anal/sorted.res"))
names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
system(paste0("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/reschr",myperm,"/*minilog | grep ncomparisons | awk '{print $14}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/reschr",myperm,"/ncomparisons.txt"))
ncomp<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/reschr",myperm,"/ncomparisons.txt"),header= FALSE )
ncomp<-sum(ncomp$V1)
dataFIN<-subset(data,data$pop==7)
dataFIN$popX[dataFIN$popX==999]<-1
combined_pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)
combined_fdr<-p.adjust(p=combined_pval, method = "fdr", n = ncomp) #notthe best approach because like this I have all pairs. This is why probably at the beginning I wasusing a smaller cutoff.
fullsign_reschr<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
combined_fdr_l[[myperm]]<-combined_fdr
combined_pval_l<-combined_pval
reschr_count[myperm]<-dim(fullsign_reschr)[1]
}
#save("combined_fdr_l",file="~/Dropbox/LDLD/removal1_combined_fdr_l.RData")
#save("combined_pval_l",file="~/Dropbox/LDLD/removal1_combined_pval_l.RData")

load("~/Dropbox/LDLD/removal1_combined_pval_l.RData")
reschr_count #15736 15717 15585

source("~/Dropbox/LDLD/analysesLDLD_header.R")
removal1_empfdr<-empiricall_fdr(combined_pval_sign,unlist(combined_pval_l),3)

nsignificant<-sum(as.numeric(removal1_empfdr<0.05))
removal1_sign<-res_sign[order(combined_pval_sign),][1:nsignificant,]
#save("removal1_sign",file="~/Dropbox/LDLD/removal1_sign_empT2.RData")
#save("removal1_empfdr",file="~/Dropbox/LDLD/removal1_empfdr.RData")
}

#infoplot
if (FALSE)
{
fullsign[,combined_pval_sign]
length(combined_pval_sign)
hist(sqrt(removal1_combined_fdr[removal1_combined_fdr<0.05]),breaks=15,col="coral",main="removal1_pval vs reshuffled",xlab="sqrt(fdr)")
hist(sqrt(combined_fdr[combined_fdr<0.05]),breaks=15,col="azure",add= TRUE )

setwd("/mnt/scratch/fabrizio/LDLD")
load("removal1_resresgen.RData") 
head(removal1_resresgen[[1]])
load("removal1_l_mutsamples.RData")
head(l_mutsamples[[1]])
setwd("~/Dropbox/LDLD")
load("removal1_logpop.RData")
head(logpop[[1]])
load("removal1_logpoptot.RData")
head(logpoptot[[1]])
load("removal1_test_bis.RData")
source("~/Dropbox/LDLD/analysesLDLD_header.R")
library(vioplot)
infoplot(i,mymaxk=5,reshufflings=test_removal1[[i]]$raw)
if (FALSE)
{
    par(mfrow=c(1,1))
    for (i in 1:12)
    {
        pdf(paste0("~/Dropbox/LDLD/ipynb/figs/coding1000g/removal1/logplot",i,".pdf"))
        infoplot(i,mymaxk=5,reshufflings=test_removal1[[i]]$raw)
        dev.off()
    }
}
}

#simple repeats
#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -AB -e "SELECT chrom, chromStart, chromEnd from simpleRepeat;" hg19 > simpleRepeats.bed

#check overlap with chips Rashmi
if (FALSE)
{
#she has a lot of strange IDs, including these illumina exm ones. Find conversion in http://support.illumina.com/downloads/humancoreexome-12v1-0_product_files.html
#downloaded as HumanCoreExome-12-v1-0-D-auxilliary-file.txt
#cat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/filtered/temp/missense.chr*.e /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/filtered/temp/synonymous.chr*.e | shuf -n 10000 | awk -v OFS='\t' '{if ($10==1){print $1,$2-1,$2,$11}}' >  /mnt/scratch/fabrizio/LDLD/randomsnps.bed

library(data.table)
badsnps<-fread("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed.bed")
randomsnps<-fread("/mnt/scratch/fabrizio/LDLD/randomsnps.bed")
MendelsError<-fread("~/Dropbox/LDLD/fromRashmi/all_mendel_errors_and_missingness.csv",sep=",",header=TRUE)
conversion_exm2dbsnp<-fread("~/Dropbox/LDLD/fromRashmi/HumanCoreExome-12-v1-0-D-auxilliary-file.txt")

if (FALSE)
{
MendelsError_converted<-MendelsError[!is.na(match(MendelsError$CHR,conversion_exm2dbsnp$Name)),]
myconversion<-conversion_exm2dbsnp[match(MendelsError_converted$CHR,conversion_exm2dbsnp$Name),]
myconversion<-cbind(MendelsError_converted,myconversion)
match(badsnps$V4,myconversion$RsID)
#conversion_exm2dbsnp[conversion_exm2dbsnp$RsID=="rs115347328",]
}

extract_FMISS<-function(mysnps)
{
names(MendelsError)
setnames(MendelsError,c("chr","Name","N_MISS","N_GENO","F_MISS"))
converted<-merge(MendelsError,conversion_exm2dbsnp,by="Name")
setnames(mysnps,c("chr","start","end","RsID"))
mysnps_errors<-merge(mysnps,converted,by="RsID")
setnames(mysnps,c("chr","start","end","Name"))
myfmiss_dbsnp<-merge(mysnps,MendelsError,by="Name")
myfmiss_dbsnp<-cbind(myfmiss_dbsnp[,],myfmiss_dbsnp$Name)
mynames<-names(myfmiss_dbsnp)
mynames[length(mynames)]<-"RsID"
setnames(myfmiss_dbsnp,mynames)
mysnps_errors<-rbind(myfmiss_dbsnp,mysnps_errors)
mysnps_errors[,chr.y:=NULL]
mynames<-names(mysnps_errors)
mynames[2]<-"chr"
setnames(mysnps_errors,mynames)
dim(mysnps_errors)
mysnps_errors<-unique(mysnps_errors)
mysnps_errors<-subset(mysnps_errors,F_MISS<0.05)
return(list(length(unique(sort(mysnps_errors$start))),mysnps_errors))
}

badsnps_errors_sign<-extract_FMISS(badsnps)
randomsnps_errors<-extract_FMISS(randomsnps)

badsnps_errors_sign[[1]]/dim(badsnps)[1]
randomsnps_errors[[1]]/dim(randomsnps)[1]

#is it possible because I have to check that present in FIN, and maybe all these are not.
#length(randomsnps_errors$F_MISS[as.numeric(randomsnps_errors$F_MISS)<0.05])/length(randomsnps_errors$F_MISS)
#[1] 0.009344424
#length(badsnps_errors_sign$F_MISS[as.numeric(badsnps_errors_sign$F_MISS)<0.05])/length(badsnps_errors_sign$F_MISS)
#[1] 0.1341463

hist(as.numeric(MendelsError_converted$F_MISS))

MendelsError_converted[!is.na(MendelsError_converted$CHR),]
match(MendelsError_converted$CHR,conversion_exm2dbsnp$Name)

rsid<-conversion_exm2dbsnp$RsID[!is.na(match(conversion_exm2dbsnp$Name,MendelsError$CHR))]
}

#BASE COMPOSITION
{
#base
library(data.table)
mylogo<-fread("/mnt/scratch/fabrizio/LDLD/logosign.bed")
badsnps<-fread("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed.bed")
mylogo2<-unique(merge(badsnps,mylogo,by='V3'))#NOT COMPLETE

temp<-mylogo2[mylogo2$V5=='A',]
temp<-temp[temp$V4.y=='T',]
temp2<-mylogo2[mylogo2$V5=='T',]
temp2<-temp2[temp2$V4.y=='A',]
AT_er<-dim(temp2)[1]+dim(temp)[1]
temp<-mylogo2[mylogo2$V5=='C',]
temp<-temp[temp$V4.y=='G',]
temp2<-mylogo2[mylogo2$V5=='G',]
temp2<-temp2[temp2$V4.y=='C',]
CG_er<-dim(temp2)[1]+dim(temp)[1]
temp<-mylogo2[mylogo2$V5=='A',]
temp<-temp[temp$V4.y=='G',]
temp2<-mylogo2[mylogo2$V5=='G',]
temp2<-temp2[temp2$V4.y=='A',]
AG_er<-dim(temp2)[1]+dim(temp)[1]
temp<-mylogo2[mylogo2$V5=='C',]
temp<-temp[temp$V4.y=='T',]
temp2<-mylogo2[mylogo2$V5=='T',]
temp2<-temp2[temp2$V4.y=='C',]
CT_er<-dim(temp2)[1]+dim(temp)[1]

temp<-mylogo[mylogo$V5=='A',]
temp<-temp[temp$V4=='T',]
temp2<-mylogo[mylogo$V5=='T',]
temp2<-temp2[temp2$V4=='A',]
AT_back<-dim(temp2)[1]+dim(temp)[1]
temp<-mylogo[mylogo$V5=='C',]
temp<-temp[temp$V4=='G',]
temp2<-mylogo[mylogo$V5=='G',]
temp2<-temp2[temp2$V4=='C',]
CG_back<-dim(temp2)[1]+dim(temp)[1]
temp<-mylogo[mylogo$V5=='A',]
temp<-temp[temp$V4=='G',]
temp2<-mylogo[mylogo$V5=='G',]
temp2<-temp2[temp2$V4=='A',]
AG_back<-dim(temp2)[1]+dim(temp)[1]
temp<-mylogo[mylogo$V5=='C',]
temp<-temp[temp$V4=='T',]
temp2<-mylogo[mylogo$V5=='T',]
temp2<-temp2[temp2$V4=='C',]
CT_back<-dim(temp2)[1]+dim(temp)[1]

mybackground<-c(AG_er,CT_er,AT_er,CG_er)/sum(c(AG_er,CT_er,AT_er,CG_er))
myerr<-c(AG_back,CT_back,AT_back,CG_back)/sum(c(AG_back,CT_back,AT_back,CG_back))
myerr<-rbind(mybackground,myerr)
colnames(myerr)<-c("AG","CT","AT","CG")
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/base_comp.pdf")
barplot(myerr,beside=TRUE,col=c("floralwhite","azure3"),ylab="count",)
legend("topleft",fill=c("floralwhite","azure3"), legend=c("biased","background"))
dev.off()
}
#------------------------------------explorative plots all pops together----------------------------v
{
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
}
#---------------------------------------------------------------------------------------------------^
#TAG 'excess of linkage'
{
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
}
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
{
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

mysnps<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/removed1.bed",stringsAsFactors =FALSE,header=FALSE)
names(mysnps)<-c("chr","start","end")
mysnps$start<-as.numeric(mysnps$start)
mysnps$end<-as.numeric(mysnps$end) 
mysnps[,1]<-sapply(mysnps[,1],function(x) paste0('chr',x))
createbedlogo(mysnps,"/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/removed1.logo")


#cat above95/coding1000g/chr*.tab | grep -v '#' | sort -Vu -k1,1 -k2,2 | awk -v OFS='\t' '{print "chr"$1,$2-1,$2}' > above95/coding1000g/chrtab.bed #background from files .tab #NB: a bit more than bed2 in $TARGET. check why!
#bedtools intersect -b above95/coding1000g/anal/snps.bed -a above95/coding1000g/chrtab.bed -v > above95/coding1000g/notsnps.bed
data<-read.table("above95/coding1000g/notsnps.bed",stringsAsFactors =FALSE,header=FALSE)
names(data)<-c("chr","start","end")
data$start<-as.numeric(data$start)
data$end<-as.numeric(data$end)
createbedlogo(data,"above95/coding1000g/logonotsign")

createbedlogo(data,"above95/coding1000g/logonotsign")

data1<-cbind(data,0,0,0,0)
write.table(myres2,file="above95/coding1000g/logonotsign.temp",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
}
#=========================CHECK IMBALANCE in .BAM FILES by using LOGO files==================
{
for MYCHR in `seq 1 22`;do
bedtools intersect -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/removed1.bed -b /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/chr${MYCHR}.vcf.gz -sorted -wo | awk -v OFS='\t' '{print $1,$2,$3,$7,$8,$9}' >> /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/removed1.bed.logo

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

setwd("/mnt/scratch/fabrizio/LDLD/")
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
}
#======================================================REMOVE CENTRAL CLUSTER
{
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
}
#======================================================trios from this
#TAG 'trios from this'
{
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
}
#======================================================================NETWORK-PLOTS
#TAG NETWORK-PLOTS
{
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
}
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
#----------------------------------------^
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
{
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
}
#epistasis imported from mathematica
{
library(vioplot)
virnorm(100,0.00001,0.000001)

}
