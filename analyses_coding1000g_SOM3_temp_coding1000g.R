if (TRUE) {
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new3"
myfolderdropbox<-"~/Dropbox/LDLD/ipynb/figs/coding1000g/march2017"
#system(paste0("cat ",myfolder,"/*minilog | grep ncomparisons | awk '{print $14}' > ",myfolder,"/ncomparisons.txt"))
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
computed_exact_pvalues<-FALSE
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
npopulations<-12
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE

setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
options(scipen=999)
setwd("/mnt/scratch/fabrizio/LDLD")
library(data.table)

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
if (TRUE)
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
#gcc ~/Dropbox/LDLD/scripts/LDLD_filter5.c -lgmp -lm -o filterbyfreq.out
#for j in "/reschr1/" "/reschr2/"; do #"/"
#for i in `seq 1 22`; do 
#./filterbyfreq.out ${myfolder}${j}chr$i.tab ${myfolder}${j}chr$i.freqlog ${myfolder}${j}chr$i.tab.mutlog #0.05 1220 & 
#done;done

#import nA
for (i in 1:22)
{
if (i==1){mutperssample<-as.numeric(read.table(paste0(myfolder,"/chr",i,".tab.mutlog")))}
else {res2<-as.numeric(read.table(paste0(myfolder,"/chr",i,".tab.mutlog")));mutperssample<-rbind(mutperssample,res2)}
}
l_mutsamples<-lapply(1:12,function(x) c())
for (i in 1:12)
{
  l_mutsamples[[i]]<-apply(mutperssample[,samples_imypopi(nspops,i)],2,sum)
}
nperm<-2
#import nA for reshuffled
l_mutsamples_reschr<-list()
for (myperm in 1:nperm){for (i in 1:22)
{
if (i==1){mutperssample<-as.numeric(read.table(paste0(myfolder,"/reschr",myperm,"/chr",i,".tab.mutlog")))}
else {res2<-as.numeric(read.table(paste0(myfolder,"/reschr",myperm,"/chr",i,".tab.mutlog")));mutperssample<-rbind(mutperssample,res2)}
}
l_mutsamples_reschr[[myperm]]<-lapply(1:12,function(x) c())
for (i in 1:12)
{
l_mutsamples_reschr[[myperm]][[i]]<-apply(mutperssample[,samples_imypopi(nspops,i)],2,sum)
}}

} 
#test if nAB inhomogeneous
{
load(paste0(myfolder,"/logpop_nAB.RData"))
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
load(paste0(myfolder,"/logpop_nAB_neg.RData"))
load(paste0(myfolder,"/logpop_nAB_nlinks.RData"))
load(paste0(myfolder,"/logpop_reschr",myperm,"_nAB.RData"))
load(paste0(myfolder,"/mylogpop_reschr",myperm,".RData"))
logpop_reschr_nAB_l<-list()
for (z in 1:nperm) { load(paste0(myfolder,"/logpop_reschr",z,"_nAB.RData")); logpop_reschr_nAB_l[[z]]<-logpop_reschr_nAB; }
logpop_reschr_nAB_neg_l<-list()
logpop_reschr_nAB_pos_l<-list()
for (z in 1:nperm) { load(paste0(myfolder,"/logpop_nAB_neg_reschr",z,"_nAB.RData")); logpop_reschr_nAB_neg_l[[z]]<-logpop_nAB_neg_reschr; }
for (z in 1:nperm) { load(paste0(myfolder,"/logpop_nAB_pos_reschr",z,"_nAB.RData")); logpop_reschr_nAB_pos_l[[z]]<-logpop_nAB_pos_reschr; }

#test_inhomogeneous_nAB_pvalues<-test_inhomogeneous_nAB(logpop_nAB,logpop_reschr_nAB_l)
#save(test_inhomogeneous_nAB_pvalues,file=paste0(myfolder,"/test_inhomogeneous_nAB_pvalues.RData"))
#test_inhomogeneous_nAB_pvalues_pos<-test_inhomogeneous_nAB(logpop_nAB_pos,logpop_reschr_nAB_pos_l)
#save(test_inhomogeneous_nAB_pvalues_pos,file=paste0(myfolder,"/test_inhomogeneous_nAB_pvalues_pos.RData"))
#test_inhomogeneous_nAB_pvalues_neg<-test_inhomogeneous_nAB(logpop_nAB_neg,logpop_reschr_nAB_neg_l)
#save(test_inhomogeneous_nAB_pvalues_neg,file=paste0(myfolder,"/test_inhomogeneous_nAB_pvalues_neg.RData"))
}

#================================================================================================^
#===============PER SAMPLE PER SNP ANALYSES: IDENTIFYING BIASED SNPs -> removal1 ================
#================================================================================================

