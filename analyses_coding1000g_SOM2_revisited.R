popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
options(scipen=999)
setwd("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/old")
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

#--------------------- creates polymorphism per individual (mutlog) file (newer)   ------------------v
{
#for i in `seq 1 22`; do 
#./filterbyfreq.out above95/coding1000g/chr$i.tab above95/coding1000g/chr$i.freqlog above95/coding1000g/chr$i.tab.mutlog 0.05 1220 & 
#done
#system("cat above95/coding1000g/chr*.tab.mutlog > above95/coding1000g/tot.mutlog")
for (i in 1:22)
{
if (i==1){mutperssample<-as.numeric(read.table(paste0("chr",i,".tab.mutlog")))}
else {res2<-as.numeric(read.table(paste0("chr",i,".tab.mutlog")));mutperssample<-rbind(mutperssample,res2)}
}
l_mutsamples<-lapply(1:12,function(x) c())
for (i in 1:12)
{
  l_mutsamples[[i]]<-apply(mutperssample[,samples_imypopi(nspops,i)],2,sum)
}
}


resres<-read.table("all.freqlog") #resres is a matrix with dim(resres) #[1] 75785  1219
#resres has 0,1,2

resresnames<-paste0(resres[,1],".",resres[,2])
load("fullsign.RData")

corr_nAvsnAB<-function(myres,logpop,l_mutsamples,snpsnames=snpsnames)
{
#it calculates the logistic correlation between the presence of alternative allele and nAB per sample
#it takes as arguments lists (each element is a population) with myres (genotypes in 0,1(het) and 2(homozygotes)), logpop(nAB per sample), l_mutsamples (number of alt allele per individual)
  pvalsnps<-c();   pvalsnps2<-c()
  for (mylines in 1:dim(myres[[1]])[1])
  {
    print(mylines)
    l_pval<-c(); l_pval2<-c()
    for (ipop in 1:length(logpop))
    {
	mygen<-unlist(myres[[ipop]][mylines,1:(dim(myres[[ipop]])[2])])
	#separates diploid genotypes into 2 haploid ones
	mygen2<-mygen
	mygen3<-mygen
	mygen2[mygen2==2]<-1
	mygen3[mygen3==1]<-0
	mygen3[mygen3==2]<-1
	#to include only heterozygotes
	mygen[mygen==2]<-0 
	if (sum(mygen)/length(mygen)>0.05 && sum(mygen)/length(mygen)<0.95) #isn't this only looking at freq of homoz??? if that was the reason I could try to replicate it 1-10-2017
	  {
	  #print(c(length(mygen),length(logpop[[ipop]][[1]]),length(l_mutsamples[[ipop]])))
	  l_pval<-append(l_pval,logreg_snps2nAB(mygen,logpop[[ipop]][[1]],l_mutsamples[[ipop]]))
	  l_pval2<-append(l_pval2,logreg_snps2nAB(c(mygen2,mygen3),c(logpop[[ipop]][[1]],logpop[[ipop]][[1]]),c(l_mutsamples[[ipop]],l_mutsamples[[ipop]])))
	  }
    }
    if (length(l_pval)>0)
      {
      pvalsnps<-c(pvalsnps,pchisq( -2*sum(log(l_pval)), df=2*length(l_pval), lower.tail=FALSE))
      pvalsnps2<-c(pvalsnps2,pchisq( -2*sum(log(l_pval2)), df=2*length(l_pval2), lower.tail=FALSE))
      } else {pvalsnps<-c(pvalsnps,2);pvalsnps2<-c(pvalsnps2,2)}
  }
res<-cbind(snpsnames,p.adjust(p=pvalsnps, method = "fdr", n = length(pvalsnps)),p.adjust(p=pvalsnps2, method = "fdr", n = length(pvalsnps2)))
return(res)
}

tempres<-corr_nAvsnAB(resres,logpop,l_mutsamples) #this is done by looking at all SNPs.



