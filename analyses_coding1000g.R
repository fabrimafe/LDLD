setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
options(scipen=999)
setwd("/mnt/scratch/fabrizio/LDLD")
namespops=c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
#index of objects
#-------------what_is_logpop?
#lists with pops and below nABf as first subelement and nAndBf as second.
#-------------what_is_gen?
#genotypes in 0,1,2 (for heterozygotes)



#======================LOAD_INITIAL_FILES====================================================
#TAG 'LOAD_INITIAL_FILES'
#============================================================================================
#===========================CREATE LOGOS================
#TAG 'CREATE LOGOS'
#=========================CHECK IMBALANCE in .BAM FILES by using LOGO files==================
#TAG 'CHECK IMBALANCE' in analyses_coding1000g_SOM.R
#======================================================REMOVE CENTRAL CLUSTER
#TAG 'REMOVE CENTRAL CLUSTER'
#======================================================trios from this
#TAG 'trios from this'
#generate BED files of snps removed.
#======================================================================NETWORK-PLOTS
#TAG NETWORK-PLOTS
#==================PER-SAMPLE-ANALYSES=====================================v
#---import files--------------------------------------------------------------------------------------v
{
info1000g<-read.table("/mnt/scratch/fabrizio/LDLD/20130606_sample_info.txt",header=TRUE,sep='\t')
samples1000g<-system("zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | head -300 | grep '#CHROM' | head -1 | awk '{for (i=10; i<=NF; i++) print $i}'",intern=TRUE)
samplesabove95<-system("cat ~/workspace/1000genomes/above95.unrelated.samples",intern=TRUE)
nspops<-as.numeric(system("cat /mnt/scratch/fabrizio/LDLD/nspops.txt",intern=TRUE))
samples_order1000g_above95<-intersect(samplesabove95,samples1000g)
infosamples<-lapply(1:12,function(x) c())
for (i in 1:12)
{
  mysamples<-samples_imypop(samples_order1000g_above95,nspops,i)
  infosamples[[i]]<-info1000g[match(mysamples,info1000g$Sample),]
}
popsamplesi<-sapply(1:12,function(x) samples_imypopi(nspops,x)) #indices of samples -> e.g. mylog2[popsamplesi[[i]]]

#for each pair of chromosome compute the nAB contribution when considering significant or all
system("cat above95/coding1000g/chr*.tabchr*.minilog | grep nAB_above_freq_threshold | grep -v and_significant > above95/coding1000g/tot.nABabovethr")
system("cat above95/coding1000g/chr*.tabchr*.minilog | grep nAB_above_freq_threshold | grep and_significant > above95/coding1000g/tot.nABabovethrsign")
mylog2<-read.table("above95/coding1000g/tot.nABabovethr") 
mylog2<-sapply(2:dim(mylog2)[2],function(x) sum(mylog2[,x]))
mylog3<-read.table("above95/coding1000g/tot.nABabovethrsign")
mylog3<-sapply(2:dim(mylog3)[2],function(x) sum(mylog3[,x]))
#

#-------------------- create whole genome log file -----------------------------
#for (i in 0:11)
#{
# system(paste0("cat above95/coding1000g/chr*.log | awk '{if ($NF==",i,"){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.pop",i,".log"))
#}
#-------------------- create whole genome log file by chr ----------------------
#system(paste0("mkdir above95/coding1000g/logs"))
#for (ipop in 0:11)
#  {
#  for (ichrA in 1:21)
#    {
#   for (ichrB in (ichrA+1):22)
#	{
#	system(paste0("cat above95/coding1000g/logs/chr*.log > above95/coding1000g/logs/chr",ichrA,".",ichrB,".pop",ipop,".log"))
#	}
#    }
#  }
#for (i in 0:11)
#{
#  system(paste0("cat above95/coding1000g/logs/chr*.log | awk '{if ($NF==",i,"){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/logs/all.pop",i,".log"))
#}
#-------------------- create whole genome freqlog file (genotype per sample) ---
#rm /mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.freqlog
#for i in `seq 1 22`; do 
#cat above95/coding1000g/chr$i.freqlog | awk -v FS='\t' -v OFS='\t' -v CHR=$i '{if (NR>1){print CHR,$0}}' >> /mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.freqlog
#done
#-------------------- creates polymorphism per individual (mutlog) file (old) --
#cat ~/workspace/1000genomes/20130606_sample_info.txt | sed s/$'\t\t'/$'\t'NA$'\t'/g > /mnt/scratch/fabrizio/LDLD/20130606_sample_info.txt
#for i in `seq 1 22`; do 
#./filter.out above95/coding1000g/chr$i.tab above95/coding1000g/chr$i.tab.mutlog 0.05 1220 & 
#done
#--------------------- creates polymorphism per individual (mutlog) file (newer)   ------------------v
#for i in `seq 1 22`; do 
#./filterbyfreq.out above95/coding1000g/chr$i.tab above95/coding1000g/chr$i.freqlog above95/coding1000g/chr$i.tab.mutlog 0.05 1220 & 
#done
#system("cat above95/coding1000g/chr*.tab.mutlog > above95/coding1000g/tot.mutlog")
for (i in 1:22)
{
if (i==1){mutperssample<-as.numeric(read.table(paste0("above95/coding1000g/chr",i,".tab.mutlog")))}
else {res2<-as.numeric(read.table(paste0("above95/coding1000g/chr",i,".tab.mutlog")));mutperssample<-rbind(mutperssample,res2)}
}
l_mutsamples<-lapply(1:12,function(x) c())
for (i in 1:12)
{
  l_mutsamples[[i]]<-apply(mutperssample[,samples_imypopi(nspops,i)],2,sum)
}
#----------------------------------------------------------------------------------------------------^
#-------logpop from log files---------(all that are printed in res files)----------------------------v
if (FALSE)
{
  logpop<-lapply(1:12,function(x) list(c(),c()))
  for (i in 0:11)
  {
    mydata<-read.table_handle.wrong.lines(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.pop",i,".log"),header=FALSE)
    logpop[[i+1]][[1]]<-sapply(3:(dim(mydata)[2]-1), function(x) sum(sapply(mydata[,x],nABf)))
    logpop[[i+1]][[2]]<-sapply(3:(dim(mydata)[2]-1), function(x) sum(sapply(mydata[,x],nAandBf)))
  }
#  save(logpop,"/mnt/scratch/fabrizio/LDLD/above95/coding1000g/logpop.RData")
  nAB_lFIN<-sapply(3:(dim(logFIN)[2]-1), function(x) sum(sapply(logFIN[,x],nABf)))#-1 because last is the population flag
  nAB_lYRI<-sapply(3:(dim(logYRI)[2]-1), function(x) sum(sapply(logYRI[,x],nABf)))
  nAB_lpop0<-sapply(3:(dim(logpop0)[2]-1), function(x) sum(sapply(logpop0[,x],nABf)))
  nAB_lpop1<-sapply(3:(dim(logpop1)[2]-1), function(x) sum(sapply(logpop1[,x],nABf)))
  nAandB_lFIN<-sapply(3:(dim(logFIN)[2]-1), function(x) sum(sapply(logFIN[,x],nAandBf)))
  nAandB_lYRI<-sapply(3:(dim(logYRI)[2]-1), function(x) sum(sapply(logYRI[,x],nAandBf)))
  nAandB_lpop0<-sapply(3:(dim(logpop0)[2]-1), function(x) sum(sapply(logpop0[,x],nAandBf)))
  nAandB_lpop1<-sapply(3:(dim(logpop1)[2]-1), function(x) sum(sapply(logpop1[,x],nAandBf)))
  length(nAandB_lpop0)#106
  length(nAandB_lpop1)#107
  length(nAandB_lFIN)#99
  FINsamples<-samples_imypop(samples_order1000g_above95,nspops,8)
  infosamplesFIN<-info1000g[match(FINsamples,info1000g$Sample),]
  YRIsamples<-samples_imypop(samples_order1000g_above95,nspops,10)
  infosamplesYRI<-info1000g[match(YRIsamples,info1000g$Sample),]
  pop0samples<-samples_imypop(samples_order1000g_above95,nspops,1)
  infosamplespop0<-info1000g[match(pop0samples,info1000g$Sample),]
  pop1samples<-samples_imypop(samples_order1000g_above95,nspops,2)
  infosamplespop1<-info1000g[match(pop1samples,info1000g$Sample),]
}
load("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/logpop.RData")
#logpop[[1]][[1]] 
#i<-1
#mydata<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.pop",i,".log"),header=FALSE)
#head(mydata)

#they do correlate!!!
if (FALSE)
{
length(l_mutsamples[[1]])
cor.test(logpop[[1]][[1]],l_mutsamples[[1]]) 
cor.test(logpop[[2]][[1]],l_mutsamples[[2]])
cor.test(logpop[[3]][[1]],l_mutsamples[[3]])
cor.test(logpop[[5]][[1]],l_mutsamples[[5]])
cor.test(logpop[[8]][[1]],l_mutsamples[[8]])
cor.test(logpop[[4]][[1]],l_mutsamples[[4]])
cor.test(logpop[[12]][[1]],l_mutsamples[[12]])
plot(logpop[[5]][[1]],l_mutsamples[[5]])


sd(l_mutsamples[[3]])/mean(l_mutsamples[[3]])
sd(l_mutsamples[[5]])/mean(l_mutsamples[[5]])
sd(l_mutsamples[[1]])/mean(l_mutsamples[[1]])
sd(l_mutsamples[[8]])/mean(l_mutsamples[[8]])
sd(l_mutsamples[[4]])/mean(l_mutsamples[[4]])
sd(l_mutsamples[[12]])/mean(l_mutsamples[[12]])
}
#------------------------------------log pop from files
#elements in contingency tables
mylogs<-list();myblocks<-list();permbychr<-list();myblocks2<-list();permbychrblocks2<-list();
for (i in 0:11)
{
print(paste("population",i))
mylogs[[i+1]]<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/logs/all.pop",i,".log")) #any pop is the same -> to fix
print("create blocks")
#myblocks[[i+1]]<-create_blocks(mylogs[[i+1]])
#myblocks2[[i+1]]<-create_blocks(mylogs[[i+1]],space_between_indep_blocks=1000000)
#print("permute by chr")
#permbychr[[i+1]]<-permute_chrAB_byblock(mylogs[[i+1]],myblocks[[i+1]],5000)
#print("permute blocks")
#permbychrblocks2[[i+1]]<-permute_chrAB_byblock(mylogs[[i+1]],myblocks2[[i+1]],10000)
#save(file="~/Dropbox/LDLD/ipynb/tempfiles/coding1000g/permbychrblocks2.RData",permbychrblocks2)
#save(file="~/Dropbox/LDLD/ipynb/tempfiles/coding1000g/permbychr.RData",permbychr)
}
load("~/Dropbox/LDLD/ipynb/tempfiles/coding1000g/permbychrblocks2.RData")
load("~/Dropbox/LDLD/ipynb/tempfiles/coding1000g/permbychr.RData")

#fixed blocks. If thr 50Mb as in Skelly 2015 I won't get any division, too many snps. Two extremes, only recomb across chromosomes. All snps.


#Other solution is using pairs of chromosome as independent blocks. This is conservative in finding differences between distributions. However, it should still have the power.
#Also, should I consider the contribution like I do now, just sum. Well, it depends on the way of permuting it. If this way it is fine.
#see it like that however: when many pairs most likely I have enough power to consider even chromosomes as single blocks.
#when lower amount of pairs, separating becomes meaningful.
if (FALSE)
{
  mylog[,1]==21
  templog<-mylog[mylog[,1]==21,]
  templog<-templog[templog[,2]==22,]
  head(templog)

  dim(mylog)
  templog<-mylog[mylog[,1]==1,]
  templog<-templog[templog[,2]==2,]
  dim(templog)

  #on single pair of chromosome these two methods return same sum as it should be
  res<-permute_block(templog,10000)
  for (i in 1:10000)
  {res[,i]<-res[,i][order(res[,i])]}
  res[106,]#391

  resorig<-sapply(5:(dim(templog)[2]-1),function (y) sum(sapply(templog[,y], function(x) nABf(x))))
  resorig[order(resorig)][106] #391

  resorig<-sapply(5:(dim(mylog)[2]-1),function (y) sum(sapply(mylog[,y], function(x) nABf(x))))
  resorig[order(resorig)][106] #22545.5
  resorig

  myres
  myblocks<-create_blocks(mylog)
  permbychr<-permute_chrAB_byblock(mylog,myblocks,100000)
  myblocks2<-create_blocks(mylog,space_between_indep_blocks=1000000)
  permbychrblocks2<-permute_chrAB_byblock(mylog,myblocks2,100000)
  #save(file="~/Dropbox/LDLD/ipynb/tempfiles/coding1000g/permbychrblocks2.RData",permbychrblocks2)
  #save(file="~/Dropbox/LDLD/ipynb/tempfiles/coding1000g/permbychr.RData",permbychr)
}
#save.image("permbychrblocks2.RData")
setwd("/mnt/scratch/fabrizio/LDLD")
#load("permbychrblocks2.RData")
#load("permbychr.RData")
#load("~/Dropbox/LDLD/ipynb/tempfiles/coding1000g/permbychrblocks2.RData")
#load("~/Dropbox/LDLD/ipynb/tempfiles/coding1000g/permbychr.RData")
if (FALSE)
{
  sum(permbychrblocks2[,1]) #428177
  sum(permbychr[,1]) #428177
  var(permbychr[,1]) #137336.6
  var(permbychrblocks2[,1]) #9129.772
  #as expected less variance in fragmenting more

  #very cool, now which is the best way to test if I can get this fragmentation by chance?

  mydata<-permbychr[,1]
  sum(permbychr[,10000]) #428177
  system.time(nulldistr1iteration<-generate_nulldistr_mixtures(permbychr,1000)) #time 448.160
  as.numeric(nulldistr1iteration[,1])
  as.numeric(nulldistr1iteration[,2])
  as.numeric(nulldistr1iteration[,3])
  as.numeric(nulldistr1iteration[,4])
  mypop[4]
  #compare value with distribution
  mypop<-partition_samples_nAB(logpop[[1]][[1]])[c(1,4,7,8)]

  sum(as.numeric(as.numeric(nulldistr1iteration[,1])>=mypop[1]))/1000 # 0.053 it seems that by taking only number of groups, both with LR and AIC I don't even get significant, although almost
  sum(as.numeric(as.numeric(nulldistr1iteration[,2])>=mypop[2]))/1000 # 0.115
  sum(as.numeric(as.numeric(nulldistr1iteration[,3])<=mypop[3]))/1000 # 0 instead by taking p-values or RL I get extremely significant.
  sum(as.numeric(as.numeric(nulldistr1iteration[,4])<=mypop[4]))/1000 #0
  #ok, so let's say that in case of this population we want to go forward in removing individuals

  system.time(permbychr<-permute_chrAB(mylog,10000)) #   user   system  elapsed 3885.680    0.410 3886.155
  #save.image("permbychr.RData")
  head(permbychr)
  dev.new()
  plot.new()
  for (i in 1:10000)
  {permbychr[,i]<-permbychr[,i][order(permbychr[,i])]}
  for (i in 1:10000)
  {permbychrblocks2[,i]<-permbychrblocks2[,i][order(permbychrblocks2[,i])]}
  library(vioplot)
  pdf("empiricaldistrperm_pop1.pdf")
  par(mfrow=c(2,2))
  vioplot(permbychr[1,],permbychr[10,],permbychr[20,],permbychr[30,],permbychr[40,],permbychr[50,],permbychr[60,],permbychr[70,],permbychr[80,],permbychr[90,],permbychr[100,],permbychr[106,],ylim=c(0,25000),col="firebrick",colMed="firebrick2")
  vioplot(permbychrblocks2[1,],permbychrblocks2[10,],permbychrblocks2[20,],permbychrblocks2[30,],permbychrblocks2[40,],permbychrblocks2[50,],permbychrblocks2[60,],permbychrblocks2[70,],permbychrblocks2[80,],permbychrblocks2[90,],permbychrblocks2[100,],permbychrblocks2[106,],add=TRUE,col="cadetblue3",colMed="cadetblue4")
  mypop<-logpop[[1]][[1]][order(logpop[[1]][[1]])]
  points(c(mypop[1],mypop[10],mypop[20],mypop[30],mypop[40],mypop[50],mypop[60],mypop[70],mypop[80],mypop[90],mypop[100],mypop[106]))
  hist(as.numeric(nulldistr1iteration[,1]),breaks=1:6,main="nblocks LR 1000 permutations",xlab="nblocks",col="gray")
  hist(as.numeric(nulldistr1iteration[,3]),breaks=50,main="LR 1000 permutations",xlab="pvalue LR",col="gray")
  hist(as.numeric(nulldistr1iteration[,4]),breaks=50,main="RL AICc 1000 permutations",xlab="RL",col="gray")
  dev.off()
  par(mfrow=c(1,1))
  pdf("empiricaldistrperm_pop1_var.pdf")
  permvar<-sapply(1:10000,function(x) var(permbychr[,x]))
  hist(log(permvar),xlim=c(log(min(permvar)),log(max(max(permvar),var(logpop[[1]][[1]]))+10)))
  abline(v=log(var(logpop[[1]][[1]])),col="red")
  dev.off()
}

}
#---------------------------------------------------------------------------------------------------^
#------------------------------------explorative plots all pops together----------------------------v
#see TAG 'explorative plots all pops together' in analyses_coding1000g_SOM.R
#---------------------------------------------------------------------------------------------------^
{
length(samples_order1000g_above95) #1217
length(mutperssample[1,]) #1218 #because I have the last 0, so 1 more
#save.image("LDLD170216.RData")
setwd("/mnt/scratch/fabrizio/LDLD")
load("LDLD170216.RData")
setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
}
#----TAG:info statistics in 1000genomes--v
#----------------------------------------^
#---histograms to illustrate separation method-----v see TAG 'histograms to illustrate separation method in analysesLDLD_header' in analyses_coding1000g_SOM.R
#------------info panels--------------------v
{
par(mfrow=c(1,1))
for (i in 1:12)
{
pdf(paste0("~/Dropbox/LDLD/ipynb/figs/coding1000g/logplot",i,".pdf"))
infoplot(i,mymaxk=5,mypopname=namespops[i])
dev.off()
}

pdf("logplot_panel.pdf")
par(mfrow=c(2,2))
for (i in 1:4)
{
infoplot(i)
}
dev.off()

#add also substructure
#first check that samples are the same in header and infosample
system("head -1 /mnt/scratch/fabrizio/LDLD/above95/coding1000g/chr22.tab > tempheader")
mytemphead<-read.csv("tempheader",head=FALSE,sep='\t')
as.character(unlist(mytemphead[,-(1:9)]))
as.vector(myinfosamples[,1])
as.character(unlist(mytemphead[,-(1:9)]))==unlist(sapply(1:12,function(x) as.data.frame(infosamples[x])[,1])) #check that all true

tbl=read.table("/mnt/scratch/fabrizio/LDLD/above95/intergenic/pop0/prunedData.3.Q")[ord,]

tbl=read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/pop0/prunedData.3.Q")[ord,]
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)

layout(matrix(c(1,1,2,2,3,3,4,4), nrow = 4, ncol = 1))
setwd("~/Dropbox/LDLD/ipynb/figs")
for (i in 1:6)
{
pdf(paste0("logplot",i,".adm.pdf"))
layout(mat=matrix(c(1,2,3),nrow = 3, ncol = 1),heights=c(2,1,1))
infoplot(i)
ord<-infoplot_ord(i)
tbl=read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/pop",i-1,"/prunedData.3.Q"))[ord,]
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
tbl=read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/intergenic/pop",i-1,"/prunedData.3.Q"))[ord,]
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()
}
plot(1)
plot(2)
plot(3)

#how to randomize properly to get confidence intervals of uniform expectation, taking into account local linkage?
#The most simple idea is uniform (or gaussian) expectation, then rank and plot.
#however, by doing this I neglect local linkage -> I need some sort of empirical distribution.
#The simple idea would be to randomize genotypes and then recalculate nAB, however this would change for sure the total number of nAB: it would still be good since I take proportions, but not really.
#Shuffle nAB blocks further than 50 Mb. I have two positions, how to do it? but A might be associated with two differen Bs.
#What do I want to randomize exactly? Hypothesis that 1 individuals cannot have all the possible linkage.

}
#--------------------------------------------------^
#-----------analyses per SAMPLE per SNP------------v
{
if (FALSE)  #from log files
{
  #one could use function: logpop2gen

  resres<-list()
  for (i in 10:10)
  {
    counter<-1
    for (ichrA in 1:21)
      {
      for (ichrB in (ichrA+1):22)
	{
	print(c(i,ichrA,ichrB))
	myn=as.numeric(system(paste0("head -12 above95/coding1000g/chr",ichrA,".tabchr",ichrB,".tab.log | awk '{if ($NF==",i,"){print NF}}'"),intern=TRUE))
	if (length(myn)<1) {myn=as.numeric(system(paste0("head -24 above95/coding1000g/chr",ichrA,".tabchr",ichrB,".tab.log | awk '{if (NR>12 && $NF==",i,"){print NF}}'"),intern=TRUE))}
	system(paste0("cat above95/coding1000g/chr",ichrA,".tabchr",ichrB,".tab.log | awk '{if ($NF==",i," && NF==",myn,"){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.pop.temp.log"))
	logsnp<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.pop.temp.log",header=FALSE)
	logsnp$chrA<-ichrA
	logsnp$chrB<-ichrB
	if (counter==1) {res<-logsnp} else {res<-rbind(res,logsnp)}
	if (ichrA==21 && ichrB==22) {resres[[i+1]]<-res}
	counter<-counter+1
	}
      }
  }
}

load("LDLD170216.RData")
source("analysesLDLD_header.R")
resres<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.freqlog")
length(resres)
rress<-list()
for (i in 1:12)
{
rress[[i]]<-resres[,samples_imypopi(nspops,i)+2] #2*i comes from chr and pos fields
}
snpsnames<-sapply(1:dim(resres)[1],function(x) paste0(resres[x,1:2],collapse="."))
resres<-rress
rm(rress)
head(resres[[1]])

tempres<-corr_nAvsnAB(resres,logpop,l_mutsamples)
save.image("resres_complete.RData")
load("resres_complete.RData")
#removed_snps<-tempres[tempres[,2]<=0.05,] #only heterozygous
removed_snps<-tempres[tempres[,3]<=0.05,] #also homozygotes dim(removed_snps) #612, 3
stringsnp<-removed_snps

removed_snps<-stringsnp2bed(removed_snps)
write.table(removed_snps,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed1.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
for (mychr in 1:22)
{
system(paste0("bedtools intersect -v -b /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed1.bed -a /mnt/scratch/fabrizio/LDLD/above95/coding1000g/chr",mychr,".tab -header | gzip > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/chr",mychr,".coding.vcf.gz"))
}

if (FALSE)
{
rm /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed.bed
touch /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed.bed
for i in `seq 1 22`; do 
echo $i
zcat /mnt/sequencedb/ucsc/goldenPath/hg19/snp142.bed.gz | grep chr$i$'\t' > /mnt/scratch/fabrizio/LDLD/temp.bed
bedtools intersect -a <(cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed1.bed | sed 's/^/chr/g') -b /mnt/scratch/fabrizio/LDLD/temp.bed -sorted -wb| awk -v OFS='\t' '{print $4,$5,$6,$7}' >>  /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed.bed
done
#bedtools intersect -a <(cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed1.bed | sed 's/^/chr/g') -b /mnt/sequencedb/ucsc/goldenPath/hg19/snp142.bed.gz -sorted | awk -v OFS='\t' '{print $4,$5,$6,$7}' >  /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed.bed
}


dim(fullsign) #from first filter #127771     24
ncomp<-read.table("above95/coding1000g/ncomparisons.txt",header=FALSE)
ncomp<-sum(ncomp$V1)
#res<-rebuild_dataLDLD(tempres[tempres[,2]<0.05,1],sign)
res<-rebuild_dataLDLD(tempres[tempres[,3]<0.05,1],fullsign)
dim(res) #57398    24
combined_fdr<-p.adjust(p=dchisq(res$X,df=2*res$popX), method = "fdr", n = ncomp) #NB: i did not correct ncomp, so they might be a bit more
length(combined_fdr) #57398
dim(res)
signleft<-res[combined_fdr<0.05,]
dim(signleft) #16578    24
myrem<-paste0(signleft[,1],".",signleft[,3],".",signleft[,2],".",signleft[,4])
length(myrem)
mysignm<-paste0(fullsign[,1],".",fullsign[,3],".",fullsign[,2],".",fullsign[,4])
mylogm<-paste0(mylog[,1],".",mylog[,3],".",mylog[,2],".",mylog[,4])
dim(rebuild_dataLDLD(myrem,fullsign,fromremovedsnps=FALSE)) #16578 as expected
dim(rebuild_dataLDLD(mysignm,fullsign,fromremovedsnps=FALSE)) #127771 as expected
dim(rebuild_dataLDLD(mylogm,mylog,fromremovedsnps=FALSE)) #127784 as expected
myres<-rebuild_dataLDLD(myrem,mylog,fromremovedsnps=FALSE)
dim(myres) #16451 some are missing #very few, annoying but they might be the few that differ in the various log files. fix that
#ok, now I should rerun analyses essentially excluding those bad snps (this last part just a quick check, actually not necessary).
#the crazy thing is that all this occurs because I have to recalculate how many comparisons for fdr, otherwise it wouldn't be needed. how can I avoid that step? I could just use version with no pvalues,
#but anyway I have to do bed intersect and stuffs. Maybe just code that does that? I don't think I can escape it.
#Actuall if I have freqs per pop..still all the pairwise. But by pop. All snps that 
#if I have file with above thresholds for each pop then I could just do bedintersect on that and take..actually best rerunning everything again.
}
#-----------------------------------------------------------------------------------------------------------------------run again ICLD on cleaned data
#nohup ./ICLD11multi.sh /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1 2 coding 1 12 1 23 nspops.txt &
#./postproc_split2.sh above95/coding1000g/removal1
#-------------------------------------------------------------ANALYSES-removal1-----------------------------------------------------------------------
{
#cat above95/coding1000g/removal1/anal/sorted.res | awk '{if (NF==23){print}}' > above95/coding1000g/removal1/anal/sorted2.res
data<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/sorted2.res")
names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
system("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/*minilog | grep ncomparisons | awk '{print $14}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/ncomparisons.txt")
ncomp<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/ncomparisons.txt",header=FALSE)
ncomp<-sum(ncomp$V1)
dataFIN<-subset(data,data$pop==7)
combined_fdr<-p.adjust(p=dchisq(dataFIN$X,df=2*dataFIN$popX), method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
length(combined_fdr)
dim(dataFIN)
dim(data)
fullsign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
#in early code called sign
rm(data,sign,dataFIN)
dim(fullsign)
str(fullsign)
rm(chrA,chrB)
signAbed<-as.data.frame(cbind(as.character(paste0('chr',fullsign$chrA)),fullsign$posA-1,fullsign$posA),stringsAsFactors =FALSE)
signBbed<-as.data.frame(cbind(as.character(paste0('chr',fullsign$chrB)),fullsign$posB-1,fullsign$posB),stringsAsFactors =FALSE)
names(signAbed)<-c("chr","start","end")
names(signBbed)<-c("chr","start","end")
signAbed$start<-as.numeric(signAbed$start)
signAbed$end<-as.numeric(signAbed$end)
signBbed$start<-as.numeric(signBbed$start)
signBbed$end<-as.numeric(signBbed$end)
#create logos
write.table(signAbed,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/snpsA.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/snpsB.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
bedtools intersect -a <(cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/snpsA.bed | sed 's/^/chr/g') -b /mnt/sequencedb/ucsc/goldenPath/hg19/snp142.bed.gz -sorted -wb | awk -v OFS='\t' '{print $4,$5,$6,$7}' >  /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/snpsA.bed
bedtools intersect -a <(cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/snpsB.bed | sed 's/^/chr/g') -b /mnt/sequencedb/ucsc/goldenPath/hg19/snp142.bed.gz -sorted -wb | awk -v OFS='\t' '{print $4,$5,$6,$7}' >  /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/snpsB.bed
#for trios
dim(fullsign[fullsign$pfisher<0.05,]) #1290
fortrios<-fullsign[fullsign$pfisher<0.05,]
#subset(fortrios,fortrios[,1]==10)[,1:4]
triosAbed<-as.data.frame(cbind(as.character(paste0('chr',fortrios$chrA)),fortrios$posA-1,fortrios$posA),stringsAsFactors =FALSE)
triosBbed<-as.data.frame(cbind(as.character(paste0('chr',fortrios$chrB)),fortrios$posB-1,fortrios$posB),stringsAsFactors =FALSE)
names(triosAbed)<-c("chr","start","end")
names(triosBbed)<-c("chr","start","end")
triosAbed$start<-as.numeric(triosAbed$start)
triosAbed$end<-as.numeric(triosAbed$end)
triosBbed$start<-as.numeric(triosBbed$start)
triosBbed$end<-as.numeric(triosBbed$end)
write.table(triosAbed,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINA.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(triosBbed,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINB.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
combined_fdr<-p.adjust(p=fullsign$pfisher, method = "fdr", n = ncomp)
length(combined_fdr[combined_fdr<0.05]) #0
combined_fdr<-p.adjust(p=fullsign$pKuli, method = "fdr", n = ncomp)
length(combined_fdr[combined_fdr<0.05]) #0
combined_fdr<-p.adjust(p=fullsign$T2, method = "fdr", n = ncomp)
length(combined_fdr[combined_fdr<0.05]) #96
fullsign[combined_fdr<0.05,]
triosT2Abed<-fullsign[combined_fdr<0.05,]
triosT2Bbed<-fullsign[combined_fdr<0.05,]
triosT2Abed<-as.data.frame(cbind(as.character(paste0('chr',triosT2Abed$chrA)),triosT2Abed$posA-1,triosT2Abed$posA),stringsAsFactors =FALSE)
triosT2Bbed<-as.data.frame(cbind(as.character(paste0('chr',triosT2Bbed$chrB)),triosT2Bbed$posB-1,triosT2Bbed$posB),stringsAsFactors =FALSE)
names(triosT2Abed)<-c("chr","start","end")
names(triosT2Bbed)<-c("chr","start","end")
triosT2Abed$start<-as.numeric(triosT2Abed$start)
triosT2Abed$end<-as.numeric(triosT2Abed$end)
triosT2Bbed$start<-as.numeric(triosT2Bbed$start)
triosT2Bbed$end<-as.numeric(triosT2Bbed$end)
write.table(triosT2Abed,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINT2A.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(triosT2Bbed,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINT2B.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
system("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINA.bed /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINB.bed | sort -Vu -k1,1 -k2,2 > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FIN.bed")
system("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINT2A.bed /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINT2B.bed | sort -Vu -k1,1 -k2,2 > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FINT2.bed")

if (FALSE)
{
rm /mnt/scratch/fabrizio/LDLD/above95/coding1000g/FIN.bed
touch /mnt/scratch/fabrizio/LDLD/above95/coding1000g/FIN.bed
cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/FIN.bed
for i in `seq 1 22`; do 
echo $i
zcat /mnt/sequencedb/ucsc/goldenPath/hg19/snp142.bed.gz | grep chr$i$'\t' > /mnt/scratch/fabrizio/LDLD/temp2.bed
bedtools intersect -a /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/anal/FIN.bed -b /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/filtered/temp/missense.chr$i.temp.e -sorted -wb  >>  /mnt/scratch/fabrizio/LDLD/above95/coding1000g/FIN.bed
done
#cat above95/coding1000g/removal1/anal/FINA.bed above95/coding1000g/removal1/anal/FINB.bed | sort -Vu -k1,1 -k2,2 | wc -l
}

#--------------------- creates polymorphism per individual (mutlog) file (newer)   ------------------v
#for i in `seq 1 22`; do 
#./filterbyfreq.out above95/coding1000g/removal1/chr$i.tab above95/coding1000g/removal1/chr$i.freqlog above95/coding1000g/removal1/chr$i.tab.mutlog 0.05 1300 & 
#done
#system("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/chr*.tab.mutlog > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/tot.mutlog")
#-------------------- create whole genome freqlog file (genotype per sample) ---------------------------------------------------
#rm /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/all.freqlog
#for i in `seq 1 22`; do 
#cat above95/coding1000g/removal1/chr$i.freqlog | awk -v FS='\t' -v OFS='\t' -v CHR=$i '{if (NR>1){print CHR,$0}}' >> /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/all.freqlog
#done

source("analysesLDLD_header.R")
#--resres is gen file for all pop together
resres<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/all.freqlog")
length(resres)
rress<-list()
for (i in 1:12)
{
rress[[i]]<-resres[,samples_imypopi(nspops,i)+2] #2*i comes from chr and pos fields
}
snpsnames<-sapply(1:dim(resres)[1],function(x) paste0(resres[x,1:2],collapse="."))
#resres<-rress #let's call resres as gen (file with the single snp genotypes
removal1_resresgen<-rress
save(removal1_resresgen,file="removal1_resresgen.RData")



for (i in 0:11)
{
  system(paste0("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/chr*.log | awk '{if ($NF==",i,"){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/all.pop",i,".log"))
}

myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1"
logpop<-pairsANDlogfiles2logpop(myfolder,mysign=fullsign,12,22)

#file.exists("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/temp.log")

for (i in 1:22)
{
if (i==1){mutperssample<-as.numeric(read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/chr",i,".tab.mutlog")))}
else {res2<-as.numeric(read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/chr",i,".tab.mutlog")));mutperssample<-rbind(mutperssample,res2)}
}
l_mutsamples<-lapply(1:12,function(x) c())
for (i in 1:12)
{
  l_mutsamples[[i]]<-apply(mutperssample[,samples_imypopi(nspops,i)],2,sum)
}
save(l_mutsamples,file="removal1_l_mutsamples.RData")

par(mfrow=c(1,1))
for (i in 1:12)
{
pdf(paste0("logplot_removal1_",i,".pdf"))
infoplot(i)
dev.off()
}

logpoptot<-pairsANDlogfiles2logpop("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1",ntot_pops=12,summarize=FALSE)
myblocks_removal1<-create_blocks(logpoptot[[1]]) #still undivided
setwd("~/Dropbox/LDLD")
save.image("removal1_0.RData")
load("removal1_0.RData")
#check every pop

system.time(test_removal1<-test_nulldistr_mixtures(12,1000,logpoptot,logpop,myblocks_removal1))
save(test_removal1,file="removal1_test.RData")
save(logpoptot,file="removal1_logpoptot.RData")
save(logpop,file="removal1_logpop.RData")
save(myblocks_removal1,file="removal1_myblocks_removal1")
load(removal1_logpop.RData);load("removal1_l_mutsamples.RData");load("removal1_resresgen.RData")
removal1_removesnps<-corr_nAvsnAB(removal1_resresgen,logpop,l_mutsamples)
save(removal1_removesnps,file="removal1_removesnps.RData")
removal1_removesnps[removal1_removesnps[,2]<=0.05,] #also homozygotes dim(removed_snps) #612, 3
removed_snps #empty! possibly because I try to clean whole data-sets, not restricting to the linked ones, so that I have a very big n for the fdr correction. #let's try restricting it.
hist()
#almost all pops are still significant for pval_sd, while for the other statistics it is much 
ipop<-5
test_removal1[[ipop]]$pval_sd
test_removal1[[ipop]]$pval_nLR
test_removal1[[ipop]]$pval_pLR
test_removal1[[ipop]]$pval_nAIC
test_removal1[[ipop]]$pval_pAIC

system.time(testfull3<-test_nulldistr_mixtures(3,1000,logpoptot,logpop,myblocks))
system.time(testfull6<-test_nulldistr_mixtures(6,1000,logpoptot,logpop,myblocks))
system.time(testfull9<-test_nulldistr_mixtures(9,1000,logpoptot,logpop,myblocks))
system.time(testfull12<-test_nulldistr_mixtures(12,1000,logpoptot,logpop,myblocks))

res[[3]]$rawmodels
res[[5]]$sd
res[[5]]$pval_sd
res[[5]]$pval_pLR
res[[5]]$pval_pAIC
plot(logpop[[5]][[1]][order(logpop[[5]][[1]])])
res[[3]]$sd
sd(logpop[[5]][[1]])
partition_samples_nAB(logpop[[5]][[1]])[c(1,4,7,8)]

sapply(dim(res[[3]]$raw)[2], function(x) sd(res[[3]]$raw[,x]))
sapply(dim(res[[3]]$raw)[2], function(x) mean(res[[3]]$raw[,x]))
sd(logpop[[3]][[1]])
mean(logpop[[3]][[1]])


dim(res[[3]]$raw)

res[[3]]$raw
partition_samples_nAB(logpop[[1]][[1]])[c(1,4,7,8)]
logpop[[1]][[1]]

sum(as.numeric(as.numeric(nulldistr1iteration[,1])>=logpop[[1]][[1]]))/1000 # 0.053 it seems that by taking only number of groups, both with LR and AIC I don't even get significant, although almost
sum(as.numeric(as.numeric(nulldistr1iteration[,2])>=logpop[[1]][[1]]))/1000 # 0.115
sum(as.numeric(as.numeric(nulldistr1iteration[,3])<=logpop[[1]][[1]]))/1000 # 0 instead by taking p-values or RL I get extremely significant.
sum(as.numeric(as.numeric(nulldistr1iteration[,4])<=logpop[[1]][[1]]))/1000 #0

sapply(nulldistr1iteration[,3],function(x) if(is.null(x)){1}else{x})
nulldistr1iteration


#very cool, now which is the best way to test if I can get this fragmentation by chance?

mydata<-permbychr[,1]
sum(permbychr[,10000]) #428177
system.time(nulldistr1iteration<-generate_nulldistr_mixtures(permbychr,1000)) #time 448.160
as.numeric(nulldistr1iteration[,1])
as.numeric(nulldistr1iteration[,2])
as.numeric(nulldistr1iteration[,3])
as.numeric(nulldistr1iteration[,4])
mypop[4]
#compare value with distribution
mypop<-partition_samples_nAB(logpop[[1]][[1]])[c(1,4,7,8)]

sum(as.numeric(as.numeric(nulldistr1iteration[,1])>=mypop[1]))/1000 # 0.053 it seems that by taking only number of groups, both with LR and AIC I don't even get significant, although almost
sum(as.numeric(as.numeric(nulldistr1iteration[,2])>=mypop[2]))/1000 # 0.115
sum(as.numeric(as.numeric(nulldistr1iteration[,3])<=mypop[3]))/1000 # 0 instead by taking p-values or RL I get extremely significant.
sum(as.numeric(as.numeric(nulldistr1iteration[,4])<=mypop[4]))/1000 #0

logpoptot[[1]]

nulldistr1iteration

partition_samples_nAB(permbychr[,3])[c(1,4,7,8)]

tempres<-partition_samples_nAB(permbychr_removal1[,10])[c(1,4,7,8)]
unlist(tempres)

permbychrblocks2<-permute_chrAB_byblock(mylog,myblocks2,100000)


#to do: permutations and test if I can go forward with removing.
#if I can, just go forward and further clean up.


tempres<-corr_nAvsnAB(resres,logpop,l_mutsamples)
save.image("resres_complete_removal1.RData")
load("resres_complete_removal1.RData")
#removed_snps<-tempres[tempres[,2]<=0.05,] #only heterozygous
removed_snps<-tempres[tempres[,3]<=0.05,] #also homozygotes dim(removed_snps) #612, 3
stringsnp<-removed_snps

removed_snps<-stringsnp2bed(removed_snps)
write.table(removed_snps,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed1.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
for (mychr in 1:22)
{
system(paste0("bedtools intersect -v -b /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removed1.bed -a /mnt/scratch/fabrizio/LDLD/above95/coding1000g/chr",mychr,".tab -header | gzip > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/chr",mychr,".coding.vcf.gz"))
}






}
#--------------------------------------------------^
#vvvvvvvvvvvvvvvvvvv-EURgroup1 -vvvvvvvvvvvvvvvvvvvvvv
#-----separate components--------------------------v
{
getwd()
i<-8
lmygroup<-c()
mygroups<-c()
nAB_l<-unlist(logpop[[i]][1])
ord<-order(nAB_l[1:length(nAB_l)])
res<-partition_samples_nAB(nAB_l[ord],maxk=6)
samples_imypopi(nspops,i)[ord][res[[2]]==1]
mygroup<-samples_imypopi(nspops,i)[ord][res[[2]]==1]
lmygroup<-c(lmygroup,length(mygroup))
mygroups<-c(mygroups,mygroup)
i<-1
nAB_l<-unlist(logpop[[i]][1])
ord<-order(nAB_l[1:length(nAB_l)])
res<-partition_samples_nAB(nAB_l[ord],maxk=6)
mygroup<-samples_imypopi(nspops,i)[ord][res[[2]]==1]
lmygroup<-c(lmygroup,length(mygroup))
mygroups<-c(mygroups,mygroup)
i<-2
nAB_l<-unlist(logpop[[i]][1])
ord<-order(nAB_l[1:length(nAB_l)])
res<-partition_samples_nAB(nAB_l[ord],maxk=6)
mygroup<-samples_imypopi(nspops,i)[ord][res[[2]]==1]
lmygroup<-c(lmygroup,length(mygroup))
mygroups<-c(mygroups,mygroup)
lmygroup #21 60 50
samples_order1000g_above95<-paste0(intersect(samplesabove95,samples1000g),collapse=',')
samples_orderEURgroup1<-paste0(intersect(samplesabove95,samples1000g)[mygroups],collapse=',')

order_samples(samples_orderEURgroup1,"/mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1","coding",chrfrom=1,chrto=22)
}
#--------------------------------------------------------------------------------------------^
#--------------analyses cleaned samples------------------------------------------------------v
{
load("LDLD170216.RData")
nspops<-c(21,60,50)
data<-read.table("above95/coding1000g/EURgroup1/anal/sorted.res")
names(data)<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
system("cat above95/coding1000g/EURgroup1/*minilog | grep ncomparisons | awk '{print $5}' > above95/coding1000g/EURgroup1/ncomparisons.txt")
ncomp<-read.table("above95/coding1000g/EURgroup1/ncomparisons.txt",header=FALSE)
ncomp<-sum(ncomp$V1)
dataFIN<-subset(data,data$pop==0)
combined_fdr<-p.adjust(p=dchisq(dataFIN$X,df=2*dataFIN$popX), method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
length(combined_fdr)
dim(dataFIN)
dim(data)
sign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
sign01<-cbind(dataFIN,combined_fdr)[combined_fdr<0.01,]
sign001<-cbind(dataFIN,combined_fdr)[combined_fdr<0.0001,]
write.table(sign,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/sign.res",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
write.table(rbind(cbind(sign[,1],sign[,3]-1,sign[,3]),cbind(sign[,2],sign[,4]-1,sign[,4])),file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/sign.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
#cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/sign.bed | sort -Vu -k1,1 -k2,2 | sed 's/^/chr/g' | grep -v V > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/EURsign.bed
dim(sign) #1284

resres<-list(data.frame(),data.frame(),data.frame()) #add chr info to log files to retrieve info for snp
for (i in 0:2)
{
  counter<-1
  for (ichrA in 1:21)
    {
    for (ichrB in (ichrA+1):22)
      {
      print(c(ichrA,ichrB))
      system(paste0("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/chr",ichrA,".tabchr",ichrB,".tab.log | awk '{if ($NF==",i,"){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/all.pop.temp.log"))
      data<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/all.pop.temp.log",header=FALSE)
      data$chrA<-ichrA
      data$chrB<-ichrB
      if (counter==1) {res<-data} else {res<-rbind(res,data)}
      if (ichrA==21 && ichrB==22) {resres[[i+1]]<-res}
      counter<-counter+1
      }
    }
}
j<-1;resres[[j]]<-resres[[j]][,c(dim(resres[[j]])[2]-1,dim(resres[[j]])[2],3:dim(resres[[j]])[2]-2)]
j<-2;resres[[j]]<-resres[[j]][,c(dim(resres[[j]])[2]-1,dim(resres[[j]])[2],3:dim(resres[[j]])[2]-2)]
j<-3;resres[[j]]<-resres[[j]][,c(dim(resres[[j]])[2]-1,dim(resres[[j]])[2],3:dim(resres[[j]])[2]-2)]


#ok, if true it would be cool. Repeat per sample analyses and see how it is now.
for (i in 0:2)
{
  system(paste0("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/chr*.log | awk '{if ($NF==",i,"){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/all.pop",i,".log"))
}

logpop<-lapply(1:3,function(x) list(c(),c())) #create logpop objects, nested lists with pop (1 level) and nAB and nAandB per sample (2 level)
for (i in 0:2)
{
  mydata<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/all.pop",i,".log"),header=FALSE)
  logpop[[i+1]][[1]]<-sapply(3:(dim(mydata)[2]-1), function(x) sum(sapply(mydata[,x],nABf)))
  logpop[[i+1]][[2]]<-sapply(3:(dim(mydata)[2]-1), function(x) sum(sapply(mydata[,x],nAandBf)))
}

infosamples<-lapply(1:3,function(x) c()) #extract info 1000g per sample
for (i in 1:3)
{
  mysamples<-samples_imypop(intersect(samplesabove95,samples1000g)[mygroups],nspops,i)
  infosamples[[i]]<-info1000g[match(mysamples,info1000g$Sample),]
}

dim(infosamples[[1]])
dim(infosamples[[2]])
dim(infosamples[[3]])

plot.new()
pdf("EURgroup1_infoplot.pdf")
par(mfrow=c(2,2))
infoplot(1)
infoplot(2)
infoplot(3)
dev.off()

dim(sign)
i<-0
mydata<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/EURgroup1/all.pop",i,".log"),header=FALSE)
paste0(sign[,1],".",sign[,3],".",sign[,2],".",sign[,4])
mydata[,1]
res<-rebuild_dataLDLD(paste0(sign[,1],".",sign[,3],".",sign[,2],".",sign[,4]),mydata,fromremovedsnps=FALSE)
head(res)

myblocks<-create_blocks(resres[[1]])
permbychr<-permute_chrAB_byblock(resres[[1]],myblocks,10000)
sum(permbychr[,1]) #12455.5
sum(logpop[[1]][[1]]) #12455.5
permvar<-sapply(1:10000,function(x) var(permbychr[,x]))
sum(as.numeric(permvar>var(logpop[[1]][[1]])))/10000 #0.0246
system.time(nulldistr1iteration<-generate_nulldistr_mixtures(permbychr,1000)) #takes long and only 21 samples
mypop<-partition_samples_nAB(logpop[[1]][[1]])[c(1,4,7,8)]
mypop<-partition_samples_nAB(logpop[[1]][[1]])
mypop[1] #funny, less divisions with AICc, an exception 2 vs 5. With 2 I would remove small cluster

permbychr<-permute_chrAB_byblock(resres[[2]],myblocks,10000)
permvar<-sapply(1:10000,function(x) var(permbychr[,x]))
sum(as.numeric(permvar>var(logpop[[2]][[1]])))/10000 #0
system.time(nulldistr1iteration<-generate_nulldistr_mixtures(permbychr,100)) #takes long and only 21 samples
mypop<-partition_samples_nAB(logpop[[2]][[1]])[c(1,4,7,8)]
mypop[1] #5 and 5
sum(as.numeric(as.numeric(nulldistr1iteration[,1])>=mypop[1]))/100 # 0.12 it seems that by taking only number of groups, both with LR and AIC I don't even get significant, although almost
sum(as.numeric(as.numeric(nulldistr1iteration[,2])>=mypop[2]))/100 # 0.18
sum(as.numeric(as.numeric(nulldistr1iteration[,3])<=mypop[3]))/100 # 0.04
sum(as.numeric(as.numeric(nulldistr1iteration[,4])<=mypop[4]))/100 # 0.05 #very nice, pvalue not so low if taking this values. Actually I have some hope of finding no differences

permbychr<-permute_chrAB_byblock(resres[[3]],myblocks,10000)
permvar<-sapply(1:10000,function(x) var(permbychr[,x]))
sum(as.numeric(permvar>var(logpop[[3]][[1]])))/10000 #0.0053
system.time(nulldistr1iteration<-generate_nulldistr_mixtures(permbychr,100)) #takes long and only 21 samples
mypop<-partition_samples_nAB(logpop[[3]][[1]])
mypop<-partition_samples_nAB(logpop[[3]][[1]])[c(1,4,7,8)] #4,4
sum(as.numeric(as.numeric(nulldistr1iteration[,1])>=mypop[1]))/100 # 0.13 it seems that by taking only number of groups, both with LR and AIC I don't even get significant, although almost
sum(as.numeric(as.numeric(nulldistr1iteration[,2])>=mypop[2]))/100 # 0.21
sum(as.numeric(as.numeric(nulldistr1iteration[,3])<=mypop[3]))/100 # 0.02
sum(as.numeric(as.numeric(nulldistr1iteration[,4])<=mypop[4]))/100 # 0.03 #very nice, pvalue not so low if taking this values. Actually I have some hope of finding no differences

for (i in 1:10000)
{permbychr[,i]<-permbychr[,i][order(permbychr[,i])]}
dev.new()
plot.new()
dim(permbychr)
library(vioplot)
pdf("empiricaldistrperm_EURpop3.pdf")
par(mfrow=c(2,1))
vioplot(permbychr[1,],permbychr[10,],permbychr[20,],permbychr[30,],permbychr[40,],permbychr[50,],ylim=c(0,1000),col="firebrick",colMed="firebrick2")
mypop<-logpop[[3]][[1]][order(logpop[[3]][[1]])]
points(c(mypop[1],mypop[10],mypop[20],mypop[30],mypop[40],mypop[50]),pch=19,col="black")
hist(log(permvar),xlim=c(log(min(permvar)),log(max(max(permvar),var(logpop[[3]][[1]]))+10)))
abline(v=log(var(logpop[[3]][[1]])),col="red")
dev.off()

save.image("LDLD170216.RData")
load("LDLD170216.RData")








#for i in /mnt/454/HGDP/genomes_Bteam/HG19Align/*.bam /mnt/sequencedb/11men/martin/hg19_1000g/BAM/*.bam ; do ./BamTable2 -r 12:159672251-159672252 $i ; done | grep 159672252
myres<-bamlogo("above95/coding1000g/logonotsign.temp",tilltheend=TRUE,startfrom=20001,filetemp="temp20k")
myres2<-as.data.frame(myres,row.names=FALSE,filetemp="temp20k")
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call","tot_samples")
write.table(myres2,file="above95/coding1000g/anal/logo2_20ktoend_notsign_temp.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE);
data1<-read.table("above95/coding1000g/anal/logo2_1to5k_notsign_temp.bed",stringsAsFactors =FALSE,header=TRUE)
data2<-read.table("above95/coding1000g/anal/logo2_5kto10k_notsign_temp.bed",stringsAsFactors =FALSE,header=TRUE)
data3<-read.table("above95/coding1000g/anal/logo2_10kto20k_notsign_temp.bed",stringsAsFactors =FALSE,header=TRUE)
names(data1)<-names(data3)
names(data2)<-names(data3)
names(data3)<-names(data3)
data<-rbind(data1[,1:12],data2[,1:12],data3[,1:12])
write.table(data,file="above95/coding1000g/anal/logo2_notsign.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)
logonotsign<-read.table("above95/coding1000g/anal/logo2_notsign.bed",header=TRUE,stringsAsFactors =FALSE)
#
logonotsign<-subset(logonotsign,logonotsign$tot_samples>0)
#NB: a large fraction of not significant are not in A and B teams! This speaks for a bias in things that are actually not significant!
dim(logonotsign)
head(logonotsign)
prHvsE(logonotsign,1)

head(data)
badata<-subset(data,is.na(data$C))
dim(data)
badata<-subset(data,data$tot_samples!=25)
head(badata)

myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=20)
myres2<-as.data.frame(myres,row.names=FALSE)
names(myres2)<-c("chr","start","end","ref","alt","seq5p","A","C","G","T","tot_second_call")
write.table(myres2,file="above95/coding1000g/anal/logo2_20.bed",row.names = FALSE,quote=FALSE,sep='\t',col.names=TRUE)
head(myres2)
library(xtable)
options(xtable.floating = FALSE)
xtable(myres2[1:20, ])

myres<-bamlogo("above95/coding1000g/logosign.bed",tilltheend=FALSE,till=20002,startfrom=20001,filetemp="temp5")
}
#======================================================================TAG_MISSENSEBIS=====================================================
#======================================================================TAG_MISSENSE========================================================

