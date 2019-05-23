#PER SAMPLE ANALYSES

setwd("~/Dropbox/LDLD")
source ("analysesLDLD_header.R")
options(scipen=999)
setwd("/mnt/scratch/fabrizio/LDLD")

#---import files--------------------------------------------------------------------------------------v
info1000g<-read.table("/mnt/scratch/fabrizio/LDLD/20130606_sample_info.txt",header= TRUE ,sep='\t')
samples1000g<-system("zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr22.phase3_shapei
t2_mvncall_integrated_v5.20130502.genotypes.vcf.gz |
head -300 | grep '#CHROM' | head -1 | awk '{for (i=10; i<=NF; i++) print $i}'",intern= TRUE )
samplesabove95<-system("cat ~/workspace/1000genomes/above95.unrelated.samples",intern= TRUE )
nspops<-as.numeric(system("cat /mnt/scratch/fabrizio/LDLD/nspops.txt",intern= TRUE ))
samples_order1000g_above95<-intersect(samplesabove95,samples1000g)
infosamples<-lapply(1:12, function (x) c())
for (i in 1:12)
{
mysamples<-samples_imypop(samples_order1000g_above95,nspops,i)
infosamples[[i]]<-info1000g[match(mysamples,info1000g$Sample),]
}
popsamplesi<-sapply(1:12, function (x) samples_imypopi(nspops,x))
#system("cat above95/coding1000g/chr*.tabchr*.minilog | grep nAB_above_freq_threshold | grep -v and_significant > above95/coding1000g/tot.nABabovethr")
#system("cat above95/coding1000g/chr*.tabchr*.minilog | grep nAB_above_freq_threshold | grep and_significant > above95/coding1000g/tot.nABabovethrsign")
#mylog2<-read.table("above95/coding1000g/tot.nABabovethr")
#mylog2<-sapply(2:dim(mylog2)[2],function(x) sum(mylog2[,x]))
#mylog3<-read.table("above95/coding1000g/tot.nABabovethrsign")
#mylog3<-sapply(2:dim(mylog3)[2],function(x) sum(mylog3[,x]))
#-------------------- create whole genome log file -----------------------------
#for (i in 0:11)
#{
# system(paste0("cat above95/coding1000g/chr*.log | awk '{if ($NF==",i,"){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/all.pop",i,".log"))
#}
#-------------------- create whole genome log file by chr ----------------------
#system(paste0("mkdir above95/coding1000g/logs"))
#for (ipop in 0:11)
# {
# for (ichrA in 1:21)
#{
# for (ichrB in (ichrA+1):22)
#{
#system(paste0("cat above95/coding1000g/logs/chr*.log > above95/coding1000g/logs/chr",ichrA,".",ichrB,".pop",ipop,".log"))
#}
#}
# }
#for (i in 0:11)
#{
# system(paste0("cat above95/coding1000g/logs/chr*.log | awk '{if ($NF==",i,"){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/logs/all.pop",i,".log"))
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
#--------------------- creates polymorphism per individual (mutlog) file (newer) ------------------v
#for i in `seq 1 22`; do
#./filterbyfreq.out above95/coding1000g/chr$i.tab above95/coding1000g/chr$i.freqlog above95/coding1000g/chr$i.tab.mutlog 0.05 1220 &
#done
#system("cat above95/coding1000g/chr*.tab.mutlog > above95/coding1000g/tot.mutlog")
for (i in 1:22)
{
if (i==1){mutperssample<-as.numeric(read.table(paste0("above95/coding1000g/chr",i,".tab.mutlog")))} else {res2<-as.numeric(read.table(paste0("above95/coding1000g/chr",i,".tab.mutlog")));mutperssample<-rbind(mutperssample,res2)}
}
l_mutsamples<-lapply(1:12, function (x) c())
for (i in 1:12)
{
l_mutsamples[[i]]<-apply(mutperssample[,samples_imypopi(nspops,i)],2,sum)
}
