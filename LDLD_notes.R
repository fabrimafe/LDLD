#reasoning about linkage
{
A<-0.1
B<-0.1
ABHW<-A*B
abHW<-(1-A)*(1-B)
ABHW-A*B
(ABHW+0.01)-A*B #to keep A and B constant, I need also to have more ab.
#more AB and ab, versus Ab and aB.
#which is really increased, AB or ab?
#it could even be that more widespread mutations as hets -> could it be that it gives negative linkage. I think only if in ML assume 1/2 weight, but actually they are most often together. I should simulate and check!

source("~/Dropbox/LDLD/analysesLDLD_inR_header.R")
#mutations that introduce hets at random: no effect on LDLD. You definitely need batches

myf<-function(){
mypop<-simulateHW_gtsin012format(100,100,1000)
null<-f_DeltaAB(gtsin012format2gts(mypop))
mypop[1,mypop[1,]==0]<-sample(c(0,1),sum(mypop[1,]==0), prob=c(0.9,0.1),replace=TRUE)
gtsin012format2gts(mypop)
temp<-f_DeltaAB(gtsin012format2gts(mypop)) #hets added on only one
mypop[2,mypop[2,]==0]<-sample(c(0,1),sum(mypop[2,]==0), prob=c(0.9,0.1),replace=TRUE)
temp2<-f_DeltaAB(gtsin012format2gts(mypop)) #hets added on both
c(null,temp,temp2)
}

res<-sapply(1:10000,function(x) myf() )
apply(res,MARGIN=1,mean)
[1] 1.035000e-05 1.124500e-05 2.373395e-05




}
#reasoning about mismatch in number of candidates for coding
{
#for coding
#I tried scan from R, and both in the old way (that I used to get many snps), or in the new way I get very bad
#overlap with CG.
#In the old way (repl1) the overlap was very good, with about 600 variants.
#With shell, intergenic has very good overlap.
#I get few variants in coding, with very good overlap.

#overall it seems that whatever I do with script works, from within R it does not. Still too few snps. Why?
#First of all, what are the files I should look at?
#with glmm on filtered populations: for which I noticed that very low power, so I get just a few snps): coding
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/removed1_pos_shell_scale.bed (short list) 
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/bad_snps_min5_header_scale/bad_snps_pos.tab p-values to see if marginally significant also enriched.

#I should actually be sure of where the files come from.
#First of all I should replace vcfs in repl2 with original ones (since overwritten)
#Then I should chekc nA files.
#I could run once this with high threshold, to be sure that I process all positions, maybe on a single chromosome. Also to check if with low threshold for glm I discard too much.
#Then I should run for min1 and min5.
#I should do that for glm, glmm all and glmm filtered, glm predicting hets.

#check nAB itself (strong outliers)
#try using log (does it change?)











#with glmm on filtered populations: for which I noticed that very low power, so I get just a few snps): intergenic
#-/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/bad_snpsremoved1_pos20.bed


#
fabrizio_mafessoni@bionc05:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g$ ls -lh repl1/min5/chr1.vcf.gz
-rw-rw-rw- 1 fabrizio_mafessoni staff 2.7M Aug 26 12:51 repl1/min5/chr1.vcf.gz
fabrizio_mafessoni@bionc05:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g$ ls -lh repl2/min5/chr1.vcf.gz
-rw-rw-rw- 1 fabrizio_mafessoni staff 2.7M Aug 26 12:50 repl2/min5/chr1.vcf.gz
fabrizio_mafessoni@bionc05:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g$ ls -lh repl3/min5/chr1.vcf.gz
-rw-rw-rw- 1 fabrizio_mafessoni staff 2.5M Aug 24 12:22 repl3/min5/chr1.vcf.gz
fabrizio_mafessoni@bionc05:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g$ ls -lh repl4/min5/chr1.vcf.gz
-rw-rw-rw- 1 fabrizio_mafessoni staff 2.6M Aug 24 12:23 repl4/min5/chr1.vcf.gz



#it does not seem that is because R crashes since no empty files
#ls -lh /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/bad_snps_min5_header_scale/*

#what are the reasons old results might have been different.
#first thing that I should make sure is that R and code are the same, to be sure of what happens in the equations.
#other possibility is that differences because:
#-glm vs glmmm
#-all populations, not only above threshold.
#I will explore both possibilities, but first let's be sure of what happens.

#glm on filtered seems to give slightly more bad snps
#glmm on everything also the same number
#it could still be glm on everything
#min1
#other thing could be of course using mut. by looking in old scripts, what I used to do is glm on filtered, but there used to be mut.
#Why this would give more power? in theory it should only correct for more mutations in individual.
#A possibility is that I truly remove a confounding factor.

#I implemented it, but it seems that pvalues do not go down. I am left with only a solution. Working on the old code and see if what changes.
#If it changes I can then start dissecting it.
#If it does not, I am fucked, I don't know what it is.



}
#implementing option parsing with popt in LDLD
{
cd ~/Dropbox/LDLD/scripts
gcc popt_example.c -o temp.out -l popt

#going through the program I noticed that I don't do alternativetominor conversion for alleleB!!!
fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal$ vim artificial.res
fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal$ cat artificial.res | awk 'BEGIN{COUNTERPOP=0;MYFLAG=0}{COUNTERPOP=COUNTERPOP+1; if ($8/$9<0.5){MYFLAG=1};if (COUNTERPOP==12){if (MYFLAG==0){print};COUNTERPOP=0;MYFLAG=0};}' 
no no no no no no no 11 20
fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal$ zcat sorted.res.gz | awk 'BEGIN{COUNTERPOP=0;MYFLAG=0}{COUNTERPOP=COUNTERPOP+1; if ($8/$9<0.5){MYFLAG=1};if (COUNTERPOP==12){if (MYFLAG==0){print};COUNTERPOP=0;MYFLAG=0};}' | less
fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal$ zcat sorted.res.gz | awk 'BEGIN{COUNTERPOP=0;MYFLAG=0}{COUNTERPOP=COUNTERPOP+1; if ($7/$9<0.5){MYFLAG=1};if (COUNTERPOP==12){if (MYFLAG==0){print};COUNTERPOP=0;MYFLAG=0};}' | less




}