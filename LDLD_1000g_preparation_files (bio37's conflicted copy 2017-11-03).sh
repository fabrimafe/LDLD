#======================LDLD shell calculations==================
#files preparation
{
#compile programs
{
#LDLD
gcc ICLDmulti13.c -o test.out -lm -lgmp -l popt
cp ICLDmulti13.c 
#filtering program that replaces all the ones I used before: can output sharing of snps, mutational load, normal filtering or freqlogs.
gcc LDLD_filter13.c -o filterbyfreq_multi.out -lgmp -lm

}

#tests ICLDmulti
{
cd ~/Dropbox/LDLD/scripts
gcc ICLDmulti10.c -o ICLDmulti.out -lgmp -lm #compile
date
./ICLDmulti.out ~/Dropbox/LDLD/ICLDtests/chr21.tab ~/Dropbox/LDLD/ICLDtests/chr22.tab ~/Dropbox/LDLD/ICLDtests/nspops.txt ~/Dropbox/LDLD/ICLDtests/file4.res ~/Dropbox/LDLD/ICLDtests/file5.log ~/Dropbox/LDLD/ICLDtests/file6.sum ~/Dropbox/LDLD/ICLDtests/file7.minilog 0.05 12 0 #example command
date
./ICLDmulti.out ~/Dropbox/LDLD/ICLDtests/chr21.tab ~/Dropbox/LDLD/ICLDtests/chr22.tab ~/Dropbox/LDLD/ICLDtests/nspops.txt ~/Dropbox/LDLD/ICLDtests/file4.res ~/Dropbox/LDLD/ICLDtests/file5.log ~/Dropbox/LDLD/ICLDtests/file6.sum ~/Dropbox/LDLD/ICLDtests/file7.minilog 0.05 12 1 
#example command
date
./ICLDmulti.out ~/Dropbox/LDLD/ICLDtests/chr21.short.tab ~/Dropbox/LDLD/ICLDtests/chr22.short.tab ~/Dropbox/LDLD/ICLDtests/nspops.txt ~/Dropbox/LDLD/ICLDtests/file4.res ~/Dropbox/LDLD/ICLDtests/file5.log ~/Dropbox/LDLD/ICLDtests/file6.sum ~/Dropbox/LDLD/ICLDtests/file7.minilog 0.05 12 1 

#implement fixed resolution
{
#what is the best way? to have homogeneous resolution for all chromosomes I should first compute resolution necessary. Then calculate it and give it ICLD, that I guess at
# this point would stay basically the same. I should simply allow for wider range.
# Now, let's check how it scales.

fabrizio_mafessoni@bionc01:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1$ zcat min1/chr21.vcf.gz | grep -v '#' | wc -l
15018
fabrizio_mafessoni@bionc01:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1$ zcat min1/chr22.vcf.gz | grep -v '#' | wc -l
30197
fabrizio_mafessoni@bionc01:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1$ ls -lh min1/anal/chr21.chr22.res
-rw-rw-rw- 1 fabrizio_mafessoni staff 19M Jul 25 10:07 min1/anal/chr21.chr22.res
fabrizio_mafessoni@bionc01:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1$ cat min1/anal/chr21.chr22.res |wc -l
85764
fabrizio_mafessoni@bionc01:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1$ zcat min5/chr21.vcf.gz | grep -v '#' | wc -l
710
fabrizio_mafessoni@bionc01:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1$ zcat min5/chr22.vcf.gz | grep -v '#' | wc -l
1346
fabrizio_mafessoni@bionc01:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1$ cat min5/anal/chr21.chr22.res |wc -l
1476

# I want a precise overview of how much it scales. That means: for number of lines, how many comparisons.
#For number of comparisons/lines how many in output.
#Actual ratio for total number of comparisons used. Notice that similar number in coding and intergenic because most is non-significant.
cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1
for myfreq in 1 5;do
for chrA in `seq 1 21`;do
echo $chrA
chrAp=$(( $chrA + 1 ))
mylinesA=$( zcat min${myfreq}/chr${chrA}.vcf.gz | grep -v '#' | wc -l )
#ls min${myfreq}/chr${chrA}.vcf.gz
for chrB in `seq $chrAp 22`;do
mylinesB=$( zcat min${myfreq}/chr${chrB}.vcf.gz | grep -v '#' | wc -l )
#ls min${myfreq}/chr${chrB}.vcf.gz
myres=$( cat min${myfreq}/anal/chr${chrA}.chr${chrB}.res | wc -l )
mycomp=$( cat min${myfreq}/chr${chrA}.tabchr${chrB}.tab.minilog | grep ncomparisons | grep -v pair | awk '{print $14}' )
echo $myfreq $chrA $chrB $mylinesA $mylinesB $myres $mycomp >> files.sizes.tab
done;done;done

#+++++++++ coding.exon min5 chr1 fucked up!!! same date as min1!  but analyses done before, so ok!
#in coding 4x more 1x than 5x! in intergenic only 2x. -> in intergenic lower threshold will not save me much if just linear with snps.
 100/1000000

mydata=read.table("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/files.sizes.tab")
pdf("~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl1/file.sizes.tab.pdf")
par(mfrow=c(2,2))
names(mydata)<-c("freq","chrA","chrB","nsnpsA","nsnpsB","res.length","ncomp")
plot(mydata$ncomp,mydata$res.length,type='n',xlab="n.comp",ylab="res.length")
points(mydata[mydata$freq==5,]$ncomp,mydata[mydata$freq==5,]$res.length,col="red",pch=19)
points(mydata[mydata$freq==1,]$ncomp,mydata[mydata$freq==1,]$res.length,col="black",pch=19)
plot(10^(-6)*mydata$nsnpsA*mydata$nsnpsB,mydata$res.length,type='n',xlab="n^2",ylab="res.length")
points(10^(-6)*mydata[mydata$freq==5,]$nsnpsA*mydata[mydata$freq==5,]$nsnpsB,mydata[mydata$freq==5,]$res.length,col="red",pch=19)
points(10^(-6)*mydata[mydata$freq==1,]$nsnpsA*mydata[mydata$freq==1,]$nsnpsB,mydata[mydata$freq==1,]$res.length,col="black",pch=19)
dev.off()

mydata5<-mydata[mydata$freq==5,]
mydata1<-mydata[mydata$freq==1,]
sum(mydata5$res.length/1000)/sum(mydata5$ncomp/1000)
sum(mydata1$res.length/1000)/sum(mydata1$ncomp/1000)
sum((10^(-6))*mydata5$res.length)/sum((10^(-6))*mydata5$nsnpsA*mydata5$nsnpsB) #0.0001960168
sum((10^(-6))*mydata1$res.length)/sum((10^(-6))*mydata1$nsnpsA*mydata1$nsnpsB) #0.0001725131
#approximatively the same for freq 1 and 5. Should be a bit more because more by chance, but less because less shared.
#So, let's try then to get a fraction that is similar to what I use for 5
mythr5=10^(-6)
sum(mythr5*data5$nsnpsA*mydata5$nsnpsB) #6260.842
#I could set a threshold that is about 10000/nspns
mythr5=
(10^(-6))*10000/sum((10^(-6))*mydata5$nsnpsA*mydata5$nsnpsB)
#heuristically I should have that 0.05/n2. I stay large because it seems otherwise I lose some significant one (n sign is not much lower than printed with res)
#---> build lookuptable
myq<-sort(outer(seq(2,10,2),sapply(17:3,function(x) 10^(-x)),'*'))[c(-1,-2)]
mys<-paste0("{",paste0(sapply(myq, function(x)
paste0("{",paste0(sapply(1:20, function(i) qchisq(1-x, df=2*i)),collapse=','),"}")),collapse=','),"}")
options(scipen=999)
paste0("{",paste0(myq,collapse=','),"}")
#I should directly feed threshold

fabrizio_mafessoni@bio74:/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl3$ ~/Dropbox/LDLD/scripts/ncomparisons.sh ./min1
5421.65
fabrizio_mafessoni@bio74:/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl3$ ~/Dropbox/LDLD/scripts/ncomparisons.sh ./min5
1295.53



}

}
#faster reorder
{
MYCHROM=1
zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | head -100 | grep CHROM
}
#generate input files
{
#overview steps:
#
#1) downloading annotation to generate bed files for filtering (annotation.bed)
#2) generate chrZ.vcf.gz files (files chrZ.vcf.gz are global by annotation but with unordered samples) filtering on annotation.bed, then by frequency. For coding and introns these
#   files include all genome, for intergenic a subsample to match amount of data in introns.
#3) reorder samples chrZ.vcf.gz to get chrZ.dest.vcf.gz (files chrZ.dest.vcf.gz are global with ordered samples)
#4) filtering by frequency from to generate chrZ.temp.vcf.gz files in /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${annotation}/min${freq}/. exons, since no further steps are required, are moved directly to final chr.vcf.gz files in destination folder
#5) for introns and intergenic I have to further subsample:
#-a) Let's get something that rather than being #100x is 10x coding. Tested: very long too run LDLD
#-b) Let's generate a subsamples.bed file, subsampled like coding, that indicate the different subsamples for replicability
#

#1) downloading annotations
{
#list of exons
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g

REAL_TAB=$(echo -e "\t")
#"s/\|/$REAL_TAB/g"
#| sed 's/n\/a/./g'

#in hgnc.bed all genes
cat hgnc.bed | sed "s/$REAL_TAB$REAL_TAB/$REAL_TAB.$REAL_TAB/g" | sed "s/$REAL_TAB$REAL_TAB/$REAL_TAB.$REAL_TAB/g" | sed 's/chr//g' | grep -v random | sort -Vu -k1,1 -k2,2  > hgnc_filter.bed
#for exons I have in field 8 number of exons and $9 their start and $10 their ends
curl  -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c |\
 awk '{n=int($8); split($9,S,/,/);split($10,E,/,/); for(i=1;i<=n;++i) {printf("%s,%s,%s,%s,%s\n",$1,$2,$3,S[i],E[i]);} }'  #from https://www.biostars.org/p/6391/

cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
#otherwise I directly go on table browser and select bed format and only coding exons for known genes. This corresponds to
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbQual=cds&fbDownBases=200&hgta_doGetBed=get+BED" | bgzip -f > coding.exons.bed.gz
#introns
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbQual=intron&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED" | bgzip -f > introns.bed.gz
#whole genes
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED" | bgzip -f > genes.bed.gz

#there is wide overlap between these categories.
bedtools subtract -a introns.bed.gz -b coding.exons.bed.gz | bgzip -f > introns.bed
bedtools subtract -a genes.bed.gz -b introns.bed.gz | bgzip -f > genes.bed

#this done after preprocessing (in december 2017): should not affect anything:
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
bedtools merge -i <( bedtools sort -i coding.exons.bed.gz ) | gzip -f > coding.exons.merged.bed.gz
}

#2) generate chrZ.vcf.gz files
{
#extract variants from 1000g with annotation        
myannotation="coding.exons"
myannotation="introns"
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
for MYCHROM in `seq 11 22`; do
echo $MYCHROM
zcat ${myannotation}.bed  | awk -v OFS='\t' -v MYCHROM=$MYCHROM '{if ($1=="chr"MYCHROM){print MYCHROM,$2,$3}}' |sort -Vu -k1,1 -k2,2n > ${myannotation}/myfilter.${MYCHROM}.bed
bedtools merge -i ${myannotation}/myfilter.${MYCHROM}.bed > mytemp.${MYCHROM}.bed; mv mytemp.${MYCHROM}.bed ${myannotation}/myfilter.${MYCHROM}.bed
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -H -R myfilter.${myannotation}/myfilter.${MYCHROM}.bed  > ${myannotation}/chr$MYCHROM.vcf
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -R ${myannotation}/myfilter.${MYCHROM}.bed  | sort -Vu -k1,1 -k2,2n >> ${myannotation}/chr${MYCHROM}.vcf
bgzip -f ${myannotation}/chr${MYCHROM}.vcf
tabix ${myannotation}/chr${MYCHROM}.vcf.gz
done

#subsample intergenic to the same number of intronic (not really) and reorder them
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
SEED=1;MYCHROM=10
while [ $MYCHROM -lt 22 ]
        do
        myannotation="intergenic"
        echo $MYCHROM
        #ARGUMENT NLINES: set here the final number of lines
        NLINES=$( zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/introns/chr${MYCHROM}.dest.vcf.gz | wc -l )
        DESTINATIONFILE=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/intergenic/chr${MYCHROM}.intergenic.seed${SEED}.vcf
        myvcfsource=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
        NLINES_SOURCE=$( zcat $myvcfsource |wc -l  )
        MYPROB=$( awk -v var1=$NLINES -v var2=$NLINES_SOURCE 'BEGIN{print var1/var2}' )
        zcat $myvcfsource | head -500 | grep '#' > $DESTINATIONFILE
        zcat $myvcfsource | awk -v MYPROB=$MYPROB '{if (rand()<MYPROB){print}}' >> $DESTINATIONFILE
        #shuf -n $NLINES <( tabix ${myvcfsource} $MYCHROM | grep -v '#' ) | sort -Vu -k1,1 -k2,2 >> $DESTINATIONFILE #shuf crashes with big files
        bgzip -f $DESTINATIONFILE
        zcat genes.bed   | awk -v OFS='\t' -v MYCHROM=$MYCHROM '{if ($1=="chr"MYCHROM){print MYCHROM,$2-5000,$3+5000}}' | awk '{if ($2<0){$2=0};print}' | sort -Vu -k1,1 -k2,2n > ${myannotation}/myfilter.${MYCHROM}.bed
        bedtools merge -i ${myannotation}/myfilter.${MYCHROM}.bed > mytemp.${MYCHROM}.bed; mv mytemp.${MYCHROM}.bed ${myannotation}/myfilter.${MYCHROM}.bed
        bedtools subtract -a ${DESTINATIONFILE}.gz -b ${myannotation}/myfilter.${MYCHROM}.bed -sorted -header | bgzip -f > /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/intergenic/chr${MYCHROM}.vcf.gz
        MYCHROM=$(( $MYCHROM + 1 ))
        done
}        

#3) reorder
{
myannotation="coding.exons"
myannotation="introns"
myannotation="intergenic"
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
for MYCHROM in `seq 1 5`; do
echo $MYCHROM 'reorder...' 
nohup Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/above95.unrelated.samples ${myannotation}/chr${MYCHROM}.vcf.gz ${myannotation}/chr${MYCHROM}.dest.vcf /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/nspops.txt 0
echo $MYCHROM 'reorder reschr...' 
#for reschr in `seq 1 2`;do #
reschr=1
   mkdir ${myannotation}/reschr${reschr}
#    nohup Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/above95.unrelated.samples ${myannotation}/chr${MYCHROM}.vcf.gz ${myannotation}/reschr${reschr}/chr${MYCHROM}.dest.vcf /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/nspops.txt 1
#done
#rm ${myannotation}/myfilter.${MYCHROM}.bed
done


#reorder (only if I want an individual reordering for already subsampled files - after step 5)
myannotation="introns"
myannotation="intergenic"
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
for MYCHROM in `seq 9 9`; do
for myrepl in `seq 2 5`; do
for myfreq in 1 5; do
echo 'chrom :' $MYCHROM 'myrep ' $myrepl 'myfreq ' $myfreq
destination_folder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/repl${myrepl}; mkdir $destination_folder
destination_folder=$destination_folder/min${myfreq}; mkdir $destination_folder
mysourcefile=${destination_folder}/chr${MYCHROM}.vcf.gz
mkdir ${destination_folder}/reschr1
cmd="Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/above95.unrelated.samples $mysourcefile ${destination_folder}/reschr1/chr${MYCHROM}.vcf /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/nspops.txt 1"
qsub -cwd -b y -l h_vmem=2G,virtual_free=2G,mem_free=2G -N "ord${myannotation}" ~/Dropbox/Vindija/scripts/qsub_runner.sh $cmd
done;done;done


}

#4)filtering by frequency per population frequency threshold
{
#if subsampling to match the variants in coding.exons myannotation_ref

#quick filtering-------------v
cd ~/Dropbox/LDLD/scripts
gcc LDLD_filter9.c -o filterbyfreq_multi.out -lgmp -lm
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
cp ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out .

#---only for exons
myannotation="coding.exons"
for myfreq in 1 5; do
for MYCHROM in `seq 1 22`; do
echo $MYCHROM
MYVCF=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/repl1/min${myfreq}/chr${MYCHROM}.vcf.gz
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out nspops.txt 0.0${myfreq} 12 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYVCF}
done;done

#---global filtered file for introns (truly global) and intergenic (subsampled)
myannotation="introns"
myannotation="intergenic"
myannotation="coding.exons"
for myfreq in 1 5; do
for MYCHROM in `seq 9 22`; do
echo "chrom: " $MYCHROM "freq: " $myfreq
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
mkdir ${myannotation}/min$myfreq
MYFREQVCF=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/min$myfreq/chr${MYCHROM}.vcf
zcat ${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out nspops.txt 0.0${myfreq} 12 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYFREQVCF}.temp.gz
done;done


}

#5a) NO #---further subsample global files for introns and intergenic matching the size of exons
{
FOLDS=10 #how many times myannotation_ref (in repl11 I used fold10)
myannotation_ref="coding.exons"
myannotation="introns"
myannotation="intergenic"
for myfreq in 1 5; do
for MYCHROM in `seq 1 22`; do
MYFREQVCF=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/min$myfreq/chr${MYCHROM}.vcf
NLINES=$( zcat /mnt/scratch/fabrizio/LDLD/above95/${myannotation_ref}.1000g/repl1/min${myfreq}/chr${MYCHROM}.vcf.gz | grep -v '#' | wc -l )
NLINES_SOURCE=$( zcat ${MYFREQVCF}.temp.gz | grep -v '#' |wc -l )
MYPROB=$( awk -v var1=$NLINES -v var2=$NLINES_SOURCE -v var3=$FOLDS 'BEGIN{res=var3*var1/var2; if (res<=1){print res} else {print 1}}' ) #if I want more lines modify this
echo $MYPROB
for myrepl in `seq 2 5`; do #myrepl=11
echo 'myrepl is ' $myrepl
destination_folder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/repl${myrepl}; mkdir $destination_folder
destination_folder=$destination_folder/min${myfreq}; mkdir $destination_folder
zcat ${MYFREQVCF}.temp.gz | awk -v MYPROB=$MYPROB -v seed=$RANDOM 'BEGIN{srand(seed);}{if (substr($1,1,1)=="#"){print} else { if (rand()<MYPROB){print}}}' | bgzip -f > ${destination_folder}/chr${MYCHROM}.vcf.gz
done;done;done
#-----------------------------^
}

#5b)
{
preparationfolder=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/
cd $preparationfolder
myannotation="introns"
myannotation="intergenic"
myannotation="coding.exons"
myannotation_ref="coding.exons"
for myfreq in 1 5; do
for MYCHROM in `seq 10 22`; do
echo $MYCHROM
MYFREQVCF=${myannotation}/min$myfreq/chr${MYCHROM}.vcf.temp.gz
NLINES=$( zcat ${myannotation_ref}/min$myfreq/chr${MYCHROM}.vcf.temp.gz | grep -v '#' | wc -l )
NLINES_SOURCE=$( zcat ${MYFREQVCF} | grep -v '#' |wc -l )
MYNSUBSAMPLES=$( awk -v var1=$NLINES -v var2=$NLINES_SOURCE 'BEGIN{print int(var2/var1)}' )
zcat $MYFREQVCF | grep -v '#' | awk -v OFS='\t' -v nsubsamples=$MYNSUBSAMPLES -v sizesubsamples=$NLINES 'BEGIN{for (i=1;i<=(nsubsamples+1);i++){ar[i]=0;}}{mypop=int(rand()*nsubsamples)+1;
while (ar[mypop]>=sizesubsamples) {mypop=mypop+1; if (mypop>(nsubsamples+1)){mypop=1};}; ar[mypop]=ar[mypop]+1;print $1,$2-1,$2,mypop};' | bgzip -f > ${myannotation}/min$myfreq/chr${MYCHROM}.subsamples.bed.gz
done;done
}
#6) prepare for LDLD (put in destination folder and reshuffle again
{
preparationfolder=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/
cd $preparationfolder
myannotation="intergenic"
myannotation="introns"
myannotation="coding.exons" #myrepl=2
for myrepl in `seq 6 10`; do
for myfreq in 1 5; do
for MYCHROM in `seq 10 22`; do
echo 'repl: ' $myrepl ' ; chrom: ' $MYCHROM
mysim=$myrepl; if [ $myannotation == "coding.exons" ]; then mysim=1;fi
MYFREQVCF=${myannotation}/min$myfreq/chr${MYCHROM}.vcf.temp.gz
MYBEDFILE=${myannotation}/min$myfreq/chr${MYCHROM}.subsamples.bed.gz
destination_folder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g
mkdir -p $destination_folder/repl${myrepl}/min${myfreq}/reschr1
bedtools intersect -a $MYFREQVCF -b <(zcat $MYBEDFILE | awk -v OFS='\t' -v repl=$mysim '{if ($4==repl){print}}' ) -sorted -header | bgzip -f > $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz
Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/above95.unrelated.samples $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz $destination_folder/repl${myrepl}/min${myfreq}/reschr1/chr${MYCHROM}.vcf /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/nspops.txt 1
done;done;done

}



#cleaning up headers
{
myfiles=$( ls ${myfolder}/chr*vcf.gz )
for myfile in $myfiles;do
echo $myfile
zcat $myfile | awk 'BEGIN{counter=0}{if (NR==1){myinit=$1} else if ((NR>1) && (substr($1,1,1)=="#") && (myinit == $1)){counter=counter+1};if ((counter==0) || (substr($1,1,1)!="#")){print}}' | bgzip -f > ${myfile}.temp.gz
mv ${myfile}.temp.gz $myfile
done


}
#==============================MAIN======================================
#generate intergenic vcf files
{
SEED=1;COUNTER=1
while [ $COUNTER -lt 23 ]
        do
        echo $COUNTER
        #ARGUMENT NLINES: set here the final number of lines
        NLINES=$( cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/chr${COUNTER}.tab | wc -l )
        NLINES=$(( 5 * $NLINES / 4 ))
        #ARGUMENT2 DESTINATIONFILE: set here destination file
        DESTINATIONFILE=/mnt/scratch/fabrizio/LDLD/above95/intergenic/chr${COUNTER}.intergenic.seed${SEED}.vcf
        #ARGUMENT3 HEADERFILE:
        myvcfsource=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${COUNTER}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
        #ARGUMENT4 SOURCEFILE: set here vcf file to subsample
        SOURCEFILE=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/filtered/above95/chr${COUNTER}.intergenic.vcf.gz
        zcat $myvcfsource | head -500 | grep '##' > $DESTINATIONFILE
        zcat $SOURCEFILE | head -500 | grep '#CHROM' >> $DESTINATIONFILE
        shuf -n $NLINES <(zcat ${SOURCEFILE} | grep -v '#' ) | sort -Vu -k1,1 -k2,2 >> $DESTINATIONFILE
        bgzip -f $DESTINATIONFILE
        COUNTER=$(( $COUNTER + 1 ))
        done

        #reorder
        nohup Rscript ~/Dropbox/LDLD/scripts/vcf2ordered_vcf.R ~/workspace/1000genomes/above95.unrelated.samples /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/intergenic intergenic.seed1.vcf &

        #reshuffle
        COUNTER=1
while [ $COUNTER -lt 3 ]
        do        
        mkdir /mnt/scratch/fabrizio/LDLD/above95/intergenic/reschr${COUNTER}
        cp /mnt/scratch/fabrizio/LDLD/above95/intergenic/*gz /mnt/scratch/fabrizio/LDLD/above95/intergenic/reschr${COUNTER} 
        nohup Rscript ~/Dropbox/LDLD/scripts/vcf2shuffled_vcf.R ~/workspace/1000genomes/above95.unrelated.samples /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/intergenic/reschr${COUNTER} intergenic.seed1.vcf &
        COUNTER=$(( $COUNTER + 1 ))
done
        
}

#verify that ordering of population is correct
{
#verify that genotypes are the same for the same individuals
{
#ARGUMENTS
#mysample: sample to compare; myvcf: destination vcf to compare; myvcfsource: original vcf;noheader: if header 0, or not 1; name_output: flag that indicate which files are;

#check intergenic
COUNTER=1
mysample=NA20504;
myvcf=/mnt/scratch/fabrizio/LDLD/above95/intergenic/chr${COUNTER}.intergenic.seed1.vcf.gz 
myvcfsource=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${COUNTER}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz 
noheader=0
name_output=intergenic.seed1
#check coding
COUNTER=1
mysample=NA20504;
myvcf=/mnt/scratch/fabrizio/LDLD/above95/coding1000g/chr1.tab
myvcfsource=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${COUNTER}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz 
noheader=1
name_output=coding1000g


if [ $noheader -gt 0 ];then
zcat $myvcfsource | head -500 | grep '##' > myheader.vcf
cat ${myvcf} >> myheader.vcf
bgzip -f myheader.vcf
myvcf=myheader.vcf.gz
fi
vcf-subset -c ${mysample} $myvcf | awk -v OFS='\t' '{print $1,$2-1,$2,$10}' | grep -v '#' | head -20 > comparisonfile.bed
bedtools intersect -a comparisonfile.bed -b <( vcf-subset -c ${mysample} <(tabix $myvcfsource -R comparisonfile.bed -h )) -sorted -wo | awk '{print $1,$2,$3,$4,$14}' > ~/Dropbox/LDLD/input_checks/${name_output}.${mysample}.chr${COUNTER}.bed
}
#verify that positions of individuals have been changed correctly
{
myvcfsource=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz 
#NA20504 is the third Tuscany individual (population 1) according to ~/workspace/1000genomes/above95.unrelated.samples
zcat /mnt/scratch/fabrizio/LDLD/above95/intergenic/chr3.intergenic.seed1.vcf.gz | grep CHROM | awk '{print $(9+3)}' #NA20504
cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/chr3.tab | grep CHROM | awk '{print $(9+3)}' #NA20504
zcat $myvcfsource | head -500 | grep CHROM | awk '{for (i=10;i<20;i++){print $i}}' #HG00096 ...

zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/coding.exons/chr20.dest.vcf.gz | grep CHROM | awk '{print $(9+3)}' | head -1
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/coding.exons/reschr1/chr20.dest.vcf.gz | grep CHROM | awk '{print $(9+3)}' | head -1
zcat $myvcfsource | head -500 | grep CHROM | awk '{for (i=10;i<20;i++){print $i}}' #HG00096 ...

}

}

}
#various additionals for further analyses
{
#creating merged vcf for intergenic to overlap with CG
{
for MYCHROM in `seq 1 22`;do
echo $MYCHROM
zcat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl1/min5/chr${MYCHROM}.vcf.gz | head -1000 | grep '#' > /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/chr${MYCHROM}.vcf
zcat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl*/min5/chr${MYCHROM}.vcf.gz | head -1000 | grep -v '#' | sort -k2,2n >> /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/chr${MYCHROM}.vcf
bgzip -f /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/chr${MYCHROM}.vcf
done
}

#generate global mutational load file for 1000g
{
#crazy snps with all hets: chr1: rs75454623
gcc LDLD_filter12.c -o filterbyfreq_multi_mutload.out -lgmp -lm
for i in `seq 1 22`;do
echo $i
zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/all.1000g.1pop.txt 0 1 > ~/Dropbox/LDLD/1000genomes/mutload/chr${i}.tab
bcftools view -Q 0.001:minor /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/all.1000g.1pop.txt 0 1 > ~/Dropbox/LDLD/1000genomes/mutload/chr${i}_maxdoubletons.tab
done

cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
myannotation="coding.exons"
for MYCHROM in `seq 1 22`;do
echo $MYCHROM
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -R ${myannotation}/myfilter.${MYCHROM}.bed -h | gzip -f > temp.gz
bcftools view -Q 0.001:minor temp.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out ~/Dropbox/LDLD/scripts/all.1000g.1pop.txt 0 1 > ~/Dropbox/LDLD/1000genomes/mutload/chr${MYCHROM}_maxdoubletons_exons.tab
done


#just to check ordering of samples
cat ~/Dropbox/LDLD/1000genomes/mutload/chr${i}_maxdoubletons.tab | grep -v '##' | awk '{if (NF>20){print}}' |
awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$2304,$2305,$2306,$2307,$2308,$2309,$2310,$2311,$2312,$2313,$2314,$2315,$2316,$2317,$2318,$2319,$2320,$2321,$2322,$2323,$2324,$2325,$2326,$2327,$2328,$2329,$2330,$2331,$2332,$2333,$2334,$2335,$2336,$2337,$2338,$2339,$2340,$2341,$2342,$2343,$2344,$2345,$2346,$2347,$2348,$2349,$2350,$2351,$2352,$2353,$2354,$2355,$2356,$2357,$2358,$2359,$2360,$2361,$2362,$2363,$2364,$2365,$2366,$2367,$2368,$2369,$2370,$2371,$2372,$2373,$2374,$2375,$2376,$2377,$2378,$2379,
$2380,$2382,$2383,$2384,$2385,$2386,$2387,$2388,$2389,$2390,$2391,$2392,$2393,$2394,$2395,$2396,$2397,$2398,$2399,$2400,$2401,$2402,$2403,$2404,$2405,$2406,$2407,$2408,$2409,$2410,$511,$512,$513,$514,$515,$516,$517,$518,$519,$520,$521,$522,$523,$524,$525,$526,$527,$528,$529,$530,$531,$532,$533,$534,$574,$575,$572,$573,$570,$571,$568,$569,$566,$567,$564,$565,$562,$563,$560,$561,$558,$559,$556,$557,$554,$555,$576,$577,$578,$579,$580,$581,$582,$583,$584,$585,$586,$587,$588,$589,$590,$591,$592,$593,$594,$595,
$596,$597,$598,$599,$600,$601,$602,$603,$604,$605,$606,$607,$608,$609,$610,$611,$612,$613,$614,$615,$616,$617,$618,$619,$620,$621,$622,$623,$848,$849,$850,$851,$852,$853,$854,$855,$856,$857,$858,$859,$860,$237,$349,$238,$239,$275,$276,$277,$278,$307,$308,$309,$331,$310,$311,$332,$333,$312,$313,$314,$315,$334,$335,$336,$337,$338,$339,$340,$341,$342,$343,$344,$345,$346,$364,$365,$350,$351,$352,$412,$413,$395,$396,$366,$367,$355,$356,$357,$358,$353,$354,$360,$361,$359,$362,$363,$414,$415,$368,$369,$370,$371,
$390,$391,$397,$398,$392,$403,$404,$393,$394,$405,$406,$399,$400,$401,$402,$407,$408,$409,$410,$411,$432,$433,$434,$435,$436,$437,$438,$439,$440,$441,$470,$471,$472,$473,$474,$475,$476,$477,$478,$479,$480,$1164,$1165,$1166,$1167,$1168,$1169,$1170,$1171,$1172,$1173,$1254,$1255,$1256,$1257,$1258,$1259,$1339,$1340,$955,$956,$957,$958,$993,$994,$995,$996,$997,$998,$999,$1002,$1003,$1004,$1005,$1006,$1007,$1008,$1009,$1016,$1017,$1018,$1019,$1020,$1022,$1023,$1024,$1025,$1026,$1027,$1028,$1029,$1030,$1042,
$1043,$1044,$1045,$1046,$1060,$1061,$1062,$1063,$1064,$1065,$1075,$1076,$1077,$1080,$1081,$1082,$1083,$1096,$1097,$1098,$1099,$1100,$1101,$1102,$1103,$1104,$1105,$1106,$1107,$1108,$1109,$1110,$1111,$1112,$1113,$1114,$1115,$1116,$1117,$1118,$1119,$1120,$1121,$1122,$1123,$1124,$1125,$1126,$1127,$1128,$1129,$1130,$1131,$1132,$1152,$1153,$1784,$1785,$1786,$1787,$1788,$1789,$1790,$1791,$1792,$1793,$1794,$1795,$1796,$1797,$1798,$1799,$1800,$1801,$1802,$1803,$1804,$1805,$1806,$1807,$1808,$1809,$1810,$1811,
$1812,$1813,$1814,$1815,$1816,$1817,$1818,$1819,$1820,$1821,$1822,$1823,$1824,$1825,$1826,$1827,$1828,$1829,$1830,$1831,$1832,$1833,$1834,$1835,$1836,$1837,$1838,$1839,$1840,$1841,$1842,$1843,$1844,$1845,$1846,$1847,$1848,$1849,$1850,$1851,$1852,$1853,$1854,$1855,$1856,$1857,$1858,$1859,$1860,$1861,$1862,$1863,$1864,$1865,$1866,$1867,$1868,$1869,$1870,$1871,$1872,$1873,$1874,$1875,$1876,$1877,$1878,$1879,$1880,$1881,$1882,$1883,$1884,$1885,$1886,$1916,$1917,$1918,$1919,$1920,$1921,$1922,$1923,$1924,
$1925,$1926,$1927,$1928,$1929,$1930,$1931,$1932,$1933,$1934,$1935,$1936,$1937,$1938,$1939,$1940,$1941,$1942,$1943,$1944,$1945,$1946,$1947,$1948,$1949,$1950,$1951,$1952,$1953,$1954,$1955,$1956,$1957,$1958,$1959,$1960,$1961,$1963,$1964,$1965,$1966,$1968,$1969,$1970,$1971,$1972,$1973,$1974,$1975,$1976,$1977,$1978,$1979,$1980,$1981,$1982,$1983,$1984,$1985,$2004,$2005,$2006,$2007,$2008,$2009,$2010,$2011,$2012,$2013,$2014,$2015,$2016,$2017,$2018,$2019,$2020,$2021,$2022,$2023,$2024,$2025,$2026,$2027,$2028,
$2029,$2030,$2031,$2032,$2033,$2034,$2035,$2036,$2037,$195,$196,$197,$198,$199,$200,$201,$202,$203,$204,$205,$206,$207,$208,$209,$210,$211,$212,$213,$214,$215,$216,$217,$218,$219,$220,$221,$222,$223,$224,$225,$226,$228,$229,$230,$231,$232,$233,$234,$235,$236,$240,$241,$242,$243,$244,$245,$246,$248,$250,$251,$252,$253,$254,$255,$256,$257,$258,$259,$260,$261,$262,$263,$264,$265,$266,$267,$268,$269,$270,$271,$272,$273,$274,$279,$280,$281,$282,$285,$286,$287,$288,$289,$290,$291,$292,$293,$294,$295,$296,$297,
$298,$299,$300,$301,$302,$303,$304,$305,$306,$65,$66,$67,$68,$69,$70,$71,$72,$73,$74,$75,$76,$77,$78,$79,$80,$81,$113,$114,$115,$116,$117,$118,$119,$120,$121,$122,$123,$124,$125,$126,$127,$128,$129,$130,$131,$132,$133,$134,$135,$136,$137,$138,$139,$140,$141,$142,$143,$144,$145,$146,$147,$148,$149,$150,$151,$152,$153,$154,$155,$156,$157,$158,$159,$160,$161,$162,$163,$164,$165,$166,$167,$168,$169,$170,$171,$172,$173,$174,$175,$176,$177,$178,$179,$180,$181,$182,$183,$184,$185,$186,$187,$188,$189,$190,$191,
$192,$193,$194,$683,$684,$685,$686,$687,$722,$688,$753,$689,$690,$693,$695,$694,$752,$737,$696,$697,$723,$738,$739,$740,$741,$751,$749,$750,$771,$772,$773,$774,$798,$901,$802,$803,$804,$822,$823,$864,$865,$893,$894,$889,$890,$877,$878,$879,$891,$880,$895,$896,$897,$898,$899,$903,$902,$900,$904,$905,$943,$944,$946,$947,$948,$949,$950,$951,$967,$952,$953,$954,$959,$960,$961,$962,$963,$964,$965,$966,$978,$972,$973,$974,$975,$976,$977,$983,$984,$985,$986,$987,$988,$989,$990,$991,$992,$1000,$1001,$1765,$1766,
$1767,$1768,$1769,$1770,$1771,$1772,$1773,$2045,$2046,$1893,$1894,$1774,$1775,$1777,$1889,$1778,$1779,$1780,$1781,$1782,$1783,$1895,$1896,$1887,$1897,$1898,$1899,$1900,$1901,$1902,$1903,$1888,$1890,$1904,$1891,$1892,$1906,$1907,$1908,$1909,$1910,$1911,$1912,$1913,$2076,$2077,$1914,$1915,$2073,$2039,$2040,$2041,$2044,$2057,$2058,$2071,$2079,$2069,$2070,$2080,$2083,$2084,$2081,$2082,$2067,$2068,$2089,$2087,$2088,$2049,$2052,$2053,$2059,$2065,$2064,$2060,$2061,$2062,$2063,$2047,$2048,$2096,$2097,$2050,
$2051,$2056,$2042,$2043,$2086,$2074,$2075,$2090,$2091,$2092,$2093,$2094,$2095,$548,$549,$550,$551,$552,$553,$649,$650,$651,$652,$653,$654,$655,$656,$657,$658,$659,$660,$661,$662,$663,$664,$665,$666,$667,$668,$669,$670,$671,$672,$673,$674,$675,$676,$677,$678,$679,$680,$681,$682,$754,$755,$756,$757,$758,$759,$761,$762,$763,$764,$765,$766,$767,$768,$769,$770,$775,$776,$777,$778,$779,$781,$782,$783,$784,$785,$786,$787,$788,$789,$790,$791,$792,$793,$794,$795,$805,$806,$807,$808,$809,$810,$811,$812,$813,$814,
$815,$816,$817,$818,$819,$820,$821,$979,$980,$981,$982,$1381,$1380,$1378,$1379,$1377,$1391,$1392,$1390,$1393,$1394,$1395,$1555,$1396,$1397,$1516,$1517,$1398,$1399,$1400,$1401,$1402,$1575,$1403,$1404,$1405,$1406,$1407,$1414,$1429,$1430,$1515,$1427,$1428,$1432,$1442,$1433,$1434,$1435,$1519,$1436,$1439,$1440,$1441,$1437,$1438,$1443,$1487,$1490,$1593,$1489,$1581,$1488,$1493,$1494,$1492,$1495,$1499,$1554,$1498,$1497,$1496,$1594,$1518,$1571,$1522,$1520,$1521,$1523,$1525,$1526,$1524,$1527,$1549,$1550,$1551,
$1552,$1557,$1556,$1553,$1573,$1574,$1570,$1572,$1576,$1580,$1592,$1595,$1597,$1596,$1598,$1617,$1618,$1653,$1607,$1615,$1616,$1661,$1662}' 

for i in `seq 1 22`;do
#MYFILE=~/Dropbox/LDLD/1000genomes/mutload/chr${i}_maxdoubletons.tab;MYDESTINATIONFILE=~/Dropbox/LDLD/1000genomes/mutload/maxdoubletons.gw.tsv
MYFILE=~/Dropbox/LDLD/1000genomes/mutload/chr${i}_maxdoubletons_exons.tab; MYDESTINATIONFILE=~/Dropbox/LDLD/1000genomes/mutload/maxdoubletons.exons.tsv
cat $MYFILE | grep -v '#' | awk '{if (NF>20){print}}' |
awk -v OFS='\t' '{print $2304,$2305,$2306,$2307,$2308,$2309,$2310,$2311,$2312,$2313,$2314,$2315,$2316,$2317,$2318,$2319,$2320,$2321,$2322,$2323,$2324,$2325,$2326,$2327,$2328,$2329,$2330,$2331,$2332,$2333,$2334,$2335,$2336,$2337,$2338,$2339,$2340,$2341,$2342,$2343,$2344,$2345,$2346,$2347,$2348,$2349,$2350,$2351,$2352,$2353,$2354,$2355,$2356,$2357,$2358,$2359,$2360,$2361,$2362,$2363,$2364,$2365,$2366,$2367,$2368,$2369,$2370,$2371,$2372,$2373,$2374,$2375,$2376,$2377,$2378,$2379,
$2380,$2382,$2383,$2384,$2385,$2386,$2387,$2388,$2389,$2390,$2391,$2392,$2393,$2394,$2395,$2396,$2397,$2398,$2399,$2400,$2401,$2402,$2403,$2404,$2405,$2406,$2407,$2408,$2409,$2410,$511,$512,$513,$514,$515,$516,$517,$518,$519,$520,$521,$522,$523,$524,$525,$526,$527,$528,$529,$530,$531,$532,$533,$534,$574,$575,$572,$573,$570,$571,$568,$569,$566,$567,$564,$565,$562,$563,$560,$561,$558,$559,$556,$557,$554,$555,$576,$577,$578,$579,$580,$581,$582,$583,$584,$585,$586,$587,$588,$589,$590,$591,$592,$593,$594,$595,
$596,$597,$598,$599,$600,$601,$602,$603,$604,$605,$606,$607,$608,$609,$610,$611,$612,$613,$614,$615,$616,$617,$618,$619,$620,$621,$622,$623,$848,$849,$850,$851,$852,$853,$854,$855,$856,$857,$858,$859,$860,$237,$349,$238,$239,$275,$276,$277,$278,$307,$308,$309,$331,$310,$311,$332,$333,$312,$313,$314,$315,$334,$335,$336,$337,$338,$339,$340,$341,$342,$343,$344,$345,$346,$364,$365,$350,$351,$352,$412,$413,$395,$396,$366,$367,$355,$356,$357,$358,$353,$354,$360,$361,$359,$362,$363,$414,$415,$368,$369,$370,$371,
$390,$391,$397,$398,$392,$403,$404,$393,$394,$405,$406,$399,$400,$401,$402,$407,$408,$409,$410,$411,$432,$433,$434,$435,$436,$437,$438,$439,$440,$441,$470,$471,$472,$473,$474,$475,$476,$477,$478,$479,$480,$1164,$1165,$1166,$1167,$1168,$1169,$1170,$1171,$1172,$1173,$1254,$1255,$1256,$1257,$1258,$1259,$1339,$1340,$955,$956,$957,$958,$993,$994,$995,$996,$997,$998,$999,$1002,$1003,$1004,$1005,$1006,$1007,$1008,$1009,$1016,$1017,$1018,$1019,$1020,$1022,$1023,$1024,$1025,$1026,$1027,$1028,$1029,$1030,$1042,
$1043,$1044,$1045,$1046,$1060,$1061,$1062,$1063,$1064,$1065,$1075,$1076,$1077,$1080,$1081,$1082,$1083,$1096,$1097,$1098,$1099,$1100,$1101,$1102,$1103,$1104,$1105,$1106,$1107,$1108,$1109,$1110,$1111,$1112,$1113,$1114,$1115,$1116,$1117,$1118,$1119,$1120,$1121,$1122,$1123,$1124,$1125,$1126,$1127,$1128,$1129,$1130,$1131,$1132,$1152,$1153,$1784,$1785,$1786,$1787,$1788,$1789,$1790,$1791,$1792,$1793,$1794,$1795,$1796,$1797,$1798,$1799,$1800,$1801,$1802,$1803,$1804,$1805,$1806,$1807,$1808,$1809,$1810,$1811,
$1812,$1813,$1814,$1815,$1816,$1817,$1818,$1819,$1820,$1821,$1822,$1823,$1824,$1825,$1826,$1827,$1828,$1829,$1830,$1831,$1832,$1833,$1834,$1835,$1836,$1837,$1838,$1839,$1840,$1841,$1842,$1843,$1844,$1845,$1846,$1847,$1848,$1849,$1850,$1851,$1852,$1853,$1854,$1855,$1856,$1857,$1858,$1859,$1860,$1861,$1862,$1863,$1864,$1865,$1866,$1867,$1868,$1869,$1870,$1871,$1872,$1873,$1874,$1875,$1876,$1877,$1878,$1879,$1880,$1881,$1882,$1883,$1884,$1885,$1886,$1916,$1917,$1918,$1919,$1920,$1921,$1922,$1923,$1924,
$1925,$1926,$1927,$1928,$1929,$1930,$1931,$1932,$1933,$1934,$1935,$1936,$1937,$1938,$1939,$1940,$1941,$1942,$1943,$1944,$1945,$1946,$1947,$1948,$1949,$1950,$1951,$1952,$1953,$1954,$1955,$1956,$1957,$1958,$1959,$1960,$1961,$1963,$1964,$1965,$1966,$1968,$1969,$1970,$1971,$1972,$1973,$1974,$1975,$1976,$1977,$1978,$1979,$1980,$1981,$1982,$1983,$1984,$1985,$2004,$2005,$2006,$2007,$2008,$2009,$2010,$2011,$2012,$2013,$2014,$2015,$2016,$2017,$2018,$2019,$2020,$2021,$2022,$2023,$2024,$2025,$2026,$2027,$2028,
$2029,$2030,$2031,$2032,$2033,$2034,$2035,$2036,$2037,$195,$196,$197,$198,$199,$200,$201,$202,$203,$204,$205,$206,$207,$208,$209,$210,$211,$212,$213,$214,$215,$216,$217,$218,$219,$220,$221,$222,$223,$224,$225,$226,$228,$229,$230,$231,$232,$233,$234,$235,$236,$240,$241,$242,$243,$244,$245,$246,$248,$250,$251,$252,$253,$254,$255,$256,$257,$258,$259,$260,$261,$262,$263,$264,$265,$266,$267,$268,$269,$270,$271,$272,$273,$274,$279,$280,$281,$282,$285,$286,$287,$288,$289,$290,$291,$292,$293,$294,$295,$296,$297,
$298,$299,$300,$301,$302,$303,$304,$305,$306,$65,$66,$67,$68,$69,$70,$71,$72,$73,$74,$75,$76,$77,$78,$79,$80,$81,$113,$114,$115,$116,$117,$118,$119,$120,$121,$122,$123,$124,$125,$126,$127,$128,$129,$130,$131,$132,$133,$134,$135,$136,$137,$138,$139,$140,$141,$142,$143,$144,$145,$146,$147,$148,$149,$150,$151,$152,$153,$154,$155,$156,$157,$158,$159,$160,$161,$162,$163,$164,$165,$166,$167,$168,$169,$170,$171,$172,$173,$174,$175,$176,$177,$178,$179,$180,$181,$182,$183,$184,$185,$186,$187,$188,$189,$190,$191,
$192,$193,$194,$683,$684,$685,$686,$687,$722,$688,$753,$689,$690,$693,$695,$694,$752,$737,$696,$697,$723,$738,$739,$740,$741,$751,$749,$750,$771,$772,$773,$774,$798,$901,$802,$803,$804,$822,$823,$864,$865,$893,$894,$889,$890,$877,$878,$879,$891,$880,$895,$896,$897,$898,$899,$903,$902,$900,$904,$905,$943,$944,$946,$947,$948,$949,$950,$951,$967,$952,$953,$954,$959,$960,$961,$962,$963,$964,$965,$966,$978,$972,$973,$974,$975,$976,$977,$983,$984,$985,$986,$987,$988,$989,$990,$991,$992,$1000,$1001,$1765,$1766,
$1767,$1768,$1769,$1770,$1771,$1772,$1773,$2045,$2046,$1893,$1894,$1774,$1775,$1777,$1889,$1778,$1779,$1780,$1781,$1782,$1783,$1895,$1896,$1887,$1897,$1898,$1899,$1900,$1901,$1902,$1903,$1888,$1890,$1904,$1891,$1892,$1906,$1907,$1908,$1909,$1910,$1911,$1912,$1913,$2076,$2077,$1914,$1915,$2073,$2039,$2040,$2041,$2044,$2057,$2058,$2071,$2079,$2069,$2070,$2080,$2083,$2084,$2081,$2082,$2067,$2068,$2089,$2087,$2088,$2049,$2052,$2053,$2059,$2065,$2064,$2060,$2061,$2062,$2063,$2047,$2048,$2096,$2097,$2050,
$2051,$2056,$2042,$2043,$2086,$2074,$2075,$2090,$2091,$2092,$2093,$2094,$2095,$548,$549,$550,$551,$552,$553,$649,$650,$651,$652,$653,$654,$655,$656,$657,$658,$659,$660,$661,$662,$663,$664,$665,$666,$667,$668,$669,$670,$671,$672,$673,$674,$675,$676,$677,$678,$679,$680,$681,$682,$754,$755,$756,$757,$758,$759,$761,$762,$763,$764,$765,$766,$767,$768,$769,$770,$775,$776,$777,$778,$779,$781,$782,$783,$784,$785,$786,$787,$788,$789,$790,$791,$792,$793,$794,$795,$805,$806,$807,$808,$809,$810,$811,$812,$813,$814,
$815,$816,$817,$818,$819,$820,$821,$979,$980,$981,$982,$1381,$1380,$1378,$1379,$1377,$1391,$1392,$1390,$1393,$1394,$1395,$1555,$1396,$1397,$1516,$1517,$1398,$1399,$1400,$1401,$1402,$1575,$1403,$1404,$1405,$1406,$1407,$1414,$1429,$1430,$1515,$1427,$1428,$1432,$1442,$1433,$1434,$1435,$1519,$1436,$1439,$1440,$1441,$1437,$1438,$1443,$1487,$1490,$1593,$1489,$1581,$1488,$1493,$1494,$1492,$1495,$1499,$1554,$1498,$1497,$1496,$1594,$1518,$1571,$1522,$1520,$1521,$1523,$1525,$1526,$1524,$1527,$1549,$1550,$1551,
$1552,$1557,$1556,$1553,$1573,$1574,$1570,$1572,$1576,$1580,$1592,$1595,$1597,$1596,$1598,$1617,$1618,$1653,$1607,$1615,$1616,$1661,$1662}' >> $MYDESTINATIONFILE
done

}


}

}
#run LDLD
{
#1000g
{
#intergenic
myfolder=/mnt/scratch/fabrizio/LDLD/above95/intergenic
nohup ./ICLD11multi.sh ${myfolder} 2 intergenic.seed1 1 12 1 23 nspops.txt &
nohup ./ICLD11multi.sh ${myfolder}/reschr1 2 intergenic.seed1 1 12 1 23 nspops.txt &
nohup ./ICLD11multi.sh ${myfolder}/reschr2 2 intergenic.seed1 1 12 1 23 nspops.txt &
/mnt/scratch/fabrizio/LDLD/postproc_split.sh ${myfolder}/reschr1
/mnt/scratch/fabrizio/LDLD/postproc_split.sh ${myfolder}/reschr2
/mnt/scratch/fabrizio/LDLD/postproc_split.sh ${myfolder}
bgzip ${myfolder}/anal/sorted.res 
bgzip ${myfolder}/reschr1/anal/sorted.res
bgzip ${myfolder}/new3/reschr2/anal/sorted.res

#coding
myfolder=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/new3
    gdriveLDLD=0B3oqoLi0PJOlbUZOVEp2UURpcHc
nohup ./ICLD11multi.sh ${myfolder} 2 coding 1 12 1 23 nspops.txt &
nohup ./ICLD11multi.sh ${myfolder}/reschr1 0 coding 1 12 1 23 nspops.txt &
nohup ./ICLD11multi.sh ${myfolder}/reschr2 0 coding 1 12 1 23 nspops.txt &
#coding with new filtering (5% overall)
myannotation=coding.exons.1000g
myannotation=introns.1000g
myannotation=intergenic.1000g
cd ~/Dropbox/LDLD/scripts
#NB: intergenic min1 myrepl=2 -> 0.00000000000001
#NB: coding min5 I copied repl1 into repl3 and repl2 into repl4, then rerun LDLD. So that I will know if results change a bit because of LDLD or filtering. It could be that filtering used to be looser, so proportionally less significant stuff.
#also I should understand why not dim(mydata) multiple than 12. about 50 in pop 8 and 9 are missing in the others...
#for myrepl in `seq 6 7`;do
myrepl=2
for myrepl in `seq 3 4`;do
myfreq=5;reschr=1
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.000001 0.0${myfreq} &
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}/reschr${reschr}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.000001 0.0${myfreq} &
myfreq=1
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}
#nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000000000001 0.0${myfreq} &
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.000000000000001 0.0${myfreq} &
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}/reschr${reschr}
#nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000000000001 0.0${myfreq} &
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.000000000000001 0.0${myfreq} &
done



cp repl3
#18-08-17
myannotation=intergenic.1000g
cd ~/Dropbox/LDLD/scripts
myrepl=6;reschr=1;
for myfreq in 1 5;do
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000001 0.0${myfreq} &
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}/reschr${reschr}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000001 0.0${myfreq} &
done
myrepl=7;reschr=1;
for myfreq in 1 5;do
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.00000001 0.0${myfreq} &
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}/reschr${reschr}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.00000001 0.0${myfreq} &
done
myrepl=8;reschr=1;myfreq=1
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000000001 0.0${myfreq} &
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}/reschr${reschr}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000000001 0.0${myfreq} &
myrepl=9;reschr=1;myfreq=1
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000000000001 0.0${myfreq} &
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}/reschr${reschr}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000000000001 0.0${myfreq} &


cd ~/Dropbox/LDLD/scripts
myannotation=introns.1000g
myrepl=11;reschr=1;
myfreq=5
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.00000001 0.0${myfreq} &
myfreq=1
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/repl${myrepl}/min${myfreq}
nohup ./ICLD11amulti.sh ${myfolder} 4 $myannotation 1 12 1 23 nspops.txt 0.0000000000001 0.0${myfreq} &

done

}

#revisions GBE
{
#simulaitons
#since filtered before, I can just run with threshold 5% and it would be the same in all of them
for MYREPL in 1 2;do
for INS in 50 100 200;do
MYFOLDERS=$( ls /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims | grep nsample${INS} )
for IFOLDER in $MYFOLDERS;do
nohup ./ICLD11multic.sh /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/${IFOLDER}/repl${MYREPL} 4 ${IFOLDER} 1 1 23 nspops${INS}.txt 0.00000001 0.01 &
nohup ./ICLD11multic.sh /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/${IFOLDER}/repl${MYREPL}/reschr1 4 ${IFOLDER} 1 1 23 nspops${INS}.txt 0.00000001 0.01 &
done;done;done

    for MYREPL in 1 2;do
    for INS in 50 100 200;do
    MYFOLDERS=$( ls /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims | grep nsample${INS} | grep sims_mu )
    for IFOLDER in $MYFOLDERS;do
    nohup ./ICLD11multic.sh /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/${IFOLDER}/repl${MYREPL}/ 4 ${IFOLDER} 1 1 23 nspops${INS}.txt 0.00000001 0.01 &
    nohup ./ICLD11multic.sh /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/${IFOLDER}/repl${MYREPL}/reschr1/ 4 ${IFOLDER} 1 1 23 nspops${INS}.txt 0.00000001 0.01 &
    done;done;done

cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims
MYFOLDERS=$( ls | grep sims )
for i in $MYFOLDERS;do
rm $i/repl*/ICLD.sims_* &
rm $i/repl*/reschr1/ICLD.sims_* &
rm $i/repl*/*res &
rm $i/repl*/*sum &
rm $i/repl*/*log &
rm $i/repl*/*tab &
rm $i/repl*/reschr1/*res &
rm $i/repl*/reschr1/*sum &
rm $i/repl*/reschr1/*log &
rm $i/repl*/reschr1/*tab &
done
    
    
}
}

#organize folders
{
for myannotation in 'intergenic' 'introns';do
for myrepl in `seq 1 5`;do
for myfreq in 1 5;do
mkdir ~/Dropbox/LDLD/ipynb/figs/${myannotation}.1000g/repl${myrepl}/
mkdir ~/Dropbox/LDLD/ipynb/figs/${myannotation}.1000g/repl${myrepl}/min${myfreq}
done;done;done
}

#merge more replicates
{
myannotation=intergenic.1000g;
myannotation=introns.1000g; 
myfreq=1
myfreq=5
destination_folder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}; cd $destination_folder
mkdir -p ${destination_folder}/merged/min${myfreq}/anal
mkdir ${destination_folder}/merged/min${myfreq}/logs
mkdir -p ${destination_folder}/merged/min${myfreq}/reschr1/anal
mkdir ${destination_folder}/merged/min${myfreq}/reschr1/logs
myfiles=$( find ${destination_folder} -name sorted.res.gz | grep -v reschr1 | grep min${myfreq} )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | bgzip -f > ${destination_folder}/merged/min${myfreq}/anal/sorted.res.gz
myfiles=$( find ${destination_folder} -name sorted.res.gz | grep reschr1 | grep min${myfreq} )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | bgzip -f > ${destination_folder}/merged/min${myfreq}/reschr1/anal/sorted.res.gz
for mypop in `seq 0 11`;do 
echo 'mypop is: ' $mypop
myfiles=$( find ${destination_folder} -name all.pop${mypop}.log.gz | grep -v merged | grep -v reschr1 | grep min${myfreq} )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | awk '{if (NR==1){myNF=NF} else if (NF==myNF){print}}' | bgzip -f > ${destination_folder}/merged/min${myfreq}/logs/all.pop${mypop}.log.gz
done
for mypop in `seq 0 11`;do 
echo 'mypop is: ' $mypop
myfiles=$( find ${destination_folder} -name all.pop${mypop}.log.gz | grep -v merged | grep reschr1 | grep min${myfreq} )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | awk '{if (NR==1){myNF=NF} else if (NF==myNF){print}}' | bgzip -f > ${destination_folder}/merged/min${myfreq}/reschr1/logs/all.pop${mypop}.log.gz
done

cp ${destination_folder}/repl1/min${myfreq}/nspops.txt ${destination_folder}/merged/min${myfreq}/
cp ${destination_folder}/repl1/min${myfreq}/reschr1/nspops.txt ${destination_folder}/merged/min${myfreq}/reschr1/
myfiles=$( find ${destination_folder} -name ncomparisons.txt | grep -v reschr1 | grep min${myfreq} )
cat $myfiles > ${destination_folder}/merged/min${myfreq}/ncomparisons.txt 
myfiles=$( find ${destination_folder} -name ncomparisons.txt | grep reschr1 | grep min${myfreq} )
cat $myfiles > ${destination_folder}/merged/min${myfreq}/reschr1/ncomparisons.txt 
myfiles=$( find ${destination_folder} -name ncomparisons_pairwise.txt | grep -v reschr1 | grep min${myfreq} )
cat $myfiles > ${destination_folder}/merged/min${myfreq}/ncomparisons_pairwise.txt 
myfiles=$( find ${destination_folder} -name all.freqlog.gz | grep -v reschr1 | grep min${myfreq} )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | bgzip -f > ${destination_folder}/merged/min${myfreq}/logs/all.freqlog.gz
myfiles=$( find ${destination_folder} -name all.freqlog.gz | grep -v reschr1 | grep min${myfreq} )
}
#--------------------------------------------------------
#--------------scan bad_snps in shell--------------------
#--------------------------------------------------------

#newer syntax
{
TARGETVCF=$1 #target vcf from which to extract biased snps
MYCHR=$2 #chromosome to process
FILERDATA=$3
MYTAG=$4 #name of the object to refer to in FILERDATA. Also used as name-tag for output
FREQFILTER=$5 # frequency filtering
PVALUE_THR=$6 #threshold for p-value (printing and glm to glmm)
DESTINATION_FOLDER=$7 # destination folder
PARSER_FILE=$8 #file in which function reorder_samples_vcf is defined to reorder samples properly
NSPOPS_FILE=$9

#15-10
MYMIN=1 #MYMIN=1
MYDEST=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min${MYMIN}
MYFILERDATA=${MYDEST}/logpop_nAB_pos20.RData
MYTAG=logpop_nAB_pos20
DESTINATION_FOLDER=${MYDEST}/bad_snps_glmm
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 13 22`;do
MYTARGETFILE=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.0${MYMIN} 0.00001 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g.sh ~/Dropbox/LDLD/scripts/nspops.txt & #for min5 0.00001
done

MYMIN=5 #MYMIN=1
MYDEST=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min${MYMIN}
MYFILERDATA=${MYDEST}/logpop_nAB_pos20.RData
MYTAG=logpop_nAB_pos20
DESTINATION_FOLDER=${MYDEST}/bad_snps_glmm
#DESTINATION_FOLDER=${MYDEST}/bad_snps_glmm_new
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 14 22`;do
MYTARGETFILE=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.0${MYMIN} 0.00001 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g.sh ~/Dropbox/LDLD/scripts/nspops.txt & #for min5 0.00001
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG
cd /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/bad_snps_glmm; cat bad_snps_logpop_nAB_pos20_glmm_all.bed | awk '{if ($8<0.000001){print}}' > bad_snps_logpop_nAB_pos20_glmm_all_p10m6.bed
cat bad_snps_logpop_nAB_pos20_glmm_all.bed | sort -g -k 8,8g | head -2000 > bad_snps_logpop_nAB_pos20_glmm_all_top.bed
cd /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/bad_snps_glmm; cat bad_snps_logpop_nAB_pos20_glmm_all.bed | awk '{if ($8<0.000001){print}}' > bad_snps_logpop_nAB_pos20_glmm_all_p10m6.bed
cat bad_snps_logpop_nAB_pos20_glmm_all.bed | sort -g -k 8,8g | head -2000 > bad_snps_logpop_nAB_pos20_glmm_all_top.bed

MYMIN=1 #MYMIN=5
MYDEST=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${MYMIN}
MYTAG=logpop_nAB_pos #logpop_nAB_pos20 logpop_nAB_pos
MYFILERDATA=${MYDEST}/${MYTAG}.RData
DESTINATION_FOLDER=${MYDEST}/bad_snps_glmm
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 8`;do
MYTARGETFILE=${MYDEST}/chr${MYCHR}.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.0${MYMIN} 0.001 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_subset1000g.sh ~/Dropbox/LDLD/scripts/nspops.txt & #for min5 0.00001
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG

#sims
for ISIM in `seq 1 4`;do
MYMIN=5 #MYMIN=1
MYDEST=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${MYMIN}
MYTAG=logpop_nAB_pos #logpop_nAB_pos20 logpop_nAB_pos
MYFILERDATA=${MYDEST}/${MYTAG}_sim${ISIM}.RData
DESTINATION_FOLDER=${MYDEST}/bad_snps_glmm_sim${ISIM}
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 5`;do
MYTARGETFILE=${MYDEST}/chr${MYCHR}.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.0${MYMIN} 0.001 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_subset1000g.sh ~/Dropbox/LDLD/scripts/nspops.txt & #for min5 0.00001
done;done


~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG



}


#this seems to suggest that the reason is logpop_nAB_pos20 file.
#I get many snps with logpop_nAB_pos20.RData file of min5, and small with min1, independently on the data and filtering I use.
#One explanations might be that with min1 I get too much linkage just by coincidence.
{
MYMIN=1 #MYMIN=1
MYDEST=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min${MYMIN}
MYFILERDATA=${MYDEST}/logpop_nAB_pos20.RData
MYTAG=logpop_nAB_pos20
DESTINATION_FOLDER=${MYDEST}/bad_snps_glmm_min5
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 11 22`;do
MYTARGETFILE=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.05 0.00001 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g.sh ~/Dropbox/LDLD/scripts/nspops.txt & #for min5 0.00001
done
MYMIN=5 #MYMIN=1
MYDEST=/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min${MYMIN}
MYFILERDATA=${MYDEST}/logpop_nAB_pos20.RData
MYTAG=logpop_nAB_pos20
DESTINATION_FOLDER=${MYDEST}/bad_snps_glmm_min1
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 3`;do
MYTARGETFILE=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.01 0.00001 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g.sh ~/Dropbox/LDLD/scripts/nspops.txt & #for min5 0.00001
done

#check difference in the two logpop_nAB_pos20 files.
load("/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/logpop_nAB_pos20.RData")
logpop_nAB_pos20_min5<-logpop_nAB_pos20
load("/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/logpop_nAB_pos20.RData")
logpop_nAB_pos20_min1<-logpop_nAB_pos20
range(unlist(logpop_nAB_pos20_min1)) # 21.5 6664.5
range(unlist(logpop_nAB_pos20_min5)) # 70.5 643.0
#ok, in min1 there are a few outliers with huge differences. Very unnatural: understand if with the same code I get the same to exclude that artifact.
} 


#one of the reasons for which I get less with min1 is because for those the first threshold is hard to overcome (as I know, pvalue glmm all are harder to achieve).
#but glmm all should be at the very least the same if not more.
MYMIN=5 #MYMIN=5
MYFILERDATA=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${MYMIN}/logpop_nAB_pos.RData
MYTAG=logpop_nAB_pos
DESTINATION_FOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${MYMIN}/bad_snps_min${MYMIN}_glmm
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 14 22`;do
MYTARGETFILE=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${MYMIN}/chr${MYCHR}.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.0${MYMIN} 0.05 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_subset1000g.sh &
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG

#low p-value (proxy with at least two pop sign) -> increase a tiny bit
MYFILERDATA=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/logpop_nAB_pos01.RData
MYTAG=logpop_nAB_pos01
DESTINATION_FOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/bad_snps_min5_glmm
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 22`;do
MYTARGETFILE=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.05 0.1 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_subset1000g.sh &
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG
}

#why lower below fdr?
#-predict hets (unlikely because difference was small). 
#I see in old code that I was doing some wrong filtering, apparently calculating the allele frequency by just looking at homoz! I should get same effect by having higher threshold for freq!! -> actually lower: maybe I should use exactly the same filtering
#-use log nAB - same results -> it is not!!!
#-sign only when in at least two populations together! -> cleaner signal -> then it should also work with lower pvalue: in fact it increases a bit, but not much.
#ref-alt?
#maybe I should run old code on new nAB!
#try to find also old nAB files!


#higher allele freq thr
MYFILERDATA=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/logpop_nAB_pos.RData
MYTAG=logpop_nAB_pos
DESTINATION_FOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/bad_snps_min5to10_glmm_higherf
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 22`;do
MYTARGETFILE=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm_log.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.1 0.1 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_subset1000g.sh &
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG

#logs
MYFILERDATA=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/logpop_nAB_pos.RData
MYTAG=logpop_nAB_pos
DESTINATION_FOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/bad_snps_min5_glmm_log
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 22`;do
MYTARGETFILE=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm_log.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.05 0.1 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_subset1000g.sh &
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG

#low p-value (proxy with at least two pop sign) -> increase a tiny bit
MYFILERDATA=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/logpop_nAB_pos01.RData
MYTAG=logpop_nAB_pos01
DESTINATION_FOLDER=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/bad_snps_min5_glmm
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 22`;do
MYTARGETFILE=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $MYTARGETFILE $MYCHR $MYFILERDATA $MYTAG 0.05 0.1 $DESTINATION_FOLDER ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_subset1000g.sh &
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG






for i in `seq 1 22`;do
bedtools intersect -a <( cat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/functional_annotation/filtered/temp/missense.chr${i}.temp | awk -v OFS='\t' '{print $1,$2-1,$2}' ) -b ./repl2/min5/removed1_pos_shell_scale.bed
done

non syn 76
syn 25
#still returning errors, but it does not stop anymore!



#intergenic
cd ~/Dropbox/dropbox_tablet/programming/plotting/
cp /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min1/empfdr_pos_hist.pdf figa.pdf
cp /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min5/empfdr_pos_hist.pdf figb.pdf
./add_legend.sh figa.pdf '(a)'
./add_legend.sh figb.pdf '(b)'
pdfjam '(a)_figa.pdf' '(b)_figb.pdf' --nup 2x1 --landscape --outfile output_ab.pdf --
mv output_ab.pdf /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/empfdr_pos_hist_min1min5.pdf
#coding
cd ~/Dropbox/dropbox_tablet/programming/plotting/
cp /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min1/empfdr_pos_hist.pdf figa.pdf
cp /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min5/empfdr_pos_hist.pdf figb.pdf
./add_legend.sh figa.pdf '(a)'
./add_legend.sh figb.pdf '(b)'
pdfjam '(a)_figa.pdf' '(b)_figb.pdf' --nup 2x1 --landscape --outfile output_ab.pdf --
mv output_ab.pdf /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/empfdr_pos_hist_min1min5.pdf


#intergenic
cd ~/Dropbox/dropbox_tablet/programming/plotting/
cp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/anal/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab_sharingCG.pdf figa.pdf
cp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/anal/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab_sharingCG.pdf figb.pdf
cp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min1/anal/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab_density.pdf figc.pdf
cp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/anal/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.bed.tab_density.pdf figd.pdf
./add_legend.sh figa.pdf '(a)'
./add_legend.sh figb.pdf '(b)'
./add_legend.sh figc.pdf '(c)'
./add_legend.sh figd.pdf '(d)'
pdfjam '(a)_figa.pdf' '(b)_figb.pdf' '(c)_figc.pdf' '(d)_figd.pdf' --nup 2x2 --landscape --outfile output_ab.pdf --
mv output_ab.pdf /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/byfreq_summary_bad_snps_logpop_nAB_pos20_glmm_all.pdf
#coding
cd ~/Dropbox/dropbox_tablet/programming/plotting/
cp /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/anal/byfreq_summary_bad_snps_logpop_nAB_pos_glmm_all.bed.tab_sharingCG.pdf figa.pdf
cp /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal/byfreq_summary_bad_snps_logpop_nAB_pos_glmm_all.bed.tab_sharingCG.pdf figb.pdf
cp /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/anal/byfreq_summary_bad_snps_logpop_nAB_pos_glmm_all.bed.tab_density.pdf figc.pdf
cp /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal/byfreq_summary_bad_snps_logpop_nAB_pos_glmm_all.bed.tab_density.pdf figd.pdf
./add_legend.sh figa.pdf '(a)'
./add_legend.sh figb.pdf '(b)'
./add_legend.sh figc.pdf '(c)'
./add_legend.sh figd.pdf '(d)'
pdfjam '(a)_figa.pdf' '(b)_figb.pdf' '(c)_figc.pdf' '(d)_figd.pdf' --nup 2x2 --landscape --outfile output_ab.pdf --
mv output_ab.pdf /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/byfreq_summary_bad_snps_logpop_nAB_pos_glmm_all.pdf




