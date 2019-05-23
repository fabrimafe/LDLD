#======================LDLD shell calculations==================
#preparation
{
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
pdf("~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/march2017/repl1/file.sizes.tab.pdf")
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
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbQual=cds&fbDownBases=200&hgta_doGetBed=get+BED" | gzip -f > coding.exons.bed.gz
#introns
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbQual=intron&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED" | gzip -f > introns.bed.gz
#whole genes
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED" | gzip -f > genes.bed.gz

#there is wide overlap between these categories.
bedtools subtract -a introns.bed.gz -b coding.exons.bed.gz | gzip -f > introns.bed
bedtools subtract -a genes.bed.gz -b introns.bed.gz | gzip -f > genes.bed


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
gzip -f ${myannotation}/chr${MYCHROM}.vcf
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
        gzip -f $DESTINATIONFILE
        zcat genes.bed   | awk -v OFS='\t' -v MYCHROM=$MYCHROM '{if ($1=="chr"MYCHROM){print MYCHROM,$2-5000,$3+5000}}' | awk '{if ($2<0){$2=0};print}' | sort -Vu -k1,1 -k2,2n > ${myannotation}/myfilter.${MYCHROM}.bed
        bedtools merge -i ${myannotation}/myfilter.${MYCHROM}.bed > mytemp.${MYCHROM}.bed; mv mytemp.${MYCHROM}.bed ${myannotation}/myfilter.${MYCHROM}.bed
        bedtools subtract -a ${DESTINATIONFILE}.gz -b ${myannotation}/myfilter.${MYCHROM}.bed -sorted -header | gzip -f > /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/intergenic/chr${MYCHROM}.vcf.gz
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
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out nspops.txt 0.0${myfreq} 12 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | gzip -f > ${MYVCF}
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
zcat ${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out nspops.txt 0.0${myfreq} 12 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | gzip -f > ${MYFREQVCF}.temp.gz
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
zcat ${MYFREQVCF}.temp.gz | awk -v MYPROB=$MYPROB -v seed=$RANDOM 'BEGIN{srand(seed);}{if (substr($1,1,1)=="#"){print} else { if (rand()<MYPROB){print}}}' | gzip -f > ${destination_folder}/chr${MYCHROM}.vcf.gz
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
for MYCHROM in `seq 1 9`; do
echo $MYCHROM
MYFREQVCF=${myannotation}/min$myfreq/chr${MYCHROM}.vcf.temp.gz
NLINES=$( zcat ${myannotation_ref}/min$myfreq/chr${MYCHROM}.vcf.temp.gz | grep -v '#' | wc -l )
NLINES_SOURCE=$( zcat ${MYFREQVCF} | grep -v '#' |wc -l )
MYNSUBSAMPLES=$( awk -v var1=$NLINES -v var2=$NLINES_SOURCE 'BEGIN{print int(var2/var1)}' )
zcat $MYFREQVCF | grep -v '#' | awk -v OFS='\t' -v nsubsamples=$MYNSUBSAMPLES -v sizesubsamples=$NLINES 'BEGIN{for (i=1;i<=(nsubsamples+1);i++){ar[i]=0;}}{mypop=int(rand()*nsubsamples)+1;
while (ar[mypop]>=sizesubsamples) {mypop=mypop+1; if (mypop>(nsubsamples+1)){mypop=1};}; ar[mypop]=ar[mypop]+1;print $1,$2-1,$2,mypop};' | gzip -f > ${myannotation}/min$myfreq/chr${MYCHROM}.subsamples.bed.gz
done;done
}
#6) prepare for LDLD (put in destination folder and reshuffle again
{
preparationfolder=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/
cd $preparationfolder
myannotation="intergenic"
myannotation="introns"
myannotation="coding.exons" #myrepl=2
for myrepl in `seq 1 10`; do
for myfreq in 1 5; do
for MYCHROM in `seq 1 22`; do
echo 'repl: ' $myrepl ' ; chrom: ' $MYCHROM
mysim=$myrepl; if [ $myannotation == "coding.exons" ]; then mysim=1;fi
MYFREQVCF=${myannotation}/min$myfreq/chr${MYCHROM}.vcf.temp.gz
MYBEDFILE=${myannotation}/min$myfreq/chr${MYCHROM}.subsamples.bed.gz
destination_folder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g
mkdir $destination_folder/repl${myrepl}
mkdir $destination_folder/repl${myrepl}/min${myfreq}
mkdir $destination_folder/repl${myrepl}/min${myfreq}/reschr1
bedtools intersect -a $MYFREQVCF -b <(zcat $MYBEDFILE | awk -v OFS='\t' -v repl=$mysim '{if ($4==repl){print}}' ) -sorted -header | gzip -f > $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz
Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/above95.unrelated.samples $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz $destination_folder/repl${myrepl}/min${myfreq}/reschr1/chr${MYCHROM}.vcf /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min5/nspops.txt 1
done;done;done

}



#cleaning up headers
{
myfiles=$( ls ${myfolder}/chr*vcf.gz )
for myfile in $myfiles;do
echo $myfile
zcat $myfile | awk 'BEGIN{counter=0}{if (NR==1){myinit=$1} else if ((NR>1) && (substr($1,1,1)=="#") && (myinit == $1)){counter=counter+1};if ((counter==0) || (substr($1,1,1)!="#")){print}}' | gzip -f > ${myfile}.temp.gz
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
        gzip -f $DESTINATIONFILE
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
gzip -f myheader.vcf
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
gzip ${myfolder}/anal/sorted.res 
gzip ${myfolder}/reschr1/anal/sorted.res
gzip ${myfolder}/new3/reschr2/anal/sorted.res

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
#for myrepl in `seq 6 7`;do
myrepl=2
for myrepl in `seq 1 3`;do
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
}

#organize folders
{
for myannotation in 'intergenic' 'introns';do
for myrepl in `seq 1 5`;do
for myfreq in 1 5;do
mkdir ~/Dropbox/LDLD/ipynb/figs/${myannotation}.1000g/march2017/repl${myrepl}/
mkdir ~/Dropbox/LDLD/ipynb/figs/${myannotation}.1000g/march2017/repl${myrepl}/min${myfreq}
done;done;done
}

#merge more intergenic
{
myannotation=intergenic.1000g; myfreq=5
destination_folder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}
mkdir ${destination_folder}/merged
mkdir ${destination_folder}/merged/min${myfreq}/
mkdir ${destination_folder}/merged/min${myfreq}/anal
mkdir ${destination_folder}/merged/min${myfreq}/logs
mkdir ${destination_folder}/merged/min${myfreq}/reschr1
mkdir ${destination_folder}/merged/min${myfreq}/reschr1/anal
mkdir ${destination_folder}/merged/min${myfreq}/reschr1/logs
myfiles=$( find ${destination_folder} -name sorted.res.gz | grep -v reschr1 | grep min5 )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | gzip -f > ${destination_folder}/merged/min${myfreq}/anal/sorted.res.gz
myfiles=$( find ${destination_folder} -name sorted.res.gz | grep reschr1 | grep min5 )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | gzip -f > ${destination_folder}/merged/min${myfreq}/reschr1/anal/sorted.res.gz
for mypop in `seq 1 2`;do 
echo 'mypop is: ' $mypop
myfiles=$( find ${destination_folder} -name all.pop${mypop}.log.gz | grep -v reschr1 | grep min5 )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | gzip -f > ${destination_folder}/merged/min${myfreq}/logs/all.pop${mypop}.log.gz
done
for mypop in `seq 7 11`;do 
echo 'mypop is: ' $mypop
myfiles=$( find ${destination_folder} -name all.pop${mypop}.log.gz | grep reschr1 | grep min5 )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | gzip -f > ${destination_folder}/merged/min${myfreq}/reschr1/logs/all.pop${mypop}.log.gz
done

cp ${destination_folder}/repl1/min${myfreq}/nspops.txt ${destination_folder}/merged/min${myfreq}/
cp ${destination_folder}/repl1/min${myfreq}/reschr1/nspops.txt ${destination_folder}/merged/min${myfreq}/reschr1/
myfiles=$( find ${destination_folder} -name ncomparisons.txt | grep -v reschr1 | grep min5 )
cat $myfiles > ${destination_folder}/merged/min${myfreq}/ncomparisons.txt 
myfiles=$( find ${destination_folder} -name ncomparisons.txt | grep reschr1 | grep min5 )
cat $myfiles > ${destination_folder}/merged/min${myfreq}/reschr1/ncomparisons.txt 
myfiles=$( find ${destination_folder} -name ncomparisons_pairwise.txt | grep -v reschr1 | grep min5 )
cat $myfiles > ${destination_folder}/merged/min${myfreq}/ncomparisons_pairwise.txt 
myfiles=$( find ${destination_folder} -name all.freqlog.gz | grep -v reschr1 | grep min5 )
zcat $myfiles | sort -k1,1n -k2,2n -k3,3n -k4,4n | gzip -f > ${destination_folder}/merged/min${myfreq}/logs/all.freqlog.gz






}

}
