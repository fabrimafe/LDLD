#==========reanalyses of data starting on 15-3-2017================
#======================LDLD shell calculations==================
{
#tests ICLDmulti
{
cd ~/Dropbox/LDLD/scripts
gcc ICLDmulti10.c -o ICLDmulti.out -lgmp -lm #compile
date
./ICLDmulti.out ~/Dropbox/LDLD/ICLDtests/chr21.tab ~/Dropbox/LDLD/ICLDtests/chr22.tab ~/Dropbox/LDLD/ICLDtests/nspops.txt ~/Dropbox/LDLD/ICLDtests/file4.res ~/Dropbox/LDLD/ICLDtests/file5.log ~/Dropbox/LDLD/ICLDtests/file6.sum ~/Dropbox/LDLD/ICLDtests/file7.minilog 0.05 12 0 #example command
date
./ICLDmulti.out ~/Dropbox/LDLD/ICLDtests/chr21.tab ~/Dropbox/LDLD/ICLDtests/chr22.tab ~/Dropbox/LDLD/ICLDtests/nspops.txt ~/Dropbox/LDLD/ICLDtests/file4.res ~/Dropbox/LDLD/ICLDtests/file5.log ~/Dropbox/LDLD/ICLDtests/file6.sum ~/Dropbox/LDLD/ICLDtests/file7.minilog 0.05 12 1 #example command
date

}

#faster reorder
{
MYCHROM=1
zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | head -100 | grep CHROM
}

#generate input files
{
#early steps
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
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbQual=cds&fbDownBases=200&hgta_doGetBed=get+BED" | gzip -f > coding.exons.bed
#introns
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbQual=intron&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED" | gzip -f > introns.bed
#whole genes
curl -s "http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=598155399_qCNZgv1Bf2iYCisjoeZlaa82NaLS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED" | gzip -f > genes.bed

myannotation="coding.exons"
myannotation="introns"
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
for MYCHROM in `seq 1 22`; do
echo $MYCHROM
zcat ${myannotation}.bed  | awk -v OFS='\t' -v MYCHROM=$MYCHROM '{if ($1=="chr"MYCHROM){print MYCHROM,$2,$3}}' |sort -Vu -k1,1 -k2,2n > ${myannotation}/myfilter.${MYCHROM}.bed
bedtools merge -i ${myannotation}/myfilter.${MYCHROM}.bed > mytemp.${MYCHROM}.bed; mv mytemp.${MYCHROM}.bed ${myannotation}/myfilter.${MYCHROM}.bed
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -H -R myfilter.${myannotation}/myfilter.${MYCHROM}.bed  > chr$MYCHROM.vcf
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -R ${myannotation}/myfilter.${MYCHROM}.bed  | sort -Vu -k1,1 -k2,2n >> chr${MYCHROM}.vcf
gzip -f chr${MYCHROM}.vcf
echo $MYCHROM 'reorder...' 
nohup Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/above95.unrelated.samples chr${MYCHROM}.vcf.gz chr${MYCHROM}.dest.vcf /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/new4/min5/nspops.txt 0
for reschr in `seq 1 1`;do
    mkdir ${myannotation}/reschr${reschr}
    nohup Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/above95.unrelated.samples chr${MYCHROM}.vcf.gz ${myannotation}/reschr${reschr}/chr${MYCHROM}.dest.vcf /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/new4/min5/nspops.txt 1
done
#rm ${myannotation}/myfilter.${MYCHROM}.bed
done
mv chr*.dest.vcf.gz  ${myannotation}

#zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHROM}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | wc -l #1103797
#zcat chr${MYCHROM}.dest.vcf.gz | wc -l #612981
#bedtools intersect -a temp.${MYCHROM}.vcf -b myfilter.${MYCHROM}.bed -sorted | wc -l
#bedtools intersect -a temp.${MYCHROM}.vcf -b myfilter.${MYCHROM}.bed -sorted -v | wc -l


#multiple flat frequency threshold
{
myannotation="coding.exons"
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
for MYCHROM in `seq 1 22`; do
echo $MYCHROM
bcftools view -q 0.05:minor ${myannotation}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/chr${MYCHROM}.min5.vcf.gz
bcftools view -q 0.01:minor -Q 0.05:minor ${myannotation}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/chr${MYCHROM}.min1max5.vcf.gz
bcftools view -q 0.005:minor -Q 0.05:minor ${myannotation}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/chr${MYCHROM}.min05max5.vcf.gz
bcftools view -Q 0.05:minor ${myannotation}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/chr${MYCHROM}.max5.vcf.gz
bcftools view -Q 0.01:minor ${myannotation}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/chr${MYCHROM}.max1.vcf.gz
mv ${myannotation}/chr${MYCHROM}.min5.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min5/
mv ${myannotation}/chr${MYCHROM}.min1max5.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1max5/
done

#reschr
for MYCHROM in `seq 11 22`; do
for reschr in 1 2; do
echo $MYCHROM
bcftools view -q 0.05:minor ${myannotation}/reschr${reschr}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/reschr${reschr}/chr${MYCHROM}.min5.vcf.gz
bcftools view -q 0.01:minor -Q 0.05:minor ${myannotation}/reschr${reschr}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/reschr${reschr}/chr${MYCHROM}.min1max5.vcf.gz
mkdir /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min5/reschr${reschr}/
mkdir /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1max5/reschr${reschr}/
mv ${myannotation}/reschr${reschr}/chr${MYCHROM}.min5.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min5/reschr${reschr}/
mv ${myannotation}/reschr${reschr}/chr${MYCHROM}.min1max5.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1max5/reschr${reschr}/
done;done
}

#flat frequency threshold 0.01, multiple locals
{
myannotation="coding.exons"
myannotation="introns"
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
mkdir /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1
mkdir /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min5
for reschr in 1 2; do
mkdir /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1/reschr${reschr}/
mkdir /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min5/reschr${reschr}/
done
for MYCHROM in `seq 1 22`; do
echo $MYCHROM
bcftools view -q 0.01:minor ${myannotation}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/chr${MYCHROM}.min1.vcf.gz
mv ${myannotation}/chr${MYCHROM}.min1.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1/
for reschr in 1 2; do
bcftools view -q 0.01:minor ${myannotation}/reschr${reschr}/chr${MYCHROM}.dest.vcf.gz | gzip -f > ${myannotation}/reschr${reschr}/chr${MYCHROM}.min1.vcf.gz
mv ${myannotation}/reschr${reschr}/chr${MYCHROM}.min1.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1/reschr${reschr}/
done;done

cp /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1/*.* /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min5/
for reschr in 1 2; do
cp /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min1/reschr${reschr}/*.* /mnt/scratch/fabrizio/LDLD/above95/${myannotation}.1000g/new4/min5/reschr${reschr}/
done

}


#reorder vcfs
if (FALSE)
{
targetsamples<-system(paste0("cat ",myorderedsamples_file),intern=TRUE) #previously called samplesabove95
sourcesamples<-system(paste0("cat ",mysourcesamples_file," | head -300 | grep '#CHROM' | head -1 | awk '{for (i=10; i<=NF; i++) print $i}'"),intern=TRUE) #previously called samples1000g
reordered_samples<-paste0(intersect(targetsamples,sourcesamples),collapse=',') #previously called samples_order1000g_above95
#nspops<-as.numeric(system("cat /mnt/scratch/fabrizio/LDLD/nspops.txt",intern=TRUE))

order_samples_vcf.gz<-function(samples_order,folder,annotation,chrfrom=1,chrto=22)
{
#samples_order: column with ordered samples: e.g. samples_order
#
        for (chr in chrfrom:chrto)
        {
        #system(paste0("vcf-subset -c ",samples_order," ",folder,"/chr",chr,".",annotation,".gz | awk '{if (NR>2){print}}' | grep -v '##' > " ,folder,"/chr",chr,".",annotation))
        system(paste0("vcf-subset -c ",samples_order," ",folder,"/chr",chr,".",annotation,".gz > " ,folder,"/chr",chr,".",annotation))
        system(paste0("gzip -f ",folder,"/chr",chr,".",annotation))
        }
}

order_samples_vcf.gz(reordered_samples,mydestination_folder,annotation,chrfrom=chrstart,chrto=chrend)

reshuffle_withinpop<-function(l_samples,l_npops)
{
res<-c()
for (i in 1:length(l_npops))
        {
        if (i==1){init<-1; end<-l_npops[1]} else {init<-init+l_npops[i-1];end<-end+l_npops[i]}
        res<-c(res,sample(l_samples[init:end],replace=FALSE))
        }
return(res)
}

}

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
#myannotation=introns.1000g
cd ~/Dropbox/LDLD/scripts
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/new4/min5
nohup ./ICLD11amulti.sh ${myfolder} 2 min1 1 12 1 23 nspops.txt 1 0.05 &
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/new4/min1
nohup ./ICLD11amulti.sh ${myfolder} 2 min1 1 12 1 23 nspops.txt 1 0.01&
for reschr in 1 2;do
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/new4/min5/reschr${reschr}
nohup ./ICLD11amulti.sh ${myfolder} 2 min1 1 12 1 23 nspops.txt 1 0.05 &
myfolder=/mnt/scratch/fabrizio/LDLD/above95/${myannotation}/new4/min1/reschr${reschr}
nohup ./ICLD11amulti.sh ${myfolder} 2 min1 1 12 1 23 nspops.txt 1 0.01 &
done

#/mnt/scratch/fabrizio/LDLD/postproc2_split.sh is identical but keeps hits that are unique in only one populations
/mnt/scratch/fabrizio/LDLD/postproc_split.sh ${myfolder}/reschr1
/mnt/scratch/fabrizio/LDLD/postproc_split.sh ${myfolder}/reschr2
/mnt/scratch/fabrizio/LDLD/postproc_split.sh ${myfolder}
/mnt/scratch/fabrizio/LDLD/postproc_split.sh 
gzip ${myfolder}/anal/sorted.res 
gzip ${myfolder}/reschr1/anal/sorted.res
gzip ${myfolder}/new3/reschr2/anal/sorted.res
}
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
#====================R initial variables=======================
{
if (FALSE) {
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic"
myfolderdropbox<-"~/Dropbox/LDLD/ipynb/figs/intergenic/"
#system(paste0("cat ",myfolder,"/*minilog | grep ncomparisons | awk '{print $14}' > ",myfolder,"/ncomparisons.txt"))
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
computed_exact_pvalues<-FALSE
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
npopulations<-12
PREPROCESSING<-FALSE
using_sign_for_nAB<-FALSE
}


if (FALSE) {
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new3"
myfolderdropbox<-"~/Dropbox/LDLD/ipynb/figs/coding1000g/march2017"
#system(paste0("cat ",myfolder,"/*minilog | grep ncomparisons | awk '{print $14}' > ",myfolder,"/ncomparisons.txt"))
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
computed_exact_pvalues<-FALSE
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
npopulations<-12
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE
}

if (TRUE) {
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/new4/min1"
myfolderdropbox<-"~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/march2017/new4/min1"
#system(paste0("cat ",myfolder,"/*minilog | grep ncomparisons | awk '{print $14}' > ",myfolder,"/ncomparisons.txt"))
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
computed_exact_pvalues<-TRUE
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
npopulations<-12
PREPROCESSING<-TRUE
using_sign_for_nAB<-TRUE
nreschr=2
}


if (TRUE) {
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/new4/min5"
myfolderdropbox<-"~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/march2017/new4/min5"
#system(paste0("cat ",myfolder,"/*minilog | grep ncomparisons | awk '{print $14}' > ",myfolder,"/ncomparisons.txt"))
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
computed_exact_pvalues<-TRUE
popnames<-c("TSI","IBS","PUR","GWD","CHB","JPT","CHS","FIN","ACB","YRI","KHV","STU")
npopulations<-12
PREPROCESSING<-TRUE
using_sign_for_nAB<-TRUE
nreschr=2
}

setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
options(scipen=999)
setwd("/mnt/scratch/fabrizio/LDLD")
library(data.table)


}
#===================SIMULATE FALSE POSITIVES==================================================
#mytab<-fread(paste0("cat ",myfolder,"/reschr1/chr",i,".tab | sed 's/|//g' "))

#======================LOAD_INITIAL_FILES====================================================
#TAG 'LOAD_INITIAL_FILES'
if (TRUE)
{
if (PREPROCESSING) {
    print("aggregating res files")
    system(paste0("mkdir ",myfolder,"/anal/"))
    system(paste0("~/Dropbox/LDLD/scripts/postproc_split.sh ",myfolder))
    print("aggregating log files for each pops")
    system(paste0("mkdir ",myfolder,"/logs"))
    for (ichrA in 1:21){ for (ichrB in (ichrA+1):22){for (ipop in 0:(npopulations-1)){system(paste0("cat ",myfolder,"/chr",ichrA,".tabchr",ichrB,".tab.log | awk -v OFS='\t' '{print ",ichrA,",",ichrB,",$0 }'> ",myfolder,"/logs/chr",ichrA,".",ichrB,".pop",ipop,".log")) }}}
    print("aggregating all log files")
    for (i in 0:(npopulations-1)) {system(paste0("cat ",myfolder,"/logs/chr*.pop",i,".log | awk '{if ($NF==",i,"){print}}' > ",myfolder,"/logs/all.pop",i,".log"))}
    for (myperm in 1:nreschr){print(paste("permutation: ",myperm))
    print("aggregating res files for permutations")
    system(paste0("mkdir ",myfolder,"/reschr",myperm,"/logs"))
    system(paste0("mkdir ",myfolder,"/reschr",myperm,"/anal"))
    system(paste0("~/Dropbox/LDLD/scripts/postproc_split.sh ",myfolder,"/reschr",myperm))
    print("aggregating log files for permutations")
    for (ichrA in 1:21){ for (ichrB in (ichrA+1):22){ for (ipop in 0:(npopulations-1)){print(paste("ichrA: ",ichrA," pop:",ipop)); system(paste0("cat ",myfolder,"/reschr",myperm,"/chr",ichrA,".tabchr",ichrB,".tab.log | awk -v OFS='\t' '{print ",ichrA,",",ichrB,",$0 }'> ",myfolder,"/reschr",myperm,"/logs/chr",ichrA,".",ichrB,".pop",ipop,".log")) }}}
    for (i in 0:(npopulations-1)) {system(paste0("cat ",myfolder,"/reschr",myperm,"/logs/chr*.pop",i,".log | awk '{if ($NF==",i,"){print}}' > ",myfolder,"/reschr",myperm,"/logs/all.pop",i,".log"))}
    }
system(paste0("cat ",myfolder,"/*minilog | grep ncomparisons | awk '{print $14}' > ",myfolder,"/ncomparisons.txt"))
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
dataFIN<-subset(mydata,mydata$pop==index_focal_pop)
ncomp<-read.table(paste0(myfolder,"/ncomparisons.txt"),header=FALSE)
ncomp<-sum(as.numeric(ncomp$V1))
combined_pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)
combined_fdr<-p.adjust(p=combined_pval, method = "fdr", n = ncomp)
length(combined_fdr)
dim(dataFIN)
dim(mydata)
mydata$pFisherpos[mydata$pFisherpos==0]<-10^(-21)#min(mydata$pFisherpos[mydata$pFisherpos!=0])
mydata$pFisherneg[mydata$pFisherneg==0]<-10^(-21)#min(mydata$pFisherneg[mydata$pFisherneg!=0])
mydata$pFisherpos[mydata$pFisherpos>1]<-1#min(mydata$pFisherpos[mydata$pFisherpos!=0])
mydata$pFisherneg[mydata$pFisherneg>1]<-1#min(mydata$pFisherneg[mydata$pFisherneg!=0])

#combined pvalues for positive and negative
min(temp$pFisherpos)
min(temp$pFisherpos)
fishersmethod<-function(z) if ( length(z)>=2) {return(sumlog(z))} else return(z)
temp<-mydata[mydata$nA/mydata$Nse>=0.01 & mydata$nB/mydata$Nse>=0.01,]
pairid<-temp[, paste(chrA,chrB,posA,posB,sep=".")]
combined_pval_pos<-aggregate(temp$pFisherpos,by=list(pairid),FUN=function(x) fishersmethod(x))
combined_pval_neg<-aggregate(temp$pFisherneg,by=list(pairid),FUN=function(x) fishersmethod(x))
combined_pval<-aggregate(temp$T2,by=list(pairid),FUN=function(x) fishersmethod(x))


fullsign<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
save(dataFIN,file=paste0(myfolder,"/dataFIN.RData"))
save(fullsign,file=paste0(myfolder,"/fullsign.RData"))
save(combined_pval,file=paste0(myfolder,"/combined_pval.RData"))
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
write.table(signAbed,file=paste0(myfolder,"/anal/snpsA.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file=paste0(myfolder,"/anal/snpsB.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);

#explorative plots to check consistency of p-values
#custom function to make color transparent
t_col <- function(color, percent = 50, name = NULL)  #color = color name, percent = % transparency,name = an optional name for the color
{
rgb.val <- col2rgb(color); # Get RGB values for named color
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100-percent)*255/100,names = name)  # Make new color using input color as base and alpha set by transparency # Save the color
invisible(t.col)  
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
if (FALSE)
{
load(paste0(myfolder,"/combined_pval.RData"))
load(paste0(myfolder,"/fullsign.RData"))
combined_pval_sign<-combined_pval
combined_fdr_sign<-fullsign$combined_fdr
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
nperm<-2
reschr_count<-rep(0,nperm)
combined_pval_l<-list();dataFIN_l<-list()
for (myperm in 1:nperm)
    {
    mydata<-fread(paste0("zcat ",myfolder,"/reschr",myperm,"/anal/sorted.res.gz"))
#    data<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/reschr",myperm,"/anal/sorted.res"))
    names23<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
    names20<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","T2","Xtot","X","pop","popX")
    if (computed_exact_pvalues) { mynames<-names23 } else { mynames<-names20 }
    #names(data)<-mynames #if data.frame
    setnames(mydata,mynames)
    #if (dim(mydata)[2]==20){names(mydata)<-names20} else if (dim(mydata)[2]==23){names(mydata)<-names23}
    system(paste0("cat ",myfolder,"/reschr",myperm,"/*minilog | grep ncomparisons | awk '{print $14}' > ",myfolder,"/reschr",myperm,"/ncomparisons.txt"))
    ncomp<-read.table(paste0(myfolder,"/reschr",myperm,"/ncomparisons.txt"),header=FALSE)
    ncomp<-sum(as.numeric(ncomp$V1))
    dataFIN<-subset(mydata,mydata$pop==index_focal_pop)
    combined_pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)    
    combined_fdr<-p.adjust(p=combined_pval, method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
    fullsign_reschr<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
    combined_pval_l[[myperm]]<-combined_pval
    dataFIN_l[[myperm]]<-dataFIN
    reschr_count[myperm]<-dim(fullsign_reschr)[1]
    }
#    save(combined_pval_l,file=paste0(myfolder,"/combined_pval_l.RData"))
#    save(dataFIN_l,file=paste0(myfolder,"/dataFIN_l.RData"))
    
empfdr<-empiricall_fdr(combined_pval_sign,unlist(combined_pval_l),nperm)
#save(empfdr,file=paste0(myfolder,"/empfdr.RData"))
load(paste0(myfolder,"/empfdr.RData"))
length(empfdr[empfdr<0.05])
plot(empfdr,type="l",xlab="index",ylab="fdr")
pdf(paste0(myfolderdropbox,"/empfdrvsT2.pdf"))
plot(combined_fdr_sign[order(combined_fdr_sign)],type="l",xlab="index",ylab="fdr",xlim=c(1,200000),col="cadetblue",lwd=2.5)
lines(empfdr,col="gold",lwd=2.5)
abline(h=0.05,col="black",lty=2,lwd=2.5)
legend(0,0.03,c("T2 fdr","T2 empirical"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
lty=c(1,1), # gives the legend appropriate symbols (lines)
lwd=c(2.5,2.5),col=c("cadetblue","gold"))
dev.off()
#NB: empirical curves obtained with #new (high res scan) and the lower resolution one are exactly the same, just that now I can explore higher fdr (that is useless anyway, apart for this plot)

pdf(paste0(myfolderdropbox,"/empfdrvsT2_hist.pdf"))
hist(log(combined_pval_sign[combined_pval_sign<0.05]+10^(-20),10),breaks=25,col="gold",main="T2 vs T2 reshuffled",xlim=c(-20,-5),xlab="log(pvalue(T2)+10^-20)")
hist(log(combined_pval_l[[1]][combined_pval_l[[1]]<0.05]+10^(-20),10),breaks=25,col="cadetblue",add=TRUE,xlim=c(-20,-5))
dev.off()

signAbed<-as.data.frame(cbind(as.character(paste0('chr',fullsign$chrA)),fullsign$posA-1,fullsign$posA),stringsAsFactors =FALSE)
signBbed<-as.data.frame(cbind(as.character(paste0('chr',fullsign$chrB)),fullsign$posB-1,fullsign$posB),stringsAsFactors =FALSE)
names(signAbed)<-c("chr","start","end")
names(signBbed)<-c("chr","start","end")
signAbed$start<-as.numeric(signAbed$start)
signAbed$end<-as.numeric(signAbed$end)
signBbed$start<-as.numeric(signBbed$start)
signBbed$end<-as.numeric(signBbed$end)
#create list of snps
write.table(signAbed,file=paste0(myfolder,"/anal/snpsA.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);
write.table(signBbed,file=paste0(myfolder,"/anal/snpsB.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE);

load(paste0(myfolder,"/dataFIN.RData"))
dataFIN_sign<-dataFIN[order(combined_pval_sign)<(sum(empfdr<0.05)+0.5),]
#save(dataFIN_sign,file=paste0(myfolder,"/dataFIN_sign.RData"))
dataFIN_top10000<-dataFIN[order(combined_pval_sign)[1:10000],]
#save(dataFIN_top10000,file=paste0(myfolder,"/dataFIN_top10000.RData"))

#calculate empirical distribution also on the basis of pKuli and prho2 (not in new3)
if (computed_exact_pvalues) {
combined_pval_sign<-combined_pval
combined_prho2pval_sign<-combined_prho2pval
combined_pKulipval_sign<-combined_pKulipval
combined_fdr_sign<-combined_fdr

combined_pvalT2_l<-list()
combined_pvalprho2_l<-list()
combined_pvalpKuli_l<-list()
#for coding low res
#{
#nperm<-4
#reschr_count<-rep(0,nperm)
#for (myperm in 2:nperm)
#}
nperm<-1
reschr_count<-rep(0,nperm)
for (myperm in 1:nperm)
    {
    data.reshuffled<-fread(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/reschr",myperm,"/anal/sorted.res"))
    names23<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
    setnames(data.reshuffled,names23)
#    system(paste0("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/*minilog | grep ncomparisons | awk '{print $14}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/reschr",myperm,"/ncomparisons.txt"))
    ncomp<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/reschr",myperm,"/ncomparisons.txt"),header=FALSE)
    ncomp<-sum(ncomp$V1)
    dataFIN<-subset(data.reshuffled,data.reshuffled$pop==7)
    combined_pval<-dchisq(dataFIN$X,df=2*dataFIN$popX)    
    combined_fdr<-p.adjust(p=combined_pval, method = "fdr", n = ncomp) #not the best approach because like this I have all pairs. This is why probably at the beginning I was using a smaller cutoff.
    fullsign_reschr<-cbind(dataFIN,combined_fdr)[combined_fdr<0.05,]
    combined_pvalT2_l[[myperm]]<-combined_pval
    reschr_count[myperm]<-dim(fullsign_reschr)[1]
    mydata<-data.reshuffled[data.reshuffled$nA>=0.05*data.reshuffled$Nse,]
    mydata<-mydata[mydata$nA<=0.95*mydata$Nse,]
    mydata<-mydata[mydata$nB>=0.05*mydata$Nse,]
    mydata<-mydata[mydata$nB<=0.95*mydata$Nse,]
    pairid<-mydata[, paste(chrA,chrB,posA,posB,sep=".")]
#    mydata[,pairid:=pairid]
    combined_prho2pval<-aggregate(mydata$prho2,by=list(pairid),FUN=function(x) 1-pchisq(-2*log(prod(x)),df=2*length(x)))
    combined_pKulipval<-aggregate(mydata$pKuli,by=list(pairid),FUN=function(x) 1-pchisq(-2*log(prod(x)),df=2*length(x)))
    combined_prho2fdr<-p.adjust(p=combined_prho2pval[,2], method = "fdr", n = ncomp) 
    combined_pKulifdr<-p.adjust(p=combined_pKulipval[,2], method = "fdr", n = ncomp) 
    combined_pvalprho2_l[[myperm]]<-combined_prho2pval
    combined_pvalpKuli_l[[myperm]]<-combined_pKulipval
    }

mypvalues<-cbind(combined_prho2pval,combined_pKulipval[,2],combined_prho2fdr,combined_pKulifdr)
sum(as.numeric(combined_prho2fdr<0.05))
sum(as.numeric(combined_pKulifdr<0.05))
#sum(as.numeric(combined_T2fdr<0.05))

#for low res
#empfdr_pKuli<-empiricall_fdr(mypvalues[,3],c(combined_pvalpKuli_l[[2]][,2],combined_pvalpKuli_l[[3]][,2],combined_pvalpKuli_l[[4]][,2]),3)
#empfdr_prho2<-empiricall_fdr(mypvalues[,2],c(combined_pvalprho2_l[[2]][,2],combined_pvalprho2_l[[3]][,2],combined_pvalprho2_l[[4]][,2]),3)
#for high res
empfdr_pKuli<-empiricall_fdr(mypvalues[,3],c(combined_pvalpKuli_l[[1]][,2]),1)
empfdr_prho2<-empiricall_fdr(mypvalues[,2],c(combined_pvalprho2_l[[1]][,2]),1)

save(empfdr_pKuli,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/empfdrpKuli.RData")
save(empfdr_prho2,file="/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/empfdrprho2.RData")

myempfdr<-c(sum(as.numeric(empfdr<0.05)),sum(as.numeric(empfdr_pKuli<0.05)),sum(as.numeric(empfdr_prho2<0.05)))
myfdr<-c(sum(as.numeric(combined_T2fdr<0.05)),sum(as.numeric(combined_pKulifdr<0.05)),sum(as.numeric(combined_prho2fdr<0.05)))
mymat<-rbind(myempfdr,myfdr)
colnames(mymat)<-c("T2","pKuli","prho2")
pdf("~/Dropbox/LDLD/ipynb/figs/coding1000g/empfdr_vs_theoreticalfdr_hist.pdf")
barplot(mymat,beside=TRUE,col=c("floralwhite","azure3"),ylab="count")
legend("topright",fill=c("floralwhite","azure3"), legend=c("empirical","theoretical"))
dev.off()
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
#----------------------------------------------------------------------------------------------------^
#=======================
#done by taking only significant links 
#=======================
#D<0 vs D>0
{
load(paste0(myfolder,"/mydata.RData"))
load(paste0(myfolder,"/dataFIN_sign.RData"))
if (!using_sign_for_nAB){load(paste0(myfolder,"/dataFIN_top10000.RData"));dataFIN_sign<-dataFIN_top10000}
#parsing the data so that I distinguish between positive and negative linkage
#notice that I could think of subdividing sites on the basis of their positive or negative linkage. Now, in that case I should calculate
#one tailed p-values (could do it in the future, not now). 
#Dstr1)Alternative is still taking T2 (since anyway I am using empirical p-values), but by subdividing into D>0 and D<0. Since most affected just a bit, take |D|>0.01 for single populations.
#Dstr2).Since this might make me lose too much power take: only those with same sign (or 0) as consensus with sum of logs of T2.
#I could also take a union of these.
#And especially consider that I am summing up over many sites, so I should not worry too much about freak cases.


myres<-dataLDLD2Dstr12_f(mydata,dataFIN_sign)
pairid<-mydata[, paste(chrA,chrB,posA,posB,sep=".")]
pairid_sign<-dataFIN_sign[, paste(chrA,chrB,posA,posB,sep=".")]
mydata<-mydata[!is.na(match(pairid,pairid_sign)),]
pairid<-mydata[, paste(chrA,chrB,posA,posB,sep=".")]
data_neg<-aggregate(mydata$D,by=list(pairid),FUN=function(x) sum(x<0))[,2]
data_pos<-aggregate(mydata$D,by=list(pairid),FUN=function(x) sum(x>0))[,2]
data_Dsum<-aggregate(mydata$D,by=list(pairid),FUN=function(x) sum(x))[,2]
data_Dsum<-aggregate(mydata$D,by=list(pairid),FUN=function(x) sum(x))
mydata<-cbind(mydata,sumlogT2=0)
mydata$sumlogT2[mydata$D<0]<- log(mydata$T2)[mydata$D<0]
mydata$sumlogT2[mydata$D>0]<- -log(mydata$T2)[mydata$D>0]
data_sumlogT2<-aggregate(mydata$sumlogT2,by=list(pairid),FUN=function(x) sum(x))
data_sumlogT2smaller0<-aggregate(mydata$sumlogT2,by=list(pairid),FUN=function(x) sum(x[x<0]+0))
data_sumlogT2greater0<-aggregate(mydata$sumlogT2,by=list(pairid),FUN=function(x) sum(x[x>0]+0))

#Dstr1:
mydata_Dstr1_neg<-myres[[1]]
mydata_Dstr1_pos<-myres[[2]]
#Dstr2:
#n pops is sum(pairid==pairid[1])
mydata_Dstr2_neg<-myres[[3]]
mydata_Dstr2_pos<-myres[[4]]
#union
#to filter unique I can set key to all columns
mydata_Dstr12neg<-myres[[5]]
mydata_Dstr12pos<-myres[[6]]

#save(mydata_Dstr12neg,file=paste0(myfolder,"/mydata_Dstr12neg.RData"))
#save(mydata_Dstr12pos,file=paste0(myfolder,"/mydata_Dstr12pos.RData"))
#save(mydata_Dstr1_pos,file=paste0(myfolder,"/mydata_Dstr1_pos.RData"))
#save(mydata_Dstr1_neg,file=paste0(myfolder,"/mydata_Dstr1_neg.RData"))


library(MASS)
require(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
mypalette<-rf(32)
#panel with patterns of linkage D>0 vs D<0
{
my.labels <- c(" (a)", " (b)", " (c)"," (d)")
my.locations <- c("topleft","topleft","topleft","topleft")
#to assign single sites to negatively or positively correlated I should do one-tailed test.
#That is certainly possible with fisher's exact test.
#Quicker option, without now having to implement it, is to take two-tailed p-value and assuming 1
#when opposite sign. This could be extrapolated exactly to one-tailed p-values, if I ever calculate
#them. Alternative is to compare sum(log(1-pvalues)), and then calculate relative likelihood to get
#a measure of how one better than the other. Or 1-pvalue*pvalue (yes cases*no cases) -> odds that one or #the other. 
histsumT2<-hist(data_sumlogT2$x,breaks=30)
histsumT2$counts<-log(histsumT2$counts)
hist1<-hist(data_neg);hist2<-hist(data_pos);histDsum<-hist(data_Dsum);
histDpos<-hist(mydata$D[mydata$D>0],breaks=30);histDneg<-hist(-mydata$D[mydata$D<0],breaks=30);
hist1log<-hist(log(subset(mydata,D<0)$T2),xlim=c(-30,0),breaks=30);
hist2log<-hist(log(subset(mydata,D>0)$T2),xlim=c(-30,0),breaks=30);
hist1rholog<-hist(log(subset(mydata,D<0)$rho2),xlim=c(-20,0),breaks=30);
hist2rholog<-hist(log(subset(mydata,D>0)$rho2),xlim=c(-20,0),breaks=30);

colneg=rgb(0,0.4,0.6,1/4)#rgb(0,0,1,1/4)
colpos=rgb(0.8,0.2,0,1/4)#rgb(1,0,0,1/4)
plot.new()
pdf(paste0(myfolderdropbox,"/hist_D0.pdf"))
par(mfrow=c(2,2))
plot(histDneg , col=colneg, main="",#main="populations with D <|> 0 per site",
xlab="|D|")  # first histogram
plot(histDpos , col=colpos, add=T)  # second
legend('topright', c("D<0","D>0"), col = c(colneg,colpos),
       border = "black",pch=15)
put.fig.letter(label=my.labels[1], location=my.locations[1], font=2)
plot(hist2log , col=colpos,xlim=c(-30,0),ylim=c(0,max(c(hist2log$counts,hist1log$counts))),main="",xlab=expression("log(T2"[i]*")"))  # first histogram
plot(hist1log , col=colneg,xlim=c(-30,0),add=T)  # second
legend('topleft', c("D<0","D>0"), col = c(colneg,colpos),
       border = "black",pch=15)
put.fig.letter(label=my.labels[2], location=my.locations[2], font=2)
plot(hist2rholog , col=colpos,xlim=c(-20,0),ylim=c(0,max(c(hist2rholog$counts,hist1rholog$counts))),main="",xlab=expression(paste("log(",rho)^2*""[i]*")"))  # first histogram
plot(hist1rholog , col=colneg,xlim=c(-20,0),add=T)  # second
legend('topleft', c("D<0","D>0"), col = c(colneg,colpos),
       border = "black",pch=15)
put.fig.letter(label=my.labels[3], location=my.locations[2], font=2)

#ggplot(mydata,aes(x=D,y=log(T2))) + stat_binhex()
#group.index <- rep(1:2, c(length(subset(mydata,D<0)$T2), length(subset(mydata,D>0)$T2)))
#sm.density.compare(c(subset(mydata,D<0)$T2,subset(mydata,D>0)$T2), group = group.index, model = "equal")
if (FALSE){ #density plots look great but have ugly kernel that make log goes beyond 0
plot.multi.dens <- function(s,mymain="",mycols=1:10,mylty=1,mylwd=2)
{
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s)) {
        junk.x = c(junk.x, density(s[[i]])$x)
        junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = mymain,col = mycols[1],lty=mylty,lwd=mylwd)
    for(i in 1:length(s)) {
        lines(density(s[[i]]), xlim = xr, ylim = yr, col = mycols[i],lty=mylty,lwd=mylwd)
    }
}

    plot.multi.dens(list(log(subset(mydata,D<0)$T2),log(subset(mydata,D>0)$T2)),mymain="log(T2) densities",mycols=c(colneg,colpos))
    legend('topleft', y = 15000, c("D<0","D>0"), col = c(colneg,colpos),
        border = "black",pch=15)
}
plot(histsumT2, main="", ylab='log(Frequency)',xlab=expression(paste("-",Sigma)*"log(T2"[i]*")"),col='gray')
put.fig.letter(label=my.labels[4], location=my.locations[4], font=2)
dev.off()
}
#figures with consistency of linkage D>0 vs D<0 across populations
{
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
mypalette<-rf(32)
mycor.test<-cor.test(data_sumlogT2smaller0$x,data_sumlogT2greater0$x,method="spearman")
#plot(data_sumlogT2smaller0$x,-data_sumlogT2greater0$x)
par(mfrow=c(1,1))
pdf(paste0(myfolderdropbox,"/npops_Dgreatersmaller0.pdf"))
plot(hist1 , col=colneg, xlim=c(0,12),main="",#main="populations with D <|> 0 per site",
xlab="n populations/site with D>0|D<0") # first histogram
plot(hist2 , col=colpos, xlim=c(0,12), add=T) # second
legend('topright', c("D<0","D>0"), col = c(colneg,colpos),
       border = "black",pch=15)
dev.off()

pdf(paste0(myfolderdropbox,"/hist_sumlogT2_D0.pdf"))
p<-ggplot(data.table(smaller=-data_sumlogT2smaller0$x,greater=data_sumlogT2greater0$x),aes(x=smaller,y=greater)) #+
h3<-p+stat_bin2d(bins=50) + scale_fill_gradientn(colours=mypalette, trans="log") +
xlab(expression(paste("-",Sigma)*"log(T2"[i]*"<0)")) +
ylab(expression(paste("-",Sigma)*"log(T2"[i]*">0)"))+
  ggtitle("concordance among population in linkage direction") +
  annotate("text", x=max(abs(data_sumlogT2smaller0$x)/2), y=max(data_sumlogT2greater0$x/2), label= paste0("spearman's rho=",mycor.test$estimate,";  p.value=",mycor.test$p.value))
h3
#stat_density2d(geom='polygon',colour='black')
#geom_density2d(aes(colour=..level..))
  #geom_density2d(aes(colour=..level..)) + 
  #scale_colour_gradient(low="green",high="red")
dev.off()
}

#Discussion with Roger:
#For nAB minor alleles go well - assumption about errors of the SC
#Calculate separatedely for positive and negative linkage and look if features are the same.
#I could even try to sum nAB of major-minor (taking both and averaging) in negatively correlated
#Notice that in per site analysis I don't care that I neglect a few sites that are positively correlated.

}
#build infoplots
{
#->log tables with only significant links
#organize_log_files(myfolder,12)
#organize_log_files(paste0(myfolder,"/reschr1"),12)
#organize_log_files(paste0(myfolder,"/reschr2"),12)

using_sign_for_nAB<-FALSE
load(paste0(myfolder,"/dataFIN_sign.RData"))
#if (!using_sign_for_nAB){load(paste0(myfolder,"/dataFIN_top10000.RData"));dataFIN_sign<-dataFIN_top10000}
load(paste0(myfolder,"/combined_pval.RData"))
load(paste0(myfolder,"/empfdr.RData"))
load(paste0(myfolder,"/dataFIN_l.RData"))
load(paste0(myfolder,"/combined_pval_l.RData"))
#calculate nAB with all minor alleles
{
mylogpop<-list();logpop_nAB<-list()
for ( ipop in 0:11)    
    {
    print(ipop)
    mydata_temp<-fread(paste0(myfolder,"/logs/all.pop",ipop,".log"))
    #mysign<-apply(dataFIN[empfdr<0.05,c(1,3,2,4),with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse="."))
    mylogpop[[ipop+1]]<-intersect_links(mydata_temp,dataFIN_sign)
    logpop_nAB[[ipop+1]]<-apply(mylogpop[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
    }
save(mylogpop,file=paste0(myfolder,"/mylogpop.RData"))
save(logpop_nAB,file=paste0(myfolder,"/logpop_nAB.RData"))

n_sign_links<-dim(dataFIN_sign)[1]
    
for ( myperm in 1:nperm){ #total nAB kept fixed
mylogpop_reschr<-list();logpop_reschr_nAB<-list()
for ( ipop in 0:11)
    {
    fraction_nlinks<-1
    mydata_temp<-fread(paste0(myfolder,"/reschr",myperm,"/logs/all.pop",ipop,".log"))
    print(ipop)
    for (i in 1:5)
        {
        sign_links_reschr<-dataFIN_l[[myperm]][order(combined_pval_l[[myperm]])<=round(n_sign_links/fraction_nlinks),c(1,2,3,4),with=FALSE]
        mylogpop_reschr[[ipop+1]]<-intersect_links(mydata_temp,sign_links_reschr)
        logpop_reschr_nAB[[ipop+1]]<-apply(mylogpop_reschr[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
        fraction_nlinks<-fraction_nlinks*(sum(logpop_reschr_nAB[[ipop+1]])/sum(logpop_nAB[[ipop+1]]))
        }
    }
save(logpop_reschr_nAB,file=paste0(myfolder,"/logpop_reschr",myperm,"_nAB.RData"))
save(mylogpop_reschr,file=paste0(myfolder,"/mylogpop_reschr",myperm,".RData"))
}
}
#calculate nAB separating D>0 and D<0
{
load(paste0(myfolder,"/mydata_Dstr12neg.RData"))
load(paste0(myfolder,"/mydata_Dstr12pos.RData"))
load(paste0(myfolder,"/mydata_Dstr1_pos.RData"))
load(paste0(myfolder,"/mydata_Dstr1_neg.RData"))

#mydata_Dstr12<-rbind(mydata_Dstr12neg,mydata_Dstr12pos)
#mydata_Dstr12<-unique(mydata_Dstr12,by=c('chrA','chrB','posA','posB')) #unique by multiple columns
#mydata_Dstr12<-mydata_Dstr12[order(chrA,chrB,posA,posB)] #sort by multiple columns

mylogpop_neg<-list();logpop_nAB_neg<-list();mylogpop_pos<-list();logpop_nAB_pos<-list();logpop_nAB<-list()
for ( ipop in 0:11)    
    {
    print(ipop)
    mydata_temp<-fread(paste0(myfolder,"/logs/all.pop",ipop,".log"))
    mysubset_neg<-mydata_Dstr12neg[mydata_Dstr12neg$pop==ipop,]
    mysubset_pos<-mydata_Dstr12pos[mydata_Dstr12pos$pop==ipop,]
    #mysign<-apply(dataFIN[empfdr<0.05,c(1,3,2,4),with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse="."))
    mylogpop_pos[[ipop+1]]<-intersect_links(mydata_temp,mysubset_pos)
    mylogpop_neg[[ipop+1]]<-intersect_links(mydata_temp,mysubset_neg)
    logpop_nAB_pos[[ipop+1]]<-apply(mylogpop_pos[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
    logpop_nAB_neg[[ipop+1]]<-apply(mylogpop_neg[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) 2-nABf(y))))
#    logpop_nAB_Dstr12[[ipop+1]]<-logpop_nAB_pos[[ipop+1]]+logpop_nAB_neg[[ipop+1]]
    }

#save(mylogpop_neg,file=paste0(myfolder,"/mylogpop_neg.RData"))
#save(mylogpop_pos,file=paste0(myfolder,"/mylogpop_pos.RData"))
#save(logpop_nAB_pos,file=paste0(myfolder,"/logpop_nAB_pos.RData"))
#save(logpop_nAB_neg,file=paste0(myfolder,"/logpop_nAB_neg.RData"))

load(paste0(myfolder,"/logpop_nAB_neg.RData"))
load(paste0(myfolder,"/logpop_nAB_pos.RData"))

n_sign_links<-dim(dataFIN[empfdr<0.05,c(1,2,3,4),with=FALSE])[1]
    
for ( myperm in 1:2){ #total nAB kept fixed
    mydata<-fread(paste0("zcat ",myfolder,"/reschr",myperm,"/anal/sorted.res.gz"))
#    data<-read.table(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding1000g/new/reschr",myperm,"/anal/sorted.res"))
    names23<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","prho2","pfisher","pKuli","T2","Xtot","X","pop","popX")
    names20<-c("chrA","chrB","posA","posB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","nAA","nBB","D","D1","rho2","T2","Xtot","X","pop","popX")
    if (computed_exact_pvalues) { mynames<-names23 } else { mynames<-names20 }
    #names(data)<-mynames #if data.frame
    setnames(mydata,mynames)
    combined_pval<-dchisq(mydata$X,df=2*mydata$popX)    
    mylogpop_neg_reschr<-list();logpop_nAB_neg_reschr<-list();mylogpop_pos_reschr<-list();logpop_nAB_pos_reschr<-list();logpop_nAB_reschr<-list()
for ( ipop in 0:11)    
    {
    print(paste("population: ",ipop))
    fraction_nlinks_neg<-1
    fraction_nlinks_pos<-1
    mydata_temp<-fread(paste0(myfolder,"/reschr",myperm,"/logs/all.pop",ipop,".log"))
    mynames_temp<-names(mydata_temp);mynames_temp[1:4]<-c("chrA","chrB","posA","posB")
    setnames(mydata_temp,mynames_temp)
    for (i in 1:5)
        {
        print(paste("iteration: ",i))
        sign_links_reschr_neg<-dataFIN_l[[myperm]][order(combined_pval_l[[myperm]])<=round(n_sign_links/fraction_nlinks_neg),c(1,2,3,4),with=FALSE]
        sign_links_reschr_pos<-dataFIN_l[[myperm]][order(combined_pval_l[[myperm]])<=round(n_sign_links/fraction_nlinks_pos),c(1,2,3,4),with=FALSE]
        #temp_neg<-intersect_links(mydata,sign_links_reschr_neg);temp_pos<-intersect_links(mydata,sign_links_reschr_pos)
        mysubset_neg<-dataLDLD2Dstr12_f(mydata,sign_links_reschr_neg,mycriterion="12")[[1]]
        mysubset_pos<-dataLDLD2Dstr12_f(mydata,sign_links_reschr_pos,mycriterion="12")[[2]]
        mysubset_neg<-mysubset_neg[mysubset_neg$pop==ipop,]
        mysubset_pos<-mysubset_pos[mysubset_pos$pop==ipop,]
        mylogpop_pos_reschr[[ipop+1]]<-intersect_links(mydata_temp,mysubset_pos)
        mylogpop_neg_reschr[[ipop+1]]<-intersect_links(mydata_temp,mysubset_neg)
        logpop_nAB_neg_reschr[[ipop+1]]<-apply(mylogpop_neg_reschr[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) 2-nABf(y))))
        logpop_nAB_pos_reschr[[ipop+1]]<-apply(mylogpop_pos_reschr[[ipop+1]][,5:(ncol(mydata_temp)-1),with=FALSE],MARGIN=2,FUN=function(x) sum(sapply(x, function(y) nABf(y))))
        fraction_nlinks_neg<-fraction_nlinks_neg*(sum(logpop_nAB_neg_reschr[[ipop+1]])/sum(logpop_nAB_neg[[ipop+1]]))
        fraction_nlinks_pos<-fraction_nlinks_pos*(sum(logpop_nAB_pos_reschr[[ipop+1]])/sum(logpop_nAB_pos[[ipop+1]]))
        }
    save(logpop_nAB_neg_reschr,file=paste0(myfolder,"/logpop_nAB_neg_reschr",myperm,"_nAB.RData"))
    save(logpop_nAB_pos_reschr,file=paste0(myfolder,"/logpop_nAB_pos_reschr",myperm,"_nAB.RData"))
    save(mylogpop_neg_reschr,file=paste0(myfolder,"/mylogpop_neg_reschr",myperm,".RData"))
    save(mylogpop_pos_reschr,file=paste0(myfolder,"/mylogpop_pos_reschr",myperm,".RData"))
    }
}

}

require(vioplot)
load(paste0(myfolder,"/logpop_nAB_pos.RData"))
load(paste0(myfolder,"/logpop_nAB_neg.RData"))
load(paste0(myfolder,"/logpop_nAB.RData"))
logpop_reschr_nAB_l<-list()
for (z in 1:nperm) { load(paste0(myfolder,"/logpop_reschr",z,"_nAB.RData")); logpop_reschr_nAB_l[[z]]<-logpop_reschr_nAB; }
logpop_reschr_nAB_neg_l<-list()
logpop_reschr_nAB_pos_l<-list()
for (z in 1:nperm) { load(paste0(myfolder,"/logpop_nAB_neg_reschr",myperm,"_nAB.RData")); logpop_reschr_nAB_neg_l[[z]]<-logpop_nAB_neg_reschr; }
for (z in 1:nperm) { load(paste0(myfolder,"/logpop_nAB_pos_reschr",myperm,"_nAB.RData")); logpop_reschr_nAB_pos_l[[z]]<-logpop_nAB_pos_reschr; }


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

#infoplots
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
save.image("~/Dropbox/LDLD/temp.RData")
#test if nAB explains differences in mutational load
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
gdriveLDLD


mytest<-cor.test(mymutsamples_rel,mylogpop_rel)

}
#test if nAB inhomogeneous
{
load(paste0(myfolder,"/logpop_nAB.RData"))
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
#from log files. One possibility is to use function logpop2gen, that goes from log files to gen files, with 0 1 and 2 for homozygotes
#Otherwise, one can use a similar procedure as attached below. However very slow. 
#Therefore I recommend using all.freqlog files generated in C++ (see above in command line comments "create whole genome freqlog file".)
#
{
#save.image("~/Dropbox/LDLD/temp.RData")
load("~/Dropbox/LDLD/temp.RData")
#load all.freqlog file #I used to call data.table with all the sites resres. Shitty name. Now I change it to nAall.
nAall<-fread(paste0("zcat ",myfolder,"/all.freqlog.gz"))
if (is.na(nAall[1,dim(nAall)[2],with=FALSE])) {nAall<-nAall[,1:(dim(nAall)[2]-1),with=FALSE]}
#nAallnames<-paste0(nAall[,1],".",nAall[,2])
nAall_l<-list()
for (i in 1:12)
{
nAall_l[[i]]<-nAall[,samples_imypopi(nspops,i)+2,with=FALSE]
}
snpsnames<-apply(nAall,MARGIN=1,function(x) paste0(x[1:2],collapse="."))

load(paste0(myfolder,"/dataFIN.RData"))
load(paste0(myfolder,"/combined_pval.RData"))
load(paste0(myfolder,"/empfdr.RData"))
load(paste0(myfolder,"/dataFIN_l.RData"))
load(paste0(myfolder,"/combined_pval_l.RData"))
load(paste0(myfolder,"/logpop_nAB.RData"))
load(paste0(myfolder,"/logpop_nAB_nlinks.RData"))
#logisticpvals_nAvsnAB_l<-corr_nAvsnAB(nAall_l,logpop_nAB,l_mutsamples)
#save(logisticpvals_nAvsnAB_l,file=paste0(myfolder,"/logisticpvals_nAvsnAB_l.RData"))
load(paste0(myfolder,"/logisticpvals_nAvsnAB_l.RData"))
load(paste0(myfolder,"/dataFIN_sign.RData"))

#in all SNPS #553
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
write.table(stringsnp2bed(as.matrix(removed_snps)),file=paste0(myfolder,"/removed1.bed"),row.names = FALSE,quote=FALSE,sep='\t',col.names=FALSE)

#write corr_nAvsnAB with random slopes
library(lme4)
#politeness.model = lmer(frequency ~ attitude + (1|subject) + (1|scenario), data=politeness)
data0<-data.table(sampleID=1:length(unlist(logpop_nAB)),population=unlist(sapply(1:length(logpop_nAB),function(x) rep(x,length(logpop_nAB[[x]])))),nAB=unlist(logpop_nAB))
data0<-rbind(data0,data0)
data0$sampleID<-as.factor(data0$sampleID)
data0$population<-as.factor(data0$population)

nsamples<-dim(nAall)[2]-2
nsites<-dim(nAall)[1]
#if small substructure, can I take that into account? if sample very different, more linkage. therefore it makes no sense to correct by distance. help comes from multiple populations!
nAall1<-nAall[,3:dim(nAall)[2],with=FALSE]
nAall2<-nAall[,3:dim(nAall)[2],with=FALSE]
for (i in names(nAall1)) { nAall1[eval(parse(text=i))==2,i]<-1;print(i) }
for (i in names(nAall2)) { nAall2[eval(parse(text=i))==1,i]<-0;print(i) }
for (i in names(nAall2)) { nAall2[eval(parse(text=i))==2,i]<-1;print(i) }
nAall1<-rbind(nAall1,nAall2)
iline<-1
data_temp<-cbind(data0,nA=c(unlist(nAall1[iline,]),unlist(nAall1[iline+nsites,])))
contr=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000000),calc.derivs=F)
data_temp$z.nAB<- scale(data_temp$nAB)

source("/home/fabrizio_mafessoni/Dropbox/dropbox_tablet/math/statistics/linearmodels/diagnostic_fcns.r")
xx.fe.re=fe.re.tab(fe.model="nA ~ z.nAB", 
re="(1|population) + (1|sampleID)",other.vars=NULL, data=as.data.frame(data_temp))
xx.fe.re$summary

my.model.random.intercepts = glmer( nA ~ z.nAB + (1|population) + (1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
fx<-function(){
my.model = glmer( nA ~ z.nAB + (1|population) + (0+z.nAB|population)+(1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
#Roger changes with this (analogous) and to drop only fixed effect but not random slope.This is like assuming that variability but that it could go in any direction
expandDoubleVerts(nA ~ z.nAB + (1+z.nAB||population)+(1|sampleID)) #nA ~ z.nAB + ((1 | population) + (0 + z.nAB | population)) + (1 | sampleID)


#is there a better proxy rather than nAB? which associated? pair for which positive. #would it solve bias towards rare? actually I should only consider positive linkage



my.model.null = glmer( nA ~ (1|population) +(1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
anova(my.model,my.model.null)[[8]]
}
system.time(fx())


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
#system("cat /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/*minilog | grep ncomparisons | awk '{print $14}' > /mnt/scratch/fabrizio/LDLD/above95/coding1000g/removal1/ncomparisons.txt")
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
