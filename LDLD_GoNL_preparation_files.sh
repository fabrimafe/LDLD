cd /mnt/restricted/GoNL/GoNL1/release5.4/03_IL_SNVs/gonl-abc_samples
zcat /mnt/restricted/GoNL/GoNL1/release5.4/03_IL_SNVs/gonl-abc_samples/gonl-abc_samples.chr1.release5.raw_SNVs.vcf.gz | head -300 | grep CHROM | awk '{for (i=10;i<=NF;i++){print $i}}' | grep -v c | grep -v d > ~/workspace/1000genomes/GoNL1.unrelated.samples
zcat /mnt/restricted/GoNL/GoNL1/release5.4/03_IL_SNVs/gonl-abc_samples/gonl-abc_samples.chr1.release5.raw_SNVs.vcf.gz | head -300 | grep CHROM | awk '{for (i=10;i<=NF;i++){print $i}}' | grep c > ~/workspace/1000genomes/GoNL1.offspring.samples 

#--------------GoNL1-------------
#2) generate chrZ.vcf.gz files
{
myannotation="coding.exons"
#myannotation="introns"
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
for MYCHROM in `seq 1 22`; do
echo $MYCHROM
zcat ${myannotation}.bed  | awk -v OFS='\t' -v MYCHROM=$MYCHROM '{if ($1=="chr"MYCHROM){print MYCHROM,$2,$3}}' |sort -Vu -k1,1 -k2,2n > ${myannotation}/myfilter.${MYCHROM}.bed
bedtools merge -i ${myannotation}/myfilter.${MYCHROM}.bed > mytemp.${MYCHROM}.bed; mv mytemp.${MYCHROM}.bed ${myannotation}/myfilter.${MYCHROM}.bed
myvcfsource=/mnt/restricted/GoNL/GoNL1/release5.4/03_IL_SNVs/gonl-abc_samples/gonl-abc_samples.chr${MYCHROM}.release5.raw_SNVs.vcf.gz
tabix $myvcfsource -H -R myfilter.${myannotation}/myfilter.${MYCHROM}.bed  > ${myannotation}/chr$MYCHROM.vcf
tabix $myvcfsource -R ${myannotation}/myfilter.${MYCHROM}.bed  | sort -Vu -k1,1 -k2,2n >> ${myannotation}/chr${MYCHROM}.vcf
bgzip -f ${myannotation}/chr${MYCHROM}.vcf
tabix ${myannotation}/chr${MYCHROM}.vcf.gz
done

#subsample intergenic to the same number of intronic (not really) and reorder them
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
myannotation=intergenic
SEED=1;MYCHROM=11
while [ $MYCHROM -lt 23 ]
        do
        echo $MYCHROM
        #ARGUMENT NLINES: set here the final number of lines
        NLINES=$( zcat coding.exons/chr${MYCHROM}.vcf.gz | wc -l)
        DESTINATIONFILE=/mnt/scratch/fabrizio/LDLD/GoNL/intergenic/chr${MYCHROM}.intergenic.seed${SEED}.vcf
        myvcfsource=/mnt/restricted/GoNL/GoNL1/release5.4/03_IL_SNVs/gonl-abc_samples/gonl-abc_samples.chr${MYCHROM}.release5.raw_SNVs.vcf.gz
        NLINES_SOURCE=$( zcat $myvcfsource |wc -l  )
        MYPROB=$( awk -v var1=$NLINES -v var2=$NLINES_SOURCE 'BEGIN{print 10*var1/var2}' )
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
for MYCHROM in `seq 11 22`; do
echo $MYCHROM 'reorder...' 
nohup Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/GoNL1.unrelated.samples ${myannotation}/chr${MYCHROM}.vcf.gz ${myannotation}/chr${MYCHROM}.dest.vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 0
echo $MYCHROM 'reorder reschr...' 
#for reschr in `seq 1 2`;do #
reschr=1
mkdir -p ${myannotation}/reschr${reschr}
nohup Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/GoNL1.unrelated.samples ${myannotation}/chr${MYCHROM}.vcf.gz ${myannotation}/reschr${reschr}/chr${MYCHROM}.dest.vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 1
#done
#rm ${myannotation}/myfilter.${MYCHROM}.bed
done

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
MYVCF=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min${myfreq}/chr${MYCHROM}.vcf.gz
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt -f 0.0${myfreq} -Q -O 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYVCF}
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 0.0${myfreq} 1 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYVCF}
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt -f 0.0${myfreq} -O 0 | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYVCF}
done;done


#---global filtered file for introns (truly global) and intergenic (subsampled)
myannotation="introns"
myannotation="intergenic"
myannotation="coding.exons"
for myfreq in 1 5; do
for MYCHROM in `seq 1 10`; do
echo "chrom: " $MYCHROM "freq: " $myfreq
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
mkdir ${myannotation}/min$myfreq
MYFREQVCF=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/min$myfreq/chr${MYCHROM}.vcf
zcat ${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt -f 0.0${myfreq} -O 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYFREQVCF}.temp.gz
done;done


}
#5b) further subsample global files for introns and intergenic matching the size of exons
{
preparationfolder=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/
cd $preparationfolder
myannotation="coding.exons"
myannotation="intergenic"
myannotation_ref="coding.exons"
for myfreq in 1 5; do
for MYCHROM in `seq 11 22`; do
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
#myannotation="coding.exons" #myrepl=2
for myrepl in `seq 1 2`; do
for myfreq in 5; do
for MYCHROM in `seq 11 22`; do
echo 'repl: ' $myrepl ' ; chrom: ' $MYCHROM
mysim=$myrepl; if [ $myannotation == "coding.exons" ]; then mysim=1;fi
MYFREQVCF=${myannotation}/min$myfreq/chr${MYCHROM}.vcf.temp.gz
MYBEDFILE=${myannotation}/min$myfreq/chr${MYCHROM}.subsamples.bed.gz
destination_folder=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/${myannotation}
mkdir -p $destination_folder/repl${myrepl}/min${myfreq}/
mkdir -p $destination_folder/repl${myrepl}/min${myfreq}/reschr1
mkdir -p $destination_folder/min${myfreq}/reschr1
bedtools intersect -a $MYFREQVCF -b <(zcat $MYBEDFILE | awk -v OFS='\t' -v repl=$mysim '{if ($4==repl){print}}' ) -sorted -header | bgzip -f > $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz
Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/GoNL1.unrelated.samples $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz $destination_folder/repl${myrepl}/min${myfreq}/reschr1/chr${MYCHROM}.vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 1
#Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/GoNL1.unrelated.samples $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz $destination_folder/repl${myrepl}/min${myfreq}/reschr1/chr${MYCHROM}.vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 1
done;done;done

}

#run LDLD
#GoNL1
{
minfreq=5
#GoNL subdividing pops
./ICLD11multi.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min${minfreq} 4 coding.exons 0 1 22 nspops_GoNL_short.txt 0.00000000001 0.0${minfreq}
./ICLD11multi.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min${minfreq}/reschr1 4 coding.exons 0 1 22 nspops_GoNL_short.txt 0.00000000001 0.0${minfreq}

for myrepl in 1 2;do
./ICLD11multi.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl}/min${minfreq} 4 intergenic 0 1 22 nspops_GoNL_short.txt 0.00000000001 0.0${minfreq}
done
for myrepl in 1 2;do
#./ICLD11multi.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl}/min${minfreq} 4 intergenic 0 1 22 /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 0.000001 0.0${minfreq}
./ICLD11multi.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl}/min${minfreq}/reschr1 4 intergenic 0 1 22 /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 0.000001 0.0${minfreq}
done

cp /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl}/min${minfreq} 


#GoNL with T2
mkdir -p /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/repl11/min${minfreq}/reschr1
cp /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min${minfreq}/reschr1/*vcf.gz /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/repl11/min${minfreq}/reschr1
cp /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min${minfreq}/*vcf.gz /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/repl11/min${minfreq}
minfreq=5
for myrepl in 1 2;do
myrepl2=$(( $myrepl + 10 ))
mkdir -p /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl2}/min${minfreq}
mkdir -p /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl2}/min${minfreq}/reschr1
cp /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl}/min${minfreq}/reschr1/*vcf.gz /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl2}/min${minfreq}/reschr1
#cp /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl}/min${minfreq}/*vcf.gz /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl2}/min${minfreq}
done

myrepl2=11
./ICLD11multic.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl2}/min${minfreq} 4 coding.exons 0 1 22 nspops_GoNL.txt 0.0000000001 0.0${minfreq}
./ICLD11multic.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl2}/min${minfreq}/reschr1 4 coding.exons 0 1 22 nspops_GoNL.txt 0.0000000001 0.0${minfreq}
for myrepl in 1 2;do
myrepl2=$(( $myrepl + 10 ))
./ICLD11multic.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl2}/min${minfreq} 4 intergenic 0 1 22 nspops_GoNL.txt 0.0000000001 0.0${minfreq}
./ICLD11multic.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl${myrepl2}/min${minfreq}/reschr1 4 intergenic 0 1 22 nspops_GoNL.txt 0.0000000001 0.0${minfreq}
done



minfreq=5
cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min${minfreq}
cat chr21.tab | head -21 > chr21.short.tab
cat chr22.tab | head -21 > chr22.short.tab
-p 0.000000000000001
./ICLDmulti.out -A chr21.short.tab -B chr22.short.tab -p 0.001 -f 0.05 -N nspops_GoNL.txt
date; ./ICLDmulti.out -A chr21.short.tab -B chr22.short.tab -p 0.00000001 -f 0.05 -N nspops_GoNL.txt; date
./ICLDmulti.out -A chr21.short.tab -B chr22.short.tab -p 0.00000001 -f 0.05 -N nspops_GoNL.txt


time ./ICLDmulti.out -A chr21.short.tab -B chr22.short.tab -p 0.00000001 -f 0.05 -N nspops_GoNL_short.txt
date; ./ICLDmulti.out -A chr21.short.tab -B chr22.short.tab -p 0.000000000000001 -f 0.05 -N nspops_GoNL.txt; wc -l chr21.short.tabchr22.short.tab.res; date;

10910347	16449050	rs150482	rs148001219	436.0	311.0	996.0	96.50	99	57	-0.0795995226	-0.0000002075	0.1079751084
10910347	17446991	rs150482	rs4819925	436.0	412.0	996.0	241.00	99	119	0.1217802939	0.0000005278	0.1886343057



./ICLD11multi.sh /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min${minfreq} 4 coding.exons 0 1 22 nspops_GoNL.txt 0.00000001 0.0${minfreq}







/mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min5$ ./LDLD.out -A chrA.vcf.gz -B chrB.vcf.gz -N nspops_GoNL.txt -f 0.01 #incanta
./ICLDmulti.out -A chr1_short.vcf -B chr2_short.vcf -N nspops.txt #va bene
./ICLDmulti.out -A chr1_short.vcf -B chr2_short.vcf -N nspops1.txt


MYCHROMA=21
MYCHROMB=22
myfolder=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min5
cp ~/Dropbox/LDLD/scripts/LDLD.out $myfolder
cd $myfolder
cp /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt $myfolder
./LDLD.out -A chr${MYCHROMA}.vcf.gz -B chr${MYCHROMB}.vcf.gz -N nspops_GoNL.txt -R test.res -L test.log -S test.sum -M test.minilog -f 0.05



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

}

#very slow, why?
{
cat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min5/chr1.tab | grep -v '#' | wc -l #5284
cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr1.tab | grep -v '#' | wc -l #9594


}

#--------------GoNL1+GoNL2-------------
#3) reorder
{
myannotation="coding.exons"
myannotation="introns"
myannotation="intergenic"
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
for MYCHROM in `seq 11 22`; do
echo $MYCHROM 'reorder...' 
nohup Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/GoNL1.offspring.samples ${myannotation}/chr${MYCHROM}.vcf.gz ${myannotation}/chr${MYCHROM}.dest.vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 0
echo $MYCHROM 'reorder reschr...' 
#for reschr in `seq 1 2`;do #
reschr=1
mkdir -p ${myannotation}/reschr${reschr}
nohup Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/GoNL1.offspring.samples  ${myannotation}/chr${MYCHROM}.vcf.gz ${myannotation}/reschr${reschr}/chr${MYCHROM}.dest.vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 1
#done
#rm ${myannotation}/myfilter.${MYCHROM}.bed
done

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
MYVCF=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/coding.exons/min${myfreq}/chr${MYCHROM}.vcf.gz
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt -f 0.0${myfreq} -Q -O 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYVCF}
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 0.0${myfreq} 1 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYVCF}
zcat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt -f 0.0${myfreq} -O 0 | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYVCF}
done;done


#---global filtered file for introns (truly global) and intergenic (subsampled)
myannotation="introns"
myannotation="intergenic"
myannotation="coding.exons"
for myfreq in 1 5; do
for MYCHROM in `seq 1 10`; do
echo "chrom: " $MYCHROM "freq: " $myfreq
cd /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g
mkdir ${myannotation}/min$myfreq
MYFREQVCF=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/${myannotation}/min$myfreq/chr${MYCHROM}.vcf
zcat ${myannotation}/chr${MYCHROM}.dest.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt -f 0.0${myfreq} -O 0  | grep -v 'Reading\|number of samples\|nsamples in pop' | bgzip -f > ${MYFREQVCF}.temp.gz
done;done


}
#5b) further subsample global files for introns and intergenic matching the size of exons
{
preparationfolder=/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/
cd $preparationfolder
myannotation="coding.exons"
myannotation="intergenic"
myannotation_ref="coding.exons"
for myfreq in 1 5; do
for MYCHROM in `seq 11 22`; do
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
#myannotation="coding.exons" #myrepl=2
for myrepl in `seq 1 2`; do
for myfreq in 5; do
for MYCHROM in `seq 11 22`; do
echo 'repl: ' $myrepl ' ; chrom: ' $MYCHROM
mysim=$myrepl; if [ $myannotation == "coding.exons" ]; then mysim=1;fi
MYFREQVCF=${myannotation}/min$myfreq/chr${MYCHROM}.vcf.temp.gz
MYBEDFILE=${myannotation}/min$myfreq/chr${MYCHROM}.subsamples.bed.gz
destination_folder=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/${myannotation}
mkdir -p $destination_folder/repl${myrepl}/min${myfreq}/
mkdir -p $destination_folder/repl${myrepl}/min${myfreq}/reschr1
mkdir -p $destination_folder/min${myfreq}/reschr1
bedtools intersect -a $MYFREQVCF -b <(zcat $MYBEDFILE | awk -v OFS='\t' -v repl=$mysim '{if ($4==repl){print}}' ) -sorted -header | bgzip -f > $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz
Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/GoNL1.unrelated.samples $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz $destination_folder/repl${myrepl}/min${myfreq}/reschr1/chr${MYCHROM}.vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 1
#Rscript ~/Dropbox/LDLD/scripts/vcf2reorderedvcf_awk.R ~/workspace/1000genomes/GoNL1.unrelated.samples $destination_folder/repl${myrepl}/min${myfreq}/chr${MYCHROM}.vcf.gz $destination_folder/repl${myrepl}/min${myfreq}/reschr1/chr${MYCHROM}.vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/nspops_GoNL.txt 1
done;done;done

}


#GoNL12_as1pop
cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5
cat *res | awk '{if (NF==18){print}}' | gzip -f > res.gz &
cat reschr1/*res | awk '{if (NF==18){print}}' | gzip -f > reschr1_res.gz &
library(data.table)
mydata<-fread("zcat res.gz")
mydatanull<-fread("zcat reschr1_res.gz")

sum(mydata$V11>0.45) #149
sum(mydatanull$V11>0.45) #93

sum(mydata$V11< -0.45) #7
sum(mydatanull$V11< -0.45) #3

sum(mydata$V11< -0.425) #1521
sum(mydatanull$V11< -0.425) #1383

sum(mydata$V12>4.5e-05) #3555
sum(mydatanull$V12>4.5e-05) #797

sum(mydata$V12>6e-05) #2755
sum(mydatanull$V12>6e-05) #387

sum(mydata$V13> 0.99) #305
sum(mydatanull$V13 > 0.99)  #228

#probably easiest way is to generate fake files with exact p-values and then parse with standard script.
zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims_pbatch0.2_prisk0.001_nsample50_perror0.1_min0.05/repl1/anal/sorted.res.gz | head
#1 2 121 1208 F F 3.0 4.0 100.0 1.50 0 0 0.0276000000 0.0000958333 0.8304039470 0.00020408163265305565 0.99999999999997823963 0.00090497737556560493 0.00000000011667122823 45.743322 45.743322 0 1 


cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5 
mkdir anal
mkdir reschr1/anal
zcat res.gz | awk '{if (NF==18){print}}' | awk '{
if ($12>6e-05) {pvals=0} else {pvals=1};
print 1,2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,pvals,pvals,pvals,pvals,pvals,$15,$16,$17,$18
}' | gzip -f > anal/sorted.red.gz
zcat reschr1_res.gz | awk '{if (NF==18){print}}' | awk '{
if ($12>6e-05) {pvals=0} else {pvals=1};
print 1,2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,pvals,pvals,pvals,pvals,pvals,$15,$16,$17,$18
}' | gzip -f > reschr1/anal/sorted.red.gz

/home/fabrizio_mafessoni/mylogs/mylogs_GoNL


#GoNL1
cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl11/min5
zcat *res.gz | awk '{if (NF==21){print}}' | gzip -f > res.gz &
zcat reschr1/*res.gz | awk '{if (NF==21){print}}' | gzip -f > reschr1_res.gz

library(data.table)
mydata<-fread("zcat *.res.gz")
mydatanull<-fread("zcat reschr1_res.gz")

sum(mydata$V11>0.45) #0
sum(mydatanull$V11>0.45) #0 

sum(mydata$V11< -0.45) #0
sum(mydatanull$V11< -0.45) #0

sum(mydata$V11< -0.425) #0
sum(mydatanull$V11< -0.425) #0

sum(mydata$V12>4.5e-05) #0
sum(mydatanull$V12>4.5e-05) #0

sum(mydata$V12>6e-05) #0
sum(mydatanull$V12>6e-05) #0 

sum(mydata$V13> 0.99) #1364
sum(mydatanull$V13 > 0.99)  #1208

sum(mydata$V13> 0.9999) #1334
sum(mydatanull$V13 > 0.9999)  #1176

sum(mydata$V13> 0.97) #7610
sum(mydatanull$V13 > 0.97)  #7594

#not informative
#sum(mydata$V19> 54) 
#sum(mydatanull$V19 > 54)  


cat chr*.tab.res | awk '{if (NF==18 && ( $12>0.000004 || $12 < -0.000004 ) ){print}}' | wc -l >> tab
cat reschr1/chr*.tab.res | awk '{if (NF==18 && ( $12>0.000004 || $12 < -0.000004 ) ){print}}' | wc -l >> tab

cat chr*.tab.res | awk '{if (NF==18 && ( $11>0.4 || $11 < -0.4 ) ){print}}' | wc -l >> tab
cat reschr1/chr*.tab.res | awk '{if (NF==18 && ( $11>0.4 || $11 < -0.4 ) ){print}}' | wc -l >> tab

cat chr*.tab.res | awk '{if (NF==18 && ( $11>0.45 || $11 < -0.45 ) ){print}}' | wc -l >> tab
cat reschr1/chr*.tab.res | awk '{if (NF==18 && ( $11>0.45 || $11 < -0.45 ) ){print}}' | wc -l >> tab

cat chr*.tab.res | awk '{if (NF==18 && $13>0.7 ){print}}' | wc -l >> tab
cat reschr1/chr*.tab.res | awk '{if (NF==18 && $13>0.7 ){print}}' | wc -l >> tab


cat reschr1/chr1.tabchr12.tab.res | awk '{if (NF==18 && ( $11>0.25 || $11 < -0.25 ) ){print}}' | wc -l



#parsing GoNL2
{
#in absence of information of who the parents are I just use kids from GoNL1
vcf_filterPASS () { 
awk '{mysub=substr($10,1,3); if ( ( substr($1,1,1)=="#") || (mysub=="1/1" || mysub=="0/1" || mysub=="1/0" || mysub=="1|1" || mysub=="0|1" || mysub=="1|0")){print}}'
}

cd /mnt/restricted/GoNL/GoNL2/vcf
MYI=1
MYFOLDERS=$( ls )
for IFOLDER in $MYFOLDERS;do
echo $MYI
IF2FOLDER=$( ls $IFOLDER)
bzcat $IFOLDER/${IF2FOLDER}/*bz2 | vcf_filterPASS | gzip -f > temp${MYI}.vcf.gz
MYI=$(( $MYI + 1 ))
done


for MYI in `seq 3 100`;do
echo $MYI
gunzip temp${MYI}.vcf.gz
bgzip temp${MYI}.vcf
done

MYI=
for IFOLDER in $MYFOLDERS;do


vcf-merge temp1.vcf.gz temp2.vcf.gz
MYI2=$(( $MYI + 1 ))
joinx vcf-merge temp${MYI}.vcf.gz temp${MYI2}.vcf.gz


echo $MYI
IF2FOLDER=$( ls $IFOLDER)
bzcat $IFOLDER/${IF2FOLDER}/*bz2 | vcf_filterPASS | gzip -f > temp${MYI}.vcf.gz

done


MYVCFS=$( ls /mnt/restricted/GoNL/GoNL2/vcf/temp2*.vcf.gz )
joinx vcf-merge $MYVCFS


}

joinx vcf-merge temp1.vcf temp2.vcf


