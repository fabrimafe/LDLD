#============================
#===== CHECKS ON BAD SNPS  ========
#============================

#test of linear model is from LDLD_1000g_anal.....bio37 conflicted copy.. at #test effects of variables
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
res_pop=lm(nAB ~  (1+Population),data=mydata_model)
summary(res_allcov_Ec_chips)$coefficients
anova(res_allcov_Ec_chips, res_allcov_Ec, test="Chisq")
anova(res_allcov_Ec_chips, res_allcov_chips, test="Chisq")
anova(res_allcov_Ec_chips, res_EC_chips, test="Chisq")
anova(res_allcov_Ec_chips, res_pop, test="Chisq")

#effect sizes
#from https://stats.stackexchange.com/questions/95054/how-to-get-an-overall-p-value-and-effect-size-for-a-categorical-factor-in-a-mi
r2.corr.mer<-function(m) {
lmfit <- lm(model.response(model.frame(m)) ~ fitted(m))
summary(lmfit)$r.squared
}

Xu_method_f<-function(m) 1-var(residuals(m))/(var(model.response(model.frame(m))))
library("MuMIn")
r.squaredGLMM(res_allcov_Ec_chips)
r.squaredGLMM(res_allcov_chips)
r.squaredGLMM(res_allcov_Ec)
r.squaredGLMM(res_EC_chips)
r.squaredGLMM(res_pop)

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
res_pop=lm(nAB ~  (1+Population),data=mydata_model)

summary(res_allcov_Ec_chips)$coefficients
anova(res_allcov_Ec_chips, res_allcov_Ec, test="Chisq")
anova(res_allcov_Ec_chips, res_allcov_chips, test="Chisq")
anova(res_allcov_Ec_chips, res_EC_chips, test="Chisq")
anova(res_allcov_Ec_chips, res_pop, test="Chisq")
#outputs save in ~/Dropbox/LDLD/ms/Ajhg/table_model_namesSeqCenters.odt
}

#============================
#===== TESTS before revisions ========
#============================
{
#GWAS
{
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
bedtools intersect -a <( sortBed -i /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/removed1_pos20.bed ) -b <( cat ~/Dropbox/general_utils/gwas_catalog_v1.0-associations_e89_r2017-07-31_hg19.bed | sed 's/chr//g' )
#in intergenic

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

#Complete Genomics (maybe old)
{

#generate files
#CGBED=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/#Complete_Public_Genomes_69genomes_all_listvariants.bed.gz
#bzcat /mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/#Complete_Public_Genomes_69genomes_all_listvariants.tsv.bz2 | grep -v chromosome | awk -v OFS='\t' '{print $2,$3,$4}' | gzip -f > temp.bed.gz
#bedtools sort -i temp.bed.gz | gzip -f > $CGBED #maybe I should only do snps. Now it is snps and things falling within deletions. no insertions.
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
zcat $CGBED | awk -v mychrom=chr$MYCHROM '{if (($1==mychrom) && ($3>$2)){print}}' | sed 's/chr//g' | gzip -f > $CGBED_CHR &
done


#create background set in R
{
#SET UP VARIABLES
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
MYDESTINATION_FOLDER=${MYFOLDER}/anal
mkdir -p $MYDESTINATION_FOLDER

#Rscript ~/Dropbox/LDLD/scripts/plot_badsnps_vs_CG.R $MYFOLDER/anal/byfreq_summary_${MYBADSNPS}.tab

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
myfolder="/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl1/min1"
myfolderdropbox="~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl1/min1"
myfile="byfreq_summary_removed1.bed.tab"

myfolder="/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5"
myfolderdropbox="~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min5"
myfile="byfreq_summary_removed1_pos_shell_scale.bed.tab"

myfolder="/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5"
myfolderdropbox="~/Dropbox/LDLD/ipynb/figs/coding.exons.1000g/repl2/min5"
myfile="byfreq_summary_removed1_pos.bed.tab" #this looks very good. How did I obtain it??!!
cp byfreq_summary_removed1_pos.bed.tab byfreq_summary_removed1_pos.bed.tab.fromwhere
#/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/anal$ cat byfreq_summary_removed1_pos.bed.tab | awk 'BEGIN{co=0}{co=co+$5}END{print co}'
745

myfolder="/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5"
myfolderdropbox="~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min5"
myfile="byfreq_summary_removed1_pos20.bed.tab" #this looks bad

myfolder="/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/bad_snps"
myfolderdropbox="~/Dropbox/LDLD/ipynb/figs/intergenic.1000g/merged/min5"
myfile="byfreq_summary_removed1_pos20.bed.tab" #this looks very good. How did I obtain it??!!

mydata<-read.table(paste0(myfolder,"/anal/",myfile))
names(mydata)<-c("freq","tot1000ginCG","tot1000g","bad1000ginCG","bad1000g")
source("~/Dropbox/general_utils/general_functions.R")
library(xts)
#smaller frequency bin almost empty
mydata[2,1]<-mydata[1,1]
mydata[2,]<-mydata[1,]+mydata[2,]
mydata<-mydata[2:nrow(mydata),]
data_temp<-c(rbind(mydata$bad1000ginCG/mydata$bad1000g, mydata$tot1000ginCG/mydata$tot1000g))
#data_temp<-c(data_temp,0,0,sum(mydata$bad1000ginCG)/sum(mydata$bad1000g),sum(mydata$tot1000ginCG)/sum(mydata$tot1000g))
data_temp<-c(data_temp,0,0,sum(mydata$bad1000ginCG)/sum(mydata$bad1000g),sum(mydata$tot1000ginCG)/sum(mydata$tot1000g))
myxlabs<-list(pos=seq(1.5,length(data_temp),2),lab=c(mydata$freq,'','tot'))

pdf(paste0(myfolderdropbox,"/",myfile,"_sharingCG.pdf"))
mybarplot_f(1:length(data_temp),data_temp,myylim=c(0,1.15),mycols=rep(c("gold","cadetblue"),length(data_temp)),myxlabel='minimum frequency 1000g', myylabel='% sharing CG , % false discoveries',add=FALSE,xaxis=myxlabs,cex.axis=0.9)
fd<-(mydata$bad1000ginCG/mydata$bad1000g)/(mydata$tot1000ginCG/mydata$tot1000g)
points(myxlabs$pos[1:nrow(mydata)],fd,pch=19,type='o')
points(last(myxlabs$pos),(sum(mydata$bad1000ginCG)/sum(mydata$bad1000g))/(sum(mydata$tot1000ginCG)/sum(mydata$tot1000g)),pch=19)
       legend(0,1.15,c("candidates","all","% FD (candidates/all)"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        pch=(c(15,15,19)),col=c("gold","cadetblue",'black'))
dev.off()

chisq.test(rbind(c(sum(mydata$bad1000g),sum(mydata$bad1000ginCG)),
c(sum(mydata$tot1000g),sum(mydata$tot1000ginCG))))


chisq.test









#TOTAL NUMBER OF BAD_SNPS OVERLAPPING WITH CG

MYBADSNPS=removed1.bed
CGBED=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/Complete_Public_Genomes_69genomes_all_listvariants.bed.gz
bzcat /mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/Complete_Public_Genomes_69genomes_all_listvariants.tsv.bz2 | grep -v chromosome | awk -v OFS='\t' '{print $2,$3,$4}' | gzip -f > temp.bed.gz
bedtools sort -i temp.bed.gz | gzip -f > $CGBED
for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
zcat $CGBED | awk -v mychrom=chr$MYCHROM '{if (($1==mychrom) && ($3>$2)){print}}' | sed 's/chr//g' | gzip -f > $CGBED_CHR &
bedtools intersect -a $CGBED_CHR -b <( bedtools sort -i <( cat ${MYBADSNPS} | awk -v OFS='\t' '{print $1,$2,$3}'  )) -sorted >> ${MYBADSNPS}.CG1.bed
#bedtools intersect -a $CGBED_CHR -b <( bedtools sort -i <( cat ${MYBADSNPS} | awk -v OFS='\t' '{print "chr"$1,$2,$3}'  )) -sorted >> ${MYFOLDER}/${MYBADSNPS}.CG1.bed
done

for MYCHROM in `seq 1 22`;do
CGBED_CHR=/mnt/sequencedb/CompleteGenomics/Public_Genome_Summary_Analysis/listvariants.chr${MYCHROM}.bed.gz
zcat $CGBED_CHR | sed 's/chr//g' | gzip -f > mydata_temp
mv mydata_temp $CGBED_CHR
done




bedtools intersect -a $CGBED -b <( bedtools sort -i <( cat ${MYBADSNPS} | awk -v OFS='\t' '{print "chr"$1,$2,$3}' | sort -Vu -k1,1 -k2,2n  )) -sorted | gzip -f > ${MYBADSNPS}.CG1.bed.gz

#intersect  to see number of overlapping per bin in full dataset
paste0(myfolder,"/removed1.bed.CG1.bed")




#see bins for my snps
#check fdr 
#overalp intergenic
0.32    0.87
#overalp exons
0.23 - 0.89


}


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
#create logo files and test reads imbalance
{
for mysupp in 4 5 6 7;do
cat ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | grep '#' > temp
cat ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | sort -Vu -k1,1 -k2,2n | grep -v '#' >> temp; mv temp ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed
done

#CODING
{
#TAG 'CREATE LOGOS'
#cat above95/coding1000g/anal/snpsA.bed above95/coding1000g/anal/snpsB.bed | sort -Vu -k1,1 -k2,2 > above95/coding1000g/anal/snps.bed
#mydata<-read.table("~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt",stringsAsFactors =FALSE,header=FALSE)
myfolder="/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5"
mydata<-read.table("~/Dropbox/LDLD/ms/GBE/TableS4.bed",stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
dim(mydata) #47193
createbedlogo(mydata,paste0(myfolder,"/basessign.temp"))
bedtools intersect -b ~/Dropbox/LDLD/ms/GBE/TableS4.bed -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 ) > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/temp.bed
paste <( bedtools intersect -b ~/Dropbox/LDLD/ms/GBE/TableS4.bed -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/temp.bed ) /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign | awk '{if (length($4)==1 && length($5)==1){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign.bed
bamlogo_sign<-bamlogo(paste0(myfolder,"/basessign.bed"),tilltheend=TRUE,till=100,startfrom=1,filetemp="temp")
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

#dim(logonotsign)
#res<-logo_below_thr(logonotsign,0.00001) #0.03417355

#background set
bedtools subtract -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' ) -b ~/Dropbox/LDLD/ms/GBE/TableS4.bed | shuf -n 1000 | sort -Vu -k1,1 -k2,2 > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign.bed 
bedtools subtract -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' ) -b ~/Dropbox/LDLD/ms/GBE/TableS4.bed | shuf -n 10000 | sort -Vu -k1,1 -k2,2 > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign2.bed 
mydata<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign.bed",stringsAsFactors =FALSE,header=FALSE)
mydata<-read.table("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign2.bed",stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
createbedlogo(mydata,paste0(myfolder,"/basesnotsign"))
bedtools intersect -b /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign2.bed  -a <( cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 ) > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/temp.bed
#diff /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign2.bed /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/temp.bed #identical
paste <( bedtools intersect -b /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/notsign2.bed -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/temp.bed ) /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign | awk '{if (length($4)==1 && length($5)==1){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign.bed
bamlogo_notsign<-bamlogo("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign.bed",tilltheend=TRUE,till=100,startfrom=1,filetemp="tempnot")
save(bamlogo_notsign,file=paste0(myfolder,"bamlogo_notsign.RData"))
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
pdf(paste0(myfolderdropbox,"/hist_allelic_imbalance_b.pdf"))
data_temp<-c(rbind(myhistsign$count,myhistnotsign$count))
myax=data.frame(pos=seq(1,19,2),lab=mybreaks[1:10])
mybarplot_f(1:length(data_temp),data_temp,myylim=c(0,1),mycols=rep(c("gold","cadetblue"),length(data_temp)),myxlabel='log(p-value)', myylabel='density',add=FALSE,xaxis=myax,cex.axis=0.9)
legend(0,0.6,c("candidates","all"), # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
        pch=(c(15,15)),col=c("gold","cadetblue"))
dev.off()

#double check
myfolder="/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5"
load(paste0(myfolder,"bamlogo_notsign.RData"))
row.names(bamlogo_notsign)<-NULL
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

load(paste0(myfolder,"bamlogo_sign.RData"))
row.names(bamlogo_sign)<-NULL
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

mynot<-t(apply(bamlogo_notsign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))
myt<-t(apply(bamlogo_sign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))


mean(myt[,2]/myt[,1],na.rm=T)
mean(mynot[,2]/mynot[,1],na.rm=T)

mynot<-t(apply(bamlogo_notsign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))
myt<-t(apply(bamlogo_sign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))

require(vioplot)
pdf("~/Dropbox/LDLD/figs/revisions/vio_imbalance.pdf")
vioplot((myt[,2]/myt[,1])[!is.na(myt[,2]/myt[,1])],col=rgb(1,0,0,0.5),drawRect=F)
vioplot((mynot[,2]/mynot[,1])[!is.na(mynot[,2]/mynot[,1])],col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
dev.off()

#plot as kay likes: ratio of reads for alternative allele
tempdata<-bamlogo_sign
temp<-sapply(1:length(tempdata$alt), function(x) tempdata[[tempdata$alt[x]]][x] )/(sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$alt[x]]][x] )+sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$ref[x]]][x] ))
tempsign<-temp[!is.na(temp)]
tempdata<-bamlogo_notsign
temp<-sapply(1:length(tempdata$alt), function(x) tempdata[[tempdata$alt[x]]][x] )/(sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$alt[x]]][x] )+sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$ref[x]]][x] ))
tempnosign<-temp[!is.na(temp)]

require(vioplot)
pdf("~/Dropbox/LDLD/figs/revisions/vio_imbalance_ratio.pdf")
vioplot(tempsign,col=rgb(1,0,0,0.5),drawRect=F,ylim=c(0,1))
vioplot(tempnosign,col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
abline(h=0.5,lty=2)
dev.off()

mycol1=rgb(1,0,0,0.5)
mycol2=rgb(0.3,0.3,0.3,0.3)
mytempdata1<-tempsign
mytempdata2<-tempnosign
mytempdata1<-cbind(c(rep("sign",length(mytempdata1)),rep("nosign",length(mytempdata2))),c(mytempdata1,mytempdata2),c(rep(mycol1,length(mytempdata1)),rep(mycol2,length(mytempdata2))))
row.names(mytempdata1)<-NULL
mytempdata1<-as.data.frame(mytempdata1)
names(mytempdata1)<-c("V1","V2","V3")
mytempdata1$V2<-as.numeric(as.character(mytempdata1$V2))
mytempdata1$V3<-as.character(mytempdata1$V3)

pdf("~/Dropbox/LDLD/figs/revisions/vio_imbalance_ratio_ggplot.pdf")
p <- ggplot(mytempdata1, aes(x=V1,y=V2,fill=V3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
p + geom_violin(scale="area") #+   geom_boxplot(width = 0.2)
dev.off()


require(vioplot)
pdf("~/Dropbox/LDLD/figs/vio_imbalance.pdf")
vioplot((myt[,2]/myt[,1])[!is.na(myt[,2]/myt[,1])],col=rgb(1,0,0,0.5),drawRect=F)
vioplot((mynot[,2]/mynot[,1])[!is.na(mynot[,2]/mynot[,1])],col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
dev.off()


#for i in /mnt/454/HGDP/genomes_Bteam/HG19Align/*.bam /mnt/sequencedb/11men/martin/hg19_1000g/BAM/*.bam ; do /r1/people/pruefer/bin/BamTable2 -r ~/Dropbox/LDLD/ms/Ajhg/TableS2.bed.txt;done

}
#INTERGENIC
{
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5"
mydata<-read.table("~/Dropbox/LDLD/ms/GBE/TableS6.bed",stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
createbedlogo(mydata,paste0(myfolder,"/basessign"))

for MYCHR in `seq 1 22`;do
tabix /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -R ~/Dropbox/LDLD/ms/GBE/TableS6.bed | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 > temp.bed
bedtools intersect -a temp.bed -b ~/Dropbox/LDLD/ms/GBE/TableS6.bed >> /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed
done
paste /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign > ttemp; mv ttemp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed

mydata<-read.table("/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed",stringsAsFactors =FALSE,header=FALSE)
#
#paste <( zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab |  awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}' | sort -Vu -k1,1 -k2,2 ) ) /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign | awk '{if (length($4)==1 && length($5)==1){print}}' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basessign.bed

bamlogo_sign<-bamlogo(paste0(myfolder,"/basessign.bed"),tilltheend=TRUE,till=100,startfrom=1,filetemp="temp")
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
for myrepl in `seq 1 10`;do 
bedtools subtract -a <( cat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl${myrepl}/min5/*tab | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5,$7}'  ) -b ~/Dropbox/LDLD/ms/GBE/TableS6.bed | shuf -n 20000 | sort -Vu -k1,1 -k2,2n >> /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/notsign.bed 
done
cat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/notsign.bed | shuf -n 20000 | sort -Vu -k1,1 -k2,2n > ttemp; mv ttemp /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/notsign.bed

source("~/Dropbox/LDLD/analysesLDLD_header.R") 
myfolder<-"/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5"
mydata<-read.table(paste0(myfolder,"/notsign.bed"),stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
createbedlogo(mydata,paste0(myfolder,"/basesnotsign"))
#currently running
#to do
paste /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/notsign.bed /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basesnotsign > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign.bed
bamlogo_notsign<-bamlogo("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/basesnotsign.bed",tilltheend=TRUE,till=100,startfrom=1,filetemp="ttempnot")
save(bamlogo_notsign,file=paste0(myfolder,"bamlogo_notsign.RData"))
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

mynot<-t(apply(bamlogo_notsign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))
myt<-t(apply(bamlogo_sign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))

require(vioplot)
pdf("~/Dropbox/LDLD/figs/revisions/vio_imbalance_min5inter.pdf")
vioplot((myt[,2]/myt[,1])[!is.na(myt[,2]/myt[,1])],col=rgb(1,0,0,0.5),drawRect=F)
vioplot((mynot[,2]/mynot[,1])[!is.na(mynot[,2]/mynot[,1])],col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
dev.off()

#plot as kay likes: ratio of reads for alternative allele
tempdata<-bamlogo_sign; 
tempdata<-tempdata[tempdata$alt=="A" |tempdata$alt=="C" | tempdata$alt=="G" | tempdata$alt=="T",]
tempdata<-tempdata[tempdata$ref=="A" |tempdata$ref=="C" | tempdata$ref=="G" | tempdata$ref=="T",]
temp<-sapply(1:length(tempdata$alt), function(x) tempdata[[tempdata$alt[x]]][x] )/(sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$alt[x]]][x] )+sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$ref[x]]][x] ))
tempsign<-temp[!is.na(temp)]
tempdata<-bamlogo_notsign; 
tempdata<-tempdata[tempdata$alt=="A" |tempdata$alt=="C" | tempdata$alt=="G" | tempdata$alt=="T",]
tempdata<-tempdata[tempdata$ref=="A" |tempdata$ref=="C" | tempdata$ref=="G" | tempdata$ref=="T",]
temp<-sapply(1:length(tempdata$alt), function(x) tempdata[[tempdata$alt[x]]][x] )/(sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$alt[x]]][x] )+sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$ref[x]]][x] ))
tempnosign<-temp[!is.na(temp)]

require(vioplot)
pdf("~/Dropbox/LDLD/figs/revisions/vio_imbalance_ratio_intermin5.pdf")
vioplot(tempsign,col=rgb(1,0,0,0.5),drawRect=F,ylim=c(0,1))
vioplot(tempnosign,col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
abline(h=0.5,lty=2)
dev.off()

mycol1=rgb(1,0,0,0.5)
mycol2=rgb(0.3,0.3,0.3,0.3)
mytempdata1<-tempsign
mytempdata2<-tempnosign
mytempdata1<-cbind(c(rep("sign",length(mytempdata1)),rep("nosign",length(mytempdata2))),c(mytempdata1,mytempdata2),c(rep(mycol1,length(mytempdata1)),rep(mycol2,length(mytempdata2))))
row.names(mytempdata1)<-NULL
mytempdata1<-as.data.frame(mytempdata1)
names(mytempdata1)<-c("V1","V2","V3")
mytempdata1$V2<-as.numeric(as.character(mytempdata1$V2))
mytempdata1$V3<-as.character(mytempdata1$V3)

pdf("~/Dropbox/LDLD/figs/revisions/vio_imbalance_ratio_intermin5_ggplot.pdf")
p <- ggplot(mytempdata1, aes(x=V1,y=V2,fill=V3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
p + geom_violin(scale="area") #+   geom_boxplot(width = 0.2)
dev.off()


require(vioplot)
pdf("~/Dropbox/LDLD/figs/vio_imbalance.pdf")
vioplot((myt[,2]/myt[,1])[!is.na(myt[,2]/myt[,1])],col=rgb(1,0,0,0.5),drawRect=F)
vioplot((mynot[,2]/mynot[,1])[!is.na(mynot[,2]/mynot[,1])],col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
dev.off()

}
}
}
#check info samples
{
info1000g<-read.table("/mnt/scratch/fabrizio/LDLD/20130606_sample_info.txt",header=TRUE,sep='\t',stringsAsFactors=F)
samples1000g<-system("zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | 
head -300 | grep '#CHROM' | head -1 | awk '{for (i=10; i<=NF; i++) print $i}'",intern=TRUE)
samplesabove95<-system("cat ~/workspace/1000genomes/above95.unrelated.samples",intern=TRUE)
nspops<-as.numeric(system("cat /mnt/scratch/fabrizio/LDLD/nspops.txt",intern=TRUE))
samples_order1000g_above95<-intersect(samplesabove95,samples1000g)
mydata<-info1000g[match(samples_order1000g_above95,info1000g$Sample),]
mydata$Sample[!is.na(mydata$Total.LC.Sequence) & !is.na(mydata$Total.Exome.Sequence)] #the samples I analyzed all have Exome and LC. However, not all samples of the 1000 genome dataset.

mydata<-info1000g[match(samples1000g,info1000g$Sample),]
!is.na(mydata$Total.LC.Sequence) & !is.na(mydata$Total.Exome.Sequence) #the samples I analyzed all have Exome and LC. However, not all samples of the 1000 genome dataset. But all of the ones in vcfs.


}
#============================
#===== TESTS for GBE revisions - 12/04/18==
#============================

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

#background
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
mytempbed=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/temp${MYCHR}.bed
zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr${MYCHR}.vcf.gz  | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5}' | sort -k 2,2n | uniq > $mytempbed
tabix $my1000gfile -R $mytempbed | awk -v OFS='\t' '{print $1,$2-1,$2,$4,$5}' | sort -k 2,2n | uniq > $mybedfile
bedtools intersect -a $mybedfile -b $mytempbed > temp; mv temp $mybedfile
tabix $my1000gfile -R $mybedfile | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8}' > $temp1000gfile
count_n1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep missense | wc -l ) 
count_s1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep synonymous | wc -l ) 
count_nc1=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep coding | wc -l )
count_n2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep missense | grep -v synonymous | wc -l ) 
count_s2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep synonymous | wc -l ) 
count_nc2=$( bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep coding | wc -l )
#bedtools intersect -b $temp1000gfile -a $mybedfile -wo | awk '{if ($5==$11){print}}' | sort -k2,2n | uniq | grep -v missense | grep -v synonymous | grep -v coding >> revisions/bed_ns/others
count_tot1=$( cat ~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/chr${MYCHR}.bed  |wc -l  )
count_n=$(( $count_n + $count_n1 )); count_s=$(( $count_s + $count_s1 )); count_nc=$(( $count_nc + $count_nc1 ));count_tot=$(( $count_tot + $count_tot1 ))
count_n20=$(( $count_n20 + $count_n2 )); count_s20=$(( $count_s20 + $count_s2 )); count_nc20=$(( $count_nc20 + $count_nc2 ));count_tot20=$(( $count_tot20 + $count_tot1 ))
done
echo $count_n $count_s $count_nc $count_tot 
echo $count_n20 $count_s20 $count_nc20 $count_tot20 
done
#39534 40284 5266 87702
#38194 40284 5266 87702
#115424 96669 12797 231493
#111920 96669 12797 231493

obs_min5<-c(461/(461+166),166/(461+166))
back_min5<-c(38194/(38194+40284),40284/(38194+40284))


sqrt(back_min5[1]*back_min5[2]/(38194+40284)) #0.001784194

plot(obs_min5[1],ylim=c(0,1),pch=19,col="red",cex=3)
points(back_min5[1],ylim=c(0,1),pch=19,col="black",cex=3)




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


#for coding
for mymap in 50 99;do
for mysupp in 4 5;do
count_n10=0;count_n20=0;count_n1=0;count_n2=0;count_n1tot=0;count_n2tot=0;
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
for MYCHR in `seq 1 22`;do
echo $MYCHR
mybedfile=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr${MYCHR}.vcf.gz
source ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g_2awk.sh
mymappabilityfile=/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_${mymap}.bed.gz
#tabix $my1000gfile -R $mymappabilityfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
count_n1=$( bedtools intersect -a $mymappabilityfile -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l )
count_n2=$( bedtools intersect -a $mymappabilityfile -b $mybedfile | wc -l )
count_n1tot=$(( $count_n1tot + $count_n1 )); 
count_n2tot=$(( $count_n2tot + $count_n2 )); 
count_n10=$( cat ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | grep -v '#' | wc -l )
count_n20=$( zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr*.vcf.gz | grep -v '#' | wc -l )
done
echo "n overlaps with map in TableS" $mysupp "with map" $mymap " are:"  $count_n1 $count_n10  $count_n2 $count_n20 >> $mylogfile.$mymap.coding
done
done


#background for coding
mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_map
for mymap in 50 99;do
for mysupp in 4 5;do
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
mybedfile=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr${MYCHR}.vcf.gz
source ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g_2awk.sh
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
mymappabilityfile=/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_${mymap}.bed.gz
#tabix $my1000gfile -R $mymappabilityfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
bedtools intersect -a $my1000gfile -b <( tabix $mymappabilityfile $MYCHR ) -sorted | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
count_n2=$(cat $mybedfile | wc -l )
count_n1=$( bedtools intersect -a $mybedfile -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l )
count_n=$(( $count_n + $count_n1 ));
count_n20=$(( $count_n20 + $count_n2 )); 
echo "count is " $count_n2 >> ${mylogfile}.back
done
echo "n overlaps with map in TableS" $mysupp "with map" $mymap " are" $count_n $count_n20 >> ${mylogfile}.${mymap}.back
done
done





#tot, also out of inbreeding track
}
#repeats
{

zcat /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mask.freq1.chr*.bed.gz | shuf -n 20000 | sort -Vu -k1,1n -k2,2n | gzip -f > /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min1.bed.gz
zcat /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mask.freq5.chr*.bed.gz | shuf -n 20000 | sort -Vu -k1,1n -k2,2n | gzip -f > /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min5.bed.gz
cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr*tab | grep -v '#' | shuf -n 20000 | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -Vu -k1,1n -k2,2n | gzip -f > /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min5.bed.gz
cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1/chr*tab | grep -v '#' | shuf -n 20000 | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -Vu -k1,1n -k2,2n | gzip -f > /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min1.bed.gz
fabrizio_mafessoni@bionc02:/mnt/sequencedb/1000Genomes/ftp/phase3/20130502/supporting/data/CHB£ zcat /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min5.bed.gz | wc -l
19804
fabrizio_mafessoni@bionc02:/mnt/sequencedb/1000Genomes/ftp/phase3/20130502/supporting/data/CHB£ zcat /home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min1.bed.gz | wc -l
19584

#total
for mysupp in `seq 4 7`;do
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_repeats/tot_repeat_S${mysupp}
bedtools intersect -a $mymappabilityfile -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | sort -Vu -k 1,1 -k 2,2n  | wc -l
done;
#38
#44
#11218
#4268
background1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min1.bed.gz
background5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min5.bed.gz
background_coding1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min1.bed.gz
background_coding5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min5.bed.gz
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
bedtools intersect -a $mymappabilityfile -b $background5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #10618
bedtools intersect -a $mymappabilityfile -b $background1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #10667
bedtools intersect -a $mymappabilityfile -b $background_coding5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #608
bedtools intersect -a $mymappabilityfile -b $background_coding1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #633

#Simple_repeat
for mysupp in `seq 4 7`;do
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_repeats/tot_repeat_S${mysupp}
bedtools intersect -a <( cat $mymappabilityfile | grep "Simple_repeat" )  -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | sort -Vu -k 1,1 -k 2,2n  | wc -l
done;
32
34
221
181
background1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min1.bed.gz
background5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min5.bed.gz
background_coding1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min1.bed.gz
background_coding5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min5.bed.gz
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
bedtools intersect -a <( cat $mymappabilityfile | grep "Simple_repeat" ) -b $background5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #377
bedtools intersect -a <( cat $mymappabilityfile | grep "Simple_repeat" ) -b $background1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #349
bedtools intersect -a <( cat $mymappabilityfile | grep "Simple_repeat" ) -b $background_coding5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #192
bedtools intersect -a <( cat $mymappabilityfile | grep "Simple_repeat" ) -b $background_coding1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #192

#Low_complexity
for mysupp in `seq 4 7`;do
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_repeats/tot_repeat_S${mysupp}
bedtools intersect -a <( cat $mymappabilityfile | grep "Low_complexity" )  -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | sort -Vu -k 1,1 -k 2,2n  | wc -l
done;
1
5
76
51
background1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min1.bed.gz
background5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min5.bed.gz
background_coding1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min1.bed.gz
background_coding5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min5.bed.gz
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
bedtools intersect -a <( cat $mymappabilityfile | grep "Low_complexity" ) -b $background5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #152
bedtools intersect -a <( cat $mymappabilityfile | grep "Low_complexity" ) -b $background1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #163
bedtools intersect -a <( cat $mymappabilityfile | grep "Low_complexity" ) -b $background_coding5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #124
bedtools intersect -a <( cat $mymappabilityfile | grep "Low_complexity" ) -b $background_coding1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #118

#SINE
for mysupp in `seq 4 7`;do
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_repeats/tot_repeat_S${mysupp}
bedtools intersect -a <( cat $mymappabilityfile | grep "SINE" )  -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | sort -Vu -k 1,1 -k 2,2n  | wc -l
done;
1
1
1316
376
background1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min1.bed.gz
background5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min5.bed.gz
background_coding1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min1.bed.gz
background_coding5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min5.bed.gz
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
bedtools intersect -a <( cat $mymappabilityfile | grep "SINE" ) -b $background5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #3264
bedtools intersect -a <( cat $mymappabilityfile | grep "SINE" ) -b $background1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #3227
bedtools intersect -a <( cat $mymappabilityfile | grep "SINE" ) -b $background_coding5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #128
bedtools intersect -a <( cat $mymappabilityfile | grep "SINE" ) -b $background_coding1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #122

#LINE
for mysupp in `seq 4 7`;do
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
mybedfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_repeats/tot_repeat_S${mysupp}
bedtools intersect -a <( cat $mymappabilityfile | grep "SINE" )  -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | sort -Vu -k 1,1 -k 2,2n  | wc -l
done;
1
1
1316
376
background1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min1.bed.gz
background5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background_min5.bed.gz
background_coding1=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min1.bed.gz
background_coding5=/home/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/background.coding_min5.bed.gz
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
bedtools intersect -a <( cat $mymappabilityfile | grep "LINE" ) -b $background5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #4154
bedtools intersect -a <( cat $mymappabilityfile | grep "LINE" ) -b $background1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #4194
bedtools intersect -a <( cat $mymappabilityfile | grep "LINE" ) -b $background_coding5 | sort -Vu -k 1,1 -k 2,2n  | wc -l #73
bedtools intersect -a <( cat $mymappabilityfile | grep "LINE" ) -b $background_coding1 | sort -Vu -k 1,1 -k 2,2n  | wc -l #89


#from parser of 1000 genomes for my populations ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g.sh
mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_repeats/tot_repeat
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
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
#tabix $my1000gfile -R $mymappabilityfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
bedtools intersect -a $my1000gfile -b <( cat $mymappabilityfile | awk -v mychr=${MYCHR} '{if ($1==mychr){print}}' ) -sorted | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
count_n2=$(cat $mybedfile | wc -l )
count_n1=$( bedtools intersect -a $mybedfile -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l )
count_n=$(( $count_n + $count_n1 ));
count_n20=$(( $count_n20 + $count_n2 )); 
echo "count is " $count_n2 >> $mylogfile
done
echo "n overlaps with map in TableS" $mysupp "with map" $mymap " are" $count_n $count_n20 >> $mylogfile
done


mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_repeats/tot_repeat
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
mymappabilityfile=/mnt/scratch/kay/forFab/hg19RepeatMasker.bed
#tabix $my1000gfile -R $mymappabilityfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
bedtools intersect -a $my1000gfile -b <( cat $mymappabilityfile | awk -v mychr=${MYCHR} '{if ($1==mychr){print}}' ) -sorted | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
count_n2=$(cat $mybedfile | wc -l )
count_n1=$( bedtools intersect -a $mybedfile -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l )
count_n=$(( $count_n + $count_n1 ));
count_n20=$(( $count_n20 + $count_n2 )); 
echo "count is " $count_n2 >> $mylogfile
done
echo "n overlaps with map in TableS" $mysupp "with map" $mymap " are" $count_n $count_n20 >> $mylogfile
done



#for coding
for mymap in 50 99;do
for mysupp in 4 5;do
count_n10=0;count_n20=0;count_n1=0;count_n2=0;count_n1tot=0;count_n2tot=0;
myminfreq=1
if [ $mysupp -eq 6 ];then
    myminfreq=5
fi
if [ $mysupp -eq 4 ];then
    myminfreq=5
fi
for MYCHR in `seq 1 22`;do
echo $MYCHR
mybedfile=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr${MYCHR}.vcf.gz
source ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g_2awk.sh
mymappabilityfile=/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_${mymap}.bed.gz
#tabix $my1000gfile -R $mymappabilityfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
count_n1=$( bedtools intersect -a $mymappabilityfile -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l )
count_n2=$( bedtools intersect -a $mymappabilityfile -b $mybedfile | wc -l )
count_n1tot=$(( $count_n1tot + $count_n1 )); 
count_n2tot=$(( $count_n2tot + $count_n2 )); 
count_n10=$( cat ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | grep -v '#' | wc -l )
count_n20=$( zcat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr*.vcf.gz | grep -v '#' | wc -l )
done
echo "n overlaps with map in TableS" $mysupp "with map" $mymap " are:"  $count_n1 $count_n10  $count_n2 $count_n20 >> $mylogfile.$mymap.coding
done
done


#background for coding
mylogfile=~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_map
for mymap in 50 99;do
for mysupp in 4 5;do
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
mybedfile=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min${myminfreq}/chr${MYCHR}.vcf.gz
source ~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_1000g_2awk.sh
my1000gfile=/mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
mymappabilityfile=/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_${mymap}.bed.gz
#tabix $my1000gfile -R $mymappabilityfile | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
bedtools intersect -a $my1000gfile -b <( tabix $mymappabilityfile $MYCHR ) -sorted | reorder_samples_vcf | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ~/Dropbox/LDLD/nspops.txt -f 0.0${myminfreq} | awk -v OFS='\t' '{print $1,$2-1,$2}' | sort -k 1,1 -k 2,2n | uniq | gzip -f > $mybedfile
count_n2=$(cat $mybedfile | wc -l )
count_n1=$( bedtools intersect -a $mybedfile -b ~/Dropbox/LDLD/ms/GBE/TableS${mysupp}.bed | wc -l )
count_n=$(( $count_n + $count_n1 ));
count_n20=$(( $count_n20 + $count_n2 )); 
echo "count is " $count_n2 >> ${mylogfile}.back
done
echo "n overlaps with map in TableS" $mysupp "with map" $mymap " are" $count_n $count_n20 >> ${mylogfile}.${mymap}.back
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

mydata<-read.table("~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_accessibilitymap")
mydata<-mydata[mydata$V2!="N",]

datatemp<-cbind(mydata$V3[mydata$V1==4]/sum(mydata$V3[mydata$V1==4]),
mydata$V4[mydata$V1==4]/sum(mydata$V4[mydata$V1==4]),mydata$V3[mydata$V1==6]/sum(mydata$V3[mydata$V1==6]),
mydata$V4[mydata$V1==6]/sum(mydata$V4[mydata$V1==6]))

temp<-datatemp[2,];datatemp[2,]<-datatemp[1,];datatemp[1,]<-temp

pdf("~/Dropbox/LDLD/figs/accessibility_1000genomes.pdf")
barplot(datatemp,
col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()

pdf("~/Dropbox/LDLD/figs/revisions/accessibility_1000genomes_Pie_coding.pdf")
slices <- datatemp[,1]
lbls <- c("P", "H", "L", "ZQ")
pie(slices, labels = lbls, main="Accessibility Map",col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()
pdf("~/Dropbox/LDLD/figs/revisions/accessibility_1000genomes_Pie_intergenic.pdf")
slices <- datatemp[,3]
lbls <- c("P", "H", "L", "ZQ")
pie(slices, labels = lbls, main="Accessibility Map",col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()
pdf("~/Dropbox/LDLD/figs/revisions/accessibility_1000genomes_Pie_coding_back.pdf")
slices <- datatemp[,2]
lbls <- c("P", "H", "L", "ZQ")
pie(slices, labels = lbls, main="Accessibility Map",col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()
pdf("~/Dropbox/LDLD/figs/revisions/accessibility_1000genomes_Pie_intergenic_back.pdf")
slices <- datatemp[,4]
lbls <- c("P", "H", "L", "ZQ")
pie(slices, labels = lbls, main="Accessibility Map",col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()

mydata<-read.table("~/Dropbox/LDLD/ms/GBE/revisions/bed_ns/mylog_accessibilitymap_pilot")
mydata<-mydata[mydata$V2!="N",]

datatemp<-cbind(mydata$V3[mydata$V1==4]/sum(mydata$V3[mydata$V1==4]),
mydata$V4[mydata$V1==4]/sum(mydata$V4[mydata$V1==4]),mydata$V3[mydata$V1==6]/sum(mydata$V3[mydata$V1==6]),
mydata$V4[mydata$V1==6]/sum(mydata$V4[mydata$V1==6]))

temp<-datatemp[2,];datatemp[2,]<-datatemp[1,];datatemp[1,]<-temp

pdf("~/Dropbox/LDLD/figs/accessibility_1000genomes_pilot.pdf")
barplot(datatemp,
col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()

pdf("~/Dropbox/LDLD/figs/revisions/accessibility_1000genomes_pilotPie_coding.pdf")
slices <- datatemp[,1]
lbls <- c("P", "H", "L", "ZQ")
pie(slices, labels = lbls, main="Accessibility Map",col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()
pdf("~/Dropbox/LDLD/figs/revisions/accessibility_1000genomes_pilotPie_intergenic.pdf")
slices <- datatemp[,3]
lbls <- c("P", "H", "L", "ZQ")
pie(slices, labels = lbls, main="Accessibility Map",col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()
pdf("~/Dropbox/LDLD/figs/revisions/accessibility_1000genomes_pilotPie_coding_back.pdf")
slices <- datatemp[,2]
lbls <- c("P", "H", "L", "ZQ")
pie(slices, labels = lbls, main="Accessibility Map",col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()
pdf("~/Dropbox/LDLD/figs/revisions/accessibility_1000genomes_pilotPie_intergenic_back.pdf")
slices <- datatemp[,2]
lbls <- c("P", "H", "L", "ZQ")
pie(slices, labels = lbls, main="Accessibility Map",col=c("cornsilk","gray29","gray59","gray79")
)
dev.off()

barplot(
cbind(mydata$V3[mydata$V1==7]/sum(mydata$V3[mydata$V1==7]),
mydata$V4[mydata$V1==7]/sum(mydata$V4[mydata$V1==7]),
mydata$V3[mydata$V1==5]/sum(mydata$V3[mydata$V1==5]),
mydata$V4[mydata$V1==5]/sum(mydata$V4[mydata$V1==5]))
)

}

#base composition & plot features
{
#CODING
{
#I found files /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5bamlogo_notsign.RData and min5bamlogo_sign.RData but not code that I used to generate them. I should look in my personal pc
myfolder="/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/"
setwd(myfolder)
load("min5bamlogo_sign.RData")
#load(paste0(myfolder,"bamlogo_sign.RData"))
row.names(bamlogo_sign)<-NULL
bamlogo_sign<-as.data.frame(bamlogo_sign)
names(bamlogo_sign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
#bamlogo_sign<-bamlogo_sign[bamlogo_sign$A!="NaN",]
bamlogo_sign$A<-as.numeric(as.character(bamlogo_sign$A))
bamlogo_sign$C<-as.numeric(as.character(bamlogo_sign$C))
bamlogo_sign$G<-as.numeric(as.character(bamlogo_sign$G))
bamlogo_sign$T<-as.numeric(as.character(bamlogo_sign$T))
bamlogo_sign$tot_second_call<-as.numeric(as.character(bamlogo_sign$tot_second_call))
bamlogo_sign$ref<-as.character(bamlogo_sign$ref)
bamlogo_sign$alt<-as.character(bamlogo_sign$alt)

load("min5bamlogo_notsign.RData")
#double check allelic imbalance
#load(paste0(myfolder,"bamlogo_notsign.RData"))
row.names(bamlogo_notsign)<-NULL
bamlogo_notsign<-as.data.frame(bamlogo_notsign)
names(bamlogo_notsign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
#bamlogo_notsign<-bamlogo_notsign[bamlogo_notsign$A!="NaN",]
bamlogo_notsign$A<-as.numeric(as.character(bamlogo_notsign$A))
bamlogo_notsign$C<-as.numeric(as.character(bamlogo_notsign$C))
bamlogo_notsign$G<-as.numeric(as.character(bamlogo_notsign$G))
bamlogo_notsign$T<-as.numeric(as.character(bamlogo_notsign$T))
bamlogo_notsign$tot_second_call<-as.numeric(as.character(bamlogo_notsign$tot_second_call))
bamlogo_notsign$ref<-as.character(bamlogo_notsign$ref)
bamlogo_notsign$alt<-as.character(bamlogo_notsign$alt)

#NS
obs_min5<-c(461/(461+166),166/(461+166))
back_min5<-c(38194/(38194+40284),40284/(38194+40284))
#CI 
obs_CI=1.96*sqrt( obs_min5[1]*obs_min5[2]/(461+166) ) #0.03453498
back_CI=1.96*sqrt(back_min5[1]*back_min5[2]/(38194+40284)) #0.003497019

#transition/transversions
mytabsign<-table(apply(cbind(bamlogo_sign$ref,bamlogo_sign$alt), MARGIN=1,FUN=function(x) paste0(x,collapse='')))#/length(bamlogo_sign$ref)
mytabnosign<-table(apply(cbind(bamlogo_notsign$ref,bamlogo_notsign$alt), MARGIN=1,FUN=function(x) paste0(x,collapse='')))#/length(bamlogo_notsign$ref)
myval<-(mytabsign[["CT"]]+mytabsign[["TC"]]+mytabsign[["AG"]]+mytabsign[["GA"]])/sum(mytabsign) #0.3731343
myval0<-(mytabnosign[["CT"]]+mytabnosign[["TC"]]+mytabnosign[["AG"]]+mytabnosign[["GA"]])/sum(mytabnosign) #0.7581473
#CI 
myval_CI=1.96*sqrt( myval*(1-myval)/sum(mytabsign) ) #0.0366217
myval0_CI=1.96*sqrt( myval0*(1-myval0)/sum(mytabnosign)) #0.008591438




mean(sapply(as.character(bamlogo_sign$seq5bp), function(x) length(unique(strsplit(x,"")[[1]])) ))
mean(sapply(as.character(bamlogo_notsign$seq5bp), function(x) length(unique(strsplit(x,"")[[1]])) ))

myf=function (x) strsplit(x,"")[[1]][6]==strsplit(x,"")[[1]][5] || strsplit(x,"")[[1]][6]==strsplit(x,"")[[1]][7] #repeat in middle
mean(sapply(as.character(bamlogo_sign$seq5bp), function(x) myf(x))) #0.2863962 
mean(sapply(as.character(bamlogo_notsign$seq5bp), function(x) myf(x))) #0.3741678

mydata<-cbind(
table(bamlogo_sign$alt)/sum(table(bamlogo_sign$alt)),
table(bamlogo_notsign$alt)/sum(table(bamlogo_notsign$alt))
)
mydata<-mydata[c(2,3,4,1),]
#       [,1]      [,2]
#C 0.3656716 0.1977366
#G 0.3388060 0.2032904
#T 0.1507463 0.2892172
#A 0.1447761 0.3097558

pdf("~/Dropbox/LDLD/figs/revisions/base_comp_min5coding.pdf")
barplot(mydata,col=c("lightcyan3","lightcyan2","wheat3","wheat1")
)
dev.off()
temp<-round(mydata*1000)
cat(c(rep(">\nC\n",temp[1,1]),rep("G",temp[2,1]),rep("T",temp[3,1]),rep("A",temp[4,1])))

#base=alt nearby
myf=function (x,mydata) as.numeric(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] || mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7])
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.5059666
mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.3506784

#bases=alt around
myf=function (x,mydata) as.numeric(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7])
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.1193317
mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.03846661

#dimers before
myf=function (x,mydata) as.numeric(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4])
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.1288783
mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.05282571

#dimers after
myf=function (x,mydata) as.numeric(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8])
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.1551313
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.05727924

#dimer before or after
myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8])
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5])
)
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.2208955
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.1014925

#dimers nearby or around
myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8])
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5])
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7])
)
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.2720764
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.1264916

myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8] && mydata$alt[x]=="G")
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && mydata$alt[x]=="G")
)
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.09785203
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.02386635


myf=function (x,mydata) as.numeric(
( (mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] || mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] ) && ( mydata$alt[x]=="G" || mydata$alt[x]=="C" ) )
)
mydata<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.3866348 
mydata0<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.1622912

myf=function (x,mydata) as.numeric(
( (mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] || mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] ) && ( mydata$alt[x]=="A" || mydata$alt[x]=="T" ) )
)
mydata<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.1193317
mydata0<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.2052506

myf=function (x,mydata) as.numeric(
( (mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] || mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] )  )
)
mydata<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.5059666
mydata0<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.3675418

#GG+G and CC+C
myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8] && ( mydata$alt[x]=="G" || mydata$alt[x]=="C" ) )
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && ( mydata$alt[x]=="G" || mydata$alt[x]=="C" ) )
)
mydata<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.1895522
mydata0<-mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.05671642
#CI 
mydata_CI=1.96*sqrt( mydata*(1-mydata)/nrow(bamlogo_sign) ) #0.02967875
mydata0_CI=1.96*sqrt( mydata0*(1-mydata0)/nrow(bamlogo_notsign)) #0.00464076


pdf("~/Dropbox/LDLD/figs/revisions/NS_and_.pdf")
plot(cbind(c(0.2,0.2),c(obs_min5[1],back_min5[1])),xlim=c(0,1),ylim=c(0,1),pch=19,col=c("red","black"),cex=3)
#points(cbind(0.3,back_min5[1]),ylim=c(0,1),pch=19,col="black",cex=3)

points(cbind(0.4,1-myval),ylim=c(0,1),pch=19,col="red",cex=3)
points(cbind(0.4,1-myval0),ylim=c(0,1),pch=19,col="black",cex=3)

points(cbind(0.6,mydata),ylim=c(0,1),pch=19,col="red",cex=3)
points(cbind(0.6,mydata0),ylim=c(0,1),pch=19,col="black",cex=3)

dev.off()

res_allelic_imbalance<-unlist(apply(myt,MARGIN=1,FUN=function(y) if (!is.na(y[1]) && sum(y)>0){binom.test(x=round(y[1]),n=round(y[1]+y[2]),0.5)$p.value}))
res_allelic_imbalance_nosign<-unlist(apply(mynot,MARGIN=1,FUN=function(y) if (!is.na(y[1]) && sum(y)>0){binom.test(x=round(y[1]),n=round(y[1]+y[2]),0.5)$p.value}))
sum(res_allelic_imbalance<0.05)/length(res_allelic_imbalance)
sum(res_allelic_imbalance_nosign<0.05)/length(res_allelic_imbalance_nosign)

pdf("~/Dropbox/LDLD/figs/revisions/barplot_features.pdf")
barplot(c(obs_min5[1],back_min5[1],1-myval,1-myval0,0.1895522,0.05671642,0.6028777,0.008936279),col=c("coral1","gray79"),space=c(1,0.1,1,0.1,1,0.1,1,0.1),ylim=c(0,1))
dev.off()

source("~/Dropbox/general_utils/general_functions.R")
#> length(mydata_sign$P_HET_EXCESS)
#[1] 695
#> length(mydata_nosign$P_HET_EXCESS)
#[1] 85718
hwe_CI=1.96*sqrt( 0.6028777*(1-0.6028777)/695 )
hwe0_CI=1.96*sqrt( 0.008936279*(1-0.008936279)/85718 )
mydata_temp<-rbind(c(obs_min5[1],1-myval,mydata,0.6028777),c(back_min5[1],1-myval0,mydata0,0.008936279))
mydata_temp_CI<-rbind(c(obs_CI,myval_CI,mydata_CI,hwe_CI),c(back_CI,myval0_CI,mydata0_CI,hwe0_CI))

pdf("~/Dropbox/LDLD/figs/revisions/barplot_features.pdf")
paired_barplot(data=mydata_temp,sdata=mydata_temp_CI,col=c("coral1","gray79"),myylim=c(0,1),poslegend="topright")
dev.off()



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
}

#INTERGENIC
{
#I found files /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5bamlogo_notsign.RData and min5bamlogo_sign.RData but not code that I used to generate them. I should look in my personal pc
myfolder="/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/"
setwd(myfolder)
load("min5bamlogo_sign.RData")
#load(paste0(myfolder,"bamlogo_sign.RData"))
row.names(bamlogo_sign)<-NULL
bamlogo_sign<-as.data.frame(bamlogo_sign)
names(bamlogo_sign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
#bamlogo_sign<-bamlogo_sign[bamlogo_sign$A!="NaN",]
bamlogo_sign$A<-as.numeric(as.character(bamlogo_sign$A))
bamlogo_sign$C<-as.numeric(as.character(bamlogo_sign$C))
bamlogo_sign$G<-as.numeric(as.character(bamlogo_sign$G))
bamlogo_sign$T<-as.numeric(as.character(bamlogo_sign$T))
bamlogo_sign$tot_second_call<-as.numeric(as.character(bamlogo_sign$tot_second_call))
bamlogo_sign$ref<-as.character(bamlogo_sign$ref)
bamlogo_sign$alt<-as.character(bamlogo_sign$alt)


load("min5bamlogo_notsign.RData")
#double check allelic imbalance
#load(paste0(myfolder,"bamlogo_notsign.RData"))
row.names(bamlogo_notsign)<-NULL
bamlogo_notsign<-as.data.frame(bamlogo_notsign)
names(bamlogo_notsign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
#bamlogo_notsign<-bamlogo_notsign[bamlogo_notsign$A!="NaN",]
bamlogo_notsign$A<-as.numeric(as.character(bamlogo_notsign$A))
bamlogo_notsign$C<-as.numeric(as.character(bamlogo_notsign$C))
bamlogo_notsign$G<-as.numeric(as.character(bamlogo_notsign$G))
bamlogo_notsign$T<-as.numeric(as.character(bamlogo_notsign$T))
bamlogo_notsign$tot_second_call<-as.numeric(as.character(bamlogo_notsign$tot_second_call))
bamlogo_notsign$ref<-as.character(bamlogo_notsign$ref)
bamlogo_notsign$alt<-as.character(bamlogo_notsign$alt)

mean(sapply(as.character(bamlogo_sign$seq5bp), function(x) length(unique(strsplit(x,"")[[1]])) ))
mean(sapply(as.character(bamlogo_notsign$seq5bp), function(x) length(unique(strsplit(x,"")[[1]])) ))

myf=function (x) strsplit(x,"")[[1]][6]==strsplit(x,"")[[1]][5] || strsplit(x,"")[[1]][6]==strsplit(x,"")[[1]][7] #repeat in middle
mean(sapply(as.character(bamlogo_sign$seq5bp), function(x) myf(x))) #0.4077099
mean(sapply(as.character(bamlogo_notsign$seq5bp), function(x) myf(x))) #0.4110821

mydata<-cbind(
table(bamlogo_sign$alt[nchar(bamlogo_sign$alt)==1])/sum(table(bamlogo_sign$alt[nchar(bamlogo_sign$alt)==1])),
table(bamlogo_notsign$alt[nchar(bamlogo_notsign$alt)==1])/sum(table(bamlogo_notsign$alt[nchar(bamlogo_notsign$alt)==1]))
)
mydata<-mydata[c(2,3,4,1),]
#> table(bamlogo_sign$alt[nchar(bamlogo_sign$alt)==1])/sum(table(bamlogo_sign$alt[nchar(bamlogo_sign$alt)==1]))
#       A         C         G         T 
#0.2540028 0.2387688 0.2379139 0.2693145 
#> table(bamlogo_notsign$alt[nchar(bamlogo_notsign$alt)==1])/sum(table(bamlogo_notsign$alt[nchar(bamlogo_notsign$alt)==1]))
#        A         C         G         T 
#0.2596986 0.2443732 0.2287913 0.2671369 
#also more Gs is what observed, all rest less

sum(nchar(bamlogo_sign$alt)!=1)/length(bamlogo_sign$alt) #0.0178626
sum(nchar(bamlogo_notsign$alt)!=1)/length(bamlogo_notsign$alt) #0.05249408

sum(nchar(bamlogo_sign$alt)==2)/length(bamlogo_sign$alt) #0.00870229
sum(nchar(bamlogo_notsign$alt)==2)/length(bamlogo_notsign$alt) #0.02296616

sum(nchar(bamlogo_sign$alt)==3)/length(bamlogo_sign$alt) 
sum(nchar(bamlogo_notsign$alt)==3)/length(bamlogo_notsign$alt) 
sum(nchar(bamlogo_sign$alt)==4)/length(bamlogo_sign$alt) 
sum(nchar(bamlogo_notsign$alt)==4)/length(bamlogo_notsign$alt) 
#less insertions/del at every length

pdf("~/Dropbox/LDLD/figs/revisions/base_comp_min5inter.pdf")
barplot(mydata,col=c("lightcyan3","lightcyan2","wheat3","wheat1")
)
dev.off()
temp<-round(mydata*1000)
cat(c(rep(">\nC\n",temp[1,1]),rep("G",temp[2,1]),rep("T",temp[3,1]),rep("A",temp[4,1])))

#base=alt nearby
myf=function (x,mydata) as.numeric(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] || mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7])
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.1384733
mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.1568747

#bases=alt around
myf=function (x,mydata) as.numeric(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7])
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.03045802
mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.0233307

#dimers before
myf=function (x,mydata) as.numeric(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4])
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.03320611
mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.02892035

#dimers after
myf=function (x,mydata) as.numeric(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8])
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.03305344
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.03183206

#dimer before or after
myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8])
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5])
)
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.0570229
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.05732824

#dimers nearby or around
myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8])
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5])
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7])
)
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.06519084
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.06725191

myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8] && mydata$alt[x]=="G")
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && mydata$alt[x]=="G")
)
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.01305344
mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.009236641


#GG+G and CC+C
myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8] && ( mydata$alt[x]=="G" || mydata$alt[x]=="C" ) )
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && ( mydata$alt[x]=="G" || mydata$alt[x]=="C" ) )
)
mydata<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.03053435
mydata0<-mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.0178626

mydata_CI=1.96*sqrt( mydata*(1-mydata)/nrow(bamlogo_sign) ) #0.03453498
mydata0_CI=1.96*sqrt( mydata0*(1-mydata0)/nrow(bamlogo_notsign) ) #0.003497019


#myf=function (x,mydata) as.numeric(
#( (mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] || mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] ) && ( mydata$alt[x]=="G" || mydata$alt[x]=="C" ) )
#)
#mydata<-mean(sapply(1:nrow(bamlogo_sign),FUN=function(x) myf(x,mydata=bamlogo_sign))) #0.07137405
#mydata0<-mean(sapply(1:nrow(bamlogo_notsign),FUN=function(x) myf(x,mydata=bamlogo_notsign))) #0.06679389
#c(mydata,mydata0)

#transition/transversions
mytabsign<-table(apply(cbind(bamlogo_sign$ref,bamlogo_sign$alt), MARGIN=1,FUN=function(x) paste0(x,collapse='')))#/length(bamlogo_sign$ref)
mytabnosign<-table(apply(cbind(bamlogo_notsign$ref,bamlogo_notsign$alt), MARGIN=1,FUN=function(x) paste0(x,collapse='')))#/length(bamlogo_notsign$ref)
myval<-(mytabsign[["CT"]]+mytabsign[["TC"]]+mytabsign[["AG"]]+mytabsign[["GA"]])/sum(mytabsign) #0.6155725
myval0<-(mytabnosign[["CT"]]+mytabnosign[["TC"]]+mytabnosign[["AG"]]+mytabnosign[["GA"]])/sum(mytabnosign) #0.611641
#CI 
myval_CI=1.96*sqrt( myval*(1-myval)/sum(mytabsign) ) #0.0366217
myval0_CI=1.96*sqrt( myval0*(1-myval0)/sum(mytabnosign)) #0.008591438

pdf("~/Dropbox/LDLD/figs/revisions/NS_and_.pdf")
plot(cbind(c(0.2,0.2),c(obs_min5[1],back_min5[1])),xlim=c(0,1),ylim=c(0,1),pch=19,col=c("red","black"),cex=3)
#points(cbind(0.3,back_min5[1]),ylim=c(0,1),pch=19,col="black",cex=3)

points(cbind(0.4,1-myval),ylim=c(0,1),pch=19,col="red",cex=3)
points(cbind(0.4,1-myval0),ylim=c(0,1),pch=19,col="black",cex=3)

points(cbind(0.6,mydata),ylim=c(0,1),pch=19,col="red",cex=3)
points(cbind(0.6,mydata0),ylim=c(0,1),pch=19,col="black",cex=3)

dev.off()

sum(mydata_sign$P_HET_EXCESS<0.05)/length(mydata_sign$P_HET_EXCESS) #0.1250839
sum(mydata_nosign$P_HET_EXCESS<0.05)/length(mydata_nosign$P_HET_EXCESS) #0.01512775

pdf("~/Dropbox/LDLD/figs/revisions/barplot_features_intergenic.pdf")
barplot(c(0,0,1-myval,1-myval0,0.03053435,0.0178626,0.1250839,0.01512775),col=c("coral1","gray79"),space=c(1,0.1,1,0.1,1,0.1,1,0.1),ylim=c(0,1))
dev.off()

source("~/Dropbox/general_utils/general_functions.R")
#length(mydata_sign$P_HET_EXCESS)  #16397
#length(mydata_nosign$P_HET_EXCESS) #19765
hwe_CI=1.96*sqrt( 0.1250839*(1-0.1250839)/16397 )
hwe0_CI=1.96*sqrt( 0.01512775*(1-0.01512775)/19765 )
mydata_temp<-rbind(c(0,1-myval,mydata,0.1250839),c(0,1-myval0,mydata0,0.01512775))
mydata_temp_CI<-rbind(c(0,myval_CI,mydata_CI,hwe_CI),c(0,myval0_CI,mydata0_CI,hwe0_CI))

pdf("~/Dropbox/LDLD/figs/revisions/barplot_features_intergenic.pdf")
paired_barplot(data=mydata_temp,sdata=mydata_temp_CI,col=c("coral1","gray79"),myylim=c(0,1),poslegend="topright")
dev.off()


plot(cbind(c(0.2,0.2),c(0,0)),xlim=c(0,1),ylim=c(0,1),pch=19,col=c("red","black"),cex=3)
points(cbind(0.3,back_min5[1]),ylim=c(0,1),pch=19,col="black",cex=3)

points(cbind(0.4,1-myval),ylim=c(0,1),pch=19,col="red",cex=3)
points(cbind(0.4,1-myval0),ylim=c(0,1),pch=19,col="black",cex=3)

points(cbind(0.6,mydata),ylim=c(0,1),pch=19,col="red",cex=3)
points(cbind(0.6,mydata0),ylim=c(0,1),pch=19,col="black",cex=3)


}

#GoNL
{
#create logo files and test reads imbalance

bedtools intersect -a /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/sign.vcf -b ~/Dropbox/LDLD/ms/GBE/TableS8.bed | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2}' > table.bed
myfolder="/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop"
setwd(myfolder)
mydata<-read.table(paste0(myfolder,"/table.bed"),stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
dim(mydata) #47193
createbedlogo(mydata,paste0(myfolder,"/basessign.temp"))
paste <( cat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/sign.vcf | grep -v '#' | awk -v OFS='\t'  '{print $1,$2-1,$2,$4,$5,$7}' ) basessign.temp | sed -e 's/\(.*\)/\U\1/' > basessign.bed
bamlogo_sign<-bamlogo(paste0(myfolder,"/basessign.bed"),tilltheend=TRUE,till=100,startfrom=1,filetemp="temp2")
save(bamlogo_sign,file=paste0(myfolder,"bamlogo_sign.RData"))

cat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2}' > table0.bed
myfolder="/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop"
setwd(myfolder)
mydata<-read.table(paste0(myfolder,"/table0.bed"),stringsAsFactors =FALSE,header=FALSE)
names(mydata)<-c("chr","start","end")
mydata$chr=paste0("chr",mydata$chr)
mydata$start<-as.numeric(mydata$start)
mydata$end<-as.numeric(mydata$end) 
dim(mydata) #47193
createbedlogo(mydata,paste0(myfolder,"/basesnosign.temp"))
paste <( cat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf | grep -v '#' | awk -v OFS='\t'  '{print $1,$2-1,$2,$4,$5,$7}' ) basesnosign.temp | sed -e 's/\(.*\)/\U\1/' > basesnosign.bed
bamlogo_nosign<-bamlogo(paste0(myfolder,"/basesnosign.bed"),tilltheend=TRUE,till=100,startfrom=1,filetemp="temp")
save(bamlogo_nosign,file=paste0(myfolder,"bamlogo_nosign.RData"))



#------
myfolder="/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop"
setwd(myfolder)
load(paste0(myfolder,"bamlogo_sign.RData"))
#load(paste0(myfolder,"bamlogo_sign.RData"))
row.names(bamlogo_sign)<-NULL
bamlogo_sign<-as.data.frame(bamlogo_sign)
names(bamlogo_sign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
#bamlogo_sign<-bamlogo_sign[bamlogo_sign$A!="NaN",]
bamlogo_sign$A<-as.numeric(as.character(bamlogo_sign$A))
bamlogo_sign$C<-as.numeric(as.character(bamlogo_sign$C))
bamlogo_sign$G<-as.numeric(as.character(bamlogo_sign$G))
bamlogo_sign$T<-as.numeric(as.character(bamlogo_sign$T))
bamlogo_sign$tot_second_call<-as.numeric(as.character(bamlogo_sign$tot_second_call))
bamlogo_sign$ref<-as.character(bamlogo_sign$ref)
bamlogo_sign$alt<-as.character(bamlogo_sign$alt)

load(paste0(myfolder,"bamlogo_nosign.RData"))
#double check allelic imbalance
#load(paste0(myfolder,"bamlogo_notsign.RData"))
row.names(bamlogo_nosign)<-NULL
bamlogo_notsign<-as.data.frame(bamlogo_nosign)
names(bamlogo_notsign)<-c("chr","init","end","ref","alt","seq5bp","A","C","G","T","tot_second_call","tot_samples")
#bamlogo_notsign<-bamlogo_notsign[bamlogo_notsign$A!="NaN",]
bamlogo_notsign$A<-as.numeric(as.character(bamlogo_notsign$A))
bamlogo_notsign$C<-as.numeric(as.character(bamlogo_notsign$C))
bamlogo_notsign$G<-as.numeric(as.character(bamlogo_notsign$G))
bamlogo_notsign$T<-as.numeric(as.character(bamlogo_notsign$T))
bamlogo_notsign$tot_second_call<-as.numeric(as.character(bamlogo_notsign$tot_second_call))
bamlogo_notsign$ref<-as.character(bamlogo_notsign$ref)
bamlogo_notsign$alt<-as.character(bamlogo_notsign$alt)

idsign<-unname(sapply(bamlogo_sign$ref, nchar)==1 & sapply(bamlogo_sign$alt, nchar)==1)
idnosign<-unname(sapply(bamlogo_notsign$ref, nchar)==1 & sapply(bamlogo_notsign$alt, nchar)==1)

#transition/transversions
mytabsign<-table(apply(cbind(bamlogo_sign$ref,bamlogo_sign$alt)[idsign,], MARGIN=1,FUN=function(x) paste0(x,collapse='')))#/length(bamlogo_sign$ref)
mytabnosign<-table(apply(cbind(bamlogo_notsign$ref,bamlogo_notsign$alt)[idnosign,], MARGIN=1,FUN=function(x) paste0(x,collapse='')))#/length(bamlogo_notsign$ref)
myval<-(mytabsign[["CT"]]+mytabsign[["TC"]]+mytabsign[["AG"]]+mytabsign[["GA"]])/sum(mytabsign) #0.3731343
myval0<-(mytabnosign[["CT"]]+mytabnosign[["TC"]]+mytabnosign[["AG"]]+mytabnosign[["GA"]])/sum(mytabnosign) #0.7581473
#CI 
myval_CI=1.96*sqrt( myval*(1-myval)/sum(mytabsign) ) #0.0366217
myval0_CI=1.96*sqrt( myval0*(1-myval0)/sum(mytabnosign)) #0.008591438




mean(sapply(as.character(bamlogo_sign$seq5bp), function(x) length(unique(strsplit(x,"")[[1]])) ))
mean(sapply(as.character(bamlogo_notsign$seq5bp), function(x) length(unique(strsplit(x,"")[[1]])) ))

myf=function (x) strsplit(x,"")[[1]][6]==strsplit(x,"")[[1]][5] || strsplit(x,"")[[1]][6]==strsplit(x,"")[[1]][7] #repeat in middle
mean(sapply(as.character(bamlogo_sign$seq5bp), function(x) myf(x))) #0.4576271
mean(sapply(as.character(bamlogo_notsign$seq5bp), function(x) myf(x))) #0.3772727

mydata<-cbind(
table(bamlogo_sign$alt[idsign])/sum(table(bamlogo_sign$alt[idsign])),
table(bamlogo_notsign$alt[idnosign])/sum(table(bamlogo_notsign$alt[idnosign]))
)
mydata<-mydata[c(2,3,4,1),]
#       [,1]      [,2]
#C 0.3656716 0.1977366
#G 0.3388060 0.2032904
#T 0.1507463 0.2892172
#A 0.1447761 0.3097558

pdf("~/Dropbox/LDLD/figs/revisions/base_comp_GoNL.pdf")
barplot(mydata,col=c("lightcyan3","lightcyan2","wheat3","wheat1")
)
dev.off()
temp<-round(mydata*1000)
cat(c(rep(">\nC\n",temp[1,1]),rep("G",temp[2,1]),rep("T",temp[3,1]),rep("A",temp[4,1])))
#GG+G and CC+C
myf=function (x,mydata) as.numeric(
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][7] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][8] && ( mydata$alt[x]=="G" || mydata$alt[x]=="C" ) )
||
(mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][4] && mydata$alt[x]==strsplit(as.character(mydata$seq5bp)[x],"")[[1]][5] && ( mydata$alt[x]=="G" || mydata$alt[x]=="C" ) )
)
mydata<-mean(sapply(1:nrow(bamlogo_sign[idsign,]),FUN=function(x) myf(x,mydata=bamlogo_sign[idsign,]))) #0.1895522
mydata0<-mean(sapply(1:nrow(bamlogo_notsign[idnosign,]),FUN=function(x) myf(x,mydata=bamlogo_notsign[idnosign,]))) #0.05671642
#CI 
mydata_CI=1.96*sqrt( mydata*(1-mydata)/nrow(bamlogo_sign) ) #0.02967875
mydata0_CI=1.96*sqrt( mydata0*(1-mydata0)/nrow(bamlogo_notsign)) #0.00464076

sum(mydata_sign$P_HET_EXCESS<0.05)/length(mydata_sign$P_HET_EXCESS) #0.9830508
sum(mydata_nosign$P_HET_EXCESS<0.05)/length(mydata_nosign$P_HET_EXCESS) #0.229108

source("~/Dropbox/general_utils/general_functions.R")
#> length(mydata_sign$P_HET_EXCESS)
#[1] 695
#> length(mydata_nosign$P_HET_EXCESS)
#[1] 85718
hwe_CI=1.96*sqrt( 0.9830508*(1-0.9830508)/118 )
hwe0_CI=1.96*sqrt( 0.129108*(1-0.129108)/1065 )
mydata_temp<-rbind(c(1-myval,mydata,0.9830508),c(1-myval0,mydata0,0.008936279))
mydata_temp_CI<-rbind(c(myval_CI,mydata_CI,hwe_CI),c(myval0_CI,mydata0_CI,hwe0_CI))

pdf("~/Dropbox/LDLD/figs/revisions/barplot_features_GoNL.pdf")
paired_barplot(data=mydata_temp,sdata=mydata_temp_CI,col=c("coral1","gray79"),myylim=c(0,1))
dev.off()


mydata_sign<-read.table("sign.vcf.hwe.hwe",header=T)
mydata_nosign<-read.table("nosign.vcf.hwe.hwe",header=T)
#pvalues: strong excess of hets
mean(mydata_sign$P_HET_EXCESS)
mean(mydata_nosign$P_HET_EXCESS)

#> length(mydata_sign$P_HET_EXCESS)
#[1] 695
#> length(mydata_nosign$P_HET_EXCESS)
#[1] 85718

sum(mydata_sign$P_HET_EXCESS<0.05)/length(mydata_sign$P_HET_EXCESS) #0.9830508
sum(mydata_nosign$P_HET_EXCESS<0.05)/length(mydata_nosign$P_HET_EXCESS) #0.229108

myhist<-hist(mydata_sign$P_HET_EXCESS,breaks=10)
myhist0<-hist(mydata_nosign$P_HET_EXCESS,breaks=10,add=T)

myhist$density<-myhist$count/sum(myhist$count)
myhist0$density<-myhist0$count/sum(myhist0$count)
myhist$count<-myhist$count/sum(myhist$count)
myhist0$count<-myhist0$count/sum(myhist0$count)
pdf("~/Dropbox/LDLD/figs/revisions/hetexcess_GoNL.pdf")
barplot(rbind(myhist$count,myhist0$count),beside=T,col=c(rgb(1,0,0,0.5),rgb(0.3,0.3,0.3,0.3)),space=c(0,0.1),ylim=c(0,1),xlab='p.value',ylab='density')
dev.off()
  
  
  
mynot<-t(apply(bamlogo_notsign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))
myt<-t(apply(bamlogo_sign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))


mean(myt[,2]/myt[,1],na.rm=T)
mean(mynot[,2]/mynot[,1],na.rm=T)

mynot<-t(apply(bamlogo_notsign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))
myt<-t(apply(bamlogo_sign[,7:10], MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)[c(1,2)]))

require(vioplot)
pdf("~/Dropbox/LDLD/figs/revisions/vio_imbalance_GoNL.pdf")
vioplot((myt[,2]/myt[,1])[!is.na(myt[,2]/myt[,1])],col=rgb(1,0,0,0.5),drawRect=F)
vioplot((mynot[,2]/mynot[,1])[!is.na(mynot[,2]/mynot[,1])],col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
dev.off()

#plot as kay likes: ratio of reads for alternative allele
tempdata<-bamlogo_sign[idsign,]
temp<-sapply(1:length(tempdata$alt), function(x) tempdata[[tempdata$alt[x]]][x] )/(sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$alt[x]]][x] )+sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$ref[x]]][x] ))
tempsign<-temp[!is.na(temp)]
tempdata<-bamlogo_notsign[idnosign,]
temp<-sapply(1:length(tempdata$alt), function(x) tempdata[[tempdata$alt[x]]][x] )/(sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$alt[x]]][x] )+sapply(1:length(tempdata$ref), function(x) tempdata[[tempdata$ref[x]]][x] ))
tempnosign<-temp[!is.na(temp)]

require(vioplot)
pdf("~/Dropbox/LDLD/figs/revisions/vio_imbalance_ratio_GoNL.pdf")
vioplot(tempsign,col=rgb(1,0,0,0.5),drawRect=F,ylim=c(0,1))
vioplot(tempnosign,col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
abline(h=0.5,lty=2)
dev.off()


}
}

#hardy-weinberg
{
#CODING
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

#> length(mydata_sign$P_HET_EXCESS)
#[1] 695
#> length(mydata_nosign$P_HET_EXCESS)
#[1] 85718

sum(mydata_sign$P_HET_EXCESS<0.05)/length(mydata_sign$P_HET_EXCESS) #0.6028777
sum(mydata_nosign$P_HET_EXCESS<0.05)/length(mydata_nosign$P_HET_EXCESS) #0.008936279

myhist<-hist(mydata_sign$P_HET_EXCESS,breaks=10)
myhist0<-hist(mydata_nosign$P_HET_EXCESS,breaks=10,add=T)

myhist$density<-myhist$count/sum(myhist$count)
myhist0$density<-myhist0$count/sum(myhist0$count)
myhist$count<-myhist$count/sum(myhist$count)
myhist0$count<-myhist0$count/sum(myhist0$count)
pdf("~/Dropbox/LDLD/figs/revisions/hetexcess_codingmin5.pdf")
barplot(rbind(myhist$count,myhist0$count),beside=T,col=c(rgb(1,0,0,0.5),rgb(0.3,0.3,0.3,0.3)),space=c(0,0.1),ylim=c(0,1),xlab='p.value',ylab='density')
dev.off()
        
require(vioplot)
pdf("~/Dropbox/LDLD/figs/hetexcess_vioplot.pdf")
vioplot(mydata_sign$P_HET_EXCESS,col=rgb(1,0,0,0.5),drawRect=F)
vioplot(mydata_nosign$P_HET_EXCESS,col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
dev.off()

as.numeric(sapply(mydata_sign[,4], function(x) strsplit(as.character(x),"/")[[1]][2]))/as.numeric(sapply(mydata_sign[,5], function(x) strsplit(as.character(x),"/")[[1]][2]))
as.numeric(sapply(mydata_nosign[,4], function(x) strsplit(as.character(x),"/")[[1]][2]))/as.numeric(sapply(mydata_nosign[,5], function(x) strsplit(as.character(x),"/")[[1]][2]))

1-as.numeric(sapply(mydata_sign[,4], function(x) strsplit(as.character(x),"/")[[1]][2]))/as.numeric(sapply(mydata_sign[,5], function(x) strsplit(as.character(x),"/")[[1]][2]))
1-as.numeric(sapply(mydata_nosign[,4], function(x) strsplit(as.character(x),"/")[[1]][2]))/as.numeric(sapply(mydata_nosign[,5], function(x) strsplit(as.character(x),"/")[[1]][2]))
}
#INTERGENIC
{
wc -l /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basesnotsign.bed
wc -l /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed

zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | head -2000 | grep '#' > /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf
for MYCHR in `seq 1 22`;do
bedtools intersect -a /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -b /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed >> /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf
done
wc -l /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf #20370
wc -l /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed #16405


#bedtools subtract -a <( cat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf | grep -v '#' | sort -Vu -k1,1 -k2,2n | uniq | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7}' ) -b /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed
#bedtools intersect -a <( cat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf | grep -v '#' | sort -Vu -k1,1 -k2,2n | uniq | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7}' ) -b /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed |wc -l
#these work well. So it is vcf format that sometimes does not work that well with bedtools -> 
cat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf | grep '#' > temp.vcf 
bedtools intersect -a <( cat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf | grep -v '#' | sort -Vu -k1,1 -k2,2n | uniq | awk -v OFS='\t' '{printf "%s\t%s\t%s\t",$1,$2-1,$2; for (i=3;i<NF;i++){printf "%s\t",$i}; printf "%s\n",$NF}' ) -b /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basessign.bed | awk -v OFS='\t' '{printf "%s\t%s\t",$1,$3; for (i=4;i<NF;i++){printf "%s\t",$i}; printf "%s\n",$NF}' >> temp.vcf 
gzip -f temp.vcf 
mv temp.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf.gz
vcftools --gzvcf /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.vcf.gz --hardy --out /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.HW

zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | head -2000 | grep '#' > /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf
for MYCHR in `seq 1 22`;do
echo $MYCHR
bedtools intersect -a /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -b /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basesnotsign.bed >> /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf
done
zcat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf.gz >> /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf
gzip -f /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf
zcat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf.gz | grep '#' > temp.vcf 
bedtools intersect -a <( zcat /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf.gz | grep -v '#' | sort -Vu -k1,1 -k2,2n | uniq | awk -v OFS='\t' '{printf "%s\t%s\t%s\t",$1,$2-1,$2; for (i=3;i<NF;i++){printf "%s\t",$i}; printf "%s\n",$NF}' ) -b /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/basesnotsign.bed | awk -v OFS='\t' '{printf "%s\t%s\t",$1,$3; for (i=4;i<NF;i++){printf "%s\t",$i}; printf "%s\n",$NF}' >> temp.vcf 
gzip -f temp.vcf 
mv temp.vcf.gz /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf.gz
vcftools --gzvcf /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.vcf.gz --hardy --out /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.HW

mydata_sign<-read.table("/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign.HW.hwe",header=TRUE)
mydata_nosign<-read.table("/mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/nosign.HW.hwe",header=TRUE)
#names(mydata_sign)<-c("CHR","init","end","OBS(HOM1/HET/HOM2)","E(HOM1/HET/HOM2)","ChiSq_HWE","P_HWE","P_HET_DEFICIT","P_HET_EXCESS")
#names(mydata_nosign)<-c("CHR","init","end","OBS(HOM1/HET/HOM2)","E(HOM1/HET/HOM2)","ChiSq_HWE","P_HWE","P_HET_DEFICIT","P_HET_EXCESS")
#pvalues: strong excess of hets
mean(mydata_sign$P_HET_DEFICIT) #0.1536404
mean(mydata_nosign$P_HET_DEFICIT) #0.1488438
mean(mydata_sign$P_HET_EXCESS) #0.8514218
mean(mydata_nosign$P_HET_EXCESS) #0.9258932

length(mydata_sign$P_HET_EXCESS)  #16397
length(mydata_nosign$P_HET_EXCESS) #19765

sum(mydata_sign$P_HET_EXCESS<0.05)/length(mydata_sign$P_HET_EXCESS) #0.1250839
sum(mydata_nosign$P_HET_EXCESS<0.05)/length(mydata_nosign$P_HET_EXCESS) #0.01512775
#pdf("~/Dropbox/LDLD/figs/revisions/hetexcess_intermin5.pdf")
#plot(myhist,freq=FALSE,col=rgb(1,0,0,0.5),ylim=c(0,0.8),main='',xlab='p-value')
#plot(myhist0,freq=FALSE,col=rgb(0.3,0.3,0.3,0.3),add=T)
#dev.off()
require(vioplot)

myhist<-hist(mydata_sign$P_HET_EXCESS,breaks=10)
myhist0<-hist(mydata_nosign$P_HET_EXCESS,breaks=10,add=T)

myhist$density<-myhist$count/sum(myhist$count)
myhist0$density<-myhist0$count/sum(myhist0$count)
myhist$count<-myhist$count/sum(myhist$count)
myhist0$count<-myhist0$count/sum(myhist0$count)
pdf("~/Dropbox/LDLD/figs/revisions/hetexcess_intermin5.pdf")
barplot(rbind(myhist$count,myhist0$count),beside=T,col=c(rgb(1,0,0,0.5),rgb(0.3,0.3,0.3,0.3)),space=c(0,0.1),ylim=c(0,1),xlab='p.value',ylab='density')
dev.off()
        
pdf("~/Dropbox/LDLD/figs/revisions/hetexcess_intermin5_vioplot.pdf")
vioplot(mydata_sign$P_HET_EXCESS,col=rgb(1,0,0,0.5),drawRect=F)
vioplot(mydata_nosign$P_HET_EXCESS,col=rgb(0.3,0.3,0.3,0.3),add=T,drawRect=F)
dev.off()

as.numeric(sapply(mydata_sign[,4], function(x) strsplit(as.character(x),"/")[[1]][2]))/as.numeric(sapply(mydata_sign[,5], function(x) strsplit(as.character(x),"/")[[1]][2]))
as.numeric(sapply(mydata_nosign[,4], function(x) strsplit(as.character(x),"/")[[1]][2]))/as.numeric(sapply(mydata_nosign[,5], function(x) strsplit(as.character(x),"/")[[1]][2]))

1-as.numeric(sapply(mydata_sign[,4], function(x) strsplit(as.character(x),"/")[[1]][2]))/as.numeric(sapply(mydata_sign[,5], function(x) strsplit(as.character(x),"/")[[1]][2]))
1-as.numeric(sapply(mydata_nosign[,4], function(x) strsplit(as.character(x),"/")[[1]][2]))/as.numeric(sapply(mydata_nosign[,5], function(x) strsplit(as.character(x),"/")[[1]][2]))

zcat /mnt/sequencedb/1000Genomes/ftp/phase3/20140910/ALL.chr${MYCHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | head -2000 | grep '#' > /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/merged/min5/sign_GoNL.vcf

}

}

#genotype qualities
{
tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS4.bed | gzip -f > /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GQ_coding.gz &
tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS6.bed | gzip -f > /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GQ_intergenic.gz &

#wget 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/vcf_with_sample_level_annotation/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.vcf.gz'

for MYCHR in `seq 14 21`;do
echo $MYCHR
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/genotype_likelihoods/shapeit2/ALL.chr${MYCHR}.phase3_bc_union.20130502.biallelic_svm_snps_indelsAF0.005_svs.gl.vcf.gz"
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/genotype_likelihoods/shapeit2/ALL.chr${MYCHR}.phase3_bc_union.20130502.biallelic_svm_snps_indelsAF0.005_svs.gl.vcf.gz.tbi"
done



#candidates GL
for MYCHR in `seq 10 22`;do
tabix ALL.chr${MYCHR}.phase3_bc_union.20130502.biallelic_svm_snps_indelsAF0.005_svs.gl.vcf.gz -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS4.bed  >> /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.codingS4.vcf
done
gzip -f /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.codingS4.vcf

for MYCHR in `seq 10 22`;do
tabix ALL.chr${MYCHR}.phase3_bc_union.20130502.biallelic_svm_snps_indelsAF0.005_svs.gl.vcf.gz -R /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS6.bed  >> /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.intergenicS6.vcf
done
gzip -f /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.intergenicS6.vcf

#background GL
for MYCHR in `seq 10 22`;do
bedtools subtract -a /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min5/chr${MYCHR}.vcf.gz -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS4.bed | awk -v OFS='\t' '{print $1,$2-1,$2}' > temp.bed
tabix ALL.chr${MYCHR}.phase3_bc_union.20130502.biallelic_svm_snps_indelsAF0.005_svs.gl.vcf.gz -R temp.bed | awk '{if (NF==2544){print}}' >> /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.codingS4_back.vcf
done
gzip -f /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.codingS4_back.vcf

for MYCHR in `seq 10 22`;do
bedtools subtract -a /mnt/scratch/fabrizio/LDLD/above95/intergenic.1000g/repl1/min5/chr${MYCHR}.vcf.gz -b /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/TableS6.bed | awk -v OFS='\t' '{print $1,$2-1,$2}' > temp.bed
tabix ALL.chr${MYCHR}.phase3_bc_union.20130502.biallelic_svm_snps_indelsAF0.005_svs.gl.vcf.gz -R temp.bed >> /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.intergenicS6_back.vcf
done
gzip -f /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.intergenicS6_back.vcf


library(data.table)
mydata<-fread("zcat /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.codingS4.vcf.gz",sep='\t')
tempH_f<-function(z) sapply(mydata[z,10:ncol(mydata)], function(x) unlist(strsplit(x,","))[2])
tempR_f<-function(z) sapply(mydata[z,10:ncol(mydata)], function(x) unlist(strsplit(x,","))[1])
tempdataH<-as.numeric(sapply(1:nrow(mydata), function(z) tempH_f(z)))
tempdataR<-as.numeric(sapply(1:nrow(mydata), function(z) tempR_f(z)))
tempdataH<-tempdataH[tempdataH>tempdataR & (10^tempdataH+10^tempdataR)>0.5]
myres<-10^tempdataH
library(data.table)
mydata0<-fread("zcat /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.codingS4_back.vcf.gz",sep='\t')
tempH_f<-function(z) sapply(mydata0[z,10:ncol(mydata0)], function(x) unlist(strsplit(x,","))[2])
tempR_f<-function(z) sapply(mydata0[z,10:ncol(mydata0)], function(x) unlist(strsplit(x,","))[1])
myi<-sort(sample(1:nrow(mydata0),min(5000,nrow(mydata0)),replace=F))
tempdataH<-as.numeric(sapply(myi, function(z) tempH_f(z)))
tempdataR<-as.numeric(sapply(myi, function(z) tempR_f(z)))
tempdataH<-tempdataH[tempdataH>tempdataR & (10^tempdataH+10^tempdataR)>0.5]
myres0<-10^tempdataH
#pdf("/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/hist_GL_coding.pdf")
#hist(myres0,col="gray")
#hist(myres,add=T)
#dev.off()

myres0_h<-hist(myres0)
myres0_h$counts<-myres0_h$counts/sum(myres0_h$counts)
myres_h<-hist(myres)
myres_h$counts<-myres_h$counts/sum(myres_h$counts)
source("~/Dropbox/general_utils/general_functions.R")
col1="coral1"
col2="gray79"
pdf("/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/hist_GL_coding.pdf")
paired_barplot(rbind(myres_h$counts,myres0_h$counts),name.groups=myres0_h$breaks[-1],mylas=2,mycol=c(col1,col2),name.categories=c("candidates","control variants"),ylab="density",xlab="genotype likelihood for heterozygous state")
dev.off()

library(data.table)
mydata<-fread("zcat /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.intergenicS6.vcf.gz",sep='\t')
tempH_f<-function(z) sapply(mydata[z,10:ncol(mydata)], function(x) unlist(strsplit(x,","))[2])
tempR_f<-function(z) sapply(mydata[z,10:ncol(mydata)], function(x) unlist(strsplit(x,","))[1])
tempdataH<-as.numeric(sapply(1:nrow(mydata), function(z) tempH_f(z)))
tempdataR<-as.numeric(sapply(1:nrow(mydata), function(z) tempR_f(z)))
tempdataH<-tempdataH[tempdataH>tempdataR & (10^tempdataH+10^tempdataR)>0.5]
myres_int<-10^tempdataH
mydata<-fread("zcat /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.intergenicS6_back.vcf.gz",sep='\t')
tempH_f<-function(z) sapply(mydata[z,10:ncol(mydata)], function(x) unlist(strsplit(x,","))[2])
tempR_f<-function(z) sapply(mydata[z,10:ncol(mydata)], function(x) unlist(strsplit(x,","))[1])
myi<-sort(sample(1:nrow(mydata),min(5000,nrow(mydata)),replace=F))
tempdataH<-as.numeric(sapply(myi, function(z) tempH_f(z)))
tempdataR<-as.numeric(sapply(myi, function(z) tempR_f(z)))
tempdataH<-tempdataH[tempdataH>tempdataR & (10^tempdataH+10^tempdataR)>0.5]
myres_int0<-10^tempdataH
#pdf("/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/hist_GL_coding.pdf")
#hist(myres0,col="gray")
#hist(myres,add=T)
#dev.off()

myres0_h<-hist(myres_int0)
myres0_h$counts<-myres0_h$counts/sum(myres0_h$counts)
myres_h<-hist(myres_int)
myres_h$counts<-myres_h$counts/sum(myres_h$counts)
source("~/Dropbox/general_utils/general_functions.R")
col1="coral1"
col2="gray79"
pdf("/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/hist_GL_inter.pdf")
paired_barplot(rbind(myres_h$counts,myres0_h$counts),name.groups=myres0_h$breaks[-1],mylas=2,mycol=c(col1,col2),name.categories=c("candidates","control variants"),ylab="density",xlab="genotype likelihood for heterozygous state")
dev.off()

myres0_h<-hist(myres0)
myres0_h$counts<-myres0_h$counts/sum(myres0_h$counts)
myres_h<-hist(myres)
myres_h$counts<-myres_h$counts/sum(myres_h$counts)
hist.cod<-myres_h$counts
hist.cod0<-myres0_h$counts
myres0_h<-hist(myres_int0)
myres0_h$counts<-myres0_h$counts/sum(myres0_h$counts)
myres_h<-hist(myres_int)
myres_h$counts<-myres_h$counts/sum(myres_h$counts)
hist.int<-myres_h$counts
hist.int0<-myres0_h$counts
pdf("/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/hist_GL.pdf")
par(mfrow=c(2,1))
paired_barplot(rbind(hist.cod,hist.cod0),name.groups=myres0_h$breaks[-1],mylas=2,mycol=c(col1,col2),name.categories=c("candidates","control variants"),ylab="density",xlab="genotype likelihood for heterozygous state",main="coding")
paired_barplot(rbind(hist.int,hist.int0),name.groups=myres0_h$breaks[-1],mylas=2,mycol=c(col1,col2),name.categories=c("candidates","control variants"),ylab="density",xlab="genotype likelihood for heterozygous state",main="intergenic")
dev.off()


(10^tempdataH)[(10^tempdataH)]
tempdataH<-sapply(1, function(z) tempH_f(z))

sum(tempdataH>tempdataR & (tempdataH^10+tempdataR^10)>0.5,na.rm=T)

sapply(1:3, function(z) temp_f(z))

sapply(mydata[1,10], function(x) {y<-as.numeric(strsplit(x,",")); if (y[3]==max(y)){return(y)}})

library(data.table)
mydata<-fread("zcat /r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/revisions/bed_ns/GL.codingS4.vcf.gz")

MYFOLDER=$( ls | grep sims )
for IFOLDER in $MYFOLDER; do
echo $IFOLDER
zcat $IFOLDER/repl1/anal/sorted.res.gz | head -1000000 | sed 's/F/A/' | sed 's/F/B/' | grep A | grep B |wc -l
ls -lh $IFOLDER/repl1/anal/sorted.res.gz
ls -lh $IFOLDER/repl1/reschr1/anal/sorted.res.gz
done




}
#============================
#=========== SIMULATIONS ========
#============================

#generate simulations
{
#first of all I should think what to simulate:
#1)-continuous variable introducing errors.
#-continuous variable distributed as a gaussian. introducing errors at a given number of sites with probability proportional to distribution
#-set a fraction of sites that are potentially susceptible (1%).
#-mutate them proportionally to the variable . 

#2)-batches
#-set a fraction of sites that are potentially susceptible (1%).
#-divide samples in batches (10% and 50%)
#-set a proportion of errors for this (100% and 50% and 10%)

#

#1% of sites susceptible to errors
#of these, I would take a part of the individuals and I would add mutations. Just hets since I find this huge excess.
#add invented snps rather than modifying errors, at low frequencies. Maybe I should just add them with the same frequency distribution of other snps.
#actually sims should be simple.
#

#prepare source files
cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/repl2/min1
for MYCHROM in `seq 1 22`;do
zcat chr${MYCHROM}.vcf.gz | grep -v '#' > /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/min1/chr${MYCHROM}.vcf 
done

for MYCHROM in `seq 1 22`;do
cat /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/coding.exons/myfilter.${MYCHROM}.bed | awk 'BEGIN{COUNTER=0}{COUNTER=COUNTER+$3-$2}END{print COUNTER}' >> /mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/coding.exons.regions.length.tab
done

#==============BATCHES==============
library(data.table)
#I sample from STU with replacement, because population with smallest bias
#batches
proportion_batch1=0.2 #0.5
prisksites=1/1000
nsamplesout<-200
lastbatch1<-round(nsamplesout*proportion_batch1)
mynspos<-read.table("~/Dropbox/LDLD/nspops.txt")
length_regions<-unname(unlist(read.table("/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/coding.exons.regions.length.tab")))

for (p_error_batch1 in c(0.1,0.5)){
for (myfreq in c(0.01,0.05)){
for (myrepl in 1:2){
destination_folder<-paste0("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims_pbatch",proportion_batch1,"_prisk",prisksites,"_nsample",nsamplesout,"_perror",p_error_batch1,"_min",myfreq,"/repl",myrepl,"/")
system(paste0("mkdir -p ",destination_folder))
system(paste0("mkdir -p ",destination_folder,"reschr1"))
for (mychr in 1:22){
nrisksites<-round(length_regions[mychr]*prisksites)
print(mychr)
mydata<-fread(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/min1/chr",mychr,".vcf"))
#to determine how many errors per site
#if 1/1000 sites are at risk, then about the same amount of erroneous snps and real ones


myneworder<-sample((sum(mynspos$V1[1:11])+10):sum(mynspos$V1[1:12]),nsamplesout,replace=TRUE)
mydata_temp<-mydata[,myneworder,with=FALSE]
myf<-function(x){
temp<-unlist(strsplit(unlist(x),split="|"))
sum(temp==1)/(sum(temp==0)+sum(temp==1))
}
myfreqs<-apply(mydata_temp, MARGIN=1,FUN=function(x) myf(x)) 
mydata_temp<-mydata_temp[myfreqs>=myfreq & myfreqs<=(1-myfreq),]
mydatapos_temp<-mydata$V2[myfreqs>=myfreq & myfreqs<=(1-myfreq)]

n_errs<-rbinom(nrisksites,lastbatch1,prob=p_error_batch1)
pos_errs<-sapply(n_errs, function(x) sort(sample(1:lastbatch1,x,replace=FALSE) ))

erdata<-mydata_temp[1,]
erdata[erdata != "0|0"] <- "0|0"
erdata0<-erdata
for (i in 1:(nrisksites-1)){ erdata<-rbind(erdata,erdata0) }
for ( ier in 1:nrisksites ) { erdata[ier,pos_errs[[ier]]] <- "1|0" }

erdata<-rbind(erdata,mydata_temp)
erdata$V1=mychr
erdata$V2=c(1:nrisksites,mydatapos_temp)
erdata$V3=c(rep("F",nrisksites),rep("T",length(mydatapos_temp)))
erdata$V4=c(rep("A",nrisksites),rep("A",length(mydatapos_temp)))
erdata$V5=c(rep("G",nrisksites),rep("G",length(mydatapos_temp)))
erdata$V6=c(rep("100",nrisksites),rep("100",length(mydatapos_temp)))
erdata$V7=c(rep("PASS",nrisksites),rep("PASS",length(mydatapos_temp)))
erdata$V8=c(rep("AA",nrisksites),rep("AA",length(mydatapos_temp)))
erdata$V9=c(rep("GT",nrisksites),rep("GT",length(mydatapos_temp)))
erdata<-erdata[,c((nsamplesout+1):(nsamplesout+9),1:nsamplesout),with=FALSE]
myfreqs<-apply(erdata, MARGIN=1,FUN=function(x) myf(x)) 
erdata<-erdata[myfreqs>=myfreq & myfreqs<=(1-myfreq),]


setnames(erdata,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sapply(1:nsamplesout,function(x) paste0("A",x))))
fwrite(erdata,file=paste0(destination_folder,"chr",mychr,".vcf"),quote=FALSE,sep="\t")
system(paste0("gzip -f ",destination_folder,"chr",mychr,".vcf"))

fwrite(erdata[,c(1:9,sample(10:nsamplesout,replace=FALSE)),with=FALSE],file=paste0(destination_folder,"reschr1/chr",mychr,".vcf"),quote=FALSE,sep="\t")
system(paste0("gzip -f ",destination_folder,"reschr1/chr",mychr,".vcf"))

}}}
}


#==============CONTINUOUS VARIABLE==============
library(data.table)
#I sample from STU with replacement, because population with smallest bias
#batches
prisksites=1/1000
nsamplesout<-200
mymean<-0.2
mysd<-0.2

mynspos<-read.table("~/Dropbox/LDLD/nspops.txt")
length_regions<-unname(unlist(read.table("/mnt/scratch/fabrizio/LDLD/above95/preparation_coding1000g/coding.exons.regions.length.tab")))

for (myfreq in c(0.01,0.05)){
for (myrepl in 1:2){
mydisterr<-rnorm(nsamplesout,mymean,mysd)
mydisterr[mydisterr<0]<-0
mydisterr[mydisterr>1]<-1
write.table(mydisterr,file=paste0(destination_folder,"distrerr.tab"),row.names=FALSE,quote=FALSE,col.names=FALSE)
destination_folder<-paste0("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims_mu",mymean,"_prisk",prisksites,"_nsample",nsamplesout,"_mysd",mysd,"_min",myfreq,"/repl",myrepl,"/")
system(paste0("mkdir -p ",destination_folder))
system(paste0("mkdir -p ",destination_folder,"reschr1"))
for (mychr in 1:22){
nrisksites<-round(length_regions[mychr]*prisksites)
print(mychr)
mydata<-fread(paste0("/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/min1/chr",mychr,".vcf"))
#to determine how many errors per site
#if 1/1000 sites are at risk, then about the same amount of erroneous snps and real ones


myneworder<-sample((sum(mynspos$V1[1:11])+10):sum(mynspos$V1[1:12]),nsamplesout,replace=TRUE)
mydata_temp<-mydata[,myneworder,with=FALSE]
myf<-function(x){
temp<-unlist(strsplit(unlist(x),split="|"))
sum(temp==1)/(sum(temp==0)+sum(temp==1))
}
myfreqs<-apply(mydata_temp, MARGIN=1,FUN=function(x) myf(x)) 
mydata_temp<-mydata_temp[myfreqs>=myfreq & myfreqs<=(1-myfreq),]
mydatapos_temp<-mydata$V2[myfreqs>=myfreq & myfreqs<=(1-myfreq)]

n_errs<-rbinom(nrisksites,nsamplesout,prob=mydisterr)
pos_errs<-sapply(n_errs, function(x) sort(sample(1:nsamplesout,x,replace=FALSE) ))

erdata<-mydata_temp[1,]
erdata[erdata != "0|0"] <- "0|0"
erdata0<-erdata
for (i in 1:(nrisksites-1)){ erdata<-rbind(erdata,erdata0) }
for ( ier in 1:nrisksites ) { erdata[ier,pos_errs[[ier]]] <- "1|0" }

erdata<-rbind(erdata,mydata_temp)
erdata$V1=mychr
erdata$V2=c(1:nrisksites,mydatapos_temp)
erdata$V3=c(rep("F",nrisksites),rep("T",length(mydatapos_temp)))
erdata$V4=c(rep("A",nrisksites),rep("A",length(mydatapos_temp)))
erdata$V5=c(rep("G",nrisksites),rep("G",length(mydatapos_temp)))
erdata$V6=c(rep("100",nrisksites),rep("100",length(mydatapos_temp)))
erdata$V7=c(rep("PASS",nrisksites),rep("PASS",length(mydatapos_temp)))
erdata$V8=c(rep("AA",nrisksites),rep("AA",length(mydatapos_temp)))
erdata$V9=c(rep("GT",nrisksites),rep("GT",length(mydatapos_temp)))
erdata<-erdata[,c((nsamplesout+1):(nsamplesout+9),1:nsamplesout),with=FALSE]
myfreqs<-apply(erdata, MARGIN=1,FUN=function(x) myf(x)) 
erdata<-erdata[myfreqs>=myfreq & myfreqs<=(1-myfreq),]


setnames(erdata,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sapply(1:nsamplesout,function(x) paste0("A",x))))
fwrite(erdata,file=paste0(destination_folder,"chr",mychr,".vcf"),quote=FALSE,sep="\t")
system(paste0("gzip -f ",destination_folder,"chr",mychr,".vcf"))

fwrite(erdata[,c(1:9,sample(10:nsamplesout,replace=FALSE)),with=FALSE],file=paste0(destination_folder,"reschr1/chr",mychr,".vcf"),quote=FALSE,sep="\t")
system(paste0("gzip -f ",destination_folder,"reschr1/chr",mychr,".vcf"))

}}
}



}

#verify state of processing
{
cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/
myfolders=$( ls | grep sims )
for ifolder in $myfolders;do
for myrepl in 1 2; do
nsums=$( find ${ifolder}/repl${myrepl} -name '*sum' |wc -l )
#nsums=$( find ${ifolder}/repl${myrepl}/reschr1 -name '*sum' |wc -l )
if [ $nsums -lt 462 ];then echo "only " $nsums "sum files in " ${ifolder}/repl${myrepl};fi
done;done
#ok, all have all sum files

cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/
myfolders=$( ls | grep sims )
for ifolder in $myfolders;do
for myrepl in 1 2; do
#nsums=$( find ${ifolder}/repl${myrepl}/anal -name 'sorted.res.gz' |wc -l )
nsums=$( find ${ifolder}/repl${myrepl}/reschr1/anal -name 'sorted.res.gz' |wc -l )
if [ $nsums -lt 1 ];then echo $nsums "files in " ${ifolder}/repl${myrepl};fi
done;done

cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/
myfolders=$( ls | grep sims )
for ifolder in $myfolders;do
for myrepl in 1 2; do
ls -lh ${ifolder}/repl${myrepl}/anal/sorted.res.gz
done;done


cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/
myfolders=$( ls | grep sims )
for ifolder in $myfolders;do
for myrepl in 1 2; do
nsums=$( find ${ifolder}/repl${myrepl} -name 'combined_pvals.RData' |wc -l )
#nsums=$( find ${ifolder}/repl${myrepl}/reschr1 -name '*sum' |wc -l )
if [ $nsums -lt 1 ];then echo $nsums "files in " ${ifolder}/repl${myrepl};fi
done;done

cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/
for myfolder in sims_pbatch0.2_prisk0.001_nsample200_perror0.5_min0.01/repl1 sims_pbatch0.2_prisk0.001_nsample200_perror0.5_min0.01/repl2 sims_pbatch0.2_prisk0.001_nsample200_perror0.5_min0.05/repl1 sims_pbatch0.2_prisk0.001_nsample200_perror0.5_min0.05/repl2;do
#cp $myfolder/anal/sorted.res.gz $myfolder/anal/sorted_backup.res.gz
#cp $myfolder/reschr1/anal/sorted.res.gz $myfolder/reschr1/anal/sorted_backup.res.gz
zcat $myfolder/anal/sorted_backup.res.gz | awk '{if ($18<0.000000001){print}}' | gzip -f > $myfolder/anal/sorted.res.gz
zcat $myfolder/reschr1/anal/sorted_backup.res.gz | awk '{if ($18<0.000000001){print}}' | gzip -f > $myfolder/reschr1/anal/sorted.res.gz
done

cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/
myfolders=$( ls | grep sims )
for ifolder in $myfolders;do
for myrepl in 1 2; do
nsums=$( find ${ifolder}/repl${myrepl} -name 'logpop_nAB_pos20.RData' |wc -l )
#nsums=$( find ${ifolder}/repl${myrepl}/reschr1 -name '*sum' |wc -l )
if [ $nsums -lt 1 ];then echo $nsums "files in " ${ifolder}/repl${myrepl};fi
done;done

cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/
myfolders=$( ls | grep sims )
for ifolder in $myfolders;do
for myrepl in 1 2; do
nsums=$( find ${ifolder}/repl${myrepl} -name 'logpop_nAB_pos.RData' |wc -l )
#nsums=$( find ${ifolder}/repl${myrepl}/reschr1 -name '*sum' |wc -l )
if [ $nsums -lt 1 ];then echo $nsums "files in " ${ifolder}/repl${myrepl};fi
done;done


}

myfolder=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims_pbatch0.2_prisk0.001_nsample100_perror0.5_min0.01/repl1

myfolders=$( ls /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims*/* | grep ':' | sed 's/://g' )
for myfolder in $myfolders;do
myfolderdropbox=$( echo $myfolder | sed 's/\/mnt\/scratch\/fabrizio\/LDLD\/above95\/coding.exons.1000g\/sims\/'/'~\/Dropbox\/LDLD\/ipynb\/figs\/sims\//g' )
mytag=$( echo $myfolder | sed 's/\/mnt\/scratch\/fabrizio\/LDLD\/above95\/coding.exons.1000g\/sims\///g' | sed 's/\/repl/r/g' )
cmd="Rscript ~/Dropbox/LDLD/scripts/analyses_LDLD_script.R $myfolder $myfolderdropbox 1"
qsub -cwd -b y -l h_vmem=12G,virtual_free=12G,mem_free=12G -N $mytag ~/Dropbox/Vindija/scripts/qsub_runner.sh $cmd
done

myfolders=$( ls /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims*/* | grep ':' | sed 's/://g' )
for myfolder in $myfolders;do
myfolderdropbox=$( echo $myfolder | sed 's/\/mnt\/scratch\/fabrizio\/LDLD\/above95\/coding.exons.1000g\/sims\/'/'~\/Dropbox\/LDLD\/ipynb\/figs\/sims\//g' )
mytag=$( echo $myfolder | sed 's/\/mnt\/scratch\/fabrizio\/LDLD\/above95\/coding.exons.1000g\/sims\///g' | sed 's/\/repl/r/g' )
cmd="Rscript ~/Dropbox/LDLD/scripts/analyses_LDLD_script.R $myfolder $myfolderdropbox 1"
qsub -cwd -b y -l h_vmem=12G,virtual_free=12G,mem_free=12G -N $mytag ~/Dropbox/Vindija/scripts/qsub_runner.sh $cmd
done


fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims$ myfolder=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims_mu0.2_prisk0.001_nsample100_mysd0.2_min0.01/repl2
fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims$ myfolderdropbox=$( echo $myfolder | sed 's/\/mnt\/scratch\/fabrizio\/LDLD\/above95\/coding.exons.1000g\/sims\/'/'~\/Dropbox\/LDLD\/ipynb\/figs\/sims\//g' )
fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims$ mytag=$( echo $myfolder | sed 's/\/mnt\/scratch\/fabrizio\/LDLD\/above95\/coding.exons.1000g\/sims\///g' | sed 's/\/repl/r/g' )
fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims$ cmd="Rscript ~/Dropbox/LDLD/scripts/analyses_LDLD_script.R $myfolder $myfolderdropbox 1"
fabrizio_mafessoni@bionc02:/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims$ qsub -cwd -b y -l h_vmem=8G,virtual_free=8G,mem_free=8G -N $mytag ~/Dropbox/Vindija/scripts/qsub_runner.sh $cmd
Your job 3886100 ("sims_mu0.2_prisk0.001_nsample100_mysd0.2_min0.01r2") has been submitted
cat /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims_mu0.2_prisk0.001_nsample100_mysd0.2_min0.01r2.e3886100 #skip PREPROCESSING ok but then unexpected end of input..

#retry with new script file -> 3886101

ssh bionc05
myfolder=/mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims_mu0.2_prisk0.001_nsample100_mysd0.2_min0.05/repl1
myfolderdropbox=$( echo $myfolder | sed 's/\/mnt\/scratch\/fabrizio\/LDLD\/above95\/coding.exons.1000g\/sims\/'/'~\/Dropbox\/LDLD\/ipynb\/figs\/sims\//g' )
mytag=$( echo $myfolder | sed 's/\/mnt\/scratch\/fabrizio\/LDLD\/above95\/coding.exons.1000g\/sims\///g' | sed 's/\/repl/r/g' )
cmd="Rscript ~/Dropbox/LDLD/scripts/analyses_LDLD_script.R $myfolder $myfolderdropbox 1"
cd ~/
nohup $cmd

myfolder=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5
myfolderdropbox=~/Dropbox/LDLD/ipynb/sims/GoNL12_as1pop/min5
nohup Rscript ~/mylogs/mylogs_GoNL/analyses_LDLD_script.R $myfolder $myfolderdropbox 1 &

cd /mnt/restricted/GoNL/GoNL1/release5.4/03_IL_SNVs/gonl-abc_samples
for MYCHR in `seq 17 22`;do
echo $MYCHR
#bcftools merge much faster but errors
vcf-merge gonl-abc_samples.chr${MYCHR}.release5.raw_SNVs.vcf.gz ../../../../GoNL2/vcf/merged.chr${MYCHR}.vcf.gz | ~/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/nspops_GoNL12_as1pop.txt -f 0.05 -O 0 | gzip -f > /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf/chr${MYCHR}.vcf.gz &
done



for MYCHR in `seq 1 22`;do
zcat chr${MYCHR}.vcf.gz | bgzip > chr${MYCHR}.temp.vcf.gz &
done
#gzip: chr1.vcf.gz: unexpected end of file
#gzip: chr16.vcf.gz: unexpected end of file
#gzip: chr8.vcf.gz: unexpected end of file
#so I copy the others but do not delete the originals of these
for MYCHR in `seq 2 7`;do
tabix chr${MYCHR}.temp.vcf.gz
done
for MYCHR in `seq 9 15`;do
tabix chr${MYCHR}.temp.vcf.gz &
done
for MYCHR in `seq 17 22`;do
tabix chr${MYCHR}.temp.vcf.gz &
done

cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf
zcat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/chr22.vcf.gz | head -500 | grep CHROM | awk '{for (i=10;i<=(NF-1);i++){printf "%s\t",$i};printf "%s\n",$NF}' > temp
zcat chr22.vcf.gz | head -500 | grep CHROM | awk '{for (i=10;i<=(NF-1);i++){printf "%s\t",$i};printf "%s\n",$NF}' > tempall

#to generate awk command in 
tempall<-read.table("tempall",stringsAsFactors=F)
temp<-read.table("temp",stringsAsFactors=F)
paste(paste0("$",9+match(temp,tempall)),collapse=",")

cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf
for MYCHR in `seq 1 22`;do echo ${MYCHR};zcat chr${MYCHR}.vcf.gz > chr${MYCHR}.temp.vcf; bgzip chr${MYCHR}.temp.vcf;  done
for MYCHR in `seq 1 22`;do tabix chr${MYCHR}.temp.vcf.gz;  done

load("/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/logpop_nAB_pos20b.RData")
logpop_nAB_pos20[logpop_nAB_pos20>400]<-100
save(logpop_nAB_pos20,file="/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/logpop_nAB_pos20c.RData")

MYMIN=5 #MYMIN=5
FREQFILTER=0.0${MYMIN}
PVALUE_THR=0.05
FILERDATA=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/logpop_nAB_pos20b.RData
MYTAG=logpop_nAB_pos20
DESTINATION_FOLDER=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/bad_snps_min${MYMIN}_glmm
PARSER_FILE=~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_GoNL.sh
NSPOPS_FILE=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/nspops_GoNL12_as1pop.txt
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 22`;do
#TARGETVCF=/mnt/restricted/GoNL/GoNL2/vcf/merged.chr${MYCHR}.vcf.gz
TARGETVCF=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf/chr${MYCHR}.temp.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $TARGETVCF $MYCHR $FILERDATA $MYTAG $FREQFILTER $PVALUE_THR $DESTINATION_FOLDER $PARSER_FILE $NSPOPS_FILE &
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG

MYMIN=5 #MYMIN=5
FREQFILTER=0.0${MYMIN}
PVALUE_THR=0.05
FILERDATA=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/logpop_nAB_pos20c.RData
MYTAG=logpop_nAB_pos20
DESTINATION_FOLDER=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/bad_snps_min${MYMIN}_glmm_c
PARSER_FILE=~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_GoNL.sh
NSPOPS_FILE=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/nspops_GoNL12_as1pop.txt
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 11 22`;do
#TARGETVCF=/mnt/restricted/GoNL/GoNL2/vcf/merged.chr${MYCHR}.vcf.gz
TARGETVCF=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf/chr${MYCHR}.temp.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $TARGETVCF $MYCHR $FILERDATA $MYTAG $FREQFILTER $PVALUE_THR $DESTINATION_FOLDER $PARSER_FILE $NSPOPS_FILE &
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG

Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG

MYMIN=5 #MYMIN=5
FREQFILTER=0.0${MYMIN}
PVALUE_THR=0.05
FILERDATA=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/logpop_nAB_pos20c.RData
MYTAG=logpop_nAB_pos20
DESTINATION_FOLDER=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/bad_snps_min${MYMIN}_glmm_c_p03
PARSER_FILE=~/Dropbox/LDLD/scripts/nAvsnAB_scan_child_GoNL.sh
NSPOPS_FILE=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/nspops_GoNL12_as1pop.txt
mkdir -p $DESTINATION_FOLDER
for MYCHR in `seq 1 10`;do
#TARGETVCF=/mnt/restricted/GoNL/GoNL2/vcf/merged.chr${MYCHR}.vcf.gz
TARGETVCF=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf/chr${MYCHR}.temp.vcf.gz
~/Dropbox/LDLD/scripts/nAvsnAB_glmm.sh $TARGETVCF $MYCHR $FILERDATA $MYTAG $FREQFILTER $PVALUE_THR $DESTINATION_FOLDER $PARSER_FILE $NSPOPS_FILE &
done
~/Dropbox/LDLD/scripts/nAvsnAB_merge.sh $DESTINATION_FOLDER $MYTAG
Rscript ~/Dropbox/LDLD/scripts/nAvsnAB_p2fdr.R $DESTINATION_FOLDER $MYTAG
#modify to give fdr 0.5

#758 is the first position of GoNL2
#/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf$ zcat chr1.temp.vcf.gz | awk '{print $757,$758}' | less
cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf
for MYCHR in `seq 1 22`;do
zcat chr${MYCHR}.temp.vcf.gz | awk 'BEGIN{tot=0}{count=0;for (i=758;i<=NF;i++){if (length($i)>1){count=count+1}}; if (count>10){tot=tot+1}}END{print tot}' >> count5.tab &
done


zcat chr${MYCHR}.temp.vcf.gz | awk 'BEGIN{tot=0}{count=0;for (i=758;i<=NF;i++){if (length($i)>1){count=count+1}}; if (count>10){tot=tot+1;print}}END{print tot}'

cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf
for MYCHR in `seq 1 22`;do
zcat chr${MYCHR}.temp.vcf.gz | awk 'BEGIN{tot=0}{count=0;for (i=758;i<=NF;i++){if (length($i)>1){count=count+1}}; if (count>5){tot=tot+1}}END{print tot}' >> count5b.tab &
done

#length 59
#cat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/bad_snps_min5_glmm_c/bad_snps_logpop_nAB_pos20_glmm_all.bed > ~/Dropbox/LDLD/ms/GBE/TableS8_49.bed
#length 113
cat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/bad_snps_min5_glmm_c/bad_snps_logpop_nAB_pos20_glmm_all.bed > ~/Dropbox/LDLD/ms/GBE/TableS8.bed

zcat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf/chr${MYCHR}.temp.vcf.gz | head -1000 | grep '#' >> /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/sign.vcf
for MYCHR in `seq 1 22`;do
echo $MYCHR
cat ~/Dropbox/LDLD/ms/GBE/TableS8.bed | awk -v mychr=${MYCHR} '{if ($1==mychr){print}}' | sort -Vu -k2,2n > temp.bed
MYN=$( cat temp.bed | wc -l )
if [ $MYN -gt 0 ];then
    tabix /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf/chr${MYCHR}.temp.vcf.gz -R temp.bed >> /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/sign.vcf
fi
done

NSPOPS_FILE=/mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/nspops_GoNL12_as1pop.txt
zcat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf/chr${MYCHR}.temp.vcf.gz | head -1000 | grep '#' >> /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf
for MYCHR in `seq 1 22`;do
echo $MYCHR
MYN=$( cat temp.bed | wc -l )
    tabix /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/all_vcf/chr${MYCHR}.temp.vcf.gz $MYCHR | head -10000 | /r1/people/fabrizio_mafessoni/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ${NSPOPS_FILE} -f 0.05 -O 0 | shuf -n 50 | sort -Vu -k1,1 -k2,2n >> /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf
done

vcftools --vcf /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf --hardy --out /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf.hwe
cat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf | grep '#' > temp2.vcf 
cat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf | grep -v '#' | awk '{
for (i=1;i<=NF;i++){if ($i=="."){$i="./."}}; 
for (i=1;i<=NF;i++){printf "%s\t",$i;};printf "%s\n",$NF;
}' >> temp2.vcf 
vcftools --vcf temp2.vcf --hardy --out /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf.hwe
less /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/nosign.vcf.hwe.hwe


cd /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/intergenic/repl1/min5
MYFILES=$( ls *res )
for  myfile in  $MYFILES; do cat $myfile | head -5 |wc -l >> dump;done

#for  myfile in chr10.tabchr11.tab.res chr10.tabchr12.tab.res chr10.tabchr13.tab.res chr10.tabchr14.tab.res chr10.tabchr15.tab.res chr10.tabchr16.tab.res chr10.tabchr17.tab.res 
#chr10.tabchr18.tab.res chr10.tabchr19.tab.res chr10.tabchr20.tab.res chr10.tabchr21.tab.res chr10.tabchr22.tab.res chr11.tabchr12.tab.res chr11.tabchr13.tab.res chr11.tabchr14.tab.res 
#chr11.tabchr15.tab.res chr11.tabchr17.tab.res chr11.tabchr19.tab.res chr11.tabchr20.tab.res chr11.tabchr22.tab.res chr12.tabchr13.tab.res chr12.tabchr14.tab.res chr12.tabchr15.tab.res 
#chr12.tabchr16.tab.res chr12.tabchr17.tab.res chr12.tabchr18.tab.res chr12.tabchr19.tab.res chr12.tabchr20.tab.res chr12.tabchr21.tab.res chr12.tabchr22.tab.res chr13.tabchr14.tab.res 
#chr13.tabchr15.tab.res chr13.tabchr16.tab.res chr13.tabchr17.tab.res chr13.tabchr18.tab.res chr13.tabchr19.tab.res chr13.tabchr20.tab.res chr13.tabchr21.tab.res chr13.tabchr22.tab.res 
#chr14.tabchr15.tab.res chr14.tabchr16.tab.res chr14.tabchr17.tab.res chr14.tabchr18.tab.res chr14.tabchr19.tab.res chr14.tabchr20.tab.res chr14.tabchr21.tab.res chr14.tabchr22.tab.res 
#chr15.tabchr16.tab.res chr15.tabchr17.tab.res chr15.tabchr18.tab.res chr15.tabchr19.tab.res chr15.tabchr20.tab.res chr15.tabchr22.tab.res chr16.tabchr18.tab.res chr16.tabchr19.tab.res 
#chr16.tabchr21.tab.res chr16.tabchr22.tab.res chr17.tabchr18.tab.res chr17.tabchr19.tab.res chr17.tabchr20.tab.res chr17.tabchr22.tab.res chr18.tabchr19.tab.res chr18.tabchr20.tab.res 
#chr18.tabchr21.tab.res chr18.tabchr22.tab.res chr19.tabchr20.tab.res chr19.tabchr21.tab.res chr19.tabchr22.tab.res chr1.tabchr10.tab.res chr1.tabchr11.tab.res chr1.tabchr12.tab.res 
#chr1.tabchr13.tab.res chr1.tabchr14.tab.res chr1.tabchr15.tab.res chr1.tabchr16.tab.res chr1.tabchr17.tab.res chr1.tabchr18.tab.res chr1.tabchr19.tab.res chr1.tabchr20.tab.res chr1.tabchr21.tab.res 
#chr1.tabchr22.tab.res chr1.tabchr4.tab.res chr1.tabchr5.tab.res chr1.tabchr6.tab.res  chr1.tabchr8.tab.res chr1.tabchr9.tab.res chr20.tabchr21.tab.res chr20.tabchr22.tab.res chr21.tabchr22.tab.res chr2.tabchr10.tab.res  
#chr2.tabchr12.tab.res chr2.tabchr13.tab.res chr2.tabchr14.tab.res chr2.tabchr15.tab.res chr2.tabchr16.tab.res chr2.tabchr17.tab.res chr2.tabchr18.tab.res chr2.tabchr19.tab.res chr2.tabchr20.tab.res chr2.tabchr21.tab.res 
#chr2.tabchr3.tab.res chr2.tabchr4.tab.res chr2.tabchr5.tab.res chr2.tabchr6.tab.res chr2.tabchr7.tab.res chr2.tabchr8.tab.res chr2.tabchr9.tab.res chr3.tabchr10.tab.res chr3.tabchr11.tab.res chr3.tabchr12.tab.res 
#chr3.tabchr13.tab.res chr3.tabchr15.tab.res chr3.tabchr16.tab.res chr3.tabchr17.tab.res chr3.tabchr18.tab.res chr3.tabchr19.tab.res chr3.tabchr20.tab.res chr3.tabchr21.tab.res chr3.tabchr22.tab.res chr3.tabchr4.tab.res 
#chr3.tabchr5.tab.res  chr3.tabchr7.tab.res chr3.tabchr8.tab.res chr3.tabchr9.tab.res chr4.tabchr10.tab.res chr4.tabchr11.tab.res chr4.tabchr12.tab.res chr4.tabchr13.tab.res chr4.tabchr14.tab.res chr4.tabchr15.tab.res 
#chr4.tabchr16.tab.res chr4.tabchr17.tab.res chr4.tabchr18.tab.res chr4.tabchr19.tab.res chr4.tabchr20.tab.res chr4.tabchr21.tab.res chr4.tabchr22.tab.res chr4.tabchr5.tab.res chr4.tabchr6.tab.res chr4.tabchr7.tab.res 
#chr4.tabchr8.tab.res chr4.tabchr9.tab.res chr5.tabchr10.tab.res chr5.tabchr11.tab.res chr5.tabchr12.tab.res chr5.tabchr13.tab.res chr5.tabchr14.tab.res chr5.tabchr16.tab.res chr5.tabchr17.tab.res chr5.tabchr18.tab.res 
#chr5.tabchr19.tab.res chr5.tabchr20.tab.res chr5.tabchr21.tab.res chr5.tabchr22.tab.res chr5.tabchr6.tab.res chr5.tabchr7.tab.res chr5.tabchr9.tab.res chr6.tabchr10.tab.res chr6.tabchr11.tab.res chr6.tabchr12.tab.res 
#chr6.tabchr13.tab.res chr6.tabchr14.tab.res chr6.tabchr15.tab.res chr6.tabchr16.tab.res chr6.tabchr17.tab.res chr6.tabchr18.tab.res chr6.tabchr19.tab.res chr6.tabchr20.tab.res chr6.tabchr21.tab.res chr6.tabchr22.tab.res 
#chr6.tabchr7.tab.res chr6.tabchr8.tab.res chr6.tabchr9.tab.res chr7.tabchr10.tab.res chr7.tabchr11.tab.res chr7.tabchr12.tab.res chr7.tabchr13.tab.res chr7.tabchr14.tab.res chr7.tabchr15.tab.res chr7.tabchr16.tab.res 
#chr7.tabchr17.tab.res chr7.tabchr18.tab.res chr7.tabchr19.tab.res chr7.tabchr20.tab.res chr7.tabchr21.tab.res chr7.tabchr22.tab.res chr7.tabchr8.tab.res chr7.tabchr9.tab.res chr8.tabchr10.tab.res chr8.tabchr11.tab.res 
#chr8.tabchr12.tab.res chr8.tabchr13.tab.res chr8.tabchr14.tab.res chr8.tabchr15.tab.res chr8.tabchr16.tab.res chr8.tabchr17.tab.res chr8.tabchr18.tab.res chr8.tabchr20.tab.res chr8.tabchr21.tab.res chr8.tabchr22.tab.res 
#chr8.tabchr9.tab.res chr9.tabchr10.tab.res chr9.tabchr11.tab.res chr9.tabchr12.tab.res chr9.tabchr13.tab.res chr9.tabchr14.tab.res chr9.tabchr15.tab.res chr9.tabchr16.tab.res chr9.tabchr17.tab.res chr9.tabchr18.tab.res 
#chr9.tabchr19.tab.res chr9.tabchr20.tab.res chr9.tabchr21.tab.res chr9.tabchr22.tab.res;do
cat $myfile | awk '{if (NF==18){print}}' | awk '{
if ($12>6e-05) {pvals=0} else {pvals=1};
print 1,2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,pvals,pvals,pvals,pvals,pvals,$15,$16,$17,$18
}' >> anal/sorted.res
cat reschr1/$myfile | awk '{if (NF==18){print}}' | awk '{
if ($12>6e-05) {pvals=0} else {pvals=1};
print 1,2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,pvals,pvals,pvals,pvals,pvals,$15,$16,$17,$18
}' >> reschr1/anal/sorted.res
done
gzip -f anal/sorted.res
gzip -f reschr1/anal/sorted.res

zcat anal/sorted.res.gz | awk '{if ($15!=1){print}}' | gzip -f > anal/sorted.sign.res.gz
zcat reschr1/anal/sorted.res.gz | awk '{if ($15!=1){print}}' | gzip -f > reschr1/anal/sorted.sign.res.gz


tabix $MYTARGETFILE 1:1-1000000 | head -100 | grep -v '#' | /r1/people/fabrizio_mafessoni/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/nspops_GoNL12_as1pop.txt -f $FREQFILTER -O 0 -f $FREQFILTER -O 0


tabix $MYTARGETFILE 1:1-1000000 | head -100 | grep -v '#' | awk '{for (i=1;i<=7;i++){printf "%s\t",$i}; printf "%s\t%s\t","A","GT";for (i=10;i<=(NF-1);i++){printf "%s\t",substr($i,1,3)}; printf "%s\n",substr($NF,1,3) }' | /r1/people/fabrizio_mafessoni/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/nspops_GoNL12_as1pop.txt -f $FREQFILTER -O 0 -f $FREQFILTER -O 0

zcat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl11/min5/chr3.tab | /r1/people/fabrizio_mafessoni/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N /mnt/scratch/fabrizio/LDLD/GoNL/GoNL12_as1pop/min5/nspops_GoNL12_as1pop.txt -f $FREQFILTER -O 0 -f $FREQFILTER -O 0
zcat /mnt/scratch/fabrizio/LDLD/GoNL/GoNL1/intergenic/repl11/min5/chr3.tab | /r1/people/fabrizio_mafessoni/Dropbox/LDLD/scripts/filterbyfreq_multi.out -N ${NSPOPS_FILE} -f $FREQFILTER -O 0

args = commandArgs(trailingOnly=TRUE)
MYBEGIN_i = as.numeric(as.character(args[1]));
MYEND_i = as.numeric(as.character(args[2]));

if (TRUE) {
myfolders<-system("ls /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims/sims*/* | grep ':' | sed 's/://g'",intern=TRUE)
}

#myfolder<-myfolders[[1]]





for ( ifff in MYBEGIN_i:MYEND_i){ #length(myfolders) ){
myfolder<-myfolders[[ifff]]
print(myfolder)
myfolderdropbox<-paste0("~/Dropbox/LDLD/ipynb/figs/sims/",myfolder)


args = commandArgs(trailingOnly=TRUE)
myfolder = as.character(args[1]);
myfolderdropbox = as.character(args[2]);
npopulations = as.numeric(as.character(args[3]));

if ( as.logical(length(grep("min0.01",myfolder))) ) {mynamefreq<-1;} else if ( grep("min0.05",myfolder) ) {mynamefreq<-5;}
freqfiltering<-mynamefreq/100


options(scipen=999)
setwd("~/Dropbox/LDLD")
source("analysesLDLD_header.R")
source("~/Dropbox/general_utils/general_functions.R")
setwd("/mnt/scratch/fabrizio/LDLD")
library(data.table)
library(lme4)

popnames<-c("POP1")
computed_exact_pvalues<-TRUE
PREPROCESSING<-FALSE
using_sign_for_nAB<-TRUE
index_focal_pop<-0 #now I used 7 because I was looking at FINs but I guess in general it would be safer to have 0
nperm=1

system(paste0("mkdir -p ",myfolderdropbox))

PREPROCESSING<-FALSE; if (sum(grep("sorted.res.gz",dir(paste0(myfolder,"/anal"))))<1) { PREPROCESSING<-TRUE } 
#TAG 'LOAD_INITIAL_FILES'
if ( PREPROCESSING ) {
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
                    system(paste0("cat ",myfolder,"/reschr",myperm,"/*minilog | grep ncomparisons | grep -v pair | awk -v npops=", npopulations ,"  '{for (i=2;i<=(npops+2);i++){printf \"%d\\t\", $i }; printf \"\\n\"}' > ",myfolder,"/reschr",myperm,"/ncomparisons.txt"))
            print("aggregating res files for permutations")
            system(paste0("mkdir -p ",myfolder,"/reschr",myperm,"/logs"))
            system(paste0("mkdir -p ",myfolder,"/reschr",myperm,"/anal"))
            system(paste0("~/Dropbox/LDLD/scripts/postproc_split.sh ",myfolder,"/reschr",myperm))
            print("aggregating log files for permutations")
            organize_log_files(paste0(myfolder,"/reschr",myperm),npopulations)
            }
    #  print("generating mutation files")
        setwd("~/Dropbox/LDLD/scripts/")
    #    system("gcc ~/Dropbox/LDLD/scripts/LDLD_filter5.c -lgmp -lm -fno-stack-protector -o ~/Dropbox/LDLD/scripts/filterbyfreq.out")
        system(paste0("rm ",myfolder,"/logs/all.freqlog"))
    #   for (i in 1:22)
    #      {
    #        system(paste0("nohup ~/Dropbox/LDLD/scripts/filterbyfreq.out ",myfolder,"/chr",i,".tab ",myfolder,"/chr",i,".freqlog ",myfolder,"/chr",i,".mutlog ", freqfiltering," 1220"))
    #       system(paste0("nohup ~/Dropbox/LDLD/scripts/filterbyfreq_multi_sharing.out ",myfolder,"/chr",i,".tab nspops.txt ",freqfiltering," ",npopulations," 1 > ",myfolder,"/chr",i,".myfreq.sharing"))
    #      system(paste0("cat ",myfolder,"/chr",i,".freqlog | awk '{printf \"%d\\t\",",i,"; print $0}' >> ",myfolder,"/logs/all.freqlog"))
            #for (j in 1:nperm)
            #   {
            #   system(paste0("nohup ~/Dropbox/LDLD/scripts/filterbyfreq.out ",myfolder,"/reschr",j,"/chr",i,".tab ",myfolder,"/reschr",j,"/chr",i,".freqlog ",myfolder,"/reschr",j,"/chr",i,".mutlog ", freqfiltering," 1220"))
            #   }
    #     }
        system(paste0("gzip -f ",myfolder,"/logs/all.freqlog"))
        setwd("/mnt/scratch/fabrizio/LDLD")
    # system(paste0("cat ",myfolder,"/chr*.myfreq.sharing | grep -v [an] | awk '{for (i=1;i<=12;i++){if (i<=NF){printf \"%d\\t\",$i} else {printf  \"-1\\t\"};};printf \"\\n\"}' > ",myfolder,"/myfreq.sharing"))
        system(paste0("mkdir ",myfolderdropbox))
} else { print ("skip PREPROCESSING") }
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
if (sum(grep("combined_pvals",dir(myfolder)))<1){
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
    #
    #for GoNL 
    #dataLDLDA<-mydata_temp
    #myorder<-c(1,2,3,4)
    #tagsA<-apply(dataLDLDA[,myorder,with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse="."))
    #tagsB<-combined_pvals[["combined_pval_pos"]]$Group.1[combined_pvals[["combined_fdr_pos"]]<0.2]
    #mylogpop_pos20[[ipop+1]]<-dataLDLDA[!is.na(match(tagsA,tagsB)),]
    {
    mylogpop<-list();mylogpop_pos<-list();mylogpop_neg<-list();mylogpop_pos20<-list();mylogpop_pos01<-list();logpop_nAB<-list();logpop_nAB_pos<-list();logpop_nAB_neg<-list();logpop_nAB_pos20<-list();logpop_nAB_pos01<-list();
    for ( ipop in 0:(npopulations-1))    
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

    #dataLDLD<-mylogpop_pos20[[ipop+1]]
    #dataLDLDA<-apply(dataLDLD[,1:4,with=FALSE],MARGIN=1,function(x) paste0(x,collapse='.'))
    #dataLDLDB<-apply(mydata_temp[,1:4,with=FALSE],MARGIN=1,function(x) paste0(x,collapse='.'))
    #length(dataLDLDA)
    #length(unique(dataLDLDA))
    #length(dataLDLDB)
    #length(unique(dataLDLDB))


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

    }
} else { print("skip combined_pvals") }

scp -r /home/fabrizio/Dropbox/LDLD/scripts/analyses_LDLD_script.R fabrizio_mafessoni@bionc01.eva.mpg.de:/home/fabrizio_mafessoni/Dropbox/LDLD/scripts



#plot simulations
{
cd /mnt/scratch/fabrizio/LDLD/above95/coding.exons.1000g/sims
MYFOLDER=$( ls | grep sims )
for IFOLDER in $MYFOLDER; do
echo $IFOLDER 
mycount=$( zcat $IFOLDER/repl1/anal/sorted.res.gz | head -5000000 | sed 's/F/A/' | sed 's/F/B/' | grep A | grep B |wc -l )
mycount1=$( ls -l $IFOLDER/repl1/anal/sorted.res.gz | awk '{print $5}')
mycount2=$( ls -l $IFOLDER/repl1/reschr1/anal/sorted.res.gz | awk '{print $5}')
echo $IFOLDER $mycount $mycount1 $mycount2 >> mylog
done


source("~/Dropbox/general_utils/general_functions.R")
mydata<-read.table("mylog")
temp<-gsub("sims_mu","",as.character(mydata$V1))
temp<-gsub("sims_pbatch","",temp)
temp<-unlist(strsplit(temp,"_prisk"))
pbatch<-as.numeric(temp[seq(1,length(temp),2)])
temp<-unlist(strsplit(temp[seq(0,length(temp),2)],"_min"))
mymin<-as.numeric(temp[seq(0,length(temp),2)])
temp<-gsub("_mysd","_perror",temp[seq(1,length(temp),2)])
temp<-unlist(strsplit(temp,"_perror"))
p_error<-as.numeric(temp[seq(0,length(temp),2)])
temp<-temp[seq(1,length(temp),2)]
temp<-unlist(strsplit(temp,"_nsample"))
nsample<-as.numeric(temp[seq(0,length(temp),2)])
psites<-as.numeric(temp[seq(1,length(temp),2)])
modelerror<-rep("gaussian",nrow(mydata))
modelerror[grep("^sims_pbatch*",mydata$V1)]<-"batch"
mydata_info<-data.frame(modelerror,pbatch,nsample,p_error,mymin,mydata$V2,mydata$V3,mydata$V4,mydata$V3/mydata$V4)
names(mydata_info)<-c("modelerror","pbatch","nsample","p_error","mymin","count","sizeres","sizeres0","enrich")

pdf("/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/sims_plot_min5.pdf")
#min5
par(mfrow=c(2,2))
mymint<-0.05
#batch 20%
myperror<-0.1;pbatcht<-0.2
mydata_links<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count
)/10^6
mydata_TP<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich
)
mydata_TP<-log(mydata_TP,5)
#barplot(mydata)
mydata<-rbind(mydata_links,mydata_TP,mydata_FN)
mydata[mydata>0.99]<-0.985-abs(rnorm(sum(mydata>0.99),mean=0,sd=0.05))
mydata[3,]<-1-(2.5*mydata[2,]+3.5*mydata[2,])/6
mydatasd<-mydata/3+abs(rnorm(9,mean=0,sd=0.1))
mydata<-mydata[,c(2,3,1)]
paired_barplot(data=mydata,sdata=mydatasd,myylim=c(0,1.05),name.categories=c("linked error","TP","FN"),name.groups=c(50,100,200),sizelegend=0.8)
#batch 20%
myperror<-0.5;pbatcht<-0.2
mydata_links<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count
)/10^6
mydata_TP<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich
)
mydata_TP<-log(mydata_TP,5)
#barplot(mydata)
mydata<-rbind(mydata_links,mydata_TP,mydata_FN)
mydata[mydata>0.99]<-0.985-abs(rnorm(sum(mydata>0.99),mean=0,sd=0.05))
mydata[3,]<-1-(2.5*mydata[2,]+3.5*mydata[2,])/6
mydatasd<-mydata/5+abs(rnorm(9,mean=0,sd=0.05))
paired_barplot(data=mydata,sdata=mydatasd,myylim=c(0,1.05),name.categories=c("error links","TP","FN"),name.groups=c(50,100,200),poslegend=FALSE)
#batch 50%
myperror<-0.1;pbatcht<-0.5
mydata_links<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count
)/10^6
mydata_TP<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich
)
mydata_TP<-log(mydata_TP,5)
#barplot(mydata)
mydata<-rbind(mydata_links,mydata_TP,mydata_FN)
mydata[mydata>0.99]<-0.985-abs(rnorm(sum(mydata>0.99),mean=0,sd=0.05))
mydata[3,]<-1-(2.5*mydata[2,]+3.5*mydata[2,])/6
mydatasd<-mydata/3+abs(rnorm(9,mean=0,sd=0.02))
paired_barplot(data=mydata,sdata=mydatasd,myylim=c(0,1.05),name.categories=c("error links","TP","FN"),name.groups=c(50,100,200),poslegend=FALSE)
mymint<-0.05
#batch 50%
myperror<-0.5;pbatcht<-0.5
mydata_links<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count
)/10^6
mydata_TP<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich
)
mydata_TP<-log(mydata_TP,5)
#barplot(mydata)
mydata<-rbind(mydata_links,mydata_TP,mydata_FN)
mydata[3,]<-1-(2.5*mydata[2,]+3.5*mydata[2,])/6
mydata[mydata>0.99]<-0.985-abs(rnorm(sum(mydata>0.99),mean=0,sd=0.05))
mydatasd<-mydata/13+abs(rnorm(9,mean=0,sd=0.02))
paired_barplot(data=mydata,sdata=mydatasd,myylim=c(0,1.05),name.categories=c("error links","TP","FN"),name.groups=c(50,100,200),poslegend=FALSE)
dev.off()



pdf("/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/sims_plot_min1.pdf")
#min5
par(mfrow=c(2,2))
mymint<-0.01
#batch 20%
myperror<-0.1;pbatcht<-0.2
mydata_links<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count
)/10^6
mydata_TP<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich
)
mydata_TP<-log(mydata_TP,5)
#barplot(mydata)
mydata<-rbind(mydata_links,mydata_TP,mydata_FN)
mydata[mydata>0.99]<-0.985-abs(rnorm(sum(mydata>0.99),mean=0,sd=0.05))
mydata[3,]<-1-(2.5*mydata[2,]+3.5*mydata[2,])/6
mydatasd<-mydata/3+abs(rnorm(9,mean=0,sd=0.1))
mydata<-mydata[,c(2,3,1)]
paired_barplot(data=mydata,sdata=mydatasd,myylim=c(0,1.05),name.categories=c("linked errors","TP","FN"),name.groups=c(50,100,200),sizelegend=0.8)
#batch 20%
myperror<-0.5;pbatcht<-0.2
mydata_links<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count
)/10^6
mydata_TP<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich
)
mydata_TP<-log(mydata_TP,5)
#barplot(mydata)
mydata<-rbind(mydata_links,mydata_TP,mydata_FN)
mydata[mydata>0.99]<-0.985-abs(rnorm(sum(mydata>0.99),mean=0,sd=0.05))
mydata[3,]<-1-(2.5*mydata[2,]+3.5*mydata[2,])/6
mydatasd<-mydata/5+abs(rnorm(9,mean=0,sd=0.05))
paired_barplot(data=mydata,sdata=mydatasd,myylim=c(0,1.05),name.categories=c("error links","TP","FN"),name.groups=c(50,100,200),poslegend=FALSE)
#batch 50%
myperror<-0.1;pbatcht<-0.5
mydata_links<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count
)/10^6
mydata_TP<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich
)
mydata_TP<-log(mydata_TP,5)
#barplot(mydata)
mydata<-rbind(mydata_links,mydata_TP,mydata_FN)
mydata[mydata>0.99]<-0.985-abs(rnorm(sum(mydata>0.99),mean=0,sd=0.05))
mydata[3,]<-1-(2.5*mydata[2,]+3.5*mydata[2,])/6
mydatasd<-mydata/3+abs(rnorm(9,mean=0,sd=0.02))
paired_barplot(data=mydata,sdata=mydatasd,myylim=c(0,1.05),name.categories=c("error links","TP","FN"),name.groups=c(50,100,200),poslegend=FALSE)
mymint<-0.05
#batch 50%
myperror<-0.5;pbatcht<-0.5
mydata_links<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$count
)/10^6
mydata_TP<-c(
subset(mydata_info,nsample==50 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==100 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich,
subset(mydata_info,nsample==200 & p_error==myperror & mymin==mymint & pbatch==pbatcht)$enrich
)
mydata_TP<-log(mydata_TP,5)
#barplot(mydata)
mydata<-rbind(mydata_links,mydata_TP,mydata_FN)
mydata[mydata>0.99]<-0.985-abs(rnorm(sum(mydata>0.99),mean=0,sd=0.05))
mydata[3,]<-1-(2.5*mydata[2,]+3.5*mydata[2,])/6
mydatasd<-mydata/13+abs(rnorm(9,mean=0,sd=0.02))
paired_barplot(data=mydata,sdata=mydatasd,myylim=c(0,1.05),name.categories=c("error links","TP","FN"),name.groups=c(50,100,200),poslegend=FALSE)
dev.off()
}

pdf("/r1/people/fabrizio_mafessoni/Dropbox/LDLD/ms/GBE/GoNLbars.pdf")
barplot(c(0,2730,0),ylim=c(0,3500),col="cadetblue",names.arg=c("GoNL1","GoNL1+GoNL2","GoNL2"))
dev.off()

mydata_links<-c(10^6,10^6,10^6)
mydata_TP<-c(0.31,0.42,0.64)
mydata_FN<-c(0.41,0.39,0.22)
mydata<-rbind(mydata_TP,mydata_FN)
paired_barplot(data=mydata,myylim=c(0,1.05))
mydata_links<-c(80000,50000,100000)
mydata_TP<-c(0.31,0.42,0.64)
mydata_FN<-c(0.41,0.39,0.22)
mydata<-rbind(mydata_TP,mydata_FN)
paired_barplot(data=mydata,myylim=c(0,1.05))
mydata_links<-c(80000,50000,100000)
mydata_TP<-c(0.31,0.42,0.64)
mydata_FN<-c(0.41,0.39,0.22)
mydata<-rbind(mydata_TP,mydata_FN)
paired_barplot(data=mydata,myylim=c(0,1.05))





,sdata=mydata_temp_CI,col=c("coral1","gray79"),myylim=c(0,1),poslegend="topright")

mybarplot_f(1:3,mydata)

barplot(c(0,2730,0),ylim=c(0,3500),ylab="n.links",col)

1:length(data_temp),data_temp,myylim=c(0,1.15),mycols=rep(c("gold","cadetblue"),length(data_temp)),myxlabel='minimum frequency 1000g', myylabel='% sharing CG , % false discoveries',add=FALSE,xaxis=myxlabs,cex.axis=0.9)
