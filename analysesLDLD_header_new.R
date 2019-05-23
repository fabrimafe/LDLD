#==========================================================REQUIRED PACKAGES==========================================
{
require(data.table)
library(stats)
library(doBy)
#library("SNPRelate")
library(Hmisc)
library(mixtools)
#library(vioplot)
}
#==========================================================GENERAL FUNCTIONS==========================================
{
lappend <- function(lst, obj) { #to append to lists
    lst[[length(lst)+1]] <- obj
    return(lst)
}
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, offset=c(0, 0), ...) {
#' @param xfrac The fraction over from the left side.
#' @param yfrac The fraction down from the top.
#' @param label The text to label with.
#' @param pos Position to pass to text()
#' @param ... Anything extra to pass to text(), e.g. cex, col.
#e.g. my.labels <- c(" (a)", " (b)", " (c)"," (d)")
#my.locations <- c("topleft","topleft","topleft","topleft")
#put.fig.letter(label=my.labels[3], location=my.locations[2], font=2)
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.015,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, ...)
}

}
#================================================PREPROCESSING FUNCTIONS (reordering individuals)=====================
{
samples_imypopi<-function(l_npops,myindex) #extract indexes of samples of population with index myindex from a l_samples file subdivided in pop with nindividuals l_npops
{
for (i in 1:myindex)
	{
	if (i==1){init<-1; end<-l_npops[1]} else {init<-init+l_npops[i-1];end<-end+l_npops[i]}
	if (i==myindex){return(init:end)}
	}
}

samples_imypop<-function(l_samples,l_npops,myindex) #extract samples of population with index myindex from a l_samples file subdivided in pop with nindividuals l_npops
{
for (i in 1:myindex)
	{
	if (i==1){init<-1; end<-l_npops[1]} else {init<-init+l_npops[i-1];end<-end+l_npops[i]}
	if (i==myindex){return(l_samples[init:end])}
	}
}

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

order_samples<-function(samples_order,folder,annotation,chrfrom=1,chrto=22)
{
#samples_order: column with ordered samples: e.g. samples_order
#
	for (chr in chrfrom:chrto)
	{
	system(paste0("vcf-subset -c ",samples_order," ",folder,"/chr",chr,".tab | awk '{if (NR>2){print}}' | gzip -9 > ",folder,"/chr",chr,".",annotation,".vcf.gz"))
	}
}

order_samples_vcf.gz<-function(samples_order,folder,annotation,chrfrom=1,chrto=22)
{
#samples_order: column with ordered samples: e.g. samples_order
#
	for (chr in chrfrom:chrto)
	{
	system(paste0("vcf-subset -c ",samples_order," ",folder,"/chr",chr,".",annotation,".gz | awk '{if (NR>2){print}}' | grep -v '##' > " ,folder,"/chr",chr,".",annotation))
	system(paste0("gzip -f ",folder,"/chr",chr,".",annotation))
	}
}

}
#==========================================================IMPORTING LDLD results=====================================
{
#cat YRI/missense/anal/*ann | sort -n -k15 | awk '{if (NF<19){for(i=1;i<=NF;i++){printf "%s ", $i}; printf "Nan\n"}else{print}}' | awk '{if (NF<19){for(i=1;i<=NF;i++){printf "%s ", $i}; printf "Nan\n"}else{print}}' | sort -u >YRI/missense/anal/res &
f_import_fdr<-function(folder,criterion_fdr)
{
  mispathfile<-paste0(folder,"/missense/anal/sorted.ann")
  synpathfile<-paste0(folder,"/synonymous/anal/sorted.ann")
  intpathfile<-paste0(folder,"/intergenic/anal/sorted.res")
  mispathfold<-paste0("/mnt/scratch/fabrizio/LDLD/counterSNPSgenome.sh ",folder,"/missense/")
  synpathfold<-paste0("/mnt/scratch/fabrizio/LDLD/counterSNPSgenome.sh ",folder,"/synonymous/")
  intpathfold<-paste0("/mnt/scratch/fabrizio/LDLD/counterSNPSgenome.sh ",folder,"/intergenic/")
  YRI.missense<-read.table(mispathfile)
  YRI.intergenic<-read.table(intpathfile)
  YRI.synonymous<-read.table(synpathfile)
  mynames<-c("chrA","chrB","snpA","snpB","nA","nB","N","nAB","nAA","nBB","D","rAB","Drandom","rABrandom","ppermD","pfisher","pKuli09","pchi","LLDA","LLDB","ENSA","HGNCA","ENTREZA","mapA","ENSB","HGNCB","ENTREZB","mapB")
  mynames_short<-c("chrA","chrB","snpA","snpB","nA","nB","N","nAB","nAA","nBB","D","rAB","Drandom","rABrandom","ppermD","pfisher","pKuli09","pchi","LLDA","LLDB")
  names(YRI.missense)<-mynames
  names(YRI.intergenic)<-mynames_short
  names(YRI.synonymous)<-mynames
  head(YRI.missense)
  #total number of comparison according to LDLDncomparisons is 530976471
  YRI.missense.ncomp<-as.numeric(system(mispathfold,intern=TRUE))
  YRI.intergenic.ncomp<-as.numeric(system(intpathfold,intern=TRUE))
  YRI.synonymous.ncomp<-as.numeric(system(synpathfold,intern=TRUE))
  print("number of snps pairs")
  print(c(YRI.missense.ncomp,YRI.intergenic.ncomp,YRI.synonymous.ncomp))
  print("number of snps pairs in files")
  print(c(length(YRI.missense$pfisher),length(YRI.synonymous$pfisher),length(YRI.intergenic$pfisher)))
  if (criterion_fdr=="pfisher")
    {
    new.missense<-YRI.missense[order(YRI.missense$pfisher),]
    new.intergenic<-YRI.intergenic[order(YRI.intergenic$pfisher),]
    new.synonymous<-YRI.synonymous[order(YRI.synonymous$pfisher),]
    new.missense.fdr<-p.adjust(p=new.missense$pfisher, method = "fdr", n = YRI.missense.ncomp)
    new.intergenic.fdr<-p.adjust(p=new.intergenic$pfisher, method = "fdr", n = YRI.intergenic.ncomp)
    new.synonymous.fdr<-p.adjust(p=new.synonymous$pfisher, method = "fdr", n = YRI.synonymous.ncomp)
    }
  else if (criterion_fdr=="pKuli09")
    {
    new.missense<-YRI.missense[order(YRI.missense$pKuli09),]
    new.intergenic<-YRI.intergenic[order(YRI.intergenic$pKuli09),]
    new.synonymous<-YRI.synonymous[order(YRI.synonymous$pKuli09),]
    new.missense.fdr<-p.adjust(p=new.missense$pKuli09, method = "fdr", n = YRI.missense.ncomp)
    new.intergenic.fdr<-p.adjust(p=new.intergenic$pKuli09, method = "fdr", n = YRI.intergenic.ncomp)
    new.synonymous.fdr<-p.adjust(p=new.synonymous$pKuli09, method = "fdr", n = YRI.synonymous.ncomp)
    }
  new.missense<-cbind(new.missense,new.missense.fdr)
  new.synonymous<-cbind(new.synonymous,new.synonymous.fdr)
  new.intergenic<-cbind(new.intergenic,new.intergenic.fdr)
  names(new.missense)<-c(mynames,"fdr")
  names(new.synonymous)<-c(mynames,"fdr")
  names(new.intergenic)<-c(mynames_short,"fdr")
  summary_data<-c(length(new.missense.fdr[new.missense.fdr<0.05]),length(new.intergenic.fdr[new.intergenic.fdr<0.05]),length(new.synonymous.fdr[new.synonymous.fdr<0.05]))
  res<-list(summary_data,new.missense[new.missense.fdr<0.05,],new.intergenic[new.intergenic.fdr<0.05,],new.synonymous[new.synonymous.fdr<0.05,])
  return(res)
}

f_import_fdr_multi<-function(folder,pop,criterion_fdr,npops)
{
  mispathfile<-paste0(folder,"/anal/",pop)
  print(mispathfile)
  mispathfold<-paste0("/mnt/scratch/fabrizio/LDLD/counterSNPSgenome.sh ",folder)
  YRI.missense<-read.table(mispathfile)
  names(YRI.missense)<-c("chrA","chrB","snpA","snpB","nA","nB","N","nAB","D","rAB","ppermD","pfisher","pKuli09","chiaggr","pop")
  YRI.missense$chiaggr<-1-pchisq(YRI.missense$chiaggr,df=npops*2)
  #total number of comparison according to LDLDncomparisons is 530976471
  YRI.missense.ncomp<-as.numeric(system(mispathfold,intern=TRUE))
  print("number of snps pairs")
  print(YRI.missense.ncomp)
  print("number of snps pairs in files")
  print(length(YRI.missense$pfisher))
  if (criterion_fdr=="pfisher")
    {
    new.missense<-YRI.missense[order(YRI.missense$pfisher),]
    new.missense.fdr<-p.adjust(p=new.missense$pfisher, method = "fdr", n = YRI.missense.ncomp)
    }
  else if (criterion_fdr=="chiaggr")
    {
    new.missense<-YRI.missense[order(YRI.missense$chiaggr),]
    new.missense.fdr<-p.adjust(p=new.missense$chiaggr, method = "fdr", n = YRI.missense.ncomp)
    }
  new.missense<-cbind(new.missense,new.missense.fdr)
  names(new.missense)<-c(names(YRI.missense),"fdr")
  summary_data<-length(new.missense.fdr[new.missense.fdr<0.05])
  res<-list(summary_data,new.missense[new.missense.fdr<0.05,])
  return(res)
}

organize_log_files<-function(myfolder,npops)
{
system(paste0("mkdir ",myfolder,"/logs"))
for (ipop in 0:(npops-1))
  {
  system(paste0("if [ -e ",myfolder,"/logs/all.pop ];then rm",myfolder,"/logs/all.pop",ipop,".log;fi "))
  for (ichrA in 1:21)
    {
   for (ichrB in (ichrA+1):22)
	{
	system(paste0("cat ",myfolder,"/chr",ichrA,".tabchr",ichrB,".tab.log | awk -v OFS='\t' '{if (NF>17 && $NF==",ipop,"){print ",ichrA,",",ichrB,",$0}}' | sort -Vu >> ",myfolder,"/logs/all.pop",ipop,".log"))
	}
    };
system(paste0("gzip -f ", myfolder,"/logs/all.pop",ipop,".log"))
  }
}

data2snpA<-function(dataLDLD){
tagsA<-apply(dataLDLD[,c("chrA","posA"),with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse="."))
return(tagsA)
}
data2snpB<-function(dataLDLD){
tagsA<-apply(dataLDLD[,c("chrB","posB"),with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse="."))
return(tagsA)
}

data2pairid<-function(dataLDLD,myorder=c(1,3,2,4)){
tagsA<-apply(dataLDLD[,myorder,with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse="."))
return(tagsA)
}

intersect_links<-function(dataLDLDA,dataLDLDB,B.is.pairid=FALSE,myorder=c(1,3,2,4)) #return dataLDLDA filtered for links present in dataLDLDB
{
tagsA<-apply(dataLDLDA[,myorder,with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse="."))
if ( B.is.pairid ) { tagsB<-dataLDLDB } else { tagsB<-apply(dataLDLDB[,myorder,with=FALSE], MARGIN=1,FUN=function(x) paste(x,collapse=".")) } 
return(dataLDLDA[!is.na(match(tagsA,tagsB)),])
}

intersect2pops<-function(data1,data2)
{
	  a1 <-data1[,1:4]
		  a2 <- data2[,1:4]
#a2 <- data1[1:100,1:4]
#sqldf('SELECT * FROM m1 EXCEPT SELECT * FROM m2') #for in m1 not in m2
		  a1Ina2 <- sqldf('SELECT * FROM a1 INTERSECT SELECT * FROM a2')
		  res<-merge(a1Ina2,data1,by=c("chrA","chrB","snpA","snpB"))
		  return(res)
}

dataLDLD2Dstr12_f<-function(mydata,dataFIN_sign,mycriterion="all"){ #function that takes a dataLDLD (raw results in data.table format), subsets them according to a data.table dataFIN_sign and return a list with data.tables used to generate nAB according to one of the possible criteria
require(data.table)
pairid<-mydata[, paste(chrA,chrB,posA,posB,sep=".")]
pairid_sign<-dataFIN_sign[, paste(chrA,chrB,posA,posB,sep=".")]
mydata<-mydata[!is.na(match(pairid,pairid_sign)),]
pairid<-mydata[, paste(chrA,chrB,posA,posB,sep=".")]
data_neg<-aggregate(mydata$D,by=list(pairid),FUN=function(x) sum(x<0))[,2]
data_pos<-aggregate(mydata$D,by=list(pairid),FUN=function(x) sum(x>0))[,2]
#data_Dsum<-aggregate(mydata$D,by=list(pairid),FUN=function(x) sum(x))[,2]
#data_Dsum<-aggregate(mydata$D,by=list(pairid),FUN=function(x) sum(x))
mydata<-cbind(mydata,sumlogT2=0)
mydata$sumlogT2[mydata$D<0]<- log(mydata$T2)[mydata$D<0]
mydata$sumlogT2[mydata$D>0]<- -log(mydata$T2)[mydata$D>0]
data_sumlogT2<-aggregate(mydata$sumlogT2,by=list(pairid),FUN=function(x) sum(x))
#data_sumlogT2smaller0<-aggregate(mydata$sumlogT2,by=list(pairid),FUN=function(x) sum(x[x<0]+0))
#data_sumlogT2greater0<-aggregate(mydata$sumlogT2,by=list(pairid),FUN=function(x) sum(x[x>0]+0))

#Dstr1:
myrho_threshold<-unname(quantile(mydata$rho2,0.75))
mydata_Dstr1<-data.table(mydata[mydata$rho2>myrho_threshold],pop=mydata$pop[mydata$rho2>myrho_threshold])
mydata_Dstr1_pos<-data.table(mydata_Dstr1[,1:4,with=FALSE][mydata_Dstr1$D>0],pop=mydata_Dstr1$pop[mydata_Dstr1$D>0])
mydata_Dstr1_neg<-data.table(mydata_Dstr1[,1:4,with=FALSE][mydata_Dstr1$D<0],pop=mydata_Dstr1$pop[mydata_Dstr1$D<0])
#Dstr2:
#n pops is sum(pairid==pairid[1])
npopulations<-sum(pairid==pairid[1])
tempindex<-(as.numeric(mydata$D>0) + as.numeric(c(sapply(data_sumlogT2$x,function(x) rep(x,npopulations)))>0))==2
mydata_Dstr2_pos<-data.table(mydata[,1:4,with=FALSE][tempindex],pop=mydata$pop[tempindex])
tempindex<-(as.numeric(mydata$D<0) + as.numeric(c(sapply(data_sumlogT2$x,function(x) rep(x,npopulations)))<0))==2
mydata_Dstr2_neg<-data.table(mydata[,1:4,with=FALSE][tempindex],pop=mydata$pop[tempindex])
#union
#to filter unique I can set key to all columns
mydata_Dstr12neg<-rbind(mydata_Dstr2_neg,mydata_Dstr1_neg) 
mydata_Dstr12pos<-rbind(mydata_Dstr2_pos,mydata_Dstr1_pos)
mydata_Dstr12neg<-unique(mydata_Dstr12neg,by=c('chrA','chrB','posA','posB','pop')) #unique by multiple columns
mydata_Dstr12pos<-unique(mydata_Dstr12pos,by=c('chrA','chrB','posA','posB','pop'))
mydata_Dstr12neg<-mydata_Dstr12neg[order(chrA,chrB,posA,posB)] #sort by multiple columns
mydata_Dstr12pos<-mydata_Dstr12pos[order(chrA,chrB,posA,posB)]
if (mycriterion=="1"){return(list(mydata_Dstr1neg,mydata_Dstr1pos))} else if (mycriterion=="2"){return(list(mydata_Dstr2neg,mydata_Dstr2pos))} else if (mycriterion=="12"){return(list(mydata_Dstr12neg,mydata_Dstr12pos))} else
{return(list(mydata_Dstr1_neg,mydata_Dstr1_pos,mydata_Dstr2_neg,mydata_Dstr2_pos,mydata_Dstr12neg,mydata_Dstr12pos))}
}


}
#====================================VCF MANIPULATION AND SUBSTRUCTURE (fst,PCA)======================================
{
f_generate_vcf<-function(dataset,target_folder)
{
  data<-cbind(dataset$chrA,dataset$snpA,dataset$snpB)
  data2<-data[firstobs(data[, 3]), ]
  data<-data2[firstobs(data2[, 2]), ][,1:2]  
#data<-unique(cbind(dataset$chrA,dataset$snpA))
  dest_file<-paste0(target_folder,"sign.vcf")
  system(paste0("head -1 ",target_folder,"chr21.tab >",dest_file))
  for (ic in 1:22)
  {
    datachr<-data[,2][data[,1]==ic]
    if (length(datachr)>0)
      {
      for (i in 1:length(datachr))
	{
	if (i==1){target_str<-paste0("awk '{if($2==",datachr[1])}
	if (i>1){target_str<-paste0(target_str,"||$2==",datachr[i])}
	}
      target_str<-paste0(target_str,"){print}}' ")
      target_command<-paste0(target_str,target_folder,"chr",ic,".tab >> ",dest_file)    
 #     print(target_command)
      system(target_command)    
      }
  }
  return(system(paste0("cat ",dest_file,"| wc -l")))
}

f_generate_random_vcf<-function(dataset,target_folder)
{
  data<-cbind(dataset$chrA,dataset$snpA,dataset$snpB)
  data2<-data[firstobs(data[, 3]), ]
  data<-data2[firstobs(data2[, 2]), ][,1:2]
#  data<-unique(cbind(dataset$chrA,dataset$snpA))
  #data<-unique(cbind(YRI.F[[2]]$chrA,YRI.F[[2]]$snpA))
  #target_folder<-"/mnt/scratch/fabrizio/LDLD/YRI/missense/"
  dest_file<-paste0(target_folder,"random.vcf")
  system(paste0("head -1 ",target_folder,"chr21.tab >",dest_file))
  for (ic in 1:22)
  {
    datachr<-data[,2][data[,1]==ic]
    if (length(datachr)>0)
      {
      target_command<-paste0("grep -v FILTER ",target_folder,"chr",ic,".tab | shuf -n ", length(datachr), " >> ",dest_file)
#      print(target_command)
      system(target_command)
      }
  }
  return(system(paste0("cat ",dest_file,"| wc -l")))
}

f_PCA<-function(vcf_file,yes.plot)
{
snpgdsVCF2GDS(vcf_file, "sign.gds",  method="biallelic.only")
sign.obj <- openfn.gds("sign.gds")
sign_pca<-snpgdsPCA(sign.obj)
sign_pca$eigenval
if (yes.plot) {plot(sign_pca$eigenvect[,1],sign_pca$eigenvect[,2] ,col="ivory4",#as.numeric(substr(sign_pca$sample, 1,3) == 'CCM')+3
pch=19,xlab="PCA1",ylab="PCA2")};
closefn.gds(sign.obj)
return((sign_pca$eigenval[1]+sign_pca$eigenval[2])/sum(sign_pca$eigenval))
}
#f_PCA(sign.vcf,TRUE)

f_PCA_random<-function(data,target_folder,nrep)
{
dest_file<-paste0(target_folder,"random.vcf")
for (i in 1:nrep)
  {
  f_generate_random_vcf(data,target_folder)
  res<-f_PCA(dest_file,FALSE)
  if (i==1){eigenvalres<-res}
  else {eigenvalres<-rbind(eigenvalres,res)}
  }
  return(eigenvalres)
}

f_PCA_multipop<-function(vcf_file,yes.plot,color)
{
snpgdsVCF2GDS(vcf_file, "sign.gds",  method="biallelic.only")
sign.obj <- openfn.gds("sign.gds")
sign_pca<-snpgdsPCA(sign.obj)
sign_pca$eigenval
if (yes.plot) {plot(sign_pca$eigenvect[,1],sign_pca$eigenvect[,2] ,col=color,#as.numeric(substr(sign_pca$sample, 1,3) == 'CCM')+3
pch=19,xlab="PCA1",ylab="PCA2")};
closefn.gds(sign.obj)
return((sign_pca$eigenval[1]+sign_pca$eigenval[2])/sum(sign_pca$eigenval))
}

f_generate_random_vcf_fromsignvcf<-function(dataset,target_folder)
{
  dest_file<-paste0(target_folder,"random.vcf")
  system(paste0("head -1 ",target_folder,"chr21.tab >",dest_file))
  for (ic in 1:22)
  {
    com2<-paste0("cat ",target_folder,"sign.vcf | awk '{if ($1==",ic,"){print}}' | wc -l")
#    print(com2)
    nperchr<-as.numeric(system(com2,intern=TRUE))
    if (nperchr>0)
      {
      target_command<-paste0("grep -v FILTER ",target_folder,"chr",ic,".tab | shuf -n ", nperchr, " >> ",dest_file)
#      print(target_command)
      system(target_command)
      }
  }
  return(system(paste0("cat ",dest_file,"| wc -l")))
}

f_PCA_random_fromvcf<-function(data,target_folder,nrep)
{
dest_file<-paste0(target_folder,"random.vcf")
for (i in 1:nrep)
  {
  f_generate_random_vcf_fromsignvcf(data,target_folder)
  res<-f_PCA(dest_file,FALSE)
  if (i==1){eigenvalres<-res}
  else {eigenvalres<-rbind(eigenvalres,res)}
  }
  return(eigenvalres)
}

randomfst<-function(data,nrep)
{
  res<-rep(0,nrep)
  for (i in 1:nrep)
  {
  f_PCA_random_fromvcf(data,"/mnt/scratch/fabrizio/LDLD/YRILWKTSIFINCHBJPT/missense/",1)
  system("head -1 YRILWKTSIFINCHBJPT/missense/random.vcf > random2.vcf")
  system("sort -n -k1 -n -k2 YRILWKTSIFINCHBJPT/missense/random.vcf >> random2.vcf")
  system("cat headercvf random2.vcf > random.vcf")
  res[i]<-as.numeric(system("vcftools --vcf random.vcf --weir-fst-pop ~/workspace/1000genomes/YRI.unrelated.samples --weir-fst-pop ~/workspace/1000genomes/LWK.unrelated.samples --weir-fst-pop ~/workspace/1000genomes/TSI.unrelated.samples --weir-fst-pop ~/workspace/1000genomes/FIN.unrelated.samples --weir-fst-pop ~/workspace/1000genomes/CHB.unrelated.samples --weir-fst-pop ~/workspace/1000genomes/JPT.unrelated.samples | grep weighted | awk '{print $7}'",intern=TRUE))
  }
return(res)
}

f_generate_vcf_1snp<-function(dataset,target_folder) #to calculate fst 1 snp at the time
{
  data<-unique(rbind(cbind(dataset$chrA,dataset$snpA),cbind(dataset$chrB,dataset$snpB)))
  dest_file<-paste0(target_folder,"sign.vcf")
  res<-rep(0,dim(data)[1])
  for (i in 1:dim(data)[1])
  {
      system(paste0("head -1 ",target_folder,"chr21.tab >",dest_file))
      target_str<-paste0("awk '{if($2==",data[i,2]," && $1==",data[i,1])
      target_str<-paste0(target_str,"){print}}' ")
      target_command<-paste0(target_str,target_folder,"chr",data[i,1],".tab >> ",dest_file)    
      system(target_command)    
      res[i]<-as.numeric(system(paste0("cat ",dest_file,"| wc -l"),intern=TRUE)) #no, not yet for fst
  }
  return(res)
}
}
#================================SFS_FUNCTIONS========================================================================
{
f_importSFS<-function(file)
{
data<-read.table(file)[[1]]
data<-sapply(data,function(x) {if(x>0.5){return(1-x)}else{return(x)}})
return(data)
}

px<-function(folder,file_nx)#nx is SFS for significant SNPs
{
  vec<-seq(0.05,0.5,0.01)
  fileSFS.mis<-paste0(folder,"/missense/SFSfile")
  fileSFS.syn<-paste0(folder,"/synonymous/SFSfile")
  fileSFS.int<-paste0(folder,"/intergenic/SFSfile")
  data_SFS.mis<-f_importSFS(fileSFS.mis)
  data_SFS.syn<-f_importSFS(fileSFS.syn)
  data_SFS.int<-f_importSFS(fileSFS.int)
  #taking nx of missense
  data<-c(file_nx[[2]]$nA/file_nx[[2]]$N,file_nx[[2]]$nB/file_nx[[2]]$N)
  data<-sapply(data,function(x) {if(x>0.5){return(1-x)}else{return(x)}})
  nx<-hist(data,breaks=vec)$density/sum(hist(data,breaks=vec)$density)
  SFSx.mis<-hist(data_SFS.mis,breaks=vec)$density/sum(hist(data_SFS.mis,breaks=vec)$density)
  SFSx.syn<-hist(data_SFS.syn,breaks=vec)$density/sum(hist(data_SFS.syn,breaks=vec)$density)
  SFSx.int<-hist(data_SFS.int,breaks=vec)$density/sum(hist(data_SFS.int,breaks=vec)$density)
  #sum(nx)  #sum(SFSx)  #
  px.mis<-nx*length(data)/(SFSx.mis*length(data_SFS.mis))
  #px.syn<-nx*length(data)/(SFSx.syn*length(data_SFS.syn))
  #px.int<-nx*length(data)/(SFSx.int*length(data_SFS.int))
  return(c(sum(px.mis*(SFSx.mis*length(data_SFS.mis))),sum(px.mis*(SFSx.syn*length(data_SFS.syn))), sum(px.mis*(SFSx.int*length(data_SFS.int)))))
  #return(c(sum(px.mis*(SFSx.mis)),sum(px.mis*(SFSx.syn)), sum(px.mis*(SFSx.int))))
}

px<-function(folder,file_nx)#nx is SFS for significant SNPs
{
  vec<-seq(0.05,0.5,0.01)
  fileSFS.mis<-paste0(folder,"/missense/SFSfile")
  fileSFS.syn<-paste0(folder,"/synonymous/SFSfile")
  fileSFS.int<-paste0(folder,"/intergenic/SFSfile")
  data_SFS.mis<-f_importSFS(fileSFS.mis)
  data_SFS.syn<-f_importSFS(fileSFS.syn)
  data_SFS.int<-f_importSFS(fileSFS.int)
  #taking nx of missense
  data<-c(file_nx[[2]]$nA/file_nx[[2]]$N,file_nx[[2]]$nB/file_nx[[2]]$N)
  data<-sapply(data,function(x) {if(x>0.5){return(1-x)}else{return(x)}})
  nx<-hist(data,breaks=vec)$density/sum(hist(data,breaks=vec)$density)
  SFSx.mis<-hist(data_SFS.mis,breaks=vec)$density/sum(hist(data_SFS.mis,breaks=vec)$density)
  SFSx.syn<-hist(data_SFS.syn,breaks=vec)$density/sum(hist(data_SFS.syn,breaks=vec)$density)
  SFSx.int<-hist(data_SFS.int,breaks=vec)$density/sum(hist(data_SFS.int,breaks=vec)$density)
  #sum(nx)  #sum(SFSx)  #
  px.mis<-nx/(SFSx.mis*length(data_SFS.mis))
  #px.syn<-nx*length(data)/(SFSx.syn*length(data_SFS.syn))
  #px.int<-nx*length(data)/(SFSx.int*length(data_SFS.int))
  return(c(sum(px.mis*(SFSx.mis*length(data_SFS.mis))),sum(px.mis*(SFSx.syn*length(data_SFS.syn))), sum(px.mis*(SFSx.int*length(data_SFS.int)))))
  #return(c(sum(px.mis*(SFSx.mis)),sum(px.mis*(SFSx.syn)), sum(px.mis*(SFSx.int))))
}

f_importSFS_folded<-function(file) 
{
data<-read.table(file)
res<-sapply(data[,3],function(x) {if(x>0.5){return(1-x)}else{return(x)}})
return(res)
}

f_importSFS_unfolded<-function(file) #new version from SFS.sh and then cat
{
data<-read.table(file)
res<-data[,3]
return(res)
}
}
#====================================LINKS PLOTTING_FUNCTIONS=========================================================
#==============================GOandTISSUESPECIFICITY_FUNCTIONS=======================================================
{
f_extractgenesHGNC<-function(dataset,col=FALSE)
{
resA<-sapply(dataset$HGNCA,function(x) {strsplit(as.character(x),";")})
resB<-sapply(dataset$HGNCB,function(x) {strsplit(as.character(x),";")})
if (col==1) {res<-sapply(resA,function(x) x[1])} else if (col==2) {res<-sapply(resB,function(x) x[1])} else {res<-c(sapply(resA,function(x) x[1]),sapply(resB,function(x) x[1]))}
return(res)
}

f_extractgenesENSG<-function(dataset,col=FALSE)
{
resA<-sapply(dataset$ENSA,function(x) {strsplit(as.character(x),";")})
resB<-sapply(dataset$ENSB,function(x) {strsplit(as.character(x),";")})
if (col==1) {res<-sapply(resA,function(x) x[1])} else if (col==2) {res<-sapply(resB,function(x) x[1])} else {res<-c(sapply(resA,function(x) x[1]),sapply(resB,function(x) x[1]))}
return(res)
}

f_extractgenes<-function(dataset)
{
resA<-sapply(dataset$HGNCA,function(x) {strsplit(as.character(x),";")})
resB<-sapply(dataset$HGNCB,function(x) {strsplit(as.character(x),";")})
res<-cbind(sapply(resA,function(x) x[1]),sapply(resB,function(x) x[1]))
return(res)
}

clean_Nan_genepairlist<-function(res) 
{
  genepairs<-c(res[1,1],res[1,2])
  for (i in 2:dim(res)[1])
    {
    if (res[i,1]!="Nan" && res[i,2]!="Nan")
      {
      genepairs<-rbind(genepairs,c(res[i,1],res[i,2]))
      }
    }
  return(genepairs)
}

random_pairs<-function(res)
{return(cbind(sample(unique(HGNCtogo$HGNC),dim(res)[1]),sample(unique(HGNCtogo$HGNC),dim(res)[1])))}

same_GO<-function(res)
{
tot<-0
for (i in 1:dim(res)[1])
  {
  resa<-HGNCtogo$GO[HGNCtogo$HGNC==res[i,1]]
  resb<-HGNCtogo$GO[HGNCtogo$HGNC==res[i,2]]
  tot=tot+length(intersect(resa,resb))
  }
return(tot)
}

gimme_GO<-function()
for (i in listgenes)
{
if(j==1){GOgene<-HGNCtogo$GO[HGNCtogo$HGNC==i]}
if(j>1){GOgene<-c(GOgene,HGNCtogo$HGNC[HGNCtogo$HGNC==i])}
j<-j+1
}

gimme_GO<-function(listIDs,IDgene2GO)
{
for (i in listIDs)
  {
  if(j==1){GOgene<-IDgene2GO[,2][IDgene2GO[,1]==i]}
  if(j>1){GOgene<-c(GOgene,IDgene2GO[,2][IDgene2GO[,1]==i])}
  j<-j+1
  }
return(GOgene)
}

load.Micha.gene.expr.data<-function()
  {i<-1;gene.expr.dat<-as.list(rep(0,length(tissues)));
  # looping through all tissues
  tissues<<-c("adipose","adrenal","blood","brain","breast","colon","heart","kidney","liver","lung","lymph","ovary","prostate","skeletal_muscle","testes","thyroid")
  for(tis in tissues){
  # loading the differential expression data
  load(paste("/r1/people/michael_dannemann/R/rdat/",tis,"se_deseq.rdat"))
  #max(gene.expr.dat[[1]]$foldChange[gene.expr.dat[[1]]$foldChange!="Inf"],na.rm=TRUE) #->26402.78  
  #  res$foldChange[res$foldChange!="Inf"]<-3000
  gene.expr.dat[[i]]<-res
  i=i+1}
  return(gene.expr.dat)
  }

pairwise.gene.corr<-function(res)
  {
  v_corr<-rep(NA,dim(res)[1])
  for (i in 1:length(v_corr))
#  for (i in 1:148)
    {
    if(res[i,1]=="Nan" || res[i,2]=="Nan" || is.na(res[i,2]) || is.na(res[i,1]) || (length(gene.expr.dat[[1]]$id[gene.expr.dat[[1]]$id==res[i,1]])==0) || (length(gene.expr.dat[[1]]$id[gene.expr.dat[[1]]$id==res[i,2]])==0)){print("missing gene")}
    else
      {	
      pair<-cbind(expr.per.gene(res[i,1]),expr.per.gene(res[i,2]))
      if (sum(!is.na(pair[,1]))>5 && sum(!is.na(pair[,2]))>5 )
	{
	v_corr[i]<-rcorr(pair, type="spearman")[[1]][1,2]
	}
      else
	v_corr[i]<-NA
      }
    }
  return(v_corr)
  }

expr.per.gene<-function(gene)
  {
  a<-sapply(1:length(tissues),function(i) gene.expr.dat[[i]]$foldChange[gene.expr.dat[[1]]$id==gene])
  b<-sapply(a,function(x) {if (x=="NaN" || x==-Inf || x==Inf){return(NA)}
  else{return(x)}})
  return(b)
  }

f_random_corr<-function(genepairs,nrep)
{
  tot<-rep(0,nrep)
  for(i in 1:nrep)
  {
    genepair.r<-random_pairs(genepairs)
    v_cor<-pairwise.gene.corr(genepair.r)
    tot[i]<-mean(v_cor[!is.na(v_cor)]^2)
  }
  return(tot)
}

isAinB<-function(A,B)
{
res<-matrix(0,nrow=length(A),ncol=2)
res[,1]<-A
for (i in 1:length(A))
  {
  res[i,2]<-length(B[B==A[i]])
  }
return(res)
}

load_tissueHAP<-function(HAPfile_raw,onlyAutosomes=TRUE,excluded_overlapping_tissue=FALSE)
{
  if (excluded_overlapping_tissue==FALSE) {excluded_overlapping_tissue<-"supercalifragidistichespiralitoso"}
  commHGNC<-paste0("cat '",HAPfile_raw,"' | grep -v ",excluded_overlapping_tissue, " | awk -v FS='\t' '{print $1}'")
  commENSG<-paste0("cat '",HAPfile_raw,"' | grep -v ",excluded_overlapping_tissue, " | awk -v FS='\t' '{print $3}'")
  commCHR<-paste0("cat '",HAPfile_raw,"' | grep -v ",excluded_overlapping_tissue, " | awk -v FS='\t' '{print $5}'")
  commENRICH<-paste0("cat '",HAPfile_raw,"' | grep -v ",excluded_overlapping_tissue, " | awk -v FS='\t' '{print $17}'")
  HGNC<-system(commHGNC,intern=TRUE)
  ENSG<-system(commENSG,intern=TRUE)
  CHR<-system(commCHR,intern=TRUE)
  ENRICH<-system(commENRICH,intern=TRUE)
  res<-cbind(ENSG,HGNC,CHR,ENRICH)
  if (onlyAutosomes) {res<-subset(res,res[,3]!="X")}
  res<-res[2:dim(res)[1],]
  return(res)
}

load_tissuespecificity_data<-function()
{
  testis<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:testis;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  fallopian<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:fallopian+tube;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  placenta<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:placenta;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  cortex<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:cerebral+cortex;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  ovary<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:ovary;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  testis_nofall<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:testis;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE,"fallopian")
  fall_notestis<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:fallopian+tube;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE,"testis")
  endo<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:endometrium;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  stomach<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:stomach;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  bonemarrow<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:bone+marrow;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  duodenum<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:duodenum;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE) 
  pancreas<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:pancreas;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  liver<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:liver;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  kidney<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:kidney;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)    
  skin<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:skin;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  colon<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:colon;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  prostate<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:prostate;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  gallbladder<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:gallbladder;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)                                
  lung<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:lung;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  muscle<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:skeletal+muscle;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)
  intestine<-load_tissueHAP("/mnt/expressions/HumanProteinAtlas/tissuespecificity/tissue_specificity_rna:small+intestine;+AND+sort_by:tissue+specific+score?format=tab",onlyAutosomes=TRUE)                                
return(list(testis,fallopian,placenta,cortex,ovary,endo,stomach,bonemarrow,duodenum,pancreas,liver,kidney,skin,colon,prostate,gallbladder,lung,muscle,intestine))
}


tissuepairing_enrichment<-function(data.pairs,testis_nofall,fall_notestis)
{
  testis_nofallA<-as.numeric(isAinB(data.pairs[,1],testis_nofall[,1])[,2])
  testis_nofallB<-as.numeric(isAinB(data.pairs[,2],testis_nofall[,1])[,2])
  fall_notestisA<-as.numeric(isAinB(data.pairs[,1],fall_notestis[,1])[,2])
  fall_notestisB<-as.numeric(isAinB(data.pairs[,2],fall_notestis[,1])[,2])
  unique_fallAtestisB<-as.numeric(fall_notestisA)+as.numeric(testis_nofallB)
  unique_fallBtestisA<-as.numeric(fall_notestisB)+as.numeric(testis_nofallA)
  indices<-c(which(unique_fallAtestisB==2),which(unique_fallBtestisA==2))[order(c(which(unique_fallAtestisB==2),which(unique_fallBtestisA==2)))]
  res<-matrix(nrow=length(indices),ncol=3,"0")
  counter<-1
  for (i in indices)
    {
    res[counter,1]<-as.character(data.pairs[i,1])
    res[counter,2]<-as.character(data.pairs[i,2])
    res[counter,3]<-as.character(data.pairs[i,3])
    counter<-counter+1
    }
  return(res)
}

tissuepairing_enrichment_all<-function(data.pairs,data.pairs2,data.pairs3,list_tissues) #gives single tissue pairings enrichment for all tissues
{
res<-matrix(nrow=length(list_tissues),ncol=3,0)
counter<-1
for (itis in list_tissues)
  {
  res1<-dim(tissuepairing_enrichment(data.pairs,itis,itis))[1]/dim(data.pairs)[1]
  res2<-dim(tissuepairing_enrichment(data.pairs2,itis,itis))[1]/dim(data.pairs2)[1]
  res3<-dim(tissuepairing_enrichment(data.pairs3,itis,itis))[1]/dim(data.pairs3)[1]
  res[counter,]<-c(res1,res2,res3)
  counter<-counter+1
  }
return(res)
}

tissuepairing_enrichment_global<-function(data.pairs,data.pairs2,data.pairs3,list_tissues) #gives single tissue pairings enrichment for all tissues
{
res1<-matrix(nrow=1,ncol=3,0)
res2<-matrix(nrow=1,ncol=3,0)
res3<-matrix(nrow=1,ncol=3,0)
counter<-1
for (itis in list_tissues)
  {
  res1<-rbind(res1,tissuepairing_enrichment(data.pairs,itis,itis))
  res2<-rbind(res2,tissuepairing_enrichment(data.pairs2,itis,itis))
  res3<-rbind(res3,tissuepairing_enrichment(data.pairs3,itis,itis))
  }
  res1<-f_clean_pairedgenes3col(res1)
  res2<-f_clean_pairedgenes3col(res2)
  res3<-f_clean_pairedgenes3col(res3)
res<-c(dim(res1)[1]/dim(data.pairs)[1],dim(res2)[1]/dim(data.pairs2)[1],dim(res3)[1]/dim(data.pairs3)[1])/length(list_tissues)
resraw<-c(dim(res1)[1],dim(res2)[1],dim(res3)[1])
return(list(rbind(resraw,res),res1,res2,res3))
}
}
#----------------------------INTERACTIONS-----------------------------------------------------------------------------
{
#STRINGdb
f_extract_pairedgenes_ens<-function(dataset)
{
resA<-sapply(dataset$ENSA,function(x) {strsplit(as.character(x),";")})
resB<-sapply(dataset$ENSB,function(x) {strsplit(as.character(x),";")})
res<-cbind(sapply(resA,function(x) x[1]),sapply(resB,function(x) x[1]))
return(res)
}

f_extract_pairedgenes_entrez<-function(dataset)
{
resA<-sapply(dataset$ENTREZA,function(x) {strsplit(as.character(x),";")})
resB<-sapply(dataset$ENTREZB,function(x) {strsplit(as.character(x),";")})
res<-cbind(sapply(resA,function(x) x[1]),sapply(resB,function(x) x[1]))
return(res)
}

f_clean_pairedgenes<-function(dataset,add0col=TRUE) #removes not unique pairs and add extra colums with 0
{
  data.cleaned2<-cbind(dataset[,1],dataset[,2],sapply(1:dim(dataset)[1],function(x) paste0(dataset[x,1],dataset[x,2])),sapply(1:dim(dataset)[1],function(x) paste0(dataset[x,2],dataset[x,1])))
  data.cleaned2<-data.cleaned2[order(data.cleaned2[,1],data.cleaned2[,2]),]
  data.cleaned3<-c(0,0)
  for (i in 1:(dim(data.cleaned2)[1]-1))
    {
    if (data.cleaned2[i,3]!=data.cleaned2[i+1,3] && data.cleaned2[i,4]!=data.cleaned2[i+1,4] ){data.cleaned3<-rbind(data.cleaned3,c(data.cleaned2[i,1],data.cleaned2[i,2]))}
    }
  data.cleaned3<-rbind(data.cleaned3,c(data.cleaned2[i+1,1],data.cleaned2[i+1,2]))
  data.cleaned3<-data.cleaned3[-1,]
#  data.cleaned3<-data.cleaned3[-dim(data.cleaned3)[1],]
  data.cleaned3<-na.omit(data.cleaned3)
  if (add0col)
  {
  data.cleaned3<-cbind(data.cleaned3,0)
  data.cleaned3<-cbind(data.cleaned3[,1],data.cleaned3[,2],data.cleaned3[,3])
  }
  else if (add0col==FALSE) {data.cleaned3<-cbind(data.cleaned3[,1],data.cleaned3[,2])}
  return(data.cleaned3)
}

f_clean_pairedgenes3col<-function(dataset) #removes not unique pairs and add extra colums with 0
{
  data.cleaned2<-cbind(dataset[,1],dataset[,2],dataset[,3],sapply(1:dim(dataset)[1],function(x) paste0(dataset[x,1],dataset[x,2])),sapply(1:dim(dataset)[1],function(x) paste0(dataset[x,2],dataset[x,1])))
  data.cleaned2<-data.cleaned2[order(data.cleaned2[,1],data.cleaned2[,2]),]
  data.cleaned3<-c(0,0,0)
  for (i in 1:(dim(data.cleaned2)[1]-1))
    {
    if (data.cleaned2[i,4]!=data.cleaned2[i+1,4] && data.cleaned2[i,5]!=data.cleaned2[i+1,5] ){data.cleaned3<-rbind(data.cleaned3,c(data.cleaned2[i,1],data.cleaned2[i,2],data.cleaned2[i,3]))}
    }
  data.cleaned3<-rbind(data.cleaned3,c(data.cleaned2[i+1,1],data.cleaned2[i+1,2],data.cleaned2[i+1,3]))
  data.cleaned3<-data.cleaned3[-1,]
#  data.cleaned3<-data.cleaned3[-dim(data.cleaned3)[1],]
  data.cleaned3<-na.omit(data.cleaned3)
  data.cleaned3<-cbind(data.cleaned3[,1],data.cleaned3[,2],data.cleaned3[,3])
  return(data.cleaned3)
}


generate_randompairs.ens<-function(npairs)
{
  ensgenes<-read.table("/mnt/expressions/STRING/ensP2G",header=FALSE)
  randensgenes<-unique(sample(unique(ensgenes$V1), 2*npairs, replace=FALSE))
  randensgenes<-data.frame(ens1=randensgenes[1:npairs],ens2=randensgenes[(npairs+1):(2*npairs)],int=rep(0,npairs))
  return(randensgenes)
}

shuffle_pairs<-function(data.pairs,add0vec=FALSE,nsamples=FALSE)
{
if (nsamples==FALSE) {nsamples<-length(data.pairs[,1])}
uniqgenes<-unique(c(data.pairs[,1],data.pairs[,2]))
res<-sample(uniqgenes,2,replace=FALSE)
for (i in 2:nsamples)
  {
  res<-rbind(res,sample(uniqgenes,2,replace=FALSE))
  }
if (add0vec) {res<-cbind(res,0)}
return(res)
}

generate_randompairs.entrez<-function(npairs)
{
  #cat /mnt/scratch/fabrizio/LDLD/entrezfromUCSC.tsv | grep -v chrY | grep -v chrX | grep -v chrM | grep -v chrUn > /mnt/scratch/fabrizio/LDLD/entrezfromUCSC_Auto.tsv
  ensgenes<-read.table("/mnt/scratch/fabrizio/LDLD/entrezfromUCSC_Auto.tsv",header=FALSE)
  ensgenes<-sapply(ensgenes[,4],function(x) return(strsplit(as.character(x),";")))
  ensgenes<-sapply(ensgenes,function(x) return(x[1]))
  ensgenes<-ensgenes[ensgenes!="n/a"]
  randensgenes<-unique(sample(unique(ensgenes), 2*npairs, replace=FALSE))
  randensgenes<-data.frame(ens1=randensgenes[1:npairs],ens2=randensgenes[(npairs+1):(2*npairs)],int=rep(0,npairs))
  return(randensgenes)
}

pairwise_interactions<-function(data.pairs,data.inter,npairs=FALSE)
{
  if (npairs==FALSE){npairs<-dim(data.pairs)[1];print(npairs)}
  for (j in 1:npairs)
  {
  if (j%%1000==0) {print(paste0("npairs is: ",j,"; n interactions so far is: ",sum(data.pairs[,3]!=0)))}
  temp1<-subset(data.inter,data.inter$V1==data.pairs[[j,1]])
  temp1<-subset(temp1,temp1$V2==data.pairs[[j,2]])$V3
  temp1 <- temp1[!is.na(temp1)]
  if (length(temp1)>0) {print(temp1);data.pairs[[j,3]]<-mean(temp1)} else
    {
    temp1<-subset(data.inter,data.inter$V2==data.pairs[[j,1]])
    temp1<-subset(temp1,temp1$V1==data.pairs[[j,2]])$V3
    temp1 <- temp1[!is.na(temp1)]
    if (length(temp1)>0) {print(temp1);data.pairs[[j,3]]<-mean(temp1)}
    }
  }
return(data.pairs)
}


pairwise_interactions_grep<-function(data.pairs,data.inter,npairs=FALSE) #NB: if I want to test the command on the shell I have to add the $ in front of the '
{
  if (npairs==FALSE){npairs<-dim(data.pairs)[1];print(npairs)}
  for (j in 1:npairs)
  {
    if (j%%500==0) {print(paste0("npairs is: ",j,"; n interactions so far is: ",sum(data.pairs[,3]!=0)))}
    comm<-paste0("grep -e '^",data.pairs[[j,1]],"\t",data.pairs[[j,2]],"\t\\|^",data.pairs[[j,2]],"\t",data.pairs[[j,1]],"\t' ",data.inter," | cut -f 3")
    temp1<-as.numeric(system(comm,intern=TRUE))
    temp1 <- temp1[!is.na(temp1)]
    if (length(temp1)>0) {print(temp1);data.pairs[[j,3]]<-mean(temp1)}
  }
return(data.pairs)
}
}
#======================================== CHECKS and POSTPROCESSING FUNCTIONS=========================================
{
postproc<-function(fromchr=1,tochr=21){
  options(scipen=999)
  counter<-0
  for (i in fromchr:tochr)
  {
  for (j in (i+1):22)
    {
    print(c(i,j))
    counter<-counter+1
    nlines<-system(paste0("cat above95/missense/anal/chr",i,".chr",j,".res | wc -l"),intern=TRUE)
    if (nlines>0){
	    datares<-read.table(paste0("above95/missense/anal/chr",i,".chr",j,".res"))
		    dataAann<-read.table(paste0("above95/missense/anal/chr",i,".chr",j,".chrA.ann"))
		    dataBann<-read.table(paste0("above95/missense/anal/chr",i,".chr",j,".chrB.ann"))
		    names(datares)<-c("chrA","chrB","snpA","snpB","dbsnpA","dbsnpB","nA","nB","Nse","nAB","D","D1","r2","bo","bo2","bo3","T2","Xtot","X","pop","popX")
		    dataAann$V1<-as.numeric(gsub("chr","",dataAann$V1))
		    dataBann$V1<-as.numeric(gsub("chr","",dataBann$V1))
		    names(dataAann)<-c("chrA","snpAm","snpA","ENSA","HGNCA","entrezA","mapA","dbsnpA")
		    names(dataBann)<-c("chrB","snpBm","snpB","ENSB","HGNCB","entrezB","mapB","dbsnpB")
		    datares[,"ENSA"]<-"0";datares[,"HGNCA"]<-"0";datares[,"entrezA"]<-"0";datares[,"mapA"]<-0;datares[,"dbsnpA"]<-0;
	    datares[,"ENSB"]<-"0";datares[,"HGNCB"]<-"0";datares[,"entrezB"]<-"0";datares[,"mapB"]<-0;datares[,"dbsnpB"]<-0;
	    if (dim(datares)[1]!=dim(dataAann)[1])
	    	{
		    for (iA in 1:dim(datares)[1])
		    {
			    for (jA in 1:dim(dataAann)[1])
			    {
				    if (datares$chrA[iA]==dataAann$chrA[jA] && datares$snpA[iA]==dataAann$snpA[jA])
				    {
					    datares$ENSA[iA]<-as.character(dataAann$ENSA[jA])
						    datares$HGNCA[iA]<-as.character(dataAann$HGNCA[jA])
						    datares$entrezA[iA]<-as.character(dataAann$entrezA[jA])
						    datares$mapA[iA]<-dataAann$mapA[jA]
						    datares$dbsnpA[iA]<-dataAann$dbsnpA[jA]
						    break
				    }
			    }
		    }
	    }
	    else
	    {
		    datares$ENSA<-as.character(dataAann$ENSA);datares$HGNCA<-as.character(dataAann$HGNCA)
			    datares$entrezA<-as.character(dataAann$entrezA);datares$mapA<-dataAann$mapA;datares$dbsnpA<-dataAann$dbsnpA
	    }
	    if (dim(datares)[1]!=dim(dataBann)[1])
	    {
		    for (iB in 1:dim(datares)[1])
		    {
			    for (jB in 1:dim(dataBann)[1])
			    {
				    if (datares$chrB[i]==dataBann$chrB[j] && datares$snpB[i]==dataBann$snpB[j])
				    {
					    datares$ENSB[iB]<-as.character(dataBann$ENSB[jB])
						    datares$HGNCB[iB]<-as.character(dataBann$HGNCB[jB])
						    datares$entrezB[iB]<-as.character(dataBann$entrezB[jB])
						    datares$mapB[iB]<-dataBann$mapB[jB]
						    datares$dbsnpB[iB]<-dataBann$dbsnpB[jB]
						    break
				    }
			    }
		    }
	    }
	    else
	    {
		    datares$ENSB<-as.character(dataBann$ENSB);datares$HGNCB<-as.character(dataBann$HGNCB)
			    datares$entrezB<-as.character(dataBann$entrezB);datares$mapB<-dataBann$mapB;datares$dbsnpB<-dataBann$dbsnpB
	    }
	    write.table(datares,file=paste0("above95/missense/anal/chr",i,".chr",j,".ann2"),row.names = FALSE,quote=FALSE);
    }}
    }
  return(datares)
  }

#return frequencies from raw data to check ICLDfreq.c
f_freqs<-function(res,npop){ #takes a raw vector of genotypes e.g. res<-c(0,0,0,0,0,0,2,0,1,0) with 1 AA and 2 Aa and vector of pop sizes e.g. nspops<-c(106,107,102,108,103,102,100,99,96,99,97,98)
	  ris<-rep(0,npop)
		  for (i in 1:npop)
		  {
			  mypop<-i-1
				  if (i>1)
				  {
					  nAa<-sum(res[(sum(nspops[1:mypop])+1):sum(nspops[1:(mypop+1)])]==2)
										nAA<-sum(res[(sum(nspops[1:mypop])+1):sum(nspops[1:(mypop+1)])]==1)
														      ris[i]<-(nAa+2*nAA)/(2*nspops[i])
				  }
				  else if (i==1)
				  {
					  nAa<-sum(res[(1:sum(nspops[1:(mypop+1)]))]==2)
						  nAA<-sum(res[(1:sum(nspops[1:(mypop+1)]))]==1)
						  ris[i]<-(nAa+2*nAA)/(2*nspops[i])
				  }
		  }
	  return(ris)
  }


}  
#===================================ARTIFACT EVIDENCE ANALYSES: removal along network, imbalanced reads, CG ==========
{
linkinCG<-function(resfile,CGbed)
{
strA<-paste0(resfile$chrA,"-",resfile$posA)
strB<-paste0(resfile$chrB,"-",resfile$posB)
strCG<-paste0(signCG$chr,"-",signCG$end)
res<-rep(0,dim(resfile)[1])
for (i in 1:dim(resfile)[1])
        {
        A1<-sum(as.numeric(strA[i]==strCG))
        B1<-sum(as.numeric(strB[i]==strCG))
        if (A1==1 && B1==1){res[i]<-1;print(i)}
        }
return(res)
}

createbedlogo<-function(data,namelogo) 
{ 
system(paste0("rm ",namelogo)) 
for ( i in 1:dim(data)[1]) 
  { 
  system(paste0("twoBitToFa /mnt/scratch/fabrizio/LDLD/hg19.2bit test.fa -seq=",data$chr[i]," -start=",data$start[i]-5," -end=",data$end[i]+5))
  nhits<-system("cat test.fa | wc -l",intern=TRUE) 
  if (nhits > 1){  system(paste0("cat test.fa | grep -v chr | head -1 >> ",namelogo))} else {system(paste0("echo NaN >> ",namelogo))} 
  if (i%%50==0){print(i)} 
  } 
} 

bamlogo<-function(mylogo,tilltheend=TRUE,till=100,startfrom=1,filetemp="temp") #mylogo e.g. "logo.bed"
{
#filetemp useful when running in parallel. Otherwise, if same filetemp it might crash or give wrong outputs
#notice that there might be Segmentation fault warnings if for some individuals BamTable2 does not find data. This is not a problem since there are others to compute statistics
datalogo<-read.table(mylogo)
head(datalogo)
if (tilltheend){upto<-dim(datalogo)[1]} else {upto<-till}
names(datalogo)<-c("chr","start","end","ref","alt","quality100g","seq5bp")
for (j in startfrom:upto)
  {
  mychr<-as.numeric(datalogo$chr[j])
  pos<-as.numeric(datalogo$end[j])
  system(paste0("for i in /mnt/454/HGDP/genomes_Bteam/HG19Align/*.bam /mnt/sequencedb/11men/martin/hg19_1000g/BAM/*.bam ; do ./BamTable2 -r ",mychr,":",pos-1,"-",pos," $i | grep ",pos, "; done > ", filetemp))  #| awk -v OFS='\t' '{res=0;for (i = 3; i <= 6; i++){if ($i>0){res=res+1}};if(res>1){print}}'; done > temp")
  nltemp<-system(paste0("cat ",filetemp," | wc -l"),intern=TRUE)
  if (nltemp>0)
    {	
    data<-read.table(filetemp)
    names(data)<-c("chr","end","V1","V2","V3","V4")
    data<-subset(data,data$end==pos)
    data<-subset(data,data$chr==mychr)
    if (dim(data)[1]==0) {print ("NAN");res<-c(0,0,0,0,0,0)}  else
      {
      data1<-subset(data,as.numeric(data$V1>0)+as.numeric(data$V2>0)+as.numeric(data$V3>0)+as.numeric(data$V4>0)>1)
      res<-c(mean(data1$V1),mean(data1$V2),mean(data1$V3),mean(data1$V4),dim(data1)[1],dim(data)[1])
      }
    }
  else
    {
    res<-c(0,0,0,0,0,0)
    }
  if (j==startfrom) {res0<-c(datalogo$chr[j],datalogo$start[j],as.character(datalogo$end[j]),as.character(datalogo$ref[j]),as.character(datalogo$alt[j]),as.character(datalogo$seq5bp[j]),res)} else
  {res0<-rbind(res0,c(datalogo$chr[j],datalogo$start[j],as.character(datalogo$end[j]),as.character(datalogo$ref[j]),as.character(datalogo$alt[j]),as.character(datalogo$seq5bp[j]),res))}
  print(c(j,dim(data)[1]))
  }
  return(res0)
}

prHvsE<-function(data,j)
#null hypothesis: heterozygous: reject it when the prob of achieving results as extreme or more extreme is unlikely. In my case 
  {
  if ((data$A[j]+data$C[j]+data$G[j]+data$T[j])==0)
  {
    res<-rep(0,4)
  }
  else
    {
    ta<-0
    tr<-0
    if (data$ref[j]=="A") {tr<-7} else if (data$ref[j]=="C"){tr<-8} else if (data$ref[j]=="G"){tr<-9} else if (data$ref[j]=="T"){tr<-10}
    if (data$alt[j]=="A") {ta<-7} else if (data$alt[j]=="C"){ta<-8} else if (data$alt[j]=="G"){ta<-9} else if (data$alt[j]=="T"){ta<-10}
    #to handle strange triallelic or indels------v
    if (nchar(data$ref[j])>1 || nchar(data$alt[j])>1 ){biallelic<-0}else{biallelic<-1}
    if (ta==0) 
      { 
      ta<-substr(data$alt[j],1,1); 
      if (ta=="A") {ta<-7} else if (ta=="C"){ta<-8} else if (ta=="G"){ta<-9} else if (ta=="T"){ta<-10}
      }
    if (tr==0) 
      { 
      tr<-substr(data$alt[j],1,1); 
      if (tr=="A") {tr<-7} else if (tr=="C"){tr<-8} else if (tr=="G"){tr<-9} else if (tr=="T"){tr<-10}
      }
    #--------------------------------------------^
    othercall<-0
    for (inonra in 7:10)
      {
      if (inonra!=ta && inonra!=tr && data[j,tr]>0 ){othercall<-othercall+data[j,tr]}
      }
    mincall<-min(data[j,tr],data[j,ta])*data$tot_second_call[j]
    maxcall<-max(data[j,tr],data[j,ta])*data$tot_second_call[j]
    othercall<-othercall*data$tot_second_call[j]
    totcalls<-mincall+maxcall
    totcallsoth<-othercall+maxcall
    imincall<-mincall
    iothercall<-othercall
    probAltEstimate<-binom.test(mincall,totcalls,p=0.5)$estimate
    if (mincall==0) {probAltNoH<-2} else { probAltNoH<-binom.test(mincall,totcalls,p=0.5)$p.value}
    if (othercall==0) {probOthNoH<-2} else { probOthNoH<-binom.test(othercall,totcallsoth,p=0.5)$p.value}
    res<-c(probAltEstimate,probAltNoH,probOthNoH,biallelic)
    }
    names(res)<-c("imbalance","pvalue_no05_alt","pvalue_no05_other","biallelic")
  return(res)
  }
#probAltH prob that by chance less extreme than this = 0.1487818

ave_LDLD<-function(chr,pos,dataLDLD)
{
dataLDLD_tA<-subset(dataLDLD,dataLDLD$chrA==chr)
dataLDLD_tA<-subset(dataLDLD_tA,dataLDLD_tA$posA==pos)
dataLDLD_tB<-subset(dataLDLD,dataLDLD$chrB==chr)
dataLDLD_tB<-subset(dataLDLD_tB,dataLDLD_tB$posB==pos)
LDLDfdr<-c(dataLDLD_tA$combined_fdr,dataLDLD_tB$combined_fdr)
return(c(length(LDLDfdr),mean(LDLDfdr),min(LDLDfdr)))
}

ave_LDLD_BED<-function(dataLDLD,mylogo,myfrom=1,myto=1)
{
  if (myto==1){myl<-dim(mylogo)[1]} else {myl<-myto}
  myf<-myfrom
  res<-cbind(rep(0,myl-myf+1),rep(0,myl-myf+1),rep(0,myl-myf+1))
  #names(res)<-c("nlinks","ave_fdr","min_fdr")
  for (i in myf:myl)
    {
    res[i-myf+1,]<-ave_LDLD(mylogo$chr[i],mylogo$end[i],dataLDLD)
    }
return(res)
}

logo_below_thr<-function(logo,max_pvalue_no05_alt)
  {
  logoinfo<-subset(logo, logo$A!=0 | logo$C!=0 | logo$T!=0 | logo$G!=0)
  logo_potentialH<-subset(logoinfo, logoinfo$tot_second_call>=0)
  logo_potentialH<-subset(logo_potentialH, logo_potentialH$imbalance>0) #informative snps
  res<-logo_potentialH[logo_potentialH$pvalue_no05_alt<max_pvalue_no05_alt,]
  print(dim(res)[1]/dim(logo_potentialH)[1])
  return(list(res,logo_potentialH))
  }

retrieve_partners<-function(chr,pos,dataLDLD)
{
  dataLDLD_tA<-subset(dataLDLD,dataLDLD$chrA==chr)
  dataLDLD_tA<-subset(dataLDLD_tA,dataLDLD_tA$posA==pos)
  dataLDLD_tB<-subset(dataLDLD,dataLDLD$chrB==chr)
  dataLDLD_tB<-subset(dataLDLD_tB,dataLDLD_tB$posB==pos)
  res<-as.data.frame(cbind(c(dataLDLD_tA$chrB,dataLDLD_tB$chrA),c(dataLDLD_tA$posB,dataLDLD_tB$posA)))
  names(res)<-c("chr","pos")
  return(res)
}

retrieve_partners_collapse<-function(mysnp,dataLDLD,generatecollapse=FALSE)
{
if (generatecollapse) 
  {
  dataLDLDA<-as.character(paste0(dataLDLD$chrA,".",dataLDLD$posA))
  dataLDLDB<-as.character(paste0(dataLDLD$chrB,".",dataLDLD$posB))
  dataLDLD<-cbind(dataLDLDA,dataLDLDB)
  }
  id1<-which(dataLDLD[,1]==mysnp)
  id2<-which(dataLDLD[,2]==mysnp)
  dataLDLD1<-dataLDLD[,2][dataLDLD[,1]==mysnp]
  dataLDLD2<-dataLDLD[,1][dataLDLD[,2]==mysnp]
  return(cbind(c(dataLDLD1,dataLDLD2),c(id1,id2)))
}

run_through_links<-function(chr,pos,dataLDLD)
  {
  mysnp<-paste0(chr,".",pos)
  dataLDLDleft<-c()
  dataLDLDlinked<-c()
  dataLDLDA<-as.character(paste0(dataLDLD$chrA,".",dataLDLD$posA))
  dataLDLDB<-as.character(paste0(dataLDLD$chrB,".",dataLDLD$posB))
  dataLDLDcollapsed<-cbind(dataLDLDA,dataLDLDB)
  myres<-retrieve_partners_collapse(mysnp,dataLDLDcollapsed,generatecollapse=FALSE)
  myres1<-myres[,1]
  myres2<-myres[,2]
  dataLDLDlinked<-append(dataLDLDlinked,myres1)
  dataLDLDcollapsed<-dataLDLDcollapsed[-as.numeric(myres2),]
  counter<-1
  while (dim(dataLDLDcollapsed)[1]>0 && length(myres1)>0)
    {
    myrest<-lapply(myres1, function(x) retrieve_partners_collapse(x,dataLDLDcollapsed,generatecollapse=FALSE))
    myres2<-unique(unlist(lapply(myrest,function(x) x[,2])))
    myres1<-unique(unlist(lapply(myrest,function(x) x[,1])))
    if (length(myres1)>0) 
      {
      dataLDLDlinked<-append(dataLDLDlinked,myres1)
      dataLDLDcollapsed<-dataLDLDcollapsed[-as.numeric(myres2),]
      }
    counter<-counter+1
    print(c(counter,dim(dataLDLDcollapsed)[1],length(dataLDLDlinked)))
    }
  return(list(dataLDLDlinked,dataLDLDcollapsed))
  }

shuffle_links<-function(nlinks,nnodes)
{
LDLD<-cbind(sample(1:nnodes,nnodes),sample(1:nnodes,nnodes))
sapply(1:nlinks,function(x) {if (LDLD[x,1]==LDLD[x,2]) {LDLD[x,2]<-sample((1:nnodes)[-LDLD[x,2]],1)}})
return(LDLD)
}

distr_k<-function()
distr_cluster_sizes<-function()

shuffle_links(10,100)
}
#---------per sample ARTIFACT_ANALYSES--------------------------------------------------------------------------------
{
#           0/0 0/0 0/1 0/0 1/1 0/1 0/1 1/1 1/1     
#           0/0 0/1 0/0 1/1 0/0 0/1 1/1 0/1 1/1    
#n11     0    0    0    0    0    0.5 1     1    2   
#n00     2    1    1    0    0    0.5 0     0    0
#n+      2     1    1   0     0    1    1    1     2
#n10     0    0    1    0    2    0.5 0     1    0
#n01     0    1    0    2    0    0.5 1     0    0
#in terms of errors n11 assumes that rare variants created as errors
#2-n11 assumes that rare variants are not seen in low quality

nABf<-function(gt) #compute nAB only for minor allele
{
	if (gt==9){res<-2} else 
	if (gt==6 || gt==8){res<-1} else 
	if (gt==5){res<-1/2} else
	{res<-0}
	return(res)
}

nAandBf<-function(gt)
{
	if (gt==7 && gt==8 && gt==9 ){res1<-2} else 
	if (gt==4 || gt==5 || gt==6){res1<-1} else {res1<-0}
	if (gt==3 && gt==6 && gt==9 ){res2<-2} else 
	if (gt==2 || gt==5 || gt==8){res2<-1} else {res2<-0}
	res<-(res1+res2)/2
	return(res)
}

nAf<-function(gt)
{
	if (gt==4 || gt==5 || gt==6){myres<-1} 
	else if (gt==7 && gt==8 && gt==9 ){myres<-2} else {myres<-0}
	return(myres)
}

nBf<-function(gt)
{
	if (gt==2 || gt==5 || gt==8){myres<-1} 
	else if (gt==3 && gt==6 && gt==9 ){myres<-2} else {myres<-0}
	return(myres)
}

pairsANDlogfiles2logpop<-function(myfolder,mysign=fullsign,ntot_pops=12,ntot_chr=22,summarize=TRUE) #with summarize=TRUE it gives only nAB per pop, if FALSE it gives full log files. in fullsign I give the snps to take
{
  logpop<-list()
  if (summarize) {for (i in 0:(ntot_pops-1)){logpop[[i+1]]<-list(0,0)}} else {for (i in 0:(ntot_pops-1)){logpop[[i+1]]<-list(0)}}
  for (i in 0:(ntot_pops-1))
  {
    for (chrA in 1:(ntot_chr-1))
      {
      for (chrB in (chrA+1):ntot_chr)
	{
	print(c(i,chrA,chrB))
	system(paste0("cat ",myfolder,"/chr",chrA,".tabchr",chrB,".tab.log | awk '{if ($NF==",i," && NF > 16){print}}'> ",myfolder,"/temp.log"))
	mydata<-read.table(paste0(myfolder,"/temp.log"),header=FALSE)
	tempsnp<-subset(mysign,mysign[,1]==chrA)
	tempsnp<-subset(tempsnp,tempsnp[,2]==chrB)
	mysignsnps<-paste0(tempsnp[,3],".",tempsnp[,4])
	mylogsnps<-paste0(mydata[,1],".",mydata[,2])
	mydata<-mydata[!is.na(match(mylogsnps,mysignsnps)),]
	if (summarize)
	  {
	  logpop[[i+1]][[1]]<-logpop[[i+1]][[1]]+sapply(3:(dim(mydata)[2]-1), function(x) sum(sapply(mydata[,x],nABf)))
	  logpop[[i+1]][[2]]<-logpop[[i+1]][[2]]+sapply(3:(dim(mydata)[2]-1), function(x) sum(sapply(mydata[,x],nAandBf)))
	  } else
	  {
	  if (logpop[[i+1]]==0){logpop[[i+1]]<-cbind(chrA,chrB,mydata)} else {logpop[[i+1]]<-rbind(logpop[[i+1]],cbind(chrA,chrB,mydata))}
	  }
	}
      }
    }
  return(logpop)
}

stringsnp2bed<-function(stringsnp){
if ( is.matrix(stringsnp)){
    for (i in 1:dim(stringsnp)[1])
    {
        temp<-as.numeric(unlist(strsplit(stringsnp[i,1],split="\\.")))
        if (i==1) {myres<-c(temp[1],temp[2]-1,temp[2])} else
        {
        myres<-rbind(myres,c(temp[1],temp[2]-1,temp[2]))
        }
    }
    return(myres)
    } else 
    {
        stringsnp.df<-as.data.frame(matrix(unlist(unname(sapply(stringsnp,function(x) strsplit(x,split='[.]')))),ncol=2,byrow=TRUE),stringsAsFactors=FALSE)
        stringsnp.df<-data.frame(chr=stringsnp.df[,1],init=as.numeric(stringsnp.df[,2])-1,end=stringsnp.df[,2])
        return(stringsnp.df)
    }
}

#multinomial goodness of fit as in http://www.r-tutor.com/elementary-statistics/goodness-fit/multinomial-goodness-fit
nAB2uniform<-function(nAB_l, pthr)
{
    keptindices<-1:length(nAB_l)
    res<-nAB_l
    while ( chisq.test(res, p=rep(1/length(res),length(res)))$p.value < pthr )
      {
#      pbinom<-sapply(res, function(x) binom.test(round(x), round(sum(res)), p = 1/length(res))$p.value)
      #removed_index<-order(pbinom)[length(pbinom)]
      removed_index<-order(res)[length(res)]
      keptindices<-keptindices[-removed_index]
      res<-res[-removed_index]
      print(c(removed_index,length(res)))
      }
return(list(keptindices,res))
}


partition_samples_nAB_samesigma<-function(mydata,maxk=25)
{
  mod1l<-sum(dnorm(mydata,mean=mean(mydata),sd=sd(mydata),log=TRUE))
  mod<-lapply(2:maxk,function(x) normalmixEM(mydata,lambda=1/x,mu=breakshist(x,mydata),sigma=10))
  #mod<-lapply(2:maxk,function(x) normalmixEM(mydata,lambda=1/x,mu=breakshist(x,mydata),sigma=rep(10,x)))
  myl<-c(mod1l,sapply(1:(maxk-1),function(i) mod[[i]]$loglik)) #loglik calculate like this (example for 2 distributions sum(log(mod2$lambda[1]*dnorm(mydata,mean=mod2$mu[1],sd=mod2$sigma[1])+mod2$lambda[2]*dnorm(mydata,mean=mod2$mu[2],sd=mod2$sigma[2]))) #this seems right formula! )
  myAICc<-sapply(1:maxk,function(i) AICc(myl[i],i,length(mydata)))#2*i,length(mydata)))
  myk<-which(myAICc==min(myAICc))[1]
  print(myk)
  mymod<-mod[[myk-1]]
  mypartitions<-apply(mymod$posterior,1,function(x) which(x==max(x)))
  return(list(mymod,mypartitions,myk,myl,myAICc))
}

partition_samples_nAB<-function(mydata,maxk=6,phtr=0.01,reorder=FALSE)
{
  #---I need ordered data, otherwise add lines below-----v
  if (reorder) { 
  #    nAB_l<-unlist(logpop[[i]][1])
  #    ord<-order(nAB_l[1:length(nAB_l)])
  #    mydata<-nAB_l[ord]
    mydata_temp<-unlist(unname(mydata))
    mydata<-mydata[order(mydata_temp)]
    }
  mod1l<-sum(dnorm(mydata,mean=mean(mydata),sd=sd(mydata),log=TRUE))
# mod<-lapply(2:maxk,function(x) normalmixEM(mydata,lambda=1/x,mu=breakshist(x,mydata),sigma=rep(10,x)))#sigma=10))
  mod<-list()
  for (x in 2:maxk)
    {
    #mod[[x-1]]<-normalmixEM(mydata,lambda=1/x,mu=breakshist(x,mydata),sigma=rep(10,x))
    mod[[x-1]]<-normalmixEM(mydata,lambda=1/x,mu=breakshist(x,mydata),arbvar=TRUE)
    }
  #mod<-lapply(2:maxk,function(x) normalmixEM(mydata,lambda=1/x,mu=breakshist(x,mydata),sigma=rep(10,x)))
  myl<-c(mod1l,sapply(1:length(mod),function(x) mod[[x]]$loglik)) #loglik calculate like this (example for 2 distributions sum(log(mod2$lambda[1]*dnorm(mydata,mean=mod2$mu[1],sd=mod2$sigma[1])+mod2$lambda[2]*dnorm(mydata,mean=mod2$mu[2],sd=mod2$sigma[2]))) #this seems right formula! )
  myAICc<-sapply(1:maxk,function(x) AICc(myl[x],2*x,length(mydata)))
  #likelihooRatioTest(myl[1],myl[i],2*i-2)
  #myLikelihoodRatioTestvsprevious<-sapply(2:maxk,function(x) likelihooRatioTest(myl[x-1],myl[x],2))
  #myLikelihoodRatioTestvs1<-sapply(2:maxk,function(x) likelihooRatioTest(myl[1],myl[x],2*x-2))  
  maxnotfound<-TRUE
  mykLR<-length(mod)+1 #maxk
  while (maxnotfound && mykLR>1)
  {
    myLikelihoodRatioTestvsprevious<-sapply(1:(mykLR),function(x) likelihooRatioTest(myl[x],myl[mykLR],(mykLR-x)*2))
    if (sum(as.numeric(myLikelihoodRatioTestvsprevious<phtr))==(mykLR-1)) {maxnotfound<-FALSE} else {mykLR<-mykLR-1}
  }
  mykAICc<-which(myAICc==min(myAICc))[1]
  mykAICc<-min(order(myAICc)[1:sum(sapply(1:length(myAICc),function(x) exp(0.5*(min(myAICc)-myAICc[x])))[order(myAICc)]>phtr)]) #rather than taking best model, I take most parsimonious with relative likelihood 0.05
  print(c(mykAICc,mykLR))
#  myk<-2
  if (mykAICc>1) {mymodAICc<-mod[[mykAICc-1]]; mypartitionsAICs<-apply(mymodAICc$posterior,1,function(x) which(x==max(x)))} else {mypartitionsAICs<-rep(1,length(mydata));mymodAICc<-1}
  if (mykLR>1) {mymod<-mod[[mykLR-1]]; mypartitionsLR<-apply(mymod$posterior,1,function(x) which(x==max(x)))} else {mypartitionsLR<-rep(1,length(mydata));mymod<-1}
  return(list(k_from_LR=mykLR,partitions_from_LR=mypartitionsLR,model_LR=mymod,k_from_AICc=mykAICc,partitions_from_AICc=mypartitionsAICs,model_fromAICc=mymodAICc,LR_pvalue=likelihooRatioTest(myl[1],myl[mykLR],(mykLR-1)*2), AICc_relativelik=exp((min(myAICc)-myAICc[1])/2)))
#  plot(nAB_l[1:length(nAB_l)][ord],pch=mypch[mypartitionsLR],col=mycols[mypartitionsAICs])
  }

#exp(1)
AICc<-function(loglik,k,n)
{
2*k-2*loglik+(2*k*(k+1))/(n-k-1)
}

likelihooRatioTest<-function(lk_null,lk_alt,df)
{
  1-pchisq(2*(lk_alt-lk_null),df=df)
}

breakshist<-function(x,mydata) sapply(0:(x+2-1),function(y) min(mydata)+y*(max(mydata)-min(mydata))/(x+2))[-x][-1]

if (FALSE) #checks to verify how it calculates the likelihood in mixtools
{
  mean(mydata)
  sd(mydata)
  sum(log(dnorm(mydata,mean=mean(mydata),sd=sd(mydata),log=FALSE)))
  sum(dnorm(mydata,mean=mean(mydata),sd=sd(mydata),log=TRUE))
  #sum(pnorm(mydata,mean=mean(mydata),sd=sd(mydata),log=TRUE))
  mod2$loglik
  mod3$loglik
  mod4$loglik
  sum(log(mod2$posterior[,1]))+sum(log(mod2$posterior[,2]))
  mod2$posterior[,1]*dnorm(mydata,mean=mod2$mu,sd=mydata$sigma)
  round(0.45)
  round(0.55)
  sum(log(mod2$posterior[,1]*dnorm(mydata,mean=mod2$mu[1],sd=mod2$sigma[1])))
  sum(log(mod2$posterior[,1]*dnorm(mydata,mean=mod2$mu[1],sd=mod2$sigma[1])+mod2$posterior[,2]*dnorm(mydata,mean=mod2$mu[2],sd=mod2$sigma[2])))
  sum(log(mod2$lambda[1]*mod2$posterior[,1]*dnorm(mydata,mean=mod2$mu[1],sd=mod2$sigma[1])+mod2$lambda[2]*mod2$posterior[,2]*dnorm(mydata,mean=mod2$mu[2],sd=mod2$sigma[2])))
  sum(log(mod2$lambda[1]*dnorm(mydata,mean=mod2$mu[1],sd=mod2$sigma[1])+mod2$lambda[2]*dnorm(mydata,mean=mod2$mu[2],sd=mod2$sigma[2]))) #this seems right formula!
  mod2$loglik
  mydata$sigma
  mydata
}

infoplot_ord<-function(i,mymaxk=6,mypthr=0.01,mynABdefault=TRUE,mynAB=0,showreal_confidence_interval=FALSE,reshufflings=0,mylog=FALSE,showmutperssample=FALSE)
{
  if (mynABdefault==1)
  {
  nAB_l<-unlist(logpop[[i]][1])
  } else {nAB_l<-mynAB}
  nAandB_l<-unlist(logpop[[i]][2])
  mres<-rep(10,length(nAB_l))
  sres<-rep(0.1,length(nAB_l))
  if (showreal_confidence_interval && reshufflings == 0)
  {
  permbychr<-permute_chrAB(mylog,10000)
  for (i in 1:10000) {permbychr[,i]<-permbychr[,i][order(permbychr[,i])]}
  } else if (reshufflings==0) {
    res<-rmultinom(1000,size=round(sum(nAB_l)),prob=rep(1/length(nAB_l),length(nAB_l)))
    mres<-sapply(1:length(nAB_l),function(x) mean(res[x,]))
    sres<-sapply(1:length(nAB_l),function(x) 1.96*sd(res[x,]))
    }
  ord<-order(nAB_l[1:length(nAB_l)])
  return(ord)
  }

infoplot<-function(i,mymaxk=6,mypthr=0.01,mynABdefault=TRUE,mynAB=0,showreal_confidence_interval=TRUE,reshufflings=0,mylog=FALSE,showmutperssample=FALSE,mypopname="",mutperssample=mymuts)
{
  if (mynABdefault)
  {
  nAB_l<-unlist(logpop[[i]][1])
  } else {nAB_l<-mynAB}
  nAandB_l<-unlist(logpop[[i]][2])
  mres<-rep(10,length(nAB_l))
  sres<-rep(0.1,length(nAB_l))
  if (showreal_confidence_interval && reshufflings==0)
  {
   reshufflings<-permbychr[[i]]
  #  permbychr<-permute_chrAB(mylog,10000)
# for (j in 1:10000) {permbychr[[i]][,j]<-permbychr[[i]][,j][order(permbychr[[i]][,j])]}
  } else if ( reshufflings==0 ) {
    res<-rmultinom(1000,size=round(sum(nAB_l)),prob=rep(1/length(nAB_l),length(nAB_l)))
    mres<-sapply(1:length(nAB_l),function(x) mean(res[x,]))
    sres<-sapply(1:length(nAB_l),function(x) 1.96*sd(res[x,]))
    } else if ( reshufflings!=0 ) {reshufflings<-reshufflings}
  upper<-mres+1.96*sres
  lower<-mres-1.96*sres
  ord<-order(nAB_l[1:length(nAB_l)])
  if (showmutperssample){mutpoints<-mutperssample[popsamplesi[[i]]][ord]}
  myinfosamples<-as.data.frame(infosamples[i])
  seqcenter<-myinfosamples$Main.project.LC.Centers[ord]
  seqcenter<-as.character(seqcenter)
  seqcenter_col<-as.character(seqcenter)
  seqcenter_col[seqcenter=='BI']<-"firebrick"
  seqcenter_col[seqcenter=='BI,MPIMG']<-"firebrick"
  #seqcenter_col[seqcenter=='BI,MPIMG']<-"firebrick1"
  seqcenter_col[seqcenter=='MPIMG']<-"coral2"
  seqcenter_col[seqcenter=='BGI']<-"cadetblue"
  seqcenter_col[seqcenter=='BCM,BGI']<-"cadetblue"
  #seqcenter_col[seqcenter=='BCM,BGI']<-"cadetblue2"
  seqcenter_col[seqcenter=='BCM']<-"aquamarine2"
  seqcenter_col[seqcenter=='ILLUMINA']<-"chartreuse3"
  seqcenter_col[seqcenter=='SC']<-"brown1"
  seqcenter_col[seqcenter=='WUGSC']<-"darkgoldenrod1"
  etplat<-as.character(myinfosamples$ET.Pilot.Platforms[ord])
  etplat[is.na(etplat)]<-0
  etplat[etplat==""]<-0
  etplat_col<-etplat
  etplat_col[etplat=="ILLUMINA"]<-"deepskyblue2"
  etplat_col[etplat=="LS454"]<-"coral2"
  #etplat[etplat=="ILLUMINA,LS454"]<-"orchid3"
  etplat_col[etplat=="ILLUMINA,LS454"]<-"deepskyblue2"
  etplat_col[etplat=="0"]<-"floralwhite"
  #seqcenter<-as.character(seqcenter)
  #length(myinfosamples$Has.Exome.LOF.Genotypes)
  myinfosamples$Has.Exome.LOF.Genotypes[is.na(myinfosamples$Has.Exome.LOF.Genotypes)]<-0
  myinfosamples$Has.Exome.LOF.Genotypes[myinfosamples$Has.Exome.LOF.Genotypes==1]<-"cornsilk4"
  myinfosamples$Has.Exome.LOF.Genotypes[myinfosamples$Has.Exome.LOF.Genotypes==0]<-"floralwhite"
  myinfosamples$Has.Axiom.Genotypes[is.na(myinfosamples$Has.Axiom.Genotypes)]<-0
  myinfosamples$Has.Axiom.Genotypes[myinfosamples$Has.Axiom.Genotypes==1]<-"cornsilk4"
  myinfosamples$Has.Axiom.Genotypes[myinfosamples$Has.Axiom.Genotypes==0]<-"floralwhite"
  myinfosamples$Has.Affy.6.0.Genotypes[is.na(myinfosamples$Has.Affy.6.0.Genotypes)]<-0
  myinfosamples$Has.Affy.6.0.Genotypes[myinfosamples$Has.Affy.6.0.Genotypes==1]<-"cornsilk4"
  myinfosamples$Has.Affy.6.0.Genotypes[myinfosamples$Has.Affy.6.0.Genotypes==0]<-"floralwhite"
  h_nABplot<-1.03
  h_covplot<-1.75
  margin<-0.03
  plot(x=1:length(nAB_l),y=mres,ylim=c(0,2.53),type = "n",xlab=paste0(mypopname," samples"),ylab="nAB                                     info                                   ",yaxs="i",yaxt="n")
  axis(2, at=c(0,0.2,0.4,0.6,0.8,0.99,h_nABplot+margin,(h_nABplot+margin+h_covplot)/2,h_covplot,(h_covplot+margin+2.03)/2,(2.03+2.28)/2,(2.28+2.53)/2),labels=c(0,0.2,0.4,0.6,0.8,1,0,0.5,1,"SC","GT","ETP"),las=2)
  #--frame-of-plot------------------------v
  abline(0,0)
  abline(h_nABplot,0)
  abline(h_covplot,0)
  #abline(1.53,0)
  abline(2.03,0)
  abline(2.28,0)
  #---------------------------------------^
  #--PLOTPOINTS--------------------------vv
  #--filtering-threshold-assuming-uniform-v
  if (FALSE)
  {
    myrese10<-nAB2uniform( nAB_l, 10^(-10))
    myrese5<-nAB2uniform( nAB_l, 10^(-50))
    myrese2<-nAB2uniform( nAB_l, 10^(-2))
    segments(length(nAB_l[myrese10[[1]]]),0,length(nAB_l[myrese10[[1]]]),1,col="red",lty=2)
    segments(length(nAB_l[myrese5[[1]]]),0,length(nAB_l[myrese5[[1]]]),1,col="red",lty=2)
    segments(length(nAB_l[myrese2[[1]]]),0,length(nAB_l[myrese2[[1]]]),1,col="red",lty=2)
  }
  #---------------------------------------^
  #-------plotting null distr
  if (reshufflings!=0)
  {
  cdat <- as.list(as.data.frame(t(as.matrix(apply(reshufflings,MARGIN=2,sort)/max(nAB_l)))))
  names(cdat)[1] <- "x"  # vioplot() needs the first element to be called 'x'
  do.call(vioplot,c(cdat,list(col="cadetblue",colMed="cadetblue",add=TRUE)))
  }
  #-------plotting points
  #  points(nAB_l[1:length(nAB_l)][ord]/max(nAB_l),pch=19)
  res<-partition_samples_nAB(nAB_l[ord], maxk=mymaxk, phtr=mypthr)
  mycols=c("black","red","black","red","black","red","black","red")
  mypch=c(19,17,19,17,19,17)
  points(1:length(nAB_l),nAB_l[1:length(nAB_l)][ord]/max(nAB_l),pch=mypch[res[[2]]],col=mycols[res[[5]]])
  #-------------------------------------^^
  covmin<-min(myinfosamples$X..Targets.Covered.to.20x.or.greater[ord])
  covmax<-max(myinfosamples$X..Targets.Covered.to.20x.or.greater[ord])
  covpoints<-(1.75-1.06)*(myinfosamples$X..Targets.Covered.to.20x.or.greater[ord]-covmin)/(covmax-covmin)+1.06
  exomin<-min(myinfosamples$Total.Exome.Sequence[ord])
  exomax<-max(myinfosamples$Total.Exome.Sequence[ord])
  exopoints<-(1.75-1.06)*(myinfosamples$Total.Exome.Sequence[ord]-exomin)/(exomax-exomin)+1.06
  is.na(myinfosamples$LC.Non.Duplicated.Aligned.Coverage)<-0
  alcmin<-min(myinfosamples$LC.Non.Duplicated.Aligned.Coverage[ord])
  alcmax<-max(myinfosamples$LC.Non.Duplicated.Aligned.Coverage[ord])
  alcpoints<-(1.75-1.06)*(myinfosamples$LC.Non.Duplicated.Aligned.Coverage[ord]-alcmin)/(alcmax-alcmin)+1.06
  lines(1:length(nAB_l),covpoints,col='aquamarine4',lwd=2)
  lines(1:length(nAB_l),exopoints,col='coral',lwd=2)
  lines(1:length(nAB_l),alcpoints,col='darkgoldenrod2',lwd=2)
  #points(nAandB_l[1:length(nAB_l)][ord]/10,pch=19, col = 'red') #wait, this is from log files. It is not good, because in log only positions that are significant
  if (showmutperssample){ points((mutpoints^2)/100,pch=19, col = 'red')}
  lines(1:length(nAB_l), mres/max(nAB_l), col = 'black')
  lines(1:length(nAB_l), upper/max(nAB_l), col = 'gray')
  lines(1:length(nAB_l), lower/max(nAB_l), col = 'gray')
  #lines(1:length(nAB_l),infosamples$Total.LC.Sequence[ord]/10000000,col='aquamarine3',lwd=2)
  lines(1:length(nAB_l),covpoints,col='aquamarine4',lwd=2)
  lines(1:length(nAB_l),exopoints,col='coral',lwd=2)
#  myinfosamples$Main.project.LC.Centers[ord]
  for (j in 1:length(nAB_l))
  {
    rect(j-0.5, 1.78, j+0.5, 2.03, col = seqcenter_col[j], border = "black") # coloured
    if (seqcenter[j]=='BI,MPIMG') {rect(j-0.5, (1.78+2.03)/2, j+0.5, 2.03, col = "coral2", border = NULL)}
    if (seqcenter[j]=='BCM,BGI') {rect(j-0.5, (1.78+2.03)/2, j+0.5, 2.03, col = "aquamarine2", border = NULL)}
    rect(j-0.5, 2.03, j+0.5, 2.03+(2.28-2.03)/3, col = myinfosamples$Has.Exome.LOF.Genotypes[ord][j], border = "black") # coloured
    rect(j-0.5, 2.03+(2.28-2.03)/3, j+0.5, 2.03+2*(2.28-2.03)/3, col = myinfosamples$Has.Axiom.Genotypes[ord][j], border = "black") # coloured
    rect(j-0.5, 2.03+2*(2.28-2.03)/3, j+0.5, 2.28, col = myinfosamples$Has.Affy.6.0.Genotypes[ord][j], border = "black") # coloured
    rect(j-0.5, 2.28, j+0.5, 2.53, col = etplat_col[j], border = NULL) # coloured
    if (etplat[j]=='BCM,BGI') {rect(j-0.5, (1.78+2.03)/2, j+0.5, 2.03, col = "coral2", border = NULL)}
  }
  #lines(1:length(nAB_l),as.numeric(seqcenter)-2000,col='cornflowerblue',lwd=2)
}

infoplot<-function(ipop,mymaxk=6,mypthr=0.01,mynABdefault=FALSE,mynAB=0,showreal_confidence_interval=TRUE,logpop_reschr_nAB_l=0,reshufflings=0,mylog=FALSE,showmutperssample=FALSE,mutpersample=FALSE,mypopname="",finite_mixture_criterion=2,myylab="nAB                              ",continuous_fields=c("X..Targets.Covered.to.20x.or.greater","Total.Exome.Sequence","LC.Non.Duplicated.Aligned.Coverage")) #version with real reshuffling of chromosomes
{
#if showreal_confidence_interval==TRUE -> logpop_reschr_nAB_l has to be supplied as a list with first level reshufflings, and second population. At least 2 reshufflings are needed to calculate the variance
  if (mynABdefault) #whether nAB info is in object logpop or other name (supplied as mynAB)
  {
  nAB_l<-unlist(logpop[[ipop]][1])
  } else {nAB_l<-mynAB[[ipop]]} 
  #nAandB_l<-unlist(logpop[[ipop]][2])
  ord<-order(nAB_l[1:length(nAB_l)])
  mres<-mean(nAB_l)
  myinfosamples<-as.data.frame(infosamples[[ipop]])
  seqcenter<-myinfosamples$Main.project.LC.Centers[ord]
  seqcenter<-as.character(seqcenter)
  seqcenter_col<-as.character(seqcenter)
  seqcenter_col[seqcenter=='BI']<-"firebrick"
  seqcenter_col[seqcenter=='BI,MPIMG']<-"firebrick"
  #seqcenter_col[seqcenter=='BI,MPIMG']<-"firebrick1"
  seqcenter_col[seqcenter=='MPIMG']<-"coral2"
  seqcenter_col[seqcenter=='BGI']<-"cadetblue"
  seqcenter_col[seqcenter=='BCM,BGI']<-"cadetblue"
  #seqcenter_col[seqcenter=='BCM,BGI']<-"cadetblue2"
  seqcenter_col[seqcenter=='BCM']<-"aquamarine2"
  seqcenter_col[seqcenter=='ILLUMINA']<-"chartreuse3"
  seqcenter_col[seqcenter=='SC']<-"brown1"
  seqcenter_col[seqcenter=='WUGSC']<-"darkgoldenrod1"
  etplat<-as.character(myinfosamples$ET.Pilot.Platforms[ord])
  etplat[is.na(etplat)]<-0
  etplat[etplat==""]<-0
  etplat_col<-etplat
  etplat_col[etplat=="ILLUMINA"]<-"deepskyblue2"
  etplat_col[etplat=="LS454"]<-"coral2"
  #etplat[etplat=="ILLUMINA,LS454"]<-"orchid3"
  etplat_col[etplat=="ILLUMINA,LS454"]<-"deepskyblue2"
  etplat_col[etplat=="0"]<-"floralwhite"
  #seqcenter<-as.character(seqcenter)
  #length(myinfosamples$Has.Exome.LOF.Genotypes)
  myinfosamples$Has.Exome.LOF.Genotypes[is.na(myinfosamples$Has.Exome.LOF.Genotypes)]<-0
  myinfosamples$Has.Exome.LOF.Genotypes[myinfosamples$Has.Exome.LOF.Genotypes==1]<-"cornsilk4"
  myinfosamples$Has.Exome.LOF.Genotypes[myinfosamples$Has.Exome.LOF.Genotypes==0]<-"floralwhite"
  myinfosamples$Has.Axiom.Genotypes[is.na(myinfosamples$Has.Axiom.Genotypes)]<-0
  myinfosamples$Has.Axiom.Genotypes[myinfosamples$Has.Axiom.Genotypes==1]<-"cornsilk4"
  myinfosamples$Has.Axiom.Genotypes[myinfosamples$Has.Axiom.Genotypes==0]<-"floralwhite"
  myinfosamples$Has.Affy.6.0.Genotypes[is.na(myinfosamples$Has.Affy.6.0.Genotypes)]<-0
  myinfosamples$Has.Affy.6.0.Genotypes[myinfosamples$Has.Affy.6.0.Genotypes==1]<-"cornsilk4"
  myinfosamples$Has.Affy.6.0.Genotypes[myinfosamples$Has.Affy.6.0.Genotypes==0]<-"floralwhite"
  h_nABplot<-1.03
  h_covplot<-1.75
  margin<-0.03
  myylab=paste0(myylab,"       info                              ")
  plot(x=1:length(nAB_l),y=rep(mres,length(nAB_l)),ylim=c(0,2.53),type = "n",xlab=paste0(mypopname," samples"),ylab=myylab,yaxs="i",yaxt="n")
  axis(2, at=c(0,0.2,0.4,0.6,0.8,0.99,h_nABplot+margin,(h_nABplot+margin+h_covplot)/2,h_covplot,(h_covplot+margin+2.03)/2,(2.03+2.28)/2,(2.28+2.53)/2),labels=c(0,0.2,0.4,0.6,0.8,1,0,0.5,1,"SC","GT","ETP"),las=2)
  #--frame-of-plot------------------------v
  abline(0,0)
  abline(h_nABplot,0)
  abline(h_covplot,0)
  #abline(1.53,0)
  abline(2.03,0)
  abline(2.28,0)
  #---------------------------------------^
  #--PLOTPOINTS--------------------------vv
  #-------plotting null distr
  if (showreal_confidence_interval) #add nAB from reshufflings of chromosomes
  {
    nperm<-length(logpop_reschr_nAB_l) 
    myreschrnAB<-sapply(1:nperm, function(x) {nAB_lr<-logpop_reschr_nAB_l[[x]][[ipop]]; ordr<-order(nAB_lr[1:length(nAB_lr)]);return(nAB_lr[ordr])})
    myreschr_mean<-apply(myreschrnAB,MARGIN=1,mean)
    myreschr_sd<-apply(myreschrnAB,MARGIN=1,sd)
    reshufflings<-sapply(1:length(myreschr_mean), function(z) rnorm(nperm,myreschr_mean[z],myreschr_sd[z])/max(nAB_l))
    cdat <- as.list(as.data.frame(reshufflings))
    names(cdat)[1] <- "x"  # vioplot() needs the first element to be called 'x'
    do.call(vioplot,c(cdat,list(col="cadetblue",colMed="cadetblue",add=TRUE)))
  }
  #-------plotting points
  if (showmutperssample){ points((mutpersample[ord])/max(mutpersample[ord]),pch=19, col = 'gray')}
  res<-partition_samples_nAB(nAB_l[ord], maxk=mymaxk, phtr=mypthr)
  mycols=c("black","red","black","red","black","red","black","red")
  mypch=c(19,17,19,17,19,17)
  if (finite_mixture_criterion==0) { mycols=c("black","black","black","black","black","black","black","black") } else if (finite_mixture_criterion==1) { mypch=c(19,19,19,19,19,19) }
  points(1:length(nAB_l),nAB_l[1:length(nAB_l)][ord]/max(nAB_l),pch=mypch[res[[5]]],col=mycols[res[[2]]])
  
  if (FALSE)
  {
  rm(res)
  res<<-partition_samples_nAB(nAB_l,maxk=maxk,pthr=mypthr)
  for (myj in 1:res[[3]])
    {
    if ((myj%%2)==0){mycolpoints<-"gold"} else {mycolpoints<-"black"}
    points((1:length(nAB_l))[res[[2]]==myj],nAB_l[1:length(nAB_l)][ord][res[[2]]==myj]/max(nAB_l),pch=19,col=mycolpoints)
    }
  }
  #-------------------------------------^^
  is.na(myinfosamples[[continuous_fields[1]]])<-0
  covmin<-min(myinfosamples[[continuous_fields[1]]][ord])
  covmax<-max(myinfosamples[[continuous_fields[1]]][ord])
  covpoints<-(1.75-1.06)*(myinfosamples[[continuous_fields[1]]][ord]-covmin)/(covmax-covmin)+1.06
  is.na(myinfosamples[[continuous_fields[2]]])<-0
  exomin<-min(myinfosamples[[continuous_fields[2]]][ord])
  exomax<-max(myinfosamples[[continuous_fields[2]]][ord])
  exopoints<-(1.75-1.06)*(myinfosamples[[continuous_fields[2]]][ord]-exomin)/(exomax-exomin)+1.06
  is.na(myinfosamples[[continuous_fields[3]]])<-0
  alcmin<-min(myinfosamples[[continuous_fields[3]]][ord])
  alcmax<-max(myinfosamples[[continuous_fields[3]]][ord])
  alcpoints<-(1.75-1.06)*(myinfosamples[[continuous_fields[3]]][ord]-alcmin)/(alcmax-alcmin)+1.06
  lines(1:length(nAB_l),covpoints,col='aquamarine4',lwd=2)
  lines(1:length(nAB_l),exopoints,col='coral',lwd=2)
  lines(1:length(nAB_l),alcpoints,col='darkgoldenrod2',lwd=2)
  #points(nAandB_l[1:length(nAB_l)][ord]/10,pch=19, col = 'red') #wait, this is from log files. It is not good, because in log only positions that are significant
#  lines(1:length(nAB_l), mres/max(nAB_l), col = 'black')
#  lines(1:length(nAB_l), upper/max(nAB_l), col = 'gray')
#  lines(1:length(nAB_l), lower/max(nAB_l), col = 'gray')
  #lines(1:length(nAB_l),infosamples$Total.LC.Sequence[ord]/10000000,col='aquamarine3',lwd=2)
  lines(1:length(nAB_l),covpoints,col='aquamarine4',lwd=2)
  lines(1:length(nAB_l),exopoints,col='coral',lwd=2)
#  myinfosamples$Main.project.LC.Centers[ord]
  for (j in 1:length(nAB_l) )#length(nAB_l))
  {
    rect(j-0.5, 1.78, j+0.5, 2.03, col = seqcenter_col[j], border = "black") # coloured
    if (seqcenter[j]=='BI,MPIMG') {rect(j-0.5, (1.78+2.03)/2, j+0.5, 2.03, col = "coral2", border = NULL)}
    if (seqcenter[j]=='BCM,BGI') {rect(j-0.5, (1.78+2.03)/2, j+0.5, 2.03, col = "aquamarine2", border = NULL)}
    rect(j-0.5, 2.03, j+0.5, 2.03+(2.28-2.03)/3, col = myinfosamples$Has.Exome.LOF.Genotypes[ord][j], border = "black") # coloured
    rect(j-0.5, 2.03+(2.28-2.03)/3, j+0.5, 2.03+2*(2.28-2.03)/3, col = myinfosamples$Has.Axiom.Genotypes[ord][j], border = "black") # coloured
    rect(j-0.5, 2.03+2*(2.28-2.03)/3, j+0.5, 2.28, col = myinfosamples$Has.Affy.6.0.Genotypes[ord][j], border = "black") # coloured
    rect(j-0.5, 2.28, j+0.5, 2.53, col = etplat_col[j], border = NULL) # coloured
    if (etplat[j]=='BCM,BGI') {rect(j-0.5, (1.78+2.03)/2, j+0.5, 2.03, col = "coral2", border = NULL)}
  }
  #lines(1:length(nAB_l),as.numeric(seqcenter)-2000,col='cornflowerblue',lwd=2)
}


logreg_snps2nAB<-function(mygen,mynAB,mymut)#calculate logistic regression to see snps that contribute the most to individuals with highest nAB
{   #roger suggests on comparing the models because even if "perfect separation" and p-value with Wald approximation is useless, still the fit is very good and likelihood ratio works perfectly
    model <- glm( mygen ~ mynAB + mymut,family=binomial(link='logit'))
    model0 <- glm( mygen ~ mymut,family=binomial(link='logit'))
    return(anova(model0, model, test = "Chisq")[[5]][2])
}

if (FALSE){
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
	if (sum(mygen)/length(mygen)>0.05 && sum(mygen)/length(mygen)<0.95)
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
}

corr_nAvsnAB<-function(nAall_l,logpop_nAB,l_mutsamples)
{
#it calculates the logistic correlation between the presence of alternative allele and nAB per sample.
#one way is having nAB as independent variable and nA as dependent. That can be a logistic, by splitting single genotypes. Other way would be a normal regression with values 0,1,2.
#version1: it takes as arguments lists (each element is a population) with myres (genotypes in 0,1(het) and 2(homozygotes)), logpop(nAB per sample), l_mutsamples (number of alt allele per individual).
#version2: dataset is going nAall, logpop_nAB,l_mutsamples
print("converting data to format")
nAall_l_temp<-nAall_l;logpop_nAB_temp<-list();l_mutsamples_temp<-list()
for (ipop in 1:length(logpop_nAB))
    {
    print(paste("population: ",ipop))
    nAall_l_temp[[ipop]]<-t(apply(nAall_l_temp[[ipop]],MARGIN=1,function(x) c(sapply(x,function(y) if (y==0){return(0)} else {return(1)}),sapply(x,function(y) if (y==2){return(1)} else {return(0)}))))
    logpop_nAB_temp[[ipop]]<-c(logpop_nAB[[ipop]],logpop_nAB[[ipop]])
    l_mutsamples_temp[[ipop]]<-c(l_mutsamples[[ipop]],l_mutsamples[[ipop]])
    }
print("calculating logistic regressions")
pvals_l<-list()
for (ipop in 1:length(logpop_nAB_temp))
    {
    print(paste("population: ",ipop))
    pvals_l[[ipop]]<-apply(nAall_l_temp[[ipop]],MARGIN=1,function(x) logreg_snps2nAB(x,logpop_nAB_temp[[ipop]],l_mutsamples_temp[[ipop]]))
    }
return(pvals_l)
}

#snps_intabfile_perpop<-c()
#for (i in 1:22)
#  {
#  snps_intabfile_perpop<-append(snps_intabfile_perpop,as.numeric(system(paste0("cat ",myfolder,"/chr",i,".freqlog | wc -l"),intern=TRUE))-1)
#  }

cleanup_snps_intabfile_perpop<-function(myfolder,myres) #create vector with number of snps that don't reach freq thr per chromosome
{
  currentchr<-1
  nextchrthr<-snps_intabfile_perpop[1]
  res<-rep(0,22)
  for (mylines in 1:dim(myres[[1]])[1])
  {
    if (mylines > nextchrthr && currentchr < 12) {currentchr<-currentchr+1;nextchrthr<-nextchrthr+snps_intabfile_perpop[currentchr+1]}
    abovethr<-0
    print(mylines)
    for (ipop in 1:12)
    {
      freq<-sum(myres[[ipop]][mylines,])/nspops[[ipop]]
      if (freq>0.05 || freq<0.95) {abovethr<-1}
    }
    if (abovethr==0) {res[currentchr]<-res[currentchr]+1}
  }
return(res)
}
#nothrinfiletab<-cleanup_snps_intabfile_perpop("above95/coding1000g",myres) #all 0..cool!

create_blocks<-function(mylog,space_between_indep_blocks=50000000)
{
res<-list()
for (mychr in 1:22)
  {
  thr<-c()
  mysnps<-unique(c(mylog[mylog[,1]==mychr,3],mylog[mylog[,2]==mychr,4]))
  mysnps<-mysnps[order(mysnps)]
  for (isnp in 1:length(mysnps))
    {
    if (isnp==1){thr<-append(thr,mysnps[1])} else
      {
      if (mysnps[isnp]>mysnps[isnp-1]+space_between_indep_blocks) { thr<-append(thr,mysnps[isnp])}
      }
    }
  res[[mychr]]<-thr
  }
return(res)
}


permute_block<-function(mydata,nperm) #starting from log file with chrA chrB posA posB gts.... pop_index generates permutations for nAB for a single block
{
tot<-sapply(5:(dim(mydata)[2]-1), function(x) sum(sapply(mydata[,x],nABf)))
tot<-sapply(1:nperm, function(x) sample(tot,replace=FALSE))
return(tot)
}


permute_chrAB<-function(mylog,nperm) #starting from genome-wide single pop log file applies permute_block considering chr pairs as blocks
{
for (chrA in 1:21)
  {
  for (chrB in (chrA+1):22)
    {
    print(c(chrA,chrB))
    templog<-mylog[mylog[,1]==chrA,]
    templog<-templog[templog[,2]==chrB,]
    if (chrA==1 && chrB==2) {myres<-permute_block(templog,nperm)} else
      {myres<-myres+permute_block(templog,nperm)}
    }
  }
return(myres)
}

permute_chrAB_byblock<-function(mylog,myblocks,nperm) #starting from genome-wide single pop log file applies permute_block considering chr pairs as blocks
{
myresnoinit<-TRUE
for (chrA in 1:21)
  {
  for (chrB in (chrA+1):22)
    {
    print(c(chrA,chrB))
    for (blockA in 1:length(myblocks[[chrA]]))
      {
      for (blockB in 1:length(myblocks[[chrB]]))
	{
	templog<-mylog[mylog[,1]==chrA,]
	templog<-templog[templog[,3]>=myblocks[[chrA]][blockA],]
	if (blockA<length(myblocks[[chrA]]))
	  {
	  templog<-templog[templog[,3]<myblocks[[chrA]][blockA+1],]
	  }
	templog<-templog[templog[,2]==chrB,]
	templog<-templog[templog[,4]>=myblocks[[chrB]][blockB],]
	if (blockB<length(myblocks[[chrB]]))
	  {
	  templog<-templog[templog[,4]<myblocks[[chrB]][blockB+1],]
	  }
	print(c(chrA,chrB,dim(templog)))
	if (dim(templog)[1]>0)
	  {
	  if (myresnoinit) {myres<-permute_block(templog,nperm); myresnoinit<-FALSE} else
	    {myres<-myres+permute_block(templog,nperm)}
	  }
	}
      }
    }
  }
return(myres)
}

generate_nulldistr_mixtures<-function(permbychr, myn,maxk=6) 
{
#generates the empirical null distribution for expected partitioning of samples. col 1 and 2 are the number of subgroups for LR and AICc, 7 and 8 LR-pvalue and RL
# permchr is the vector of permutations, myn is how many samples to take
#suited from many permutations, for example when one reshuffle chromosomes by blocks
for (x in 1:myn)
  {
  t_res<-0
  mymaxk<-maxk
  while (length(t_res)==1) #this internal loop prevent errors from non-convergence and excessive repetitions when convergence will never happen
    {
    if (mymaxk>1)
      {
      try(t_res<-unlist(partition_samples_nAB(permbychr[,x],maxk=mymaxk)[c(1,4,7,8)])) 
      mymaxk<-mymaxk-1
      } else {t_res<-c(1,1,1,1)}
    }
    if (x==1) {res<-t_res} else {res<-rbind(res,t_res)}
  }
return(res)
}

test_inhomogeneous_nAB<-function(logpop_nAB,logpop_reschr_nAB_l) 
#test that nAB is distributed non uniformly across samples by comparing with a random null datasets (e.g. shuffled chromosomes)
#suited for small number of reshuffling (use a t_test)
#arguments are 
#logpop_nAB: (list with vectors of nAB for samples, where each element of list is a population
#logpop_reschr_nAB_l: (list of lists like logpop_nAB, with first level that is the number of reshuffling
#it returns a p-value based on log p-values of finite mixtures, and one based on standard deviation, that seems more powerful
{
res<-list()
res$nAB_mean<-sapply(1:length(logpop_reschr_nAB_l), function(x) sapply(1:length(logpop_reschr_nAB_l[[x]]), function(y) mean(logpop_reschr_nAB_l[[x]][[y]])))
res$nAB_sd_random<-sapply(1:length(logpop_reschr_nAB_l), function(x) sapply(1:length(logpop_reschr_nAB_l[[x]]), function(y) sd(logpop_reschr_nAB_l[[x]][[y]])))
for (i in 1:10){res$nAB_sd_random<-t(apply(res$nAB_sd_random,MARGIN=1,function(x) if (max(x)-min(x)>0){return(x)} else {return(rnorm(length(x),x,10^(-8)))}))}
res$nAB_sd<-sapply(1:length(logpop_nAB), function(y) sd(logpop_nAB[[y]]))
res$nAB_sd_t.test<-sapply(1:length(logpop_nAB), function(y) t.test(res$nAB_sd_random[y,],mu=res$nAB_sd[y],alternative="less")$p.value)
mixtures_data<-lapply(1:length(logpop_nAB),function(ipop) partition_samples_nAB(logpop_nAB[[ipop]],reorder=TRUE))
mixtures_reschr<-lapply(1:length(logpop_reschr_nAB_l), function(y) lapply(1:length(logpop_reschr_nAB_l[[1]]),function(ipop) partition_samples_nAB(logpop_reschr_nAB_l[[y]][[ipop]],reorder=TRUE)))
res$AICc_relativelik_reschr<-sapply(1:length(mixtures_reschr),function(x) sapply(1:length(mixtures_reschr[[1]]),function(ipop) mixtures_reschr[[x]][[ipop]]$AICc_relativelik))
res$AICc_relativelik_data<-sapply(1:length(mixtures_data),function(ipop) mixtures_data[[ipop]]$AICc_relativelik)
#to avoid errors in cases in which both p-values are identical.
for (i in 1:10){res$AICc_relativelik_reschr<-t(apply(res$AICc_relativelik_reschr,MARGIN=1,function(x) if (max(x)-min(x)>0){return(x)} else {return(rnorm(length(x),x,10^(-10)))}))}

#Since fraction t-test does not work-> log
res$AICc_relativelik_t.test<-sapply(1:length(logpop_nAB), function(y) t.test(log(res$AICc_relativelik_reschr[y,]+10^(-10)),mu=log(res$AICc_relativelik_data[y]+10^(-10)),alternative="greater")$p.value)
return(res)
}




rebuild_dataLDLD<-function(removed_snps,dataLDLD,fromremovedsnps=TRUE) #when fromremovedsnps=FALSE give vector with both snps . separated
{    #extract links or any pairwise information from a single vector of snps either to be kept or removed
if (fromremovedsnps)
    {
    dataLDLDA<-as.character(paste0(dataLDLD[,1],".",dataLDLD[,3]))
    dataLDLD<-dataLDLD[is.na(match(dataLDLDA,removed_snps)),]
    dataLDLDB<-as.character(paste0(dataLDLD[,2],".",dataLDLD[,4]))
    dataLDLD<-dataLDLD[is.na(match(dataLDLDB,removed_snps)),]
    } else
    {
    dataLDLDA<-as.character(paste0(dataLDLD[,1],".",dataLDLD[,3],".",dataLDLD[,2],".",dataLDLD[,4]))
    dataLDLD<-dataLDLD[!is.na(match(dataLDLDA,removed_snps)),]
    }
    return(dataLDLD)
}

test_nulldistr_mixtures<-function(npops,nreps,logpoptot,logpop,myblocks,popstart=1) 
{ #return list with all the permutations as well all the pvalues. It needs blocks generated by permute_chrAB_byblock and nABbypop file with all the genotypes per link, all for all populations.
  #notice that logpoptot is the raw log with all genotypes, while nABbypop is the output of pairsANDlogfiles2logpop with only the nAB per individual as [[ipop]][[1]]
  myres<-list()
  for (ipop in popstart:npops)
    {
    permbychr_removal1<-permute_chrAB_byblock(logpoptot[[ipop]],myblocks,nreps)
    permbychr_removal1<-sapply (1:dim(permbychr_removal1)[2],function(x) permbychr_removal1[,x][order(permbychr_removal1[,x])])
    res<-list()
    res$raw<-permbychr_removal1
    res$mean_byindiv<-sapply (1:dim(permbychr_removal1)[1],function(x) mean(permbychr_removal1[x,]))
    res$sd_byindiv<-sapply (1:dim(permbychr_removal1)[1],function(x) sd(permbychr_removal1[x,]))
    res$sd<-sapply (1:dim(permbychr_removal1)[2],function(x) sd(permbychr_removal1[,x]))
    res$pval_sd<-sum(as.numeric(res$sd>=sd(logpop[[ipop]][[1]])))/length(res$sd)
    res$rawmodels<-generate_nulldistr_mixtures(permbychr_removal1,nreps)
    fitmodel<-0
    mymaxk<-6
    while (length(fitmodel)==1) #this internal loop prevent errors from non-convergence and excessive repetitions when convergence will never happen
      {
      if (mymaxk>1)
	{
	try(fitmodel<-partition_samples_nAB(logpop[[ipop]][[1]],maxk=mymaxk)[c(1,4,7,8)])
	mymaxk<-mymaxk-1
	} else {fitmodel<-c(1,1,1,1)}
      }
    res$pval_nLR<-sum(as.numeric(as.numeric(res$rawmodels[,1])>=unlist(fitmodel)[1]))/nreps
    res$pval_nAIC<-sum(as.numeric(as.numeric(res$rawmodels[,2])>=unlist(fitmodel)[2]))/nreps
    res$pval_pLR<-sum(as.numeric(as.numeric(res$rawmodels[,3])<=unlist(fitmodel)[3]))/nreps
    res$pval_pAIC<-sum(as.numeric(as.numeric(res$rawmodels[,4])<=unlist(fitmodel)[4]))/nreps
    myres[[ipop]]<-res
    }
  return(myres)
}

logpop2gen<-function(mylog)#go from log files to gen files, with 0 1 and 2 for homozygotes
  {
    mysnps1<-paste(mylog[,dim(mylog)[2]-1],mylog[,1],sep=".")
    mysnps2<-paste(mylog[,dim(mylog)[2]],mylog[,2],sep=".")
    mysnps1u<-unique(mysnps1)
    mysnps2u<-setdiff(unique(mysnps2),mysnps1)
    for (isnp1 in 1:length(mysnps1u))
      {
	print(isnp1)
	temp<-subset(mylog,mysnps1==mysnps1u[isnp1])[1,]
	mynA<-t(sapply(3:(length(temp)-3), function(x) sapply(temp[x],nAf)))
	if (isnp1==1) {mydata<-as.data.frame(list(mynA,mypos=mysnps1u[isnp1]))} else
	  {
	  mydata<-rbind(mydata,as.data.frame(list(mynA,mypos=mysnps1u[isnp1])))
	  }
      }
    for (isnp2 in 1:length(mysnps2u))
      {
	print(isnp2)
	temp<-subset(mylog,mysnps2==mysnps2u[isnp2])[1,]
	mynA<-t(sapply(3:(length(temp)-3), function(x) sapply(temp[x],nBf)))
	mydata<-rbind(mydata,as.data.frame(list(mynA,mypos=mysnps2u[isnp2])))
      }
    return(mydata)
  }

#slower version
empiricall_fdr<-function(distr1,distrnull,howmanynull) #howmanynull=number of null permutations of set
{
    res<-rep(1,length(distr1))
    distr1ordered<-distr1[order(distr1)]
    for ( i in 1:length(distr1))
      {
      res[i]<-sum(as.numeric(distrnull<=distr1ordered[i]))/(i*howmanynull)
      }
return(res)
}
#faster version with ecdf. distr1 is distr of observed; distrnull is a list of pvalues from empirical null; howmanynull is the number of how many null simulations
empiricall_fdr<-function(distr1,distrnull,howmanynull) {
    null_fdr_f<-ecdf(distrnull);
    obs_fdr_f<-ecdf(distr1);
    fdr_f<-function(y) (length(distrnull)*null_fdr_f(y)/howmanynull)/(length(distr1)*obs_fdr_f(y))
    myfdr<-sapply(distr1,function(y) fdr_f(y))
    return(myfdr)
}
distr1<-rep(0.01,100)
distrnull<-c(0.01,rep(0.02,99))




}
#===================================SNPS discovery==================================================
{

require(metap)
#Roger changes with this (analogous) and to drop only fixed effect but not random slope.This is like assuming that variability but that it could go in any direction
#expandDoubleVerts(nA ~ z.nAB + (1+z.nAB||population)+(1|sampleID)) #nA ~ z.nAB + ((1 | population) + (0 + z.nAB | population)) + (1 | sampleID)
#is there a better proxy rather than nAB? which associated? pair for which positive. #would it solve bias towards rare? actually I should only consider positive linkage
fx<-function(data_temp){
my.model = glmer( nA ~ nAB + (1|population) + (0+nAB|population)+(1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
my.model.null = glmer( nA ~ (1|population) +(1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
anova(my.model,my.model.null)[['Pr(>Chisq)']][[2]]
}
fx1<-function(data_temp){
my.model = glmer( nA ~ nAB + (1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
my.model.null = glmer( nA ~ (1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
anova(my.model,my.model.null)[['Pr(>Chisq)']][[2]]
}


fishersmethod<-function(z) if ( length(z)>=2) {return(sumlog(z)$p)} else return(z)

#this is by restricting only on populations where I have power enough because of freq. However in principle this should also capture overall

myfx_parse<-function(iline,myfreq=0.05,mypvalue_threshold=0.05) { 
#used if dataset provided from general analyses in block of data0 and nAall
data_temp0<-cbind(nA=t(nAall[iline,]),data0)
setnames(data_temp0,c('nA','sampleID','population','nAB','nAmut'))
return(myfx(data_temp0,myfreq,mypvalue_threshold))
}
myfx<-function(data_temp0,myfreq=0.05,mypvalue_threshold=0.05){
data_temp<-data_temp0[unlist(unname(aggregate(data_temp0$nA,by=list(data_temp0$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
if ( nrow(data_temp)>0 ){
glmres<-preliminary_check_glm_nAvsnAB(data_temp$nA,data_temp$nAB,data_temp$population)
if (glmres<mypvalue_threshold ) { if (length(unique(data_temp$population))>1) {return(c(fx(data_temp),fx(data_temp0),glmres))} else {return(c(fx1(data_temp),fx(data_temp0),glmres))} }}
}

myfx_glm<-function(iline,myfreq=0.05) {
data_temp<-cbind(nA=t(nAall[iline]),data0)
setnames(data_temp,c('nA','sampleID','population','nAB'))
data_temp<-data_temp[unlist(unname(aggregate(data_temp$nA,by=list(data_temp$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
return(preliminary_check_glm_nAvsnAB(data_temp$nA,data_temp$nAB,data_temp$population))
}
myfx_glmm<-function(iline,myfreq=0.05) {
data_temp<-cbind(nA=t(nAall[iline]),data0)
setnames(data_temp,c('nA','sampleID','population','nAB'))
data_temp<-data_temp[unlist(unname(aggregate(data_temp$nA,by=list(data_temp$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
if (length(unique(data_temp$population))>1) {return(fx(data_temp))} else {return(fx1(data_temp))}
}
myfx_glmm_nofilter<-function(iline) {
data_temp<-cbind(nA=t(nAall[iline]),data0)
setnames(data_temp,c('nA','sampleID','population','nAB'))
return(fx(data_temp))
}
preliminary_check_glm_nAvsnAB<-function(mynA,mynAB,population){
res<-sapply(unique(population), function(i) logreg_snps2nAB(mynA[population==i],mynAB[population==i]))
res<-res[!is.na(res)]
return(fishersmethod(res))
}
logreg_snps2nAB<-function(mygen,mynAB,myfreq=0.05)#calculate logistic regression to see snps that contribute the most to individuals with highest nAB
{   #roger suggests on comparing the models because even if "perfect separation" and p-value with Wald approximation is useless, still the fit is very good and likelihood ratio works perfectly
    if ( sum(mygen)/length(mygen) >= myfreq & sum(mygen)/length(mygen) < 0.95 ) {
    model <- glm( mygen ~ mynAB,family=binomial(link='logit'))
    model0 <- glm( mygen ~ 1,family=binomial(link='logit'))
    return(anova(model0, model, test = "Chisq")[[5]][2])
    } else return(NA)
}    
logreg_snps2nAB<-function(mynA,mynAB)#calculate logistic regression to see snps that contribute the most to individuals with highest nAB
{   #roger suggests on comparing the models because even if "perfect separation" and p-value with Wald approximation is useless, still the fit is very good and likelihood ratio works perfectly
    model <- glm( mynA ~ mynAB,family=binomial(link='logit'))
    model0 <- glm( mynA ~ 1,family=binomial(link='logit'))
    pvalue<-anova(model0, model, test = "Chisq")[[5]][2]
    if (is.na(pvalue)){return(1)} else {return(pvalue)}
    }    

    
myfx_parse_mut<-function(iline,myfreq=0.05,mypvalue_threshold=0.05) { 
#used if dataset provided from general analyses in block of data0 and nAall
data_temp0<-cbind(nA=t(nAall[iline,]),data0)
setnames(data_temp0,c('nA','sampleID','population','nAB','nAmut'))
return(myfx_mut(data_temp0,myfreq,mypvalue_threshold))
}
myfx_mut<-function(data_temp0,myfreq=0.05,mypvalue_threshold=1){
data_temp<-data_temp0[unlist(unname(aggregate(data_temp0$nA,by=list(data_temp0$population),function(x) if (sum(x)/length(x)>=myfreq & sum(x)/length(x)<=(1-myfreq)){return(rep(TRUE,length(x)))}else return(rep(FALSE,length(x))))[,2])),]
if ( nrow(data_temp)>0 ){
glmres<-preliminary_check_glm_nAvsnAB_mut(data_temp$nA,data_temp$nAB,data_temp$population,data_temp$nAmut)
if (glmres<mypvalue_threshold ) { if (length(unique(data_temp$population))>1) {return(c(fx_mut(data_temp),fx_mut(data_temp0),glmres))} else {return(c(fx1_mut(data_temp),fx_mut(data_temp0),glmres))} }}
}
preliminary_check_glm_nAvsnAB_mut<-function(mynA,mynAB,population,mynAmut){
res<-sapply(unique(population), function(i) logreg_snps2nAB_mut(mynA[population==i],mynAB[population==i],mynAmut[population==i]))
res<-res[!is.na(res)]
return(fishersmethod(res))
}
logreg_snps2nAB_mut<-function(mynA,mynAB,mynAmut)#calculate logistic regression to see snps that contribute the most to individuals with highest nAB
{   #roger suggests on comparing the models because even if "perfect separation" and p-value with Wald approximation is useless, still the fit is very good and likelihood ratio works perfectly
    model <- glm( mynA ~ mynAB+mynAmut,family=binomial(link='logit'))
    model0 <- glm( mynA ~ mynAmut,family=binomial(link='logit'))
    pvalue<-anova(model0, model, test = "Chisq")[[5]][2]
    if (is.na(pvalue)){return(1)} else {return(pvalue)}
    }    
fx_mut<-function(data_temp){
my.model = glmer( nA ~ nAB + nAmut + (1|population) + (0+nAmut|population)+(0+nAB|population)+(1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
my.model.null = glmer( nA ~ nAmut + (0+nAmut|population) + (1|population) +(1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
anova(my.model,my.model.null)[['Pr(>Chisq)']][[2]]
}
fx1_mut<-function(data_temp){
my.model = glmer( nA ~ nAB + nAmut + (1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
my.model.null = glmer( nA ~ nAmut + (1|sampleID), data=data_temp,family=binomial(link='logit'),control=contr)
anova(my.model,my.model.null)[['Pr(>Chisq)']][[2]]
}


}
