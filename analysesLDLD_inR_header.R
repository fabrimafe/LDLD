#==================================================
#=========== p-values and fdr =====================
#==================================================
#===========fisher's method
if (FALSE)
{
  pvalues<-c(0.2,0.4,0.5)
  df<-2*length(pvalues)
  pchisq( -2*sum(log(pvalues)), df, lower.tail=FALSE)
  pvalues<-c(0.2)
  df<-2*length(pvalues)
  pchisq( -2*sum(log(pvalues)), df, lower.tail=FALSE)
  pvalues<-c(0.05)
  df<-2*length(pvalues)
  pchisq( -2*sum(log(pvalues)), df, lower.tail=FALSE)
  #but from my calculation don't coincide. Why?
  pchisq(27.082154, 2, lower.tail=FALSE) #0.000001315785 $#actually they do!
  #1-0.999995
  #0.000001 #in theory it should never print 999 case, but there might be discrepancies between T2 and pfisher. In this cases actually either Kuli or pfisher less conservative than T2
}
#==========================


#==================================================
#=========== PEGAS + APE ==========================
#==================================================
library("pegas")
x<-c("A-A","A-a","a-a")
gts2loci<-function(gts)
{
loc1<-c(rep("A-A",gts[[1]]+gts[[2]]+gts[[3]]),rep("A-T",gts[[4]]+gts[[5]]+gts[[6]]),rep("T-T",gts[[7]]+gts[[8]]+gts[[9]]))
loc2<-c(rep("A-A",gts[[1]]),rep("A-T",gts[[2]]),rep("T-T",gts[[3]]),rep("A-A",gts[[4]]),rep("A-T",gts[[5]]),rep("T-T",gts[[6]]),rep("A-A",gts[[7]]),rep("A-T",gts[[8]]),rep("T-T",gts[[9]]))
#return(cbind(as.loci(loc1),as.loci(loc2)))
return(as.loci(as.data.frame(cbind(loc1,loc2)),allele.sep="-"))
}

#verified that LD2 from Pegas and my formula give same res
#f_DeltaAB(1:9) #-0.02666667
#LD2(gts2loci(1:9)) -0.0266667




#===================================================
#====== HARDY WEINBERG PERMUTATIONS OF GENOTYPES #==
#===================================================
#function that simulate HW genotypes, and calculate f_DeltaAB. This version seems broken, while HWE2 works
HWE<-function(nA,nB,Ni)
  {
  s<-sample(c(rep(1,nA),rep(0,2*Ni-nA)),Ni,replace=FALSE)
  sA<-s+sample(c(rep(1,nA-sum(s)),rep(0,2*Ni-(nA-sum(s)))),Ni,replace=FALSE)
  s<-sample(c(rep(1,nB),rep(0,2*Ni-nB)),Ni,replace=FALSE)
  sB<-s+sample(c(rep(1,nB-sum(s)),rep(0,2*Ni-(nB-sum(s)))),Ni,replace=FALSE)
  s<-rep(0,9)
  for (i in 1:Ni) #generate HWE genotypes
    {
    if ( sA[[i]]== 0 && sB[[i]] == 0 ) s[1]<-s[1]+1
    else if ( sA[[i]]==0 && sB[[i]] ==1 ) s[2]<-s[2]+1
    else if ( sA[[i]]==0 && sB[[i]] ==2 ) s[3]<-s[3]+1
    else if ( sA[[i]]==1 && sB[[i]] ==0 ) s[4]<-s[4]+1
    else if ( sA[[i]]==1 && sB[[i]] ==1 ) s[5]<-s[5]+1
    else if ( sA[[i]]==1 && sB[[i]] ==2 ) s[6]<-s[6]+1
    else if ( sA[[i]]==2 && sB[[i]] ==0 ) s[7]<-s[7]+1
    else if ( sA[[i]]==2 && sB[[i]] ==1 ) s[8]<-s[8]+1
    else if ( sA[[i]]==2 && sB[[i]] ==2 ) s[9]<-s[9]+1
    }
  f_DeltaAB(s)
  }

HWE2<-function(nA,nB,Ni)
  {
  f_DeltaAB(gtsin012format2gts(simulateHWgts(nA,nB,Ni)))
  }
  
simulateHW_gtsin012format<-function(nA,nB,Ni)
  {
  s<-sample(c(rep(1,nA),rep(0,2*Ni-nA)),2*Ni,replace=FALSE)
  sA<-s[1:Ni]+s[(Ni+1):(2*Ni)]
  s<-sample(c(rep(1,nB),rep(0,2*Ni-nB)),2*Ni,replace=FALSE)
  sB<-s[1:Ni]+s[(Ni+1):(2*Ni)]
  rbind(sA,sB)
  }

gtsin012format2gts<-function(s)
{
  sA<-s[1,]
  sB<-s[2,]
  s<-rep(0,9)
  for (i in 1:length(sA)) #generate HWE genotypes
    {
	 if ( sA[[i]]==0 && sB[[i]] ==0 ) s[1]<-s[1]+1
    else if ( sA[[i]]==0 && sB[[i]] ==1 ) s[2]<-s[2]+1
    else if ( sA[[i]]==0 && sB[[i]] ==2 ) s[3]<-s[3]+1
    else if ( sA[[i]]==1 && sB[[i]] ==0 ) s[4]<-s[4]+1
    else if ( sA[[i]]==1 && sB[[i]] ==1 ) s[5]<-s[5]+1
    else if ( sA[[i]]==1 && sB[[i]] ==2 ) s[6]<-s[6]+1
    else if ( sA[[i]]==2 && sB[[i]] ==0 ) s[7]<-s[7]+1
    else if ( sA[[i]]==2 && sB[[i]] ==1 ) s[8]<-s[8]+1
    else if ( sA[[i]]==2 && sB[[i]] ==2 ) s[9]<-s[9]+1
    }
s
}
  

  
f_DeltaAB<-function(s) #delta for unlinked loci as in LD2 of pegas
{
  nAhat <- 2 * ( s[7] + s[8] + s[9] ) + s[4] + s[5] + s[6]; #nA1
  nBhat <- 2 * ( s[3] + s[6] + s[9] ) + s[2] + s[5] + s[8]; #nB1
  Nse <- 2 * (s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9]);
  nABhat = (2 * s[9]) + s[8] + s[6] + (s[5] / 2); #nA1B1
  if ( Nse > 0 ) { DeltaAB <- 2* nABhat / Nse- 2 * #the first 2 is because at pag 126 of weir98 DeltaABhat is like this, it uses number of indiv. I don't get honestly why, but it is true, it is correct like this!
(nAhat * nBhat) / (Nse*Nse)} else {DeltaAB <- 0}
  DeltaAB
}

f_rhoAB<-function(s)
{
  nAhat <- 2 * ( s[7] + s[8] + s[9] ) + s[4] + s[5] + s[6]; #nA1
  nBhat <- 2 * ( s[3] + s[6] + s[9] ) + s[2] + s[5] + s[8]; #nB1
  Nse <- 2 * (s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9]);
  nABhat = (2 * s[9]) + s[8] + s[6] + (s[5] / 2); #nA1B1
  if ( Nse > 0 ) { DeltaAB <- 2* nABhat / Nse- 2 * #the first 2 is because at pag 126 of weir98 DeltaABhat is like this, it uses number of indiv. I don't get honestly why, but it is true, it is correct like this!
(nAhat * nBhat) / (Nse*Nse)} else {DeltaAB <- 0}
  return(c(DeltaAB,DeltaAB/sqrt(nAhat*(Nse-nAhat)*nBhat*(Nse-nBhat)/(Nse^4))))
}

f_rhoAB_nAB<-function(nA,nB,nAB,Nse) {
DeltaAB<-f_DeltaAB_nAB(nA,nB,nAB,Nse)
return(DeltaAB/sqrt(nA*nB*(Nse-nA)*(Nse-nB)*(1/(Nse*Nse*Nse*Nse))))
}

#FLAG_RAB_DIVIDEDBY2
#f_rhoAB(c(10,0,0,0,0,0,0,0,10)) #NB: I should divide by 2 in order to have proper r because D is given by gametic+genotypic. In fact 0.5*0.5-0.5=0.25 but in my case 0.5! #++
#In a way I think DAB as in weir is a bit unnatural, and more natural when it is divided by 2. -> I will have it like that.

frAB<-function(DeltaAB,nA,nB,Nse)
{
  RAB<-0;
  if ( Nse > 0 && nA>0 && nB>0 && nB<Nse && nA<Nse) { RAB<- DeltaAB/sqrt(nA*nB*(Nse-nA)*(Nse-nB)*(1/(Nse*Nse*Nse*Nse)));} else {RAB <- 0}
  return(RAB)
}

f_DeltaAB_nAB<-function(nA,nB,nAB,Nse)
{
  if ( Nse > 0 ) { 
    #DeltaAB <- 2* nAB / Nse- 2 * #the first 2 is because at pag 126 of weir98 DeltaABhat is like this, it uses number of indiv. I don't get honestly why, but it is true, it is correct like this!
    #But all values by 2..Actually I think it is wrong (17-09-15)
    DeltaAB <- nAB / Nse-(nA * nB) / (Nse*Nse)} else {DeltaAB <- 0}
  DeltaAB
}

HWD<-function(s,row)
{
  if (row==FALSE)
  {
  #pBB-pB*pB
  res<-(s[3]+s[6]+s[9])/sum(s)-((s[2]+s[5]+s[8])+2*(s[3]+s[6]+s[9]))/((2*sum(s))^2)
  } 
  else if (row==TRUE)
  {
  #pAA-pA*pA
  res<-sum(s[7:9])/sum(s)-(sum(s[4:6])+2*sum(s[7:9]))/((2*sum(s))^2)
  }
return(res)
}

f_rhoABWeir<-function(s)
{
  nAhat <- 2 * ( s[7] + s[8] + s[9] ) + s[4] + s[5] + s[6]; #nA1
  nBhat <- 2 * ( s[3] + s[6] + s[9] ) + s[2] + s[5] + s[8]; #nB1
  Nse <- 2 * (s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9]);
  nABhat = (2 * s[9]) + s[8] + s[6] + (s[5] / 2); #nA1B1
  if ( Nse > 0 ) { DeltaAB <- 2* nABhat / Nse- 2 *(nAhat * nBhat) / (Nse*Nse)} else {DeltaAB <- 0}
  #pA(1-pA)+pAA-pA^2=pAA+pA-pA*pA
  denA=2*( s[7] + s[8] + s[9] )/Nse+(nAhat/Nse)-2*(nAhat/Nse)^2
  denB=2*( s[3] + s[6] + s[9] )/Nse+(nBhat/Nse)-2*(nBhat/Nse)^2
  return(c(DeltaAB,(DeltaAB^2)/(denA*denB)))
}

a2A<-function(s)
{
  temp<-s[1:3]
  s[1:3]<-s[7:9]
  s[7:9]<-temp
  return(s)
}
b2B<-function(s)
{
  temp<-c(s[1],s[4],s[7])
  s[1]<-s[3];s[4]<-s[6];s[7]<-s[9]
  s[3]<-temp[1];s[6]<-temp[2];s[9]<-temp[3];
  return(s)
}

Rc2<-function(s)
{
f_rhoABWeir(s)[[2]]+f_rhoABWeir(a2A(s))[[2]]+f_rhoABWeir(b2B(s))[[2]]+f_rhoABWeir(a2A(b2B(s)))[[2]]
}

#f_rhoAB_nAB(13,13,0,20)
#f_rhoAB_nAB(13,13,0,20)
#f_rhoAB(c(1,1,1,1,20,1,1,1,1))
#f_rhoAB(c(10,0,0,0,0,0,0,0,10))#I get 2 for LD. This is wrong, it should be 1.

#check f_DeltaAB
#mean(sapply(rep(1,100000),function(x) HWE(10,10,100)))
#f_DeltaAB_nAB(32,32,32,100)
#f_DeltaAB_nAB(10,10,1,100)
#f_DeltaAB_nAB(11,10,1,100)
#f_DeltaAB_nAB(9,10,1,100)
#f_DeltaAB(c(1,1,1,1,20,1,1,1,1))
#f_DeltaAB(c(1,1,1,1,0,1,1,1,1))
#f_DeltaAB(c(1,1,1,1,0,1,1,1,0))
#f_DeltaAB(c(1,1,1,1,0,1,1,1,20))
#f_DeltaAB(rep(1,9))
#f_DeltaAB(c(rep(1,4),100,rep(1,4)))
#f_DeltaAB(1:9)
#f_DeltaAB(10:90)
#f_DeltaAB(c(72, 29, 5, 0, 1, 0, 0, 0, 0))
#f_DeltaAB(c(22,35,21,0,0, 1, 0, 0, 0, 0))
#mean(sapply(rep(1,100000),function(x) HWE(4,2,10)))
#HWE(2,4,10)
#sum(s==1)
#HWE(10,0,10)
#DeltaAB

#check of frAB
#double frAB(double DeltaAB , double nA, double nB, double Nse)
#  {
#  double RAB=0;
#  //printf("RAB: %f %f %f \n",nA,nB,Nse);
#  if ( Nse > 0 && nA>0 && nB>0 && nB<Nse && nA<Nse) { RAB= DeltaAB/sqrt(nA*nB*(Nse-nA)*(Nse-nB)*(1/(Nse*Nse*Nse*Nse)));} else {RAB = 0;};
#    return RAB;
#  }




#READING FILES
#fileName="/mnt/scratch/fabrizio/LDLD/analyses/YRI/chr11.tabchr10.tab.res"
#res <- readLines(fileName)
#res<-strsplit(res,split="\t")
#res2<-lapply(res,function(x) as.numeric(c(x[3],x[4],x[5],x[6],x[7],x[8])))
#sapply(res2[1:100], function(x) HWE(x[1],x[2],x[3]))
#sapply(res2[1:100], function(x) HWE(x[1],x[2],x[3]))
#sapply(res2[1:100], function(x) x[5])
#res[1:10][[1]][3]

#con=file(fileName,open="r")
#line=readLines(con) 
#long=length(line)
#for (i in 1:long){
#    linn=readLines(con,1)
#    print(linn)
#}
#close(con)

#GENERATE LOOKUP TABLE FOR p-values by permutations # -> way to slow
lookupbuild<-function(nA,nB,Nt,Nrep)
  {
  a<-sort(sapply(rep(1,Nrep),function(x) HWE(nA,nB,Nt)))
  s<-rep(0,1+min(nA,nB)) #nAB is at max max(nA,nB)
  for (i in 0:min(nA,nB))
    {
    if (f_DeltaAB_nAB(nA,nB,i,Nt) != 0) s[i+1]<-(sum(a<= -abs(f_DeltaAB_nAB(nA,nB,i,Nt))) + sum(a >= abs(f_DeltaAB_nAB(nA,nB,i,Nt))))/Nrep #pvalue
    else s[i+1]<-1
    }
#    print(a)
#  print(a> abs(f_DeltaAB_nAB(nA,nB,1,Nt)))
  s
  }
#check
#system.time(a<-lookupbuild(10,10,100,1000000))
#system.time(lookupbuild(10,10,200,100000))
#lookupbuild(10,10,100,10000000)
#lookupbuild(10,10,200,100000)
#f_DeltaAB_nAB(10,10,0,200)
#sort.list(sapply(rep(1,1000),function(x) HWE(4,2,10)))
#sort(sapply(rep(1,1000),function(x) HWE(4,2,10)))
#max(1,10)
#sort(sapply(rep(1,1000),function(x) HWE(10,10,100)))

#==============================================================
#=======  P-values for linkage - permutations and fisher #=====
#==============================================================

#==============chisquare test as described in Weir pag. 127=====
#definition of piAhat and DAhat given at pages 119 and 95, respectively.
#chi2AB=n DeltahatAB^2/((piA+DAhat)*(piB+DBhat))

DAhat<-function(pA,pAA) pAA-pA^2
#f_DeltaAB(c(72, 29, 5, 0, 1, 0, 0, 0, 0))
#f_DeltaAB(c(90, 20, 0, 20, 0, 0, 0, 0, 0))
#f_DeltaAB(c(90, 20, 0, 0, 25, 0, 0, 0, 2))
piAhat<-function(pA) pA*(1-pA)
chi2AB<-function(s)
  {
  nAhat <- 2 * ( s[7] + s[8] + s[9] ) + s[4] + s[5] + s[6]; #nA1
  nBhat <- 2 * ( s[3] + s[6] + s[9] ) + s[2] + s[5] + s[8]; #nB1
  Nse <- 2 * (s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9]);
  nABhat = (2 * s[9]) + s[8] + s[6] + (s[5] / 2); #nA1B1
  (Nse/2)*f_DeltaAB_nAB(nAhat,nBhat,nABhat,Nse)^2/((piAhat(nAhat/Nse)+DAhat(nAhat/Nse, 2*( s[7] + s[8] + s[9] )/Nse))*(piAhat(nBhat/Nse)+DAhat(nBhat/Nse, 2*( s[3] + s[6] + s[9] )/Nse)))
  }

#chi2AB(c(72, 29, 5, 0, 1, 0, 0, 0, 0))
#chi2AB(c(90, 20, 0, 20, 0, 0, 0, 0, 0))
#pchisq(chi2AB(c(90, 20, 0, 20, 0, 0, 0, 0, 0)),4)
#chi2AB(c(90, 20, 0, 20, 30, 0, 0, 0, 0))
#pchisq(chi2AB(c(90, 20, 0, 20, 30, 0, 0, 0, 0)),4)
#chi2AB(c(40, 20, 0, 10, 30, 0, 0, 0, 0))
#pchisq(chi2AB(c(40, 20, 0, 10, 50, 0, 0, 0, 0)),4)
#pchisq(chi2AB(c(40, 20, 0, 10, 10, 0, 0, 0, 10)),4)

#==============Fisher's exact test===================
Weir98Fisher<-function(gts) #with big numbers it returns errors because factorial are really huge
  {
  #numerator
  naafact<-factorial(gts[1]+gts[2]+gts[3])
  nAafact<-factorial(gts[4]+gts[5]+gts[6])
  nAAfact<-factorial(gts[7]+gts[8]+gts[9])
  nbbfact<-factorial(gts[1]+gts[4]+gts[7])
  nBbfact<-factorial(gts[2]+gts[5]+gts[8])
  nBBfact<-factorial(gts[3]+gts[6]+gts[9])
  num<-nBBfact*nBbfact*nbbfact*nAAfact*nAafact*naafact
  #denominator
  print(num)
  print(log(nBBfact,10)+log(nBbfact,10)*log(nbbfact,10)*log(nAAfact,10)*log(nAafact,10)*log(naafact,10))
  print(10^log(num,10))
  den<-1
  for (i in 1:9)
    {
    den<-factorial(gts[i])*den
    }
  (num/(den*factorial(gts[1]+gts[2]+gts[3]+gts[4]+gts[5]+gts[6]+gts[7]+gts[8]+gts[9])))
  }

Weir98Fisher_log<-function(gts) #turnaround to avoid arbitrary precision issues by using logarithm transform
  {
  #numerator
  naafact<-factorial(gts[1]+gts[2]+gts[3])
  nAafact<-factorial(gts[4]+gts[5]+gts[6])
  nAAfact<-factorial(gts[7]+gts[8]+gts[9])
  nbbfact<-factorial(gts[1]+gts[4]+gts[7])
  nBbfact<-factorial(gts[2]+gts[5]+gts[8])
  nBBfact<-factorial(gts[3]+gts[6]+gts[9])
  num<-log(nBBfact,10)+log(nBbfact,10)+log(nbbfact,10)+log(nAAfact,10)+log(nAafact,10)+log(naafact,10)
  #denominator
  den<-1
  for (i in 1:9)
    {
    den<-factorial(gts[i])*den
    }
  10^(num-log(den,10)-log(factorial(gts[1]+gts[2]+gts[3]+gts[4]+gts[5]+gts[6]+gts[7]+gts[8]+gts[9]),10))
  }
#check
#f_DeltaAB(c(72, 29, 5, 0, 1, 0, 0, 0, 0))
#f_DeltaAB(c(22,35,21,0,0, 1, 0, 0, 0, 0))
#Weir98Fisher(c(72, 29, 5, 0, 1, 0, 0, 0, 0))
#Weir98Fisher(c(22,35,21,0,0, 1, 0, 0, 0, 0))
#Weir98Fisher_log(c(22,35,21,0,0, 1, 0, 0, 0, 0))
#Weir98Fisher(c(72, 0, 0, 29, 1, 0, 5, 0, 0))
#Weir98Fisher_log(c(72, 0, 0, 29, 1, 0, 5, 0, 0))
#Weir98Fisher_log(c(0, 0, 0, 29, 1, 0, 5, 0, 10))
nABhat_f<-function(s) (2 * s[9]) + s[8] + s[6] + (s[5] / 2);
nAhat_f<-function(s) 2 * ( s[7] + s[8] + s[9] ) + s[4] + s[5] + s[6]
nBhat_f<-function(s) 2 * ( s[3] + s[6] + s[9] ) + s[2] + s[5] + s[8]
Nse_f <-function(s) 2 * (s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9])

#generator of full permutations with all possible genotypes on the basis of allele freq: what one needs for HWE - to be checked better
#clearly not one needs for LDE to use Daub/Weir's formula since it always sum up to more than 1!
if (FALSE){
exactFishtest<-function(s)
  {
  nAhat <- nAhat_f(s)
  nBhat <- nBhat_f(s)
  Nse <- Nse_f(s)
  gtst<-rep(0,9)
  res<-c()
  for (row3 in 0:(nAhat%/%2)) #diploid A
    {
    for (col3 in 0:(nBhat%/%2)) #diploid B
      {
      for (gt9 in 0:min(row3,col3)) #both diploids
	{
	for (gt5 in 0:min(nAhat-2*row3,nBhat-2*col3,Nse/2-gt9)) #both heteroz
	  {
	  for (gt8 in 0:min(row3-gt9,nBhat-2*col3-gt5,Nse/2-gt9-gt5)) #A dipl, Bh
	    {
	    for (gt6 in 0:min(col3-gt9,nAhat-2*row3-gt5,Nse/2-gt9-gt5-gt8)) #A h, Bh
	      {
	      gtst[2]<-nBhat-2*col3-gt5-gt8
	      gtst[3]<-col3-gt6-gt9
	      gtst[4]<-nAhat-2*row3-gt5-gt6
	      gtst[5]<-gt5
	      gtst[6]<-gt6
	      gtst[7]<-row3-gt8-gt9
	      gtst[8]<-gt8
	      gtst[9]<-gt9
	      gtst[1]<-Nse/2-sum(gtst[2:9])
	    if (sum(gtst<0)==0)
	      {
	      #check that assortment consistent
	      if (nAhat_f(gtst)==nAhat && nBhat_f(gtst)==nBhat && Nse_f(gtst)==Nse)
		{
	      print(gtst)
	      res<-c(res,Weir98Fisher_log(gtst))
	     }}}   } 
	  }
	}
      }
    }
  res
 }
}
#check
#exactFishtest(c(1, 0, 0, 0, 0, 0, 0, 0, 1))  
#exactFishtest(c(0, 0, 0, 0, 0, 0, 0, 0, 10))
#exactFishtest(c(10, 0, 0, 0, 0, 0, 0, 0, 0))
#exactFishtest(c(9, 1, 0, 0, 0, 0, 0, 0, 0))
#exactFishtest(c(9, 2, 0, 0, 0, 0, 0, 0, 0))
#exactFishtest(c(10, 0, 0, 0, 0, 0, 0, 0, 1))
#exactFishtest(c(2, 0, 0, 0, 0, 0, 0, 0, 2))
#system.time(a<-exactFishtest(c(90, 0, 0, 0, 0, 0, 0, 0, 10)))
#system.time(a<-exactFishtest(c(140, 0, 2, 8, 1, 10, 2, 0, 1)))
#system.time(a<-exactFishtest(c(72, 29, 5, 0, 1, 0, 0, 0, 0)))
#sum(a)

#generator of arrangements for LDE: maybe it is shuffling keeping fixed the genotypes (amount of diploids AA and BB)
#in fact it always sums up to 1
#install.packages("Rmpfr")
#install(Rmpfr)
library("Rmpfr")
#library("/tmp/RtmpQcW13M/downloaded_packages/Rmpfr")

#gamma(as(1000,"mpfr"))
exactFishtest<-function(s)
  {
  nAhat <- nAhat_f(s)
  nBhat <- nBhat_f(s)
  Nse <- Nse_f(s)
  gtst<-rep(0,9)
  res<-c()
  row3<-s[7]+s[8]+s[9]
  col3<-s[3]+s[6]+s[9]
      for (gt9 in 0:min(row3,col3)) #both diploids
	{
	for (gt5 in 0:min(nAhat-2*row3,nBhat-2*col3,Nse/2-gt9)) #both heteroz
	  {
	  for (gt8 in 0:min(row3-gt9,nBhat-2*col3-gt5,Nse/2-gt9-gt5)) #A dipl, Bh
	    {
	    for (gt6 in 0:min(col3-gt9,nAhat-2*row3-gt5,Nse/2-gt9-gt5-gt8)) #A h, Bh
	      {
	      gtst[2]<-nBhat-2*col3-gt5-gt8
	      gtst[3]<-col3-gt6-gt9
	      gtst[4]<-nAhat-2*row3-gt5-gt6
	      gtst[5]<-gt5
	      gtst[6]<-gt6
	      gtst[7]<-row3-gt8-gt9
	      gtst[8]<-gt8
	      gtst[9]<-gt9
	      gtst[1]<-Nse/2-sum(gtst[2:9])
	      }
	    if (sum(gtst<0)==0)
	      {
	      #check that assortment consistent
	      if (nAhat_f(gtst)==nAhat && nBhat_f(gtst)==nBhat && Nse_f(gtst)==Nse)
		{
	      print(gtst)
	      res<-c(res,Weir98Fisher_log(gtst))
		}
	      }
	    }
	  }
	}
  print(sum(res))	  
  res
 }
  
#exactFishtest(c(1, 0, 0, 0, 0, 0, 0, 0, 1))  
#exactFishtest(c(0, 0, 0, 0, 0, 0, 0, 0, 10))
#exactFishtest(c(10, 0, 0, 0, 0, 0, 0, 0, 0))
#exactFishtest(c(9, 1, 0, 0, 0, 0, 0, 0, 0))
#exactFishtest(c(9, 2, 0, 0, 0, 0, 0, 0, 0))
#exactFishtest(c(10, 0, 0, 0, 0, 0, 0, 0, 1))
#exactFishtest(c(2, 0, 0, 0, 0, 0, 0, 0, 2))
#all good, sum to 1
#system.time(exactFishtest(c(90, 0, 0, 0, 0, 0, 0, 0, 10)))
#system.time(exactFishtest(c(70, 0, 0, 0, 0, 0, 0, 0, 10)))
#system.time(exactFishtest(c(70, 5, 0, 0, 0, 0, 0, 0, 10)))
#system.time(exactFishtest(c(70, 5, 5, 0, 0, 0, 0, 0, 10)))
#system.time(exactFishtest(c(70, 0, 0, 5, 0, 0, 0, 0, 0)))
#system.time(exactFishtest(c(70, 0, 0, 5, 5, 0, 0, 0, 0)))
#system.time(exactFishtest(c(70, 0, 0, 5, 5, 0, 5, 0, 0)))
#when other elements starts being wrong
#system.time(exactFishtest(c(70, 5, 5, 5, 0, 0, 0, 0, 10)))
#system.time(exactFishtest(c(70, 0, 0, 5, 0, 0, 0, 0, 10)))
#system.time(exactFishtest(c(70, 0, 0, 5, 5, 5, 0, 5, 0)))


#let's make it simpler
#sum1
#exactFishtest(c(1, 0, 0, 0, 0, 0, 0, 0, 1))
#exactFishtest(c(1, 1, 0, 0, 0, 0, 0, 0, 1))
#exactFishtest(c(1, 1, 1, 0, 0, 0, 0, 0, 1))
#exactFishtest(c(1, 1, 1, 0, 0, 0, 0, 1, 1))
#exactFishtest(c(1, 1, 1, 0, 0, 0, 1, 1, 1))
#exactFishtest(c(1, 0, 0, 0, 1, 0, 0, 0, 0))

#sum<1
#exactFishtest(c(1, 1, 1, 0, 1, 0, 1, 1, 1))
#exactFishtest(c(1, 1, 1, 1, 0, 0, 1, 1, 1))
#exactFishtest(c(1, 1, 1, 0, 0, 1, 1, 1, 1))
#exactFishtest(c(1, 1, 1, 0, 0, 1, 0, 0, 0))
#exactFishtest(c(1, 1, 0, 0, 0, 1, 0, 0, 0))
#exactFishtest(c(1, 0, 0, 0, 0, 1, 0, 0, 0))

exactDtest<-function(s) #corrected, works all the times
  {
  myDeltam<--abs(f_DeltaAB(s))
  myDeltap<-abs(f_DeltaAB(s))
  nAhat <- nAhat_f(s)
  nBhat <- nBhat_f(s)
  Nse <- Nse_f(s)
  gtst<-rep(0,9)
  res<-c()
  pvalue<-0
  row3<-s[7]+s[8]+s[9]
  col3<-s[3]+s[6]+s[9]
  row2<-s[4]+s[5]+s[6]
  col2<-s[2]+s[5]+s[8]
  if (myDeltap==0) {pvalue<-1}
  else
    {
    for (gt9 in 0:min(row3,col3)) #both diploids
      {
      for (gt5 in 0:min(nAhat-2*row3,nBhat-2*col3,Nse/2-gt9)) #both heteroz
	{
	for (gt8 in 0:min(row3-gt9,col2-gt5,Nse/2-gt9-gt5)) #A dipl, Bh
	  {
	  for (gt6 in 0:min(col3-gt9,row2-gt5,Nse/2-gt9-gt5-gt8)) #A h, Bh
	    {
	    gtst[2]<-col2-gt5-gt8
	    gtst[3]<-col3-gt6-gt9
	    gtst[4]<-row2-gt5-gt6
	    gtst[5]<-gt5
	    gtst[6]<-gt6
	    gtst[7]<-row3-gt8-gt9
	    gtst[8]<-gt8
	    gtst[9]<-gt9
	    gtst[1]<-Nse/2-sum(gtst[2:9])
	    if (sum(gtst<0)==0)
	      {
#	      print(gtst)
	      res<-c(res,Weir98Fisher_log(gtst))
	      if (f_DeltaAB(gtst)<= myDeltam || f_DeltaAB(gtst)>= myDeltap)
		{
		pvalue<-pvalue+Weir98Fisher_log(gtst)
		}
	      }
	    }
	  }
	}
    }}
  print(sum(res))	  
  pvalue
 }
#a<-c(2,43,46,5,2,6,1,2,0)
#nAhat_f(a)
#nBhat_f(a)
#Nse_f(a)
#f_DeltaAB(a)
#exactDtest(c(2,43,46,5,2,6,1,2,0))

#a<-c(72,19,0,12,4,0,0,0,0);exactDtest(a)

#let's make it simpler
#sum1
#exactDtest(c(1, 0, 0, 0, 0, 0, 0, 0, 1))
#exactDtest(c(1, 1, 0, 0, 0, 0, 0, 0, 1))
#exactDtest(c(1, 1, 1, 0, 0, 0, 0, 0, 1))
#exactDtest(c(1, 1, 1, 0, 0, 0, 0, 1, 1))
#exactDtest(c(1, 1, 1, 0, 0, 0, 1, 1, 1))
#exactDtest(c(1, 0, 0, 0, 1, 0, 0, 0, 0))

#sum<1
#exactDtest(c(1, 1, 1, 0, 1, 0, 1, 1, 1))
#exactDtest(c(1, 1, 1, 1, 0, 0, 1, 1, 1))
#exactDtest(c(1, 1, 1, 0, 0, 1, 1, 1, 1))
#exactDtest(c(1, 1, 1, 0, 0, 1, 0, 0, 0))
#exactDtest(c(1, 1, 0, 0, 0, 1, 0, 0, 0))
#exactDtest(c(1, 0, 0, 0, 0, 1, 0, 0, 0))

#system.time(a<-exactDtest(c(140, 0, 2, 8, 1, 10, 2, 0, 1)))
#system.time(a<-exactDtest(c(72, 29, 5, 0, 1, 0, 0, 0, 0)))

#exactDtest(c(140, 0, 2, 8, 1, 10, 2, 0, 1))
#exactDtest(c(2, 2, 2, 2, 2, 2, 2, 2, 2))
#exactDtest(c(2, 2, 2, 2, 2, 2, 2, 2, 5))
#exactDtest(c(0, 2, 2, 2, 2, 2, 2, 2, 5))
#exactDtest(c(0, 2, 2, 2, 2, 2, 2, 2, 10))
#exactDtest(c(1, 2, 2, 2, 2, 2, 2, 2, 2))
#exactDtest(c(0, 2, 2, 2, 2, 2, 2, 2, 3))
#exactDtest(c(0, 2, 2, 2, 2, 2, 2, 2, 4))

#exactDtest(c(0, 2, 2, 2, 2, 2, 2, 2, 2))
#exactDtest(c(0, 2, 2, 2, 2, 2, 2, 2, 2))
#exactDtest(c(200, 2, 2, 2, 2, 2, 2, 2, 2))
#exactDtest(c(0, 2, 2, 2, 2, 2, 2, 2, 2))

#now I test if I can use D for p-value like kulinskaya2009
pexactpermrho<-function(s) #
  {
  myDeltam<--abs(f_rhoAB(s))
  myDeltap<-abs(f_rhoAB(s))
  nAhat <- nAhat_f(s)
  nBhat <- nBhat_f(s)
  Nse <- Nse_f(s)
  gtst<-rep(0,9)
  res<-c()
  pvalue<-0
  row3<-s[7]+s[8]+s[9]
  col3<-s[3]+s[6]+s[9]
  row2<-s[4]+s[5]+s[6]
  col2<-s[2]+s[5]+s[8]
  if (myDeltap==0) {pvalue<-1}
  else
    {
    for (gt9 in 0:min(row3,col3)) #both diploids
      {
      for (gt5 in 0:min(nAhat-2*row3,nBhat-2*col3,Nse/2-gt9)) #both heteroz
	{
	for (gt8 in 0:min(row3-gt9,col2-gt5,Nse/2-gt9-gt5)) #A dipl, Bh
	  {
	  for (gt6 in 0:min(col3-gt9,row2-gt5,Nse/2-gt9-gt5-gt8)) #A h, Bh
	    {
	    gtst[2]<-col2-gt5-gt8
	    gtst[3]<-col3-gt6-gt9
	    gtst[4]<-row2-gt5-gt6
	    gtst[5]<-gt5
	    gtst[6]<-gt6
	    gtst[7]<-row3-gt8-gt9
	    gtst[8]<-gt8
	    gtst[9]<-gt9
	    gtst[1]<-Nse/2-sum(gtst[2:9])
	    if (sum(gtst<0)==0)
	      {
#	      print(gtst)
	      res<-c(res,Weir98Fisher_log(gtst))
	      if (f_rhoAB(gtst)<= myDeltam || f_rhoAB(gtst)>= myDeltap)
		{
		pvalue<-pvalue+Weir98Fisher_log(gtst)
		}
	      }
	    }
	  }
	}
 }   }
  print(sum(res))	  
  pvalue
 }
#exactDtest(c(140, 0, 2, 8, 1, 10, 2, 0, 1))
#pexactpermrho(c(140, 0, 2, 8, 1, 10, 2, 0, 1))
#exactDtest(c(72, 29, 5, 0, 1, 0, 0, 0, 0))
#pexactpermrho(c(72, 29, 5, 0, 1, 0, 0, 0, 0))
#exactDtest(c(22, 29, 5, 13, 4, 7, 2, 4, 6))
#pexactpermrho(c(22, 29, 5, 13, 4, 7, 2, 4, 6))


#pexactpermrho(c(0, 0, 0, 0, 9, 10, 2, 0, 1))








#NB: problem with big factorials!
#install package rmpfr
#http://stackoverflow.com/questions/22668381/error-while-installing-rmpfr-on-ubuntu
#require("Rmpfr")
#getGroupMembers("Math")
#factorialMpfr(200)

Weir98Fisher_log<-function(gts) #turnaround to avoid arbitrary precision issues by using logarithm transform
  {
  #numerator
  naafact<-factorialMpfr(gts[1]+gts[2]+gts[3])
  nAafact<-factorial(gts[4]+gts[5]+gts[6])
  nAAfact<-factorial(gts[7]+gts[8]+gts[9])
  nbbfact<-factorialMpfr(gts[1]+gts[4]+gts[7])
  nBbfact<-factorial(gts[2]+gts[5]+gts[8])
  nBBfact<-factorial(gts[3]+gts[6]+gts[9])
  num<-log(nBBfact,10)+log(nBbfact,10)+log10(nbbfact)+log(nAAfact,10)+log(nAafact,10)+log10(naafact)
  #denominator
  den<-1
  for (i in 1:9)
    {
    den<-factorialMpfr(gts[i])*den
    }
  10^(num-log10(den)-log10(factorialMpfr(gts[1]+gts[2]+gts[3]+gts[4]+gts[5]+gts[6]+gts[7]+gts[8]+gts[9])))
  }
  
 #a<-factorialMpfr(250)
 #a*200 
  
  #Weir98Fisher_log(c(200, 2, 2, 2, 2, 2, 2, 2, 2))
  
  #log10(a)
  

#generate look up table for chisquare

#a<-seq(0.001,1,0.001)
#b<-seq(1.02,50,0.02)
#aa<-sapply(a,function(x) 1-pchisq(x, df=1))
#bb<-sapply(b,function(x) 1-pchisq(x, df=1))
#write(paste(aa,collapse=", "),"lookup1")
#write(paste(bb,collapse=", "),"lookup1")


#============ check that genotypes parsed correctly =================================
#I do
#./go.out YRI/missense/file1 YRI/missense/file2 file3 file4 file5 file6 0.05 > temp
#cat temp | grep iline | sed s/": "/"."/ |sed s/"\t"/"<-c("/ | sed 's/\t/,/g' | sed 's/[,]$/)/g' > temp2
#cat chr21.tab | grep 15011942 | awk '{for (i = 10; i <= NF; i++) {if (i<NF) {print $i "\t"} else {print $i "\t"}}}' | sed ':a;N;$!ba;s/\n//g' | sed 's/\t/,,/g' | sed 's/0.0/0/g' | sed 's/0.1/2/g' | sed 's/1.0/2/g' | sed 's/1.1/1/g' | sed 's/,,/,/g' | sed 's/^,//g' | sed 's/,$//g'
#ok 1007432
#ok 1686040
#ok 3759784
#12009956
ilines1.0<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
ilines1.1<-c(0,0,0,2,0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,2,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0)
ilines1.2<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0)
ilines1.3<-c(0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,2,0,2,0,0,0,0,2,0,2,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,0,2,0,0,0,0,0,2,0,2,0,2,0,0,2,0,2,0,0)
ilines1.4<-c(1,2,1,1,2,1,1,1,1,2,1,1,1,1,1,1,1,1,2,1,1,1,1,2,2,2,1,1,1,1,1,2,1,2,1,1,1,1,1,1,1,1,1,2,2,2,1,1,2,2,1,1,1,1,0,1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,2,1,1,1,2,2,1,1,1,1,1,1,2,1,1,1,1,2,1,2,1,1,1,1,2,1,1,1,1,1,2)
ilines2.0<-c(2,1,1,1,1,2,1,1,1,1,1,1,1,1,2,1,1,1,1,2,2,2,1,1,1,1,1,2,1,2,1,1,1,1,1,1,1,1,1,2,2,2,1,1,2,2,1,1,1,1,0,1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,2,1,1,1,2,2,1,1,1,1,1,1,2,1,1,1,1,2,1,2,1,1,1,1,2,1,1,1,1,1,2,0,0,0,0)
ilines2.1<-c(0,0,0,2,0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,2,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,2,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0)
ilines2.2<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0)
ilines2.3<-c(0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,2,0,2,0,0,0,0,2,0,2,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,0,2,0,0,0,0,0,2,0,2,0,2,0,0,2,0,2,0,0)
ICLD15chrA<-c(0,0,0,0,0,0,2,0,0,0,0,2,2,0,2,0,1,2,0,0,1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,2,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0)
ICLD15chrB<-c(2,2,2,1,0,0,1,2,0,0,2,2,1,2,2,2,1,2,0,2,1,2,1,2,0,0,2,2,2,2,0,0,2,2,0,0,2,2,2,1,0,0,2,2,1,0,2,0,2,2,2,0,1,0,2,2,2,0)


buildtable_begin<-function(ilines1,ilines2){
  gts<-c(0,0,0,0,0,0,0,0,0)
  parray3<-ilines1
  for ( ibt in 1:length(ilines1) )
  {
    if ( ilines1[ibt] == 0 && ilines2[ibt] == 0 )
    {
      parray3[ibt] = 0
      gts[1]=gts[1]+1
    }
    else if ( ilines1[ibt] == 1 && ilines2[ibt] == 1 )
    {
      parray3[ibt] = 9
      gts[9]=gts[9]+1
    }
    else if ( ilines1[ibt] == 2 && ilines2[ibt] == 2 )
    {
      parray3[ibt] = 5
      gts[5]=gts[5]+1
    }
    else if ( ilines1[ibt] == 3 || ilines2[ibt] == 3 )
    {
      parray3[ibt] = 30 
    }
    else if ( ilines1[ibt] == 4 || ilines2[ibt] == 4 )
    {
      parray3[ibt] = 40 
    }
    else if ( ilines1[ibt] == 9 || ilines2[ibt] == 9 )
    {
      parray3[ibt] = 90 
    }
    else if ( ilines1[ibt] == 0 && ilines2[ibt] == 2 )
    {
      parray3[ibt] = 2 
      gts[2]=gts[2]+1
    }
    else if ( ilines1[ibt] == 0 && ilines2[ibt] == 1 )
    {
      parray3[ibt] = 3
      gts[3]=gts[3]+1
    }
    else if ( ilines1[ibt] == 2 && ilines2[ibt] == 0 )
    {
      parray3[ibt] = 4 
      gts[4]=gts[4]+1
    }
    else if ( ilines1[ibt] == 2 && ilines2[ibt] == 1 )
    {
      parray3[ibt] = 6
      gts[6]=gts[6]+1
    }
    else if ( ilines1[ibt] == 1 && ilines2[ibt] == 0 )
    {
      parray3[ibt] = 7 
      gts[7]=gts[7]+1
    }
    else if ( ilines1[ibt] == 1 && ilines2[ibt] == 2 )
    {
      parray3[ibt] = 8
      gts[8]=gts[8]+1
    }
  }
return(gts)
}
#buildtable_begin(ilines1.1,ilines2.1)
#buildtable_begin(ilines1.1,ilines2.2)
#buildtable_begin(ilines1.2,ilines2.1)
#buildtable_begin(ilines1.2,ilines2.2)
#f_rhoAB(buildtable_begin(ilines1.1,ilines2.1))

#buildtable_begin(ilines1.3,ilines2.3) #ok, genotypes parsed ok

#buildtable_begin(ICLD15chrA,ICLD15chrB)

#==================================================================
#============ simulating haplotype linkage ========================
#==================================================================


generateGenome<-function( freq, N) #NB, here with Genome I can think both of a single indiv or a single base for a whole pop!
  {
  pop1<-rbinom(2*N,1,freq)
  pop1<-rbind(pop1[1:N],pop1[(N+1):(2*N)])
  return(pop1)
  }	
compressGenome<-function(genome)
  {
  cgen<-sapply(1:length(genome[1,]), function(x) {
  if (genome[1,x]==0 && genome[2,x]==0) {return(0)}
  else if (genome[1,x]==1 && genome[2,x]==1) {return(1.1)}
  else if (genome[1,x]==0 && genome[2,x]==1) {return(0.1)}
  else if (genome[1,x]==1 && genome[2,x]==0) {return(1.0)}
  })
  return(cgen)
  }
mutatepos<-function(state,mu) {if (rbinom(1,1,mu)==1) {return(abs(state-1))} else {return(state)}}
mutatepos_c<-function(state,mu) {if (rbinom(1,1,mu)==1) {return(sample(c(0,1.1,0.1,1.0),1))} else {return(state)}}

mutate<-function(genome,mu) #I can see this both as mutating a single individual genome or a whole pop at a single pos
  {
  for (i in 1:length(genome))
    {genome[i]<-mutatepos_c(genome[i],mu)}
  return(genome)
  }
genome<-compressGenome(generateGenome( 0.3, 20))
genome
mutate(genome,0.1)

init_DNAstretch<-function(N, freq, mu , nsnps)
  {
  d<-compressGenome(generateGenome( freq, N))
  a<-d
  for (i in 1:(nsnps-1))
    {
    d<-mutate(d,mu)
    a<-rbind(a,d)
    }
  rownames(a)<-1:nsnps
  return(t(a))
  }
#pop<-init_DNAstretch(20,0.3,0.1,7)

#random population
#t(sapply(1:20,function(x)compressGenome(generateGenome( 0.3, 10))))

#
#generateGenome( 0.3, 10)
#pop[,1][pop[,1]==1.1]


#round(3/2)






