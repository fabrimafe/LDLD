#define functions
{
#genotypes are in order AB, Ab, aB and ab
#full simulations
{
nex.gen.r<-function(mygenotypes,rho){ 
irec<-sample(1:Ne,rbinom(1,Ne,rho)) #recombining pairs
genotype_mat<-cbind(mygenotypes[1:Ne],mygenotypes[(Ne+1):(2*Ne)])
genotype_mat[irec,]<-t(apply(genotype_mat[irec,],MARGIN=1,FUN=function(x) if (x[1]==1 && x[2]==4){return(c(2,3))} else if (x[1]==2 && x[2]==3){return(c(1,4))} else return(x) ))
c(genotype_mat)
}
nex.gen.s<-function(mygenotypes,s_vec){ 
sample(1:4,2*Ne,table_f(mygenotypes)*s_vec,replace=TRUE)
}
table_f<-function(vec){c(sum(vec==1),sum(vec==2),sum(vec==3),sum(vec==4))}
simf<-function(genotypes_start,s_vec,ngenerations=10000,Ne=Ne){
genotypes<-genotypes_start
mylog<-t(sapply(1:ngenerations, function(x)(rep(0,length(s_vec)+1))))
for (igen in 1:ngenerations)
    {
    genotypes<-nex.gen.s(nex.gen.r(genotypes,0.5),s_vec)
    mylog[igen,1:4]<- table_f(genotypes)*0.5/Ne
    mylog[igen,5]<-mylog[igen,1]-(mylog[igen,1]+mylog[igen,2])*(mylog[igen,1]+mylog[igen,3])
    }
return(mylog)
}
}
#count simulations
{
mygenotypes<-c(5000,5000,6000,4000)
nex.gen.r<-function(mygenotypes,rho){ 
rec_gen1<-rbinom(1,mygenotypes[1],rho) #number of AB recombining
rec_gen1<-rhyper(1,mygenotypes[4],sum(mygenotypes[1:3]),rec_gen1) #number of AB recombining with ab
rec_gen2<-rbinom(1,mygenotypes[2],rho) #number of AB recombining
rec_gen2<-rhyper(1,mygenotypes[3],sum(mygenotypes[c(1,2,4)]),rec_gen2) #number of AB recombining with ab
mygenotypes+c(-rec_gen1+rec_gen2,rec_gen1-rec_gen2,rec_gen1-rec_gen2,-rec_gen1+rec_gen2)
}
nex.gen.s<-function(mygenotypes,s_vec){ 
rmultinom(1,sum(mygenotypes),prob=mygenotypes*s_vec)
}
simf<-function(mygenotypes,s_vec,ngenerations=10000){
Ne<-sum(mygenotypes)
mylog<-t(sapply(1:ngenerations, function(x)(rep(0,5))))
for (igen in 1:ngenerations)
    {
    mygenotypes<-c(nex.gen.s(nex.gen.r(mygenotypes,0.5),s_vec))
    mylog[igen,1:4]<-mygenotypes 
    mylog[igen,5]<-mygenotypes[1]-(mygenotypes[1]+mygenotypes[2])*(mygenotypes[1]+mygenotypes[3])/Ne
    }
return(mylog/Ne)
}

Ne<-10000
pA<-0.05
pB<-0.05
freqAB<-pA*pB
freqAb<-pA*(1-pB)
freqaB<-(1-pA)*pB
freqab<-(1-pA)*(1-pB)
mygenotypes<-c(rmultinom(1,2*Ne,prob=c(freqAB,freqAb,freqaB,freqab)))
nsims<-1000
s_vec<-c(1.1,1,1,1)
mysims_vstrong<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1,1-0.1,1-0.1,1)
mysims_vstrong_balancing<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1.01,1,1,1)
mysims_strong<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1,1-0.01,1-0.01,1)
mysims_strong_balancing<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1.0004,1,1,1)
mysims_weak<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1,1-0.0004,1-0.0004,1)
mysims_weak_balancing<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
#making 3D array from list of matrices
Ne<-10000
pA<-0.2
pB<-0.2
freqAB<-pA*pB
freqAb<-pA*(1-pB)
freqaB<-(1-pA)*pB
freqab<-(1-pA)*(1-pB)
mygenotypes<-c(rmultinom(1,2*Ne,prob=c(freqAB,freqAb,freqaB,freqab)))
nsims<-1000
s_vec<-c(1.1,1,1,1)
mysims_vstrong20<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1,1-0.1,1-0.1,1)
mysims_vstrong_balancing20<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1.01,1,1,1)
mysims_strong20<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1,1-0.01,1-0.01,1)
mysims_strong_balancing20<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1.0004,1,1,1)
mysims_weak20<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
s_vec<-c(1,1-0.0004,1-0.0004,1)
mysims_weak_balancing20<-lapply(1:nsims,function(x) simf(mygenotypes,s_vec,ngenerations=10000))
save.image('~/Dropbox/LDLD/figs/sims.RData')



arr <- array( unlist(mysims_weak_balancing20) , c(10000,5,nsims) )
myres<-apply( arr , 1:2 , mean ) #rowMeans( arr , dims = 2 ) #does mean over the first two dimensions
mysd<-apply( arr , 1:2 , sd )


arr <- array( unlist(mysims_strong) , c(10000,5,nsims) )
myres<-apply( arr , 1:2 , mean ) #rowMeans( arr , dims = 2 ) #does mean over the first two dimensions
mysd<-apply( arr , 1:2 , sd )

load('~/Dropbox/LDLD/figs/sims.RData')
for ( mycondition in c('mysims_weak_balancing20','mysims_strong_balancing20','mysims_weak20','mysims_strong20'))
{
arr <- array( unlist(get(mycondition)) , c(10000,5,nsims) )
myres<-apply( arr , 1:2 , mean ) #rowMeans( arr , dims = 2 ) #does mean over the first two dimensions
mylow<-apply( arr , 1:2 , function(x) quantile(x,0.05) )
myhigh<-apply( arr , 1:2 , function(x) quantile(x,0.95) )
png(paste0('~/Dropbox/LDLD/figs/',mycondition,'.png'))
par(mar=c(5,4,4,5)+.1)
plot(myres[,1],ylim=c(0,1),type='n',xlab="generation",ylab='allele frequency')
mtext("linkage disequilibrium (D)", side=4, line=3)
abline(h=0.5,lwd=1,lty=2)
#axis(4,at=c(0,0.25,0.5,0.75,1),labels=c(expression(-10^'-3'),expression(paste(-5,'*',10^'-3')),0,expression(paste(5,'*',10^'-3')),expression(10^'-3')))
axis(4,at=c(0,0.25,0.5,0.75,1),labels=c(expression(-0.1),expression(-0.05),0,expression(0.05),expression(0.1)))
#loAB <- loess(myres[,1]~c(1:10000))
#loAb <- loess(myres[,2]~c(1:10000))
#loaB <- loess(myres[,3]~c(1:10000))
loD <- loess(myres[,5]~c(1:10000))
loDm <- loess(mylow[,5]~c(1:10000))
loDp <- loess(myhigh[,5]~c(1:10000))
#lines(predict(loAB), col='black', lwd=2)
#lines(predict(loAb), col='blue', lwd=2)
#lines(predict(loaB), col='forestgreen', lwd=2)
polygon(c(1:10000,10000:1),c(predict(loDm),rev(predict(loDp)))*10+0.5,col='coral1',border='coral1')
lines(predict(loD)*10+0.5, col='red', lwd=2)
lines(myres[,1], col='black', lwd=2)
lines(myres[,2], col='blue', lwd=2)
lines(myres[,3], col='forestgreen', lwd=2)
#lines(myres[,5]*1000+0.5, col='red', lwd=2)
legend(0.3,1,c('AB','Ab','aB','D'),col=c('black','blue','forestgreen','red'),lty=rep(1,4),lwd=2)
dev.off()
}


#lines((myres[,5]*100+0.5),col='red')

#HWE
Ne<-10000
pA<-0.05
pB<-0.05
freqAB<-pA*pB
freqAb<-pA*(1-pB)
freqaB<-(1-pA)*pB
freqab<-(1-pA)*(1-pB)

#freqAb<-0.2
#freqaB<-0.2
#freqAB<-0.2
#freqab<-1-freqAB-freqaB-freqAb
s_vec<-c(1.01,1,1,1)


genotypes<-sample(1:4,2*Ne,prob=c(freqAB,freqAb,freqaB,freqab),replace=TRUE)


s_vec<-c(1.01,1,1,1)
system.time(simf(genotypes,s_vec,1000,Ne))



irec<-sample(1:Ne,rbinom(1,Ne,rho)) #recombining pairs
genotype_mat<-cbind(genotypes[1:Ne],genotypes[(Ne+1):(2*Ne)])
genotype_mat[irec,]<-t(apply(genotype_mat[irec,],MARGIN=1,FUN=function(x) if (x[1]==1 && x[2]==4){return(c(2,3))} else if (x[1]==2 && x[2]==3){return(c(1,4))} else return(x) ))
c(genotype_mat)
}

nex.gen.r(genotypes,0.5)

table(genotypes)*s_vec

nex.gen.s
genotypes