library(mixtools)
set.seed(100)
y <- matrix(rpois(70, 6), 10, 7)
y <- 1:10
cuts <- c(2, 5, 7)
out1 <- makemultdata(y, cuts = cuts) #does not do any EM, simply split pop for later giving initial conditions to EM
out1 #return simply x$ = input i.e. y, while y is how many samples per split group

# The sulfur content of the coal seams in Texas
A <- c(1.51, 1.92, 1.08, 2.04, 2.14, 1.76, 1.17)
B <- c(1.69, 0.64, .9, 1.41, 1.01, .84, 1.28, 1.59)
C <- c(1.56, 1.22, 1.32, 1.39, 1.33, 1.54, 1.04, 2.25, 1.49)
D <- c(1.3, .75, 1.26, .69, .62, .9, 1.2, .32)
E <- c(.73, .8, .9, 1.24, .82, .72, .57, 1.18, .54, 1.3)
length(A)
length(B)
length(C)
length(D)
length(E)
dis.coal <- makemultdata(A, B, C, D, E,cuts = median(c(A, B, C, D, E)))
em.out <- multmixEM(dis.coal)
em.out[1:4]
dis.coal <- makemultdata(c(A, B, C, D, E),cuts = median(c(A, B, C, D, E)))
em.out <- multmixEM(dis.coal)
summary(em.out )
summary.mixEM(em.out)
em.out
print.mvnpEM(em.out)

compCDF(dis.coal$x, em.out$posterior, xlab="Sulfur", ylab="", main="empirical CDFs")
density.spEM(em.out)
plot(em.out)
em.out

median(c(A, B, C, D, E))
dis.coal <- makemultdata(c(A, B, C, D, E),cuts = 1.1)
em.out <- multmixEM(dis.coal)
em.out 
plot(c(A, B, C, D, E))



data(faithful)
attach(faithful)
waiting
hist(waiting, main="Time between Old Faithful eruptions",xlab="Minutes", ylab="", cex.main=1.5, cex.lab=1.5, cex.axis=1.4)
wait1 <- normalmixEM(waiting, lambda = .5, mu = c(55, 80), sigma = 5)
wait1 

mydata<-nAB_l[order(nAB_l)]
plot(mydata)
hist(mydata)

mod2<-normalmixEM(mydata,lambda=.5,mu=c(2400,2900),sigma=10)
plot(mod2,density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8)
mod3<-normalmixEM(mydata,lambda=.5,mu=c(2400,2700,2900),sigma=10)
plot(mod3,density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8)

myn<-100
res<-rep(0,myn)
myp<-0.3
for (j in 1:2)
  {
    for (i in 1:100000)
    {
      if (i==1){res<-rep(0,myn)}
      res<-res+rbinom(p=myp,n=myn,size=2)
    }
  if (j==1){obj<-hist(res);rangeres<-max(res)-min(res);mybreaks<-seq(min(res)-rangeres/5,max(res)+rangeres/5,round((max(res)-min(res))/20))}
  obj<-hist(res,breaks=mybreaks)
}
obj
range(res)
range(mybreaks)

totnAB<-sum(res)
res<-rmultinom(p=rep(1/myn,myn),size=totnAB,n=1)

hist(res)

