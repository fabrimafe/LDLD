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
if (glmres<mypvalue_threshold ) { if (length(unique(data_temp$population))>1) {return(c(fx(data_temp),fx(data_temp0),glmres))} else {return(c(fx1(data_temp),fx(data_temp0),glmres))} } else return(1)}
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
