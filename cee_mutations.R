pp<-function(probs,i,j)
{
probs[i]/(1-)
}
probs<-c(0.3,0.3,0.3,0.1)
res<-c(rep("A",300),rep("G",300),rep("C",300),rep("T",100))

300*pp(table(res)/length(res),1,2)
300*pp(table(res)/length(res),1,3)
100*pp(table(res)/length(res),1,4)
table(res)

length(res[res=="G"])+length(res[res=="A"])*pp(table(res)/length(res),1,2)+length(res[res=="C"])*pp(table(res)/length(res),3,2)+length(res[res=="T"])*pp(table(res)/length(res),4,2)-
(length(res[res=="G"])*pp(table(res)/length(res),2,1)+length(res[res=="G"])*pp(table(res)/length(res),2,3)+length(res[res=="G"])*pp(table(res)/length(res),2,4))


300*pp(table(res)/length(res),1,3)
100*pp(table(res)/length(res),1,4)



