library(data.table)
library(ggplot2)
library(cowplot)
library(pROC)

###My functions###
#' Calculates the log Bayes factor of each S against the alternative with the highest lML
#' 
#' This function is not needed when the number of alternative models is only two, in this case, a simple cyclic shift would work (see example below)
#' x[order(-lML),.(bf=lML-shift(lML,type="cyclic"),trueS,model.alignment.stemCells,alignment.stemCells),by=simulationCondition]
#' 
#' @param x dataTable with columns simulationCondition, lML, simulationCondition, and model.alignment.stemCells
#' @value copy of x with a lBF column added. bf is the log Bayes factor of the given model against the best for the same simulationCondition group without itself
getlBF=function(x){
  if(!identical(key(x),c("simulationCondition","lML")))
  {
    setkey(x,simulationCondition,lML) #this is sorting lML in ascending order
  }
  lBF=unlist(lapply(x[,unique(simulationCondition)],FUN=function(condition,theData){
    out=vector(mode="numeric",length=theData[condition,.N])
    for (iS in 1:length(out)){
      out[iS]=theData[condition,][iS,lML]-theData[condition,][-iS,last(lML)]
    }
    out
  },theData=x))
  outTable=copy(x)
  outTable[,`:=`(lBF=lBF)]
}
###

###Main###
#Data parsing
theData=fread("/Users/diego/Documents/Ciencia/Postdoc/projects/flipFlop/simulationStudy/simStudy/analysis/currentSummary.csv")
theData[,`:=`(simulationCondition=paste(sep="_",alignment.stemCells,flipflop.lambda,flipflop.mu,flipflop.gamma,errorModel.kappaScale,nleaves))]
setkey(theData,simulationCondition,lML)
theData[,`:=`(trueS=alignment.stemCells==model.alignment.stemCells)]

#Checking for a good ESS threshold trying to maximize the precision at different logBFs 
resultsData=data.table(minEss=as.numeric(NULL),minlBF=as.numeric(NULL),ppv=as.numeric(NULL),n=as.numeric(NULL))
for(minEss in seq(from=50,to=500,length=100)){
  validConditions=theData[alignment.stemCells!=20 & minESSb>=minEss & maxRhat<=1.05,.N,by=simulationCondition][N==2,] #I need to remove the samples with alignment.stemCells!=20 because for those we do not have the result of the true model!
  bfData=getlBF(theData[validConditions,])
  for(minlBF in seq(from=5, to=10, by=1)){
    ppv=bfData[lBF>minlBF,mean(trueS)]
    n=bfData[lBF>minlBF,.N]
    resultsData=rbind(resultsData,data.table(minEss=minEss,minlBF=minlBF,ppv=ppv,n=n))
  }
}

resultsData[,`:=`(minlBFF=factor(minlBF))]
ggplot(resultsData,aes(x=minEss,y=ppv,color=minlBFF))+
  geom_point() +
  scale_color_brewer(type = "seq" , palette=3, name="Min log(BF)") +
  scale_y_continuous(name="Positive predictive value (precision)") +
  scale_x_continuous(name="Min ESS") +
  theme_cowplot()

#The classic >=200 is looking good here

#Trying to find the best logBF threshold and making the ROC curve
resultsData2=data.table(minEss=as.numeric(NULL),
                        minlBF=as.numeric(NULL),
                        p=as.numeric(NULL),
                        n=as.numeric(NULL),
                        pp=as.numeric(NULL),
                        pn=as.numeric(NULL),
                        tp=as.numeric(NULL),
                        fp=as.numeric(NULL),
                        tn=as.numeric(NULL),
                        fn=as.numeric(NULL),
                        nTrueInPP=as.numeric(NULL),
                        nTrueResolved=as.numeric(NULL))
minEss=200
validConditions=theData[alignment.stemCells!=20 & minESSb>=minEss & maxRhat<=1.05,.N,by=simulationCondition][N==2,] #I need to remove the samples with alignment.stemCells!=20 because for those we do not have the result of the true model!
bfData=getlBF(theData[validConditions,])

#ROC curve
roc_curve=roc(bfData$trueS,bfData$lBF,plot=T,print.auc=T)
my.coords <- coords(roc=roc_curve, x = "all", transpose = FALSE)
setDT(my.coords)
my.coords[,`:=`(optSS=sqrt(((1-specificity)^2)+((1-sensitivity)^2)))]
selectedThresholdSS=my.coords[order(-optSS,sensitivity),last(threshold)] ##Selected threshold if we weight specificity and sensitivity equally

#tops=theData[validConditions,][order(-lML),.(bf=lML-shift(lML,fill=-Inf,type="lead"),trueS,model.alignment.stemCells,alignment.stemCells),by=simulationCondition][,.SD[1],by=simulationCondition]
for(minlBF in seq(from=bfData[,min(lBF)]-1, to=bfData[,max(lBF)]+1, length=10000)){
  p=bfData[trueS==T,.N]
  n=bfData[trueS==F,.N]
  pp=bfData[lBF>=minlBF,.N]
  pn=bfData[lBF<minlBF,.N]
  tp=bfData[lBF>=minlBF & trueS==T,.N]
  fp=bfData[lBF>=minlBF & trueS==F,.N]
  tn=bfData[lBF<minlBF & trueS==F,.N]
  fn=bfData[lBF<minlBF & trueS==T,.N]
  nTrueInPP=bfData[lBF>=minlBF & trueS==T,.N,by=simulationCondition][,.N]#How many simulation conditions have the true in the positive
  nTrueResolved=bfData[lBF>=minlBF,.(.N,containTrue=sum(trueS)>=1),by=simulationCondition][N==1 & containTrue,.N]#How many simulation conditions have only the true in the positive
  resultsData2=rbind(resultsData2,data.table(minEss=minEss,minlBF=minlBF,p=p,n=n,pp=pp,pn=pn,tp=tp,fp=fp,tn=tn,fn=fn,nTrueInPP=nTrueInPP,nTrueResolved=nTrueResolved))
}
resultsData2[,`:=`(sensitivity=tp/p,specificity=tn/n,ppv=tp/pp,npv=tn/pn,mk=(tp/pp)+(tn/pn)-1,propTrueInPP=nTrueInPP/bfData[,length(unique(simulationCondition))],propTrueResolved=nTrueResolved/bfData[,length(unique(simulationCondition))])]

ggplot(resultsData2,aes(x=ppv,y=npv,color=minlBF))+
  geom_point() +
  scale_y_continuous(name="Negative predictive value") +
  scale_x_continuous(name="Positive predictive value") +
  scale_color_distiller(name="lBF",type="seq",palette = 3,direction = 1) +
  theme_cowplot()

ggplot(resultsData2,aes(x=pn/(p+n),y=npv,color=minlBF))+
  geom_point() +
  scale_y_continuous(name="Negative predictive value",lim=c(0,1)) +
  scale_x_continuous(name="Proportion Predicted negative") +
  scale_color_distiller(name="lBF",type="seq",palette = 3,direction = 1) +
  theme_cowplot()

selectedThresholdMkNpv=resultsData2[order(-mk,-npv),][!is.nan(mk),first(minlBF)]
selectedThresholdNpv=resultsData2[order(-npv),][npv!=1,first(minlBF)]
#resultsData2[,`:=`(opt=sqrt((1-npv)^2+(1-(pn/(p+n)))^2))]


##Proportion of conditions with a single true answer vs with the true answer in the "accepted answers" category
##Something between 0 and -12 seems ideal, if we weight them equaly
ggplot(resultsData2,aes(x=propTrueInPP,y=propTrueResolved,color=minlBF))+
  geom_point() +
  scale_y_continuous(name="Proportion conditions correctly fully resolved",lim=c(0,1)) +
  scale_x_continuous(name="Proportion of conditions with the true answer accepted") +
  scale_color_distiller(name="lBF",type="seq",palette = 3,direction = 1,lim=c(-12,0)) +
  theme_cowplot()

##Something between 0 and -12 seems ideal, if we weight them equally
resultsData2[,`:=`(optPtipPtr=sqrt(((1-propTrueInPP)^2)+((1-propTrueResolved)^2)))]
selectedThresholdPtipPtr=resultsData2[order(optPtipPtr),first(minlBF)]

selectedThresholdSS
selectedThresholdMkNpv
selectedThresholdNpv
selectedThresholdPtipPtr



