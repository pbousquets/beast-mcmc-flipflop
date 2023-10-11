suppressPackageStartupMessages({ 
  library(rwty)
  library(rstan)
  library(LaplacesDemon)
  library(HDInterval)
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(phangorn)
})

#My functions
#####

#' Modified version of analyze.rwty
#' 
my.analyze.rwty <- function(chains, burnin=0, window.size=20, treespace.points = 100, n.clades = 20,
                         min.freq = 0.0, fill.color = NA, filename = NA, 
                         overwrite=FALSE, facet=TRUE, free_y=FALSE, autocorr.intervals=100, ess.reps = 20,
                         treedist = 'PD', params = NA, max.sampling.interval = NA, 
                         plotParams=F,
                         ...){
  
  chains <- check.chains(chains)
  
  N <- length(chains[[1]]$trees)
  
  rwty:::rwty.params.check(chains, N, burnin, window.size, treespace.points, min.freq, filename, overwrite)
  
  plots=NULL
  
  if(plotParams){
    # check to see if ptables exist, make related plots
    if(all(unlist(lapply(chains, function(x) length(x$ptable[,1])))) > 0){ 
      # plot parameters for all chains
      parameter.plots <- makeplot.all.params(chains, burnin = burnin, facet=facet, strip = 1)
      parameter.correlations <- makeplot.pairs(chains, burnin = burnin, params = params, treedist = treedist, strip = 1)
      names(parameter.correlations) <- paste0(names(parameter.correlations), ".correlations")
    }
    else{
      parameter.plots <- makeplot.topology(chains, burnin = burnin, facet = facet)
      parameter.correlations <- NA
    }
    plots <- c(parameter.plots,
               parameter.correlations)
  }
  
  # plot autocorrelation
  if(N < 200){
    autocorr.plot <- NULL
  } else {
    autocorr.plot <- makeplot.autocorr(chains, burnin = burnin, autocorr.intervals = autocorr.intervals, facet = facet, max.sampling.interval = max.sampling.interval) 
  }
  
  # plot sliding window sf plots
  splitfreq.sliding <- makeplot.splitfreqs.sliding(chains, burnin=burnin, n.clades = n.clades, window.size = window.size, facet = facet)
  acsf.sliding <- makeplot.acsf.sliding(chains, burnin=burnin, window.size = window.size, facet = facet)
  
  # plot cumulative sf plots
  splitfreq.cumulative <- makeplot.splitfreqs.cumulative(chains, burnin=burnin, n.clades = n.clades, window.size = window.size, facet = facet)
  acsf.cumulative <- makeplot.acsf.cumulative(chains, burnin=burnin, window.size = window.size, facet = facet)
  
  # plot treespace for all chains
  treespace.plots <- makeplot.treespace(chains, n.points = treespace.points, burnin = burnin, fill.color = fill.color)
  
  # Add citations for all packages
  citations <- list(
    citation('rwty'),
    citation('ape'),
    citation('phangorn'),
    citation('ggplot2'),
    citation('coda'),
    citation('viridis'),
    citation('ggdendro'),
    citation('GGally'),
    citation('plyr'),
    citation('reshape2'),
    citation('stats')
  )
  
  
  plots <- c(plots,
             autocorr.plot,
             splitfreq.sliding,
             acsf.sliding,
             splitfreq.cumulative,
             acsf.cumulative,
             treespace.plots)
  
  
  # plot multichain plots when appropriate
  if(length(chains) > 1){
    asdsf.plot <- makeplot.asdsf(chains, burnin = burnin, window.size = window.size, min.freq = min.freq)
    
    splitfreq.matrix.plots <- makeplot.splitfreq.matrix(chains, burnin = burnin)
    
    plots <- c(plots, asdsf.plot, splitfreq.matrix.plots)
  }
  
  plots[["citations"]] <- citations
  
  # Print all to pdf if filename provided
  if(!is.na(filename)){
    print(sprintf("Saving plots to file: %s", filename))
    pdf(file=filename, width = 10, height = 7, ...)
    print(plots)
    dev.off()
  }
  
  return(plots)
  
}


#' Uses external BEAST programs to estimate the MCC tree after combining 
#' several tree traces
#' 
#' @param inputFiles vector with the .trees files to parse
#' @param outDir directory to write the intermediate and final MCC tree
#' @param burninStates number of states to remove from each tree trace
#' @param traceFileName name (or suffix) of the intermediate combined tree trace
#' @param mccFileName name (or suffix) of the final MCC tree file
#' @param appendDirname logical to indicate if the fileNames are suffixes or
#' full names
#' @param sumFun Name of the function to use to summarize node heights
#' @param removeTrace Removes the combined tree trace after having estimated the MCC tree succesfully
#' @return the MCC as an ape::phylo
estimateMCCtree=function(inputFiles,outDir,burninStates,traceFileName="combinedPruned.trees",mccFileName="MCC.tree",appendDirname=T,sumFun=c("mean","median"),removeTrace=T){
  sumFun=match.arg(sumFun)
  
  #Prepare output names
  if(appendDirname){
    traceFileName=paste0(basename(outDir),"_",traceFileName)
    mccFileName=paste0(basename(outDir),"_",mccFileName)
  }
  
  fullTraceFileName=paste(sep="/",outDir,traceFileName)
  fullMccFileName=paste(sep="/",outDir,mccFileName)
  
  if(!file.exists(fullMccFileName)){
    #Test if required BEAST executables are available
    if(system2('which','logcombiner',stdout = F, stderr = F)!=0){
      stop("Required logcombiner executable not in PATH")
    }
    if(system2('which','treeannotator',stdout = F, stderr = F)!=0){
      stop("Required treeannotator executable not in PATH")
    }
    
    #Run logcombiner
    logCombinerArgs=c("-trees","-burnin",burninStates,inputFiles,fullTraceFileName)
    logCombinerCode=system2('logcombiner',logCombinerArgs,stdout = F, stderr = F)
    
    #Run treeannotator
    treeAnnotatorArgs=c("-heights",sumFun,fullTraceFileName,fullMccFileName)
    treeAnnotatorCode=system2('treeannotator',treeAnnotatorArgs,stdout = F, stderr = F)

    #Read MCC tree and return it
    if(logCombinerCode == 0 & treeAnnotatorCode == 0){
      if(removeTrace){
        file.remove(fullTraceFileName)
      }
      return(read.nexus(fullMccFileName))
    } else {
      stop("Problem generating the MCC tree estimate")
      return (NA)
    }
  } else {
    message("Re-using previously-calculated MCC tree")
    return(read.nexus(fullMccFileName))
  }
}

#' Summary function for point estimates
#' 
#' @param x data
#' @param type character indicating if we want to use the median, mean, or mode
#' @value point estimate
sumFun=function(x,type=c("mean","median","mode")){
  type=match.arg(type)
  if(type=="mean"){
    return(mean(x))
  } else if (type=="median"){
    return(median(x))
  } else if (type=="mode") {
    return(mode(x))
  }
}

#' Tests if a value is within a closed interval
#' 
#' @param x value to test
#' @param interval named numeric vector with a lower and upper interval
#' @value logic value indicating if x is within interval
withinInterval=function(x,interval){
  if(all(names(interval)!=c("lower","upper"))){
    warning("Interval is not an interval, returning NA")
    return(NA)
  }else{
    if(x>=interval[1] & x<=interval[2]){
      return(T)
    }
  }
  return(F)
}

#' Method-of-moments estimator of AICM
#' Posterior simulation-based analogue of Akaike's information criterion (AIC)
#' through Markov chain Monte Carlo
#' Raftery et al. 2007
#' Translated from BEAST's logMarginalLikelihoodAICM function
#' 
#' @param loglikelihoods vector of posterior sample of loglikelihoods
#' @retun AICM
AICM=function(loglikelihoods){
  return(2*var(loglikelihoods)-2*mean(loglikelihoods))
}

#' Generates a data.table merging data from multiple MCMC chains
#' 
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param burninP Proportion of MCMC samples to discard as burnin
#' @return data.table with same columns as chain$ptable + a chain factor
getDataTable=function(chains,burninP){
  returnTable=data.table(chains[[1]]$ptable,chain=1)
  for(ichain in 2:length(chains)){
    returnTable=rbind(returnTable,data.table(chains[[ichain]]$ptable,chain=ichain))
  }
  returnTable[,`:=`(chain=as.factor(chain),burnin=ifelse(state<max(state)*burninP,T,F))]
  return(returnTable)
}

#' Generates a data.table with distances to the true tree from multiple MCMC tree chains
#' 
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param burninP Proportion of MCMC samples to discard as burnin
#' @param trueTree Phylo true tree
#' @param distance Distance to use to compare the tree samples against the true tree
#' @param normalize normalize the RF or wRF distance (phangorn)
#' @param rooted take bipartitions for rooted trees into account (phangorn)
#' @param check.labels compares labels of the trees (phangorn)
#' @return data.table with a state column, distance, and a chain factor
getDistanceToTrueTreeData=function(chains,burninP,trueTree,distance=c("RF","wRF"),normalize=T,rooted=T,check.labels=T){
  distance=match.arg(distance)
  calculateDistance=function(tree1,tree2,distance,normalize=T,rooted=T,check.labels=T){
    switch(distance,
           "RF" = RF.dist(tree1=tree1,tree2=tree2,normalize=normalize,rooted=rooted,check.labels = check.labels),
           "wRF" = wRF.dist(tree1=tree1,tree2=tree2,normalize=normalize,rooted=rooted,check.labels = check.labels),
           NA)
  }
  returnTable=data.table(state=chains[[1]]$ptable$state,
                         chain=1,
                         distance=sapply(chains[[1]]$trees,
                                         FUN=function(x,trueTree){calculateDistance(trueTree,x,distance,normalize=normalize,rooted=rooted,check.labels=check.labels)},
                                         trueTree=trueTree,
                                         simplify = T,
                                         USE.NAMES = F)
  )
  for(ichain in 2:length(chains)){
    returnTable=rbind(returnTable,data.table(state=chains[[ichain]]$ptable$state,
                                             chain=ichain,
                                             distance=sapply(chains[[ichain]]$trees,
                                                             FUN=function(x,trueTree){calculateDistance(trueTree,x,distance,normalize=normalize,rooted=rooted,check.labels=check.labels)},
                                                             trueTree=trueTree,
                                                             simplify = T,
                                                             USE.NAMES = F)
    )
    )
  }
  setnames(returnTable,old="distance",new=ifelse(normalize==T,paste0("normalized.",distance),paste0("absolute.",distance)))
  returnTable[,`:=`(chain=as.factor(chain),burnin=ifelse(state<max(state)*burninP,T,F))]
  return(returnTable)
}

#' Merges a list of rwty.chains into one, after proper burnin
#' 
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param burninP Proportion of MCMC samples to discard as burnin
#' @return rwty.chain with the chains concatenated with the burnin trees and all paramter values (ptable) removed.
mergeTreeTrace=function(chains,burninP){
  returnChain=chains[[1]]
  returnChain$ptable=NULL
  last=length(returnChain$trees)
  first=round(burninP*last)
  returnChain$trees=returnChain$trees[first:last]
  for(ichain in 2:length(chains)){
    returnChain$trees=c(returnChain$trees,chains[[ichain]]$trees[first:last])
  }
  return(returnChain)
}

#' Generates a data.table with accuracy statistics
#' 
#' @param theData MCMC data in data.table format after getDataTable
#' @param params Named vector with true (simulated) parameter values
#' @param fun Function to calculate the parameter estimate
#' @param credMass Mass for the hight posterior density interval
#' @return data.table with a row per parameter and columns for relative error and value within 0.9HPDI
getAccuracyTableFromDataTable=function(theData,params,fun=sumFun,credMass=0.95){
  returnTable=data.table(
    t(
      sapply(names(params),FUN=function(paramName,params){
        trueVal=params[paramName]
        names(trueVal)=NULL
        c("withinHDI"=withinInterval(trueVal,hdi(theData[[paramName]],credMass = credMass)),
          "rr"=(sumFun(theData[[paramName]])-trueVal)/trueVal)
      },simplify = T,params=params)
    ),
    keep.rownames = T
  )
  setnames(returnTable,old="rn",new="param")
  returnTable[,`:=`(withinHDI=as.logical(withinHDI))]
  return(returnTable)
}

#' Generates a table with summary statistics for all parameters
#'
#' @param theData MCMC data in data.table format after getDataTable
#' @param params List of parameter names to analyze
#' @param credMass Mass for the hight posterior density interval
#' @return data.table with a row per parameter and colums for mean, median, hdilower and hdiupper
getSummaryTableFromDataTable=function(theData,params,credMass=0.95){
  returnTable=data.table(
    t(theData[,lapply(.SD,function(x,credMass=0.95){
      theHDI=hdi(x,credMass=credMass)
      c("mean"=mean(x,na.rm=T),
        "median"=median(x,na.rm=T),
        "HDILower"=theHDI[1],
        "HDIupper"=theHDI[2])},
      credMass=0.95),.SDcols=params]
    ),
    keep.rownames = T
  )
  setnames(returnTable,c("param","mean","median","HDILower","HDIUpper"))
  return(returnTable)
}

#' Generates a data.table with convergence statistics for all parameters
#' 
#' @param theData MCMC data in data.table format after getDataTable
#' @param params List of parameter names to analyze
#' @return data.table with a row per parameter and colums for Rhat, ESSB and ESST
getConvergenceTableFromDataTable=function(theData,params){
  returnTable=data.table(
    t(
      sapply(
        params,FUN=function(param){
          stanMatrix=as.matrix(dcast(theData[,c("state",param,"chain"),with=F],state~chain,value.var = param)[,-1])
          c("Rhat"=Rhat(stanMatrix),"essB"=ess_bulk(stanMatrix),"essT"=ess_tail(stanMatrix))
        },simplify = T
      )
    ),
    keep.rownames = T
  )
  setnames(returnTable,old="rn",new="param")
  return(returnTable)
}

#' Generates a table with convergence statistics for all parameters
#'
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param params List of parameter names to analyze
#' @param burnin Number of MCMC samples to discard as burnin
#' @return data.table with a row per parameter and colums for Rhat, ESSB and ESST
getConvergenceTableFromChains=function(chains,params,burnin){
  returnTable=data.table(
    t(
      sapply(
        params,FUN=function(param){
          stanTable=data.table(chains[[1]]$ptable[-c(1:burnin),param])
          for(ichain in 2:length(chains)){
            colname=paste0("V",ichain)
            stanTable[,(colname):=(data.table(chains[[ichain]]$ptable[-c(1:burnin),param]))]
          }
          stanMatrix=as.matrix(stanTable)
          c("Rhat"=Rhat(stanMatrix),"essB"=ess_bulk(stanMatrix),"essT"=ess_tail(stanMatrix))
        },simplify = T
      )
    ),
    keep.rownames = T
  )
  setnames(returnTable,old="rn",new="param")
  return(returnTable)
}

#####

#Functions to generate Prior plots. WARNING: HARDCODED
#####

#' Prior density of a parameter, given its value and name
#' 
#' @param x value
#' @param parameter name of the parameter
#' @return density of x for the parameter given
dPrior=function(x,parameter=c("flipflop.lambda","flipflop.mu","flipflop.gamma","errorModel.kappaScale","errorModel.deltaOffset","errorModel.etaOffset","constant.popSize","luca_branch")){
  parameter=match.arg(parameter)
  if(parameter=="flipflop.lambda"){
    sd=1
    return(dhalfnorm(x,scale = sqrt(pi/2)/sd))
  } else if(parameter=="flipflop.mu" | parameter=="flipflop.gamma"){
    sd=0.05
    return(dhalfnorm(x,scale = sqrt(pi/2)/sd))
  } else if (parameter=="errorModel.kappaScale"){
    meanlog=4.56
    sdlog=0.3
    return(dlnorm(x,meanlog=meanlog,sdlog=sdlog))
  } else if (parameter=="errorModel.etaOffset"){
    shape1=95
    shape2=5
    return(dbeta(x,shape1,shape2))
  } else if (parameter=="errorModel.deltaOffset"){
    shape1=5
    shape2=95
    return(dbeta(x,shape1,shape2))
  } else if (parameter=="constant.popSize"){
    return(1/x)
  } else if (parameter=="luca_branch"){
    min=0
    max=50
    return(dunif(x,min=min,max=max))
  } else {
    return(NA)
  }
}

#' Domain of a Prior giving start and end cumulative probabilities
#' Used to overlay the prior distribution on a posterior sample
#' 
#' @param n number of samples of points of prior's domain to return
#' @param parameter name of the parameter
#' @param minP lower cumulative probability to start the domain
#' @param maxP higher cumulative probability to stop the domain
#' @return data.table with an x column with the domain of the prior
domainPrior=function(n=1000,parameter=c("flipflop.lambda","flipflop.mu","flipflop.gamma","errorModel.kappaScale","errorModel.deltaOffset","errorModel.etaOffset","constant.popSize","luca_branch"), minP=0.05, maxP=0.95){
  parameter=match.arg(parameter)
  from=0
  to=0
  
  if(parameter=="flipflop.lambda"){
    sd=1
    from=qhalfnorm(minP,sqrt(pi/2)/sd)
    to=qhalfnorm(maxP,sqrt(pi/2)/sd)
  } else if(parameter=="flipflop.mu" | parameter=="flipflop.gamma"){
    sd=0.05
    from=qhalfnorm(minP,sqrt(pi/2)/sd)
    to=qhalfnorm(maxP,sqrt(pi/2)/sd)
  } else if (parameter=="errorModel.kappaScale"){
    meanlog=4.56
    sdlog=0.3
    from=qlnorm(minP,meanlog=meanlog,sdlog=sdlog)
    to=qlnorm(maxP,meanlog=meanlog,sdlog=sdlog)
  } else if (parameter=="errorModel.etaOffset"){
    shape1=95
    shape2=5
    from=qbeta(minP,shape1,shape2)
    to=qbeta(maxP,shape1,shape2)
  } else if (parameter=="errorModel.deltaOffset"){
    shape1=5
    shape2=95
    from=qbeta(minP,shape1,shape2)
    to=qbeta(maxP,shape1,shape2)
  } else if (parameter=="constant.popSize"){
    #Quantiles not available, we'll go from 0.1 to 100
    from=1
    to=100
  } else if (parameter=="luca_branch"){
    min=0
    max=50
    from=qunif(minP,min=min,max=max)
    to=qunif(maxP,min=min,max=max)
  } else {
    return(NA)
  }
  return(data.table(x=seq(length=n,from=from,to=to)))
}
#####

#Plotting functions
#####
#' Plots posterior sample, prior, true, and estimated values for a parameter in pdf
#' 
#' @param thisData data.table, with burnin removed if needed
#' @param trueParameterValues named vector of true/simulated values
#' @param outDir output directory
#' @param outName output name
#' @param ndigits number of digits to round numbers printed in the plot
#' @param width width in inches
#' @param height height in inches
#' @param ndigits number of digits to round the estimated and simulated values labeled on the plot
#' @return void
writePlotEstimatedContinuousParameters=function(thisData,trueParameterValues,outDir,outName,ndigits=3,width=7,height=7){

  pdf(file=paste0(outDir,"/",outName),width = width, height = height)

  for(param in names(trueParameterValues)){
    theplot=ggplot(data=thisData,aes(x=.data[[param]],group=chain,fill=chain,color=chain)) +
      geom_density(alpha=0.3,linewidth=1) +
      #stat_function(fun=dPrior,args=list(parameter=param),color="red",group=1) +
      geom_area(data=domainPrior(parameter=param),stat='function',aes(x=x,color="Prior",fill="Prior"),fun=dPrior,args=list(parameter=param),alpha=0.3,inherit.aes = F) +
      geom_vline(aes(xintercept=trueParameterValues[param],linetype="Simulated"),linewidth=1,color="red") +
      stat_summary(aes(xintercept = after_stat(x),y=0,linetype="Estimated",fill=NULL,group=NULL),fun=sumFun,geom="vline",orientation="y",color="black") +
      stat_summary(fun=sumFun,geom="text",aes(x=.data[[param]],label=round(after_stat(x),digits=ndigits),y=0),orientation="y",vjust=1.5,hjust=-0.5,inherit.aes=F)+
      annotate(geom = "text",x=trueParameterValues[param],y=0,label=round(trueParameterValues[param],digits=ndigits),color="red",vjust=1.5,hjust=1.5) +
      scale_linetype(name="") +
      scale_color_brewer(name="Chain",type = "qual",palette = 2) +
      scale_fill_brewer(name="Chain", type = "qual",palette = 2) +
      scale_x_continuous(name="Parameter value") +
      scale_y_continuous(name="Density")+
      labs(title=param)+
      theme_cowplot()
    print(theplot)
    #save_plot(filename = paste0(outDir,"/",param,format),theplot,base_height = 6)
  }
  dev.off()
}

#' Plots posterior sample and estimated values for a parameter
#' 
#' @param thisData data.table, with burnin removed if needed
#' @param parameterNames vector with the names of the parameters to print
#' @param outDir output directory
#' @param outName output name
#' @param ndigits number of digits to round numbers printed in the plot
#' @param width width in inches
#' @param height height in inches
#' @param ndigits number of digits to round the estimated and simulated values labeled on the plot
#' @return void
writePlotContinuousParameters=function(thisData,parameterNames,outDir,outName,ndigits=3,width=7,height=7){
  pdf(file=paste0(outDir,"/",outName),width = width, height = height)
  for(param in parameterNames){
    theplot=ggplot(data=thisData,aes(x=.data[[param]],group=chain,fill=chain,color=chain)) +
      geom_density(alpha=0.3,linewidth=1) +
      stat_summary(aes(xintercept = after_stat(x),y=0,linetype="Estimated",fill=NULL,group=NULL),fun=sumFun,geom="vline",orientation="y",color="black") +
      stat_summary(fun=sumFun,geom="text",aes(x=.data[[param]],label=round(after_stat(x),digits=ndigits),y=0),orientation="y",vjust=1.5,hjust=-0.5,inherit.aes=F)+
      scale_linetype(name="") +
      scale_color_brewer(name="Chain",type = "qual",palette = 2) +
      scale_fill_brewer(name="Chain", type = "qual",palette = 2) +
      scale_x_continuous(name="Parameter value") +
      scale_y_continuous(name="Density")+
      labs(title=param)+
      theme_cowplot()
    print(theplot)
    #save_plot(filename = paste0(outDir,"/",param,format),theplot,base_height = 6)
  }
  dev.off()
  
}
#####

#DEBUG
#baseDir="~/projects/flipFlop/simulationStudy/simStudy"
#args=c(0.1,8,paste(sep="/",baseDir,"analysis"),paste(sep="/",baseDir,"333/sim_3_0.1_0.001_0.001_100_9_3cells.trees"),paste(sep="/",baseDir,"666/sim_3_0.1_0.001_0.001_100_9_3cells.trees"),paste(sep="/",baseDir,"999/sim_3_0.1_0.001_0.001_100_9_3cells.trees"))
#args=c(0.1,8,paste(sep="/",baseDir,"analysis"),paste(sep="/",baseDir,"333/sim_10_1.0_0.01_0.01_100_9_10cells.trees"),paste(sep="/",baseDir,"666/sim_10_1.0_0.01_0.01_100_9_10cells.trees"),paste(sep="/",baseDir,"999/sim_10_1.0_0.01_0.01_100_9_10cells.trees"))

#####

#Parsing command-line arguments and checking input
#####
args <- commandArgs(trailingOnly = TRUE)

outDir=""
if (length(args) <=4) {
  print("Usage: script burnin_proportion n_cores baseOutdir sample1.tree ... sampleN.tree")
  q()
} else {
  print(paste0("Running the script in the trace files (",paste(collapse =",",args[-c(1:3)]),") with a burnin of ",args[1]," proportion of the trees and ",args[2]," processors"))
}

rwty.processors=as.numeric(args[2])
setDTthreads(rwty.processors)
burninP=as.numeric(args[1])
baseDir=args[3]
files=args[-c(1:3)]

if(!all(file.exists(files))){
  stop("ERROR: Not all input files exist. ",paste(files,collapse=", "))
}
#####

#Preparing output
#####
outDir=paste(sep="/",baseDir,gsub(x = basename(files[1]),pattern = ".trees",replacement = ""))
dir.create(outDir,recursive=T)
print(paste0("Saving outputs in ",outDir))
#####

#Hardcoded output config
#####
summaryOutputSuffix="summary.csv"
detailedOutputSuffix="detailed.csv"
plotName.estimatedContinuousParameters="plotsEstimatedContinuousParameters.pdf"
plotName.otherContinuousParameters="plotsOtherContinuousParameters.pdf"
plotName.treeAccuracyParameters="plotsTreeAccuracyParameters.pdf"
plotName.rwty="plotsRWTY.pdf"
plotName.prefix.correlations="plotCorrelations."
#####

#Parsing Simulation parameterization
#####
#WARNING: there was a small prior parameterization error in the xml files, luca_branch should be a U[0,45] instead of [0,50]
#TODO check units of constant.popsize

#Simulation parameters that are fixed for all replicates. WARNING: HARDCODED
fixedArguments=c("errorModel.deltaOffset"=0.05,"errorModel.etaOffset"=0.93,"constant.popSize"=5,"luca_branch"=13.34)
#Variable simulation parameters: Parsed from filename. WARNING: HARDCODED
simulationArguments=as.numeric(regmatches(files[1],regexec('sim_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)cells.trees',files[1]))[[1]][-1])
names(simulationArguments)=c("alignment.stemCells","flipflop.lambda","flipflop.mu","flipflop.gamma","errorModel.kappaScale","nleaves","model.alignment.stemCells")
#All simulation parameters
simulationArguments=c(simulationArguments,fixedArguments)
#True values for all continuous variables to be estimated. WARNING: HARDCODED
trueContinuousParameters=simulationArguments[!names(simulationArguments)%in%c("alignment.stemCells","nleaves","model.alignment.stemCells")]
#True tree. WARNING: HARDCODED
trueTreeFile=paste0(dirname(baseDir),"/data/","test_",simulationArguments["nleaves"],"samples.tree")
#####

#Parsing input traces
#####
#Constant parameters that should be eliminated from the traces. WARNING: HARDCODED.
constantParams=c("luca_height","alignment.stemCells","clock.rate","treeModel.rootHeight") #treeModel.rootHeight is not constant, but it is a linear transformation of luca_branch
##TODO TEST IF THE FILES EXIST, OR OTHERWISE TERMINATE
sink('/dev/null') #shut up!
chains=invisible(lapply(files,FUN=function(x,constantParams){
  chain=load.trees(x,format="BEAST",logfile = gsub(x,pattern = ".trees",replacement = ".log"));
  chain$ptable=chain$ptable[,!names(chain$ptable)%in%constantParams];
  chain}
  ,constantParams=constantParams))
sink()
nSamples=unique(lapply(chains,FUN=function(x){length(x$trees)}))
if(length(nSamples) != 1){
  stop("ERROR: chains are not of the same length")
}
#####

#Data reorganization
#####
thisData=getDataTable(chains,burninP)
burninStates=thisData[,max(state)]*burninP
burninSamples=thisData[chain==1 & burnin==T,.N]
#####

#Assessing convergence of continuous variables
#####
paramsForConvergence=names(thisData)
paramsForConvergence=paramsForConvergence[!paramsForConvergence%in%c("state","chain","burnin")]
continuousConvergence=getConvergenceTableFromDataTable(thisData[burnin==F,],paramsForConvergence)
#####

#Assessing estimation accuracy of continuous variables
#####
continuousAccuracy=getAccuracyTableFromDataTable(thisData[burnin==F,],trueContinuousParameters)
#####

#Get marginal likelihood estimates from the posterior sample (EVIL, we should be doing PS or SS)
#####
thisLML=LML(LL=thisData[burnin==F,likelihood],method="HME")$LML
thisAICM=AICM(thisData[burnin==F,likelihood]) #The lower the better
#####

#Plotting continuous parameters
#####
#Overlap of prior, and posterior sample, including point estimate and true value
writePlotEstimatedContinuousParameters(thisData[burnin==F],trueContinuousParameters,outDir,plotName.estimatedContinuousParameters)

#Posterior sample for parameters without true values
writePlotContinuousParameters(thisData[burnin==F],paramsForConvergence[!paramsForConvergence%in%names(trueContinuousParameters)],outDir,plotName.otherContinuousParameters)
#####

#Assess tree convergence
#####

#RWTY approxESS and pseudoESS
#topologyApproxESS=rwty::topological.approx.ess(chains,burnin=burninSamples,treedist = 'RF')
#topologyPseudoESS=rwty::topological.pseudo.ess(chains,burnin=burninSamples,treedist = 'RF')

#RWTY approxESS and pseudoESS after merging the traces to be able to compare better with the continuous parameters
mergedchain=mergeTreeTrace(chains,burninP)
topologyApproxESSMerged=rwty::topological.approx.ess(mergedchain,burnin=0,treedist = 'RF')
topologyPseudoESSMerged=rwty::topological.pseudo.ess(mergedchain,burnin=0,treedist = 'RF')
finalTopologyApproxESS=ifelse(topologyApproxESSMerged$operator=="=",topologyApproxESSMerged$approx.ess,paste0(topologyApproxESSMerged$operator,topologyApproxESSMerged$approx.ess))
finalTopologyPseudoESS=mean(topologyPseudoESSMerged$Chain.1)

#Comparison of sum vs merge. We can see they are very similar, at least in this case
#sum(topologyApproxESSMerged$approx.ess)-topologyApproxESSMerged$approx.ess
#mean(apply(topologyPseudoESS,FUN=sum,MARGIN = 1))-mean(topologyPseudoESSMerged$Chain.1)

#Equivalent to pseudoESS but with the true tree as focal point, to then calculate Rhat and ESSB
#pseudoESS is just calculating the ESS of a continuous variable that summarizes the tree, in this case, the distance to a focal tree
#RWTY paper recommends calculating several and using the mean, their estimate is very close to the approxESS but is not really an ESS
#Here we know the true tree, so I am using it as focal tree. I am not 100% sure this is the best way to proceed, but I am also using the approxESS as gold standard
trueTree=read.tree(trueTreeFile)
theTreeTopologyData=getDistanceToTrueTreeData(chains,burninP,trueTree,"RF",rooted=F) #RWTY unroots the trees automatically, so I am doing the same here
theTreeData=suppressMessages(getDistanceToTrueTreeData(chains,burninP,trueTree,"wRF",rooted=F))
theTreeData=merge(theTreeTopologyData,theTreeData,by = c("state","chain","burnin"))
treeConvergence=getConvergenceTableFromDataTable(theTreeData[burnin==F,],c("normalized.RF","normalized.wRF"))

#Assess tree accuracy
#####
#I need to re-estimate the MCC tree because I did it independently for each replicate
estimatedTree=estimateMCCtree(files,outDir,burninStates,sumFun="mean")
#Tree topology: RF distance
thisRFdist=RF.dist(trueTree,estimatedTree,normalize=T,rooted = T)
#Tree topology with branch lengths: RFL distance with k = 1 (Robinson and Foulds 1979, see https://academic.oup.com/sysbio/article/64/2/205/1630737) 
thiswRFdist=wRF.dist(trueTree,estimatedTree,normalize = T,rooted = T)
#Plot
writePlotContinuousParameters(theTreeData[burnin==F,],c("normalized.RF","normalized.wRF"),outDir,plotName.treeAccuracyParameters)

#Table
treeAccuracy=getAccuracyTableFromDataTable(theTreeData[burnin==F,],c("normalized.RF"=0,"normalized.wRF"=0))
treeAccuracy[rr==Inf,`:=`(rr=NA)] #We can't actually calculate the RR because the true value is 0. We are keeping this as NA to have the same format in all parameters
#####

#RWTY topology and correlation outputs
#####
chainsWithWRF=chains
for(ichain in 1:length(chainsWithWRF)){
  chainsWithWRF[[ichain]]$ptable$normalized.wRF=theTreeData[chain==ichain,normalized.wRF]
}

#Topology-related plots if they make sense
if(!any(theTreeData[burnin==F,.(variance=var(normalized.RF)),by=chain]$variance==0))
{
  rwtyPlots=my.analyze.rwty(chains, burnin=burninSamples, fill.color = 'posterior',filename = paste0(outDir,"/",plotName.rwty))
}
#continuous parameters and topology (if variable)
pairsPlot=makeplot.pairs(chainsWithWRF,burnin=burninSamples,focal.tree = trueTree,params = c(paramsForConvergence[!paramsForConvergence%in%"coalescent"],"normalized.wRF"))

for (plotname in names(pairsPlot)){
  jpeg(filename = paste0(outDir,"/",plotName.prefix.correlations,plotname,".jpeg"),width = 1200, height = 1200,type = "quartz")
  print(pairsPlot[[plotname]])
  dev.off()
}

#####

#Write detailed output
#####

#Continuous parameters
continuousSummary=getSummaryTableFromDataTable(thisData[burnin==F,],paramsForConvergence)
setkey(continuousSummary,param)
setkey(continuousAccuracy,param)
setkey(continuousConvergence,param)
continuousResults=merge(continuousAccuracy,continuousConvergence,all=T)
continuousResults=merge(continuousResults,continuousSummary,all=T)

#Tree parameters
treeSummary=getSummaryTableFromDataTable(theTreeData[burnin==F,],c("normalized.RF","normalized.wRF"))
setkey(treeSummary,param)
setkey(treeAccuracy,param)
setkey(treeConvergence,param)
treeResults=merge(treeAccuracy,treeConvergence,all=T)
treeResults=merge(treeResults,treeSummary,all=T)

#All
allDetailedResults=rbind(continuousResults,treeResults)
allDetailedResults=cbind(data.table(t(simulationArguments)),allDetailedResults)##Adding the simulation arguments as columns to make an id later after merging the data from all simulations

detailedOutputFileName=paste0(outDir,"/",basename(outDir),"_",detailedOutputSuffix)
write.csv(allDetailedResults,file = detailedOutputFileName,quote = F,row.names = F)
#####

#Write summary
#####
#TODO Eliminate the popSize filter, only here because I need to transform units!
summaryOutputFileName=paste0(outDir,"/",basename(outDir),"_",summaryOutputSuffix)

#I could add this here too, but I do have this information on the detailed output, so this is probably not a good idea
#meansForSummary=data.table(t(allDetailedResults[,.(mean)]))
#setnames(meansForSummary,c(allDetailedResults[,paste0("mean.",param)]))
#

summaryData=data.table(t(simulationArguments),RFdist=thisRFdist,wRFdist=thiswRFdist,approxTopologyESS=finalTopologyApproxESS,pseudoTopologyESS=finalTopologyPseudoESS,lML=thisLML,AICm=thisAICM,maxRR=allDetailedResults[param!="constant.popSize",][!is.na(rr)][order(-abs(rr)),][1,rr],allWithinHDI=allDetailedResults[,all(withinHDI,na.rm=T)],maxRhat=allDetailedResults[,max(Rhat,na.rm=T)],minESSb=allDetailedResults[,min(essB,na.rm=T)])
write.csv(summaryData,file = summaryOutputFileName,quote = F,row.names = F)
#####