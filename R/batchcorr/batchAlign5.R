#' Get a peak table from an XCMS object
#'
#' peakTab will extract a peak table from an XCMS object with samples as rows and features as columns. Features are named by 'm/z@rt'
#' @param XS an XCMS object
#' @return a peak table with samples as rows and features as columns, features (colnames) named by 'm/z@rt'
#' @export
peakTab=function(XS) {
  pTab=peakTable(XS)
  id=paste('RP',pTab$mz,'@',pTab$rt,sep='')
  X=t(pTab[,-1:-13])
  colnames(X)=id
  return(X)
}

#' BA: Grab batch relevant data for a subtype of samples
#'
#' grabAlign will bring out sampletype:'grp' from batch:'batch' of an XCMS object
#' @param XS an XSMS object
#' @param batch a batch identifier
#' @param grp a sample type identifier
#' @return a peak table with samples as rows and features as columns, features (colnames) named by 'm/z@rt'
#' @export
grabAlign=function(XS,batch,grp) {
  incl=(XS@phenoData[,1]==batch & XS@phenoData[,2]==grp)
  peakTab=peakTab(XS)
  QC=peakTab[incl,]
  return(QC)
}

#' BA: Aggregate presense/missingness per batch
#'
#' batchFlag will aggregate presentness/missingness per feature at batch level within the combination: batch x sample type
#' @param PTnofill a peaktable with missing values (ie without either hard filling or imputation)
#' @param meta a matrix of metadata with batch in column 1 and sample type in column 2
#' @param peakInfo a matrix with m/z (col1) and rt (col2) of features (rows)
#' @param NAhard proportion of NAs within batch for feature to be considered missing
#' @param NAsoft experimental: Use at own risk
#' @param quantileSoft experimental: Use at own risk
#' @return meta: metadata aggretated by batch x sample type
#' @return flagHard: a table of presentness/missingness per batch x sample type
#' @return mz: m/z values of the features
#' @return rt: rt values of the features
#' @return flagAll: experimental: Use at own risk
#' @return flagSoft: experimental: Use at own risk
#' @return flagType: experimental: Use at own risk
#' @export
batchFlag=function(PTnofill,meta,peakInfo,NAhard=0.8,NAsoft=0.5,quantileSoft=0.1) {
  quant=quantile(PTnofill,quantileSoft,na.rm=TRUE)
  batch=meta[,1]
  uniqBatch=unique(batch)
  grp=meta[,2]
  uniqGrp=unique(grp)
  n=length(uniqBatch)*length(uniqGrp)
  flagAll=flagType=flagHard=flagSoft=matrix(nrow=n,ncol=ncol(PTnofill))
  colnames(flagAll)=colnames(flagType)=colnames(PTnofill)
  batchMeta=matrix(nrow=n,ncol=2)
  colnames(batchMeta)=c('batch','grp')
  i=0
  for (b in uniqBatch) {
    for (g in uniqGrp) {
      i=i+1
      PTsub=PTnofill[batch==b & grp==g,]
      NAs=apply(PTsub,2,function(x) sum(is.na(x))/length(x))
      intMean=apply(PTsub,2,function(x) mean(x,na.rm=TRUE))
      flagAll[i,]=ifelse(NAs>=NAhard,1,ifelse(NAs>=NAsoft & intMean<quant,1,0))
      flagType[i,]=ifelse(NAs>=NAhard,'hard',ifelse(NAs>=NAsoft & intMean<quant,'soft',''))
      flagHard[i,]=ifelse(NAs>=NAhard,1,0)
      flagSoft[i,]=flagAll[i,]-flagHard[i,]
      batchMeta[i,]=c(b,g)
    }
  }
  # Convert flags to !NA flags
  flagAll=1-flagAll
  flagHard=1-flagHard
  flagSoft=1-flagSoft
  return(batchFlag=list(meta=batchMeta,flagAll=flagAll,flagType=flagType,flagHard=flagHard,flagSoft=flagSoft,mz=peakInfo[,1],rt=peakInfo[,2]))
}

#' BA: Find alignment candidates
#'
#' align will find candidates for alignment among features that have systematic missingness "within box" of maximum m/z & rt differences. Used internally.
#' @param flags a table of presentness/missingness per batch x sample type
#' @param mz m/z values of the features
#' @param rt rt values of the features
#' @param mzdiff maximum distance in m/z to be considered for alignment
#' @param rtdiff maximum distance in rt to be considered for alignment
#' @return events: all possible alignments
#' @return features: all features involved in possible alignments
#' @return clusters: groups of features linked through mutual events
#' @export
align=function(flags,mz,rt,mzdiff=0.005,rtdiff=10) {
  misIndex=which(colSums(flags)!=nrow(flags) & colSums(flags)!=0) # Take out misalignment candidates
  event=list()
  c=0
  for (i in 1:(length(misIndex)-1)) {
    m=misIndex[i] # Take out next available misalignment candidate
    misPlus=misIndex[misIndex>m]# Look at the ones with higher index
    rtPlus=rt[misPlus] # Their rts
    mzPlus=mz[misPlus] # and mzs
    inBox=which(abs(rtPlus-rt[m])<rtdiff & abs(mzPlus-mz[m])<mzdiff) # Which are in the mz*rt box?
    if (length(inBox)>0) {
      for (j in 1:length(inBox)) {
        n=misPlus[inBox[j]]
        dotProd=sum(flags[,m]*flags[,n]) # Check overlapping features within box
        if (dotProd==0) { # For systematic misalignment to occur, these must be missed between features, ie dotProd=0
          c=c+1
          event[[c]]=c(c,m,n,dotProd,sum(flags[,m]),sum(flags[,n]),sum(flags[,m])+sum(flags[,n])) # Take out data
        }
      }
    }
  }
  if (length(event)==0) {
    return(list())
  }
  event=matrix(unlist(event),ncol=7,byrow=T)  # Rearrange
  event=as.data.frame(event)
  colnames(event)=c('event','m','n','dotProd','mSum','nSum','newSum')
  uniqFeat=unique(c(event$m,event$n))
  allMN=c(event$m,event$n)
  uniqAll=unique(allMN)
  duplAll=allMN[duplicated(c(allMN))]
  uniqM=unique(event$m)
  duplM=event$m[duplicated(event$m)]
  uniqN=unique(event$n)
  duplN=event$n[duplicated(event$n)]
  duplMN=event$n[event$n%in%event$m]
  evTemp=event
  features=as.data.frame(matrix(nrow=1,ncol=8+nrow(flags)))
  colnames(features)=c('cluster','featureIndex','mz','rt','to','from','MN','NN',paste('bF',1:nrow(flags),sep=''))
  ## Cluster together same m's
  for (j in 1:length(uniqM)) {
    uF=uniqM[j]
    mNum=sum(event$m==uF)
    nNum=sum(event$n==uF)
    featIndex=c(j,uF,mz[uF],rt[uF],0,0,0,0,flags[,uF])
    features=rbind(features,featIndex)
    for (mIndex in 1:mNum) {
      eventNum=which(event$m==uF)[mIndex]
      nM=event$n[eventNum]
      evTemp[eventNum,]=NA
      MN=ifelse(nM%in%duplMN,1,0)
      NN=ifelse(nM%in%duplN,1,0)
      featIndex=c(j,nM,mz[nM],rt[nM],0,1,MN,NN,flags[,nM])
      features=rbind(features,featIndex)
    }
  }
  features=features[-1,]
  ## cluster together same n's
  for (j in duplN) {
    clu=unique(features[features$featureIndex%in%j,1])
    for (k in 2:length(clu)) {
      features[features[,1]%in%clu[k],1]=clu[1]
    }
  }
  features[features$featureIndex%in%duplN,6]=features[features$featureIndex%in%duplN,6]+1  # Add flag for combination of features
  ## cluster together common m&n's
  for (j in duplMN) {
    clu=unique(features[features$featureIndex%in%j,1])
    for (k in 2:length(clu)) {
      features[features[,1]%in%clu[k],1]=clu[1]
    }
  }
  features[features$featureIndex%in%duplMN,5]=features[features$featureIndex%in%duplMN,5]+1  # Add flag for combination of features
  features=features[!duplicated(features$featureIndex),]
  features$cluster=as.numeric(as.factor(features$cluster))
  clusters=clusterMatrix(features)

  return(list(events=event,features=features,clusters=clusters))
}

#' BA: Sort features (alignment candidates) into clusters
#'
#' clusterMatrix will sort features into clusters. Used internally.
#' @param features all features involved in possible alignments
#' @return clusters: groups of features linked through mutual events
#' @export
clusterMatrix=function(features) {
  nClust=length(unique(features$cluster))
  clustMat=matrix(nrow=nClust,ncol=3+length(grep('bF',colnames(features))))
  colnames(clustMat)=c('cluster','nFeat','dotProd',paste('bF',1:length(grep('bF',colnames(features))),sep=''))
  for (c in 1:nClust) {
    clustSub=features[features$cluster==c,]
    vectMult=as.numeric(clustSub[1,grep('bF',colnames(features))])
    for (f in 2:nrow(clustSub)) {
      dotProd=sum(vectMult*clustSub[f,grep('bF',colnames(features))])
      vectMult=vectMult+as.numeric(clustSub[f,grep('bF',colnames(features))])
    }
    clustMat[c,]=c(c,nrow(clustSub),dotProd,vectMult)
  }
  return(as.data.frame(clustMat))
}

#' BA: Split overcrowded clusters into subclusters
#'
#' clustSplit will split clusters by iteratively removing maximum distances until "single possible candidates" are identified. The function works recursively of severely overcrowded clusters. Used internally.
#' @param clustFlags a table of presentness/missingness per batch x sample type for the features within cluster
#' @param mz m/z values of the features
#' @param rt rt values of the features
#' @return A vector with subcluster identifiers
#' @export
clustSplit=function(clustFlags,mz,rt) {
  clustList=rep(1,ncol(clustFlags))
  dotFlag=matrix(0,ncol=ncol(clustFlags),nrow=ncol(clustFlags))
  for (i in 1:ncol(clustFlags)) {
    for (j in 1:ncol(clustFlags)) {
      dotFlag[i,j]=ifelse(sum(clustFlags[,i]*clustFlags[,j])==0,1,0)
    }
  }
  rtDist=as.matrix(dist(rt,diag=TRUE,upper=TRUE))
  mzDist=as.matrix(dist(mz,diag=TRUE,upper=TRUE))
  heurFlag=dotFlag*mzDist*rtDist
  heurFlag=ifelse(heurFlag==0,NA,heurFlag)
  hFlag=heurFlag
  nDist=apply(hFlag,1,function(x) sum(!is.na(x)))
  while (min(nDist)>1) {
    hFlag[which(hFlag==max(hFlag,na.rm=TRUE))]=NA
    nDist=apply(hFlag,1,function(x) sum(!is.na(x)))
  }
  cand=which(nDist==1)
  coCand=numeric()
  heur=numeric()
  for (n in cand) {
    coCand=c(coCand,which.min(hFlag[,n]))
    heur=c(heur,min(hFlag[,n],na.rm=TRUE))
  }
  f1=cand[which.min(heur)]
  f2=coCand[which.min(heur)]
  clustList=rep(2,ncol(clustFlags))
  clustList[cand[which.min(heur)]]=clustList[coCand[which.min(heur)]]=1
  hFlagNew=heurFlag[!1:nrow(heurFlag)%in%c(f1,f2),!1:nrow(heurFlag)%in%c(f1,f2)]
  if (is.matrix(hFlagNew)) {
    nDistNew=apply(hFlagNew,1,function(x) sum(!is.na(x)))
    if (sum(nDistNew)==0) {
      clustList[clustList==2]=seq(1:sum(clustList==2))+1
    } else {
      subFeats=as.data.frame(cbind(clustList,t(clustFlags)))
      colnames(subFeats)[1]='cluster'
      subClustMat=clusterMatrix(subFeats)
      if (subClustMat$dotProd[2]!=0) {
        ### clustSplit(new matrix)
        cat('\n Warning: 2nd round of clustSplit highly experimental \n')
        newSplit=clustSplit(clustFlags[,clustList==2],mz[clustList==2],rt[clustList==2])
        clustList[clustList==2]=newSplit+1
      }
    }
  }
  return(clustList)
}

#' BA: Plot the identified clusters
#'
#' plotClust will plot clusters identified based on sample type. Used internally.
#' @param batchflag a table of presentness/missingness per batch x sample type for the features within cluster
#' @param grpFlag a vector corresponding to the rows of batchFlag used for sample type
#' @param cluster extracted cluster to be plotted
#' @param text cluster identifier
#' @param color vector with colors
#' @param mzwidth plot span of m/z
#' @param rtwidth plot span of rt
#' @export
plotClust=function(batchflag,grpFlag,cluster,text,color=2,mzwidth=0.02,rtwidth=100) {
  bF=batchflag
  mzspan=c(mean(cluster$mz)-mzwidth/2,mean(cluster$mz)+mzwidth/2)
  rtspan=c(mean(cluster$rt)-rtwidth/2,mean(cluster$rt)+rtwidth/2)
  plot(cluster$rt,cluster$mz,col=color,pch=16,xlim=rtspan,ylim=mzspan,xlab='Retention time (s)',ylab='m/z')
  text(rtspan[1],mzspan[2]-mzwidth/20,text,cex=1.3,pos=4)
  bFs=apply(cluster[,grep('bF',colnames(cluster))],1,function(x) paste(x,collapse =" "))
  text(cluster$rt,cluster$mz,cluster$featureIndex,pos=3)
  text(cluster$rt,cluster$mz,bFs,pos=1)
  whichPlot=abs(bF$mz-mean(cluster$mz))<mzwidth/2 & abs(bF$rt-mean(cluster$rt))<rtwidth/2
  mzPlot=bF$mz[whichPlot]
  rtPlot=bF$rt[whichPlot]
  bfPlot=apply(bF$flagHard[grpFlag,whichPlot],2,function(x) paste(x,collapse=' '))
  points(rtPlot,mzPlot,col=1,pch=1)
  text(rtPlot,mzPlot,which(whichPlot),pos=3)
  text(rtPlot,mzPlot,bfPlot,pos=1)
}

#' BA: Find alignment candidates
#'
#' alignIndex will find features systematically misaligned between batches using "sample type" information.
#' @param batchflag a table of presentness/missingness per batch x sample type for the features within cluster
#' @param flagType This is experimental. Use at your own risk.
#' @param grpType sample type to be used to find batch alignments
#' @param mzdiff maximum distance in m/z to be considered for alignment
#' @param rtdiff maximum distance in rt to be considered for alignment
#' @param report boolean for plotting results of owercrowded cluster splitting for visual examination
#' @param reportName string with name of report file
#' @return An object (list) consisting of the following features:
#' @return grpType: sample type used to find batch alignments
#' @return shift: list with features to combine
#' @return events: all possible alignment candidates (events)
#' @return features: features after splitting overcrowded clusters
#' @return clusters: clusters after splitting overcrowded clusters
#' @return oldFeatures: features before splitting overcrowded clusters
#' @return oldClusters: clusters before splitting overcrowded clusters
#' @export
alignIndex=function(batchflag,flagType=c('Hard','Soft','All'),grpType='Q',mzdiff=0.005,rtdiff=10,report=TRUE,reportName='splits') {
  bF=batchflag
  if (missing(flagType)) flagType='Hard'
  if (flagType=='Hard') grpSub=bF$flagHard[bF$meta[,2]==grpType,] # Take out matrix based on group type (such as Q for QC or R for Ref)
  if (flagType=='Soft') grpSub=bF$flagSoft[bF$meta[,2]==grpType,] # Take out matrix based on group type (such as Q for QC or R for Ref)
  if (flagType=='All') grpSub=bF$flagAll[bF$meta[,2]==grpType,] # Take out matrix based on group type (such as Q for QC or R for Ref)
  mz=bF$mz
  rt=bF$rt
  a2=a1=align(grpSub,mz,rt,mzdiff=mzdiff,rtdiff=rtdiff)
  a1Clust=a1$clusters
  splits=which(a1Clust$dotProd!=0)
  if (report==TRUE) pdf(file=paste(reportName,'.pdf',sep=''))
  for (s in splits) {
    cluster=a1$features[a1$features$cluster==s,]
    # print(cluster)
    mzSp=cluster$mz
    rtSp=cluster$rt
    flagSp=t(cluster[9:ncol(a1$features)])
    alignSplit=clustSplit(flagSp,mzSp,rtSp)
    text=paste('Original cluster',s)
    if (report==TRUE) plotClust(bF,bF$meta[,2]==grpType,cluster,text=text,color=alignSplit+1)
    newClust=ifelse(alignSplit==1,s,max(a2$features$cluster)+alignSplit-1)
    a2$features$cluster[a2$features$cluster==s]=newClust
  }
  if (report==TRUE) dev.off()
  ### Remove duplicates
  a2$features=a2$features[a2$features$cluster%in%unique(a2$features$cluster[duplicated(a2$features$cluster)]),]
  a2$features=a2$features[order(a2$features$cluster),]
  a2$features$cluster=as.numeric(as.factor(a2$features$cluster))
  a2$clusters=clusterMatrix(a2$features)
  ### Make index string for alignment
  nClust=nrow(a2$clusters)
  shiftList=numeric(ncol(bF$flagHard))
  shiftGrp=character(ncol(bF$flagHard))
  for (c in 1:nClust) {
    feats=a2$features$featureIndex[a2$features$cluster==c]
    shiftList[feats]=feats[1]
    shiftGrp[feats]=as.character(grpType)
    shift=data.frame(list=shiftList,grp=shiftGrp)
  }
  return(list(grpType=grpType,shift=shift,events=a2$events,features=a2$features,clusters=a2$clusters,oldFeatures=a1$features,oldClusters=a1$clusters))
}

#' BA: Plot clusters of aligned features
#'
#' plotAlign will plot clusters of aligned features
#' @param batchflag a table of presentness/missingness per batch x sample type for the features within cluster
#' @param alignindex An object (list) consisting of alignment information
#' @param clust which cluster(s) to plot. If missing, plots all clusters.
#' @param plotType Whether to plot internally ('plot'), to pdf ('pdf') or png ('png') file.
#' @param reportName name of plotfile
#' @param mzwidth plot span of m/z
#' @param rtwidth plot span of rt
#' @export
plotAlign=function(batchflag,alignindex,clust,plotType=c('plot','pdf','png'),reportName='clustPlot',mzwidth=0.02,rtwidth=100) {
  if (missing(plotType)) plotType='plot'
  aI=alignindex
  bF=batchflag
  if (plotType=='pdf') pdf(file=paste(reportName,'.pdf',sep=''))
  if (plotType=='png') png(filename=paste(reportName,'.png',sep=''))
  if (missing(clust)) {
    clustPlots=1:dim(aI$clusters)[1]
  } else clustPlots=clust
  # cat(clustPlots)
  if (length(clustPlots)==2) par(mfrow=c(2,1))
  if (length(clustPlots)>2) par(mfrow=c(3,2))
  for (c in clustPlots) {
    cluster=aI$features[aI$features$cluster==c,]
    # print(cluster)
    mz=cluster$mz
    rt=cluster$rt
    flags=t(cluster[9:ncol(aI$features)])
    text=paste('Cluster',c)
    plotClust(batchflag=bF,grpFlag=bF$meta[,2]==aI$grpType,cluster=cluster,text=text,color=2)
  }
  if (plotType=='pdf' | plotType=='png') dev.off()
}

#' BA: Aggregate clustering information from two sample types
#'
#' aggregateIndex will aggregate alignment indices IF present for both sample types
#' @param aI1 an alignIndex object for sample type 1
#' @param aI2 an alignIndex object for sample type 2
#' @return An object (list) with aggregated alignment information
#' @export
aggregateIndex=function(aI1,aI2) {
  list1=aI1$shift$list
  list2=aI2$shift$list
  lists=rbind(list1,list2)
  grp1=as.character(aI1$shift$grp)
  grp2=as.character(aI2$shift$grp)
  grps=rbind(grp1,grp2)
  ## Find conflicts
  conflict=apply(lists,2,function(x) ifelse(x[1]==x[2],F,ifelse(any(x==0),F,T)))
  confClust=unique(as.numeric(lists[,conflict]))
  lists[lists%in%confClust]=0
  listNew=ifelse(list1==list2,list1,0)
  grpNew=ifelse(list1==list2,paste(grp1,grp2,sep=''),'')
  shift=data.frame(list=listNew,grp=grpNew)
  return(list(shift=shift,shift1=aI1$shift,shift2=aI2$shift))
}

#' BA: Alignment of peaktable based on alignIndex and batchFlag data
#'
#' batchAlign will use the 'batchFlag' and 'alignIndex' information to align sample peaks that are systematically misaligned across batches.
#' @param batchflag a table of presentness/missingness per batch x sample type for the features within cluster
#' @param alignindex An object (list) consisting of alignment information
#' @param peaktable_filled a peaktable without missing values (ie after hard filling or imputation)
#' @param meta a matrix of metadata with batch in column 1 and sample type in column 2
#' @return An object (list) consisting of the following features:
#' @return PTalign: A peaktable with batch aligned peak areas
#' @return boolAveragedAlign: boolean vector of features where alignment has been made using feature averaging (i.e. where batches are missing within features). Length: same as final number of features
#' @return PTfill: peaktable without missing values (indata)
#' @return boolKeep: boolean vector of features kept after alignment (the rest of the combined features are deleted)
#' @return boolAveragedFill: boolean vector of features where alignment has been made using feature averaging (i.e. where batches are missing within features). Length: same as orgiginal number of features
#' @return aI: alignIndex object (indata)
#' @return meta: metadata (indata)
#' @export
batchAlign=function(batchflag,alignindex,peaktable_filled,meta) {
  bF=batchflag
  flags=bF$flagHard
  batchFlags=as.matrix(aggregate(flags,list(bF$meta[,1]),'sum')[,-1])
  batch=bF$meta[,1]
  uniqBatch=unique(batch)
  aI=alignindex
  shift=aI$shift$list
  uniqSh=unique(shift[shift!=0])
  grp=aI$shift$grp
  PTfill=peaktable_filled
  boolKeep=rep(TRUE,ncol(PTfill))
  boolAveraged=!boolKeep
  for (sh in uniqSh) {
    newFeat=numeric(nrow(PTfill))
    feats=which(shift==sh)
    subVect=flags[,feats]
    if (is.matrix(subVect)) {
      vectAdd=subVect[,1]
      subBatchVect=batchFlags[,feats]
      batchVectAdd=subBatchVect[,1]
      subFeats=PTfill[,feats]
      aveFeat=rowSums(t(colSums(subVect)*t(subFeats)))/sum(subVect)
      for (f2 in 2:length(feats)) {
        vect2=subVect[,f2]
        dotProd=sum(vectAdd*vect2)
        vectAdd=vectAdd+vect2
        batchVect2=subBatchVect[,f2]
        batchDotProd=sum(batchVectAdd*batchVect2)
        batchVectAdd=batchVectAdd+batchVect2
      }
      if (dotProd==0 & batchDotProd==0) {
        for (f in feats) {
          boolFeat=rep(FALSE,nrow(PTfill))
          bFlag=which(batchFlags[,f]>0)
          for (bFl in bFlag) {
            batchID=uniqBatch[bFl]
            boolFeat=boolFeat|ifelse(meta[,1]==batchID,T,F)
          }
          newFeat[boolFeat]=PTfill[boolFeat,f]
        }
        ### How to deal with non-present batches???
        zeros=which(rowSums(subBatchVect)==0)
        if (length(zeros)>0) {
          boolAveraged[feats]=TRUE
          boolFeat=rep(FALSE,nrow(PTfill))
          for (z in zeros) {
            batchID=uniqBatch[z]
            boolFeat=boolFeat|ifelse(meta[,1]==batchID,T,F)
          }
          newFeat[boolFeat]=aveFeat[boolFeat]
        }
        PTfill[,feats[1]]=newFeat
        PTfill[,feats[-1]]=0
        boolKeep[feats[-1]]=FALSE
      }
    }
  }
  PTfill=PTfill[,boolKeep]
  return(list(PTalign=PTfill, boolAveragedAlign=boolAveraged[boolKeep], PTfill=peaktable_filled, boolKeep=boolKeep, boolAveragedFill=boolAveraged, aI=alignindex, meta=meta))
}
