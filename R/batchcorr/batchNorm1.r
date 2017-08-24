#' Extract m/z and rt from peak table
#'
#' Extract features from peak table and report their m/z and rt values
#' @param PT a peak table with variables as columns
#' @param start character from which to start the read of peakInfo (from PT colnames)
#' @return a matrix with m/z and rt of features as columns
#' @export
peakInfo=function(PT,start=3) {
  peakInfo=matrix(unlist(strsplit(colnames(PT),'@')),ncol=2,byrow=TRUE)
  peakInfo[,1]=substr(peakInfo[,1],start,max(nchar(peakInfo[,1])))
  peakInfo=matrix(as.numeric(peakInfo),ncol=2)
  colnames(peakInfo)=c('mz','rt')
  rownames(peakInfo)=paste('feature',1:nrow(peakInfo),sep='_')
  return(peakInfo)
}

#' Limit variables/features to those common between batches
#'
#' Extract features from peak table and report their m/z and rt values
#' @param batchFeats list with feature names per batch
#' @param batchLimit lower limit of number of batches in which the feature needs to be present to be included in the final list (defaults to length of `batchFeats` list)
#' @return a vector with the features present in at least `batchLimit` batches
#' @export
featComb=function(batchFeats,batchLimit) {
  nBatch=length(batchFeats)
  if (missing(batchLimit)) batchLimit=nBatch
  vars=batchFeats[[1]]
  for (b in 2:nBatch) {
    vars=c(vars,batchFeats[[b]])
  }
  varTab=table(vars) #Tabularize number of occurrences of features
  finalVars=t(as.data.frame(varTab[varTab>=batchLimit])) # Bring out the features present >= batchLimit times NB! Not correctly sorted!
  featInfo=peakInfo(finalVars) # Bring out mz and rt to sort them properly
  finalVars=colnames(finalVars)  # take names of features only
  finalVars=finalVars[order(featInfo[,1],featInfo[,2])]  # Sort them according to mz and rt
  return(finalVars)
}


#' Extract features from multiple batch data
#'
#' Extract features present in all batches and combine them into a master peaktable.
#' @param batchObjs a list with batch objects (from within batch drift correction)
#' @param batchLimit lower limit of number of batches in which the feature needs to be present to be included in the final list (defaults to length of `batchFeats` list)
#' @param finalVars a vector with the features to bring out from batch PTs. NB! All features need be present in all batches (defaults to all features present in all batches)
#' @return a combined, but NOT normalised, peaktable - limited to common features
#' @export
batchComb=function(batchObjs,batchLimit,finalFeats) {
  nBatch=length(batchObjs)
  if (missing(batchLimit)) batchLimit=nBatch
  if (missing(finalFeats)) {
    batchFeats=list()
    for (b in 1:nBatch) {
      batchFeats[[b]]=batchObjs[[b]]$finalVars
    }
    finalFeats=featComb(batchFeats,batchLimit)
  }
  PTComb=subset(batchObjs[[1]]$TestFeatsFinal,select=finalFeats)
  for (b in 2:nBatch) {
    PTComb=rbind(PTComb,subset(batchObjs[[b]]$TestFeatsFinal,select=finalFeats))
  }
  return(PTComb)
}

#' BN: Info on reference samples aggregated on batch level
#'
#' Reference samples are aggregated on batch level
#' @param PT a multi-batch master peak table
#' @param meta metadata with batch (col1) and sample type (col2)
#' @param grpType a sample type identifier for reference samples
#' @param CVlimit CV criterion to pass for Ref samples per batch
#' @return an object containing:
#' @return CV: boolean per batch & feature in CV<limit
#' @return aveInt: average reference intensity per batch & feature
#' @export
refOut=function(PT,meta,grpType='R',CVlimit=0.3) {
  batch=meta[,1]
  grp=meta[,2]
  uniqBatch=unique(batch)
  CVMat=matrix(nrow=length(uniqBatch),ncol=ncol(PT))
  rownames(CVMat)=uniqBatch
  aveIntMat=matrix(nrow=length(uniqBatch),ncol=ncol(PT))
  rownames(aveIntMat)=uniqBatch
  for (b in 1:length(uniqBatch)) {
    bat=uniqBatch[b]
    PTbatch=PT[batch==bat & grp==grpType,]
    CVMat[b,]=ifelse(cv(PTbatch)<=CVlimit,TRUE,FALSE)
    aveIntMat[b,]=apply(PTbatch,2,mean)
  }
  return(list(CV=CVMat,aveInt=aveIntMat))
}

#' BN: Between batch normalisation by Ref samples
#'
#' Batches are normalised by Ref samples if passing heuristic distance criterion.
#' @param PT a multi-batch master peak table
#' @param meta metadata with batch (col1) and sample type (col2)
#' @param refs a refOut object
#' @param FCLimit Fold-change criterion for intensity (in relation to average intensity FC between batches)
#' @return an object containing:
#' @return PTRef: Reference sample-normalised multi-batch peak table
#' @return refCorr: Boolean matrix with info on which batches were normalised by reference samples
#' @return PTOrg (indata peaktable)
#' @export
refCorr=function (PT, meta, refs, FCLimit = 5){
  batch = meta[, 1]
  PTcorr = PT
  cvm = refs$CV
  nBatch = nrow(cvm)
  uniqBatch = rownames(cvm)
  nFeat = ncol(cvm)
  refCorrMat = matrix(FALSE, nrow = nBatch, ncol = nFeat)
  cvFlags = apply(cvm, 2, function(x) sum(x, na.rm = TRUE))
  aveInt = refs$aveInt
  meanIntRat = matrix(1, nrow = nBatch, ncol = nBatch)
  for (b in 1:(nBatch - 1)) {
    for (bb in (b + 1):nBatch) {
      meanIntRat[bb, b] = mean(aveInt[bb, ])/mean(aveInt[b,
                                                         ])
      meanIntRat[b, bb] = 1/meanIntRat[bb, b]
    }
  }
  whichFeatsCV = which(cvFlags > 1)
  lenCV = length(whichFeatsCV)
  for (lc in 1:lenCV) {
    feat = whichFeatsCV[lc]
    featIntRat = matrix(1, nrow = nBatch, ncol = nBatch)
    for (b in 1:(nBatch - 1)) {
      for (bb in (b + 1):nBatch) {
        featIntRat[bb, b] = aveInt[bb, feat]/aveInt[b,feat]
        featIntRat[b, bb] = 1/featIntRat[bb, b]
      }
    }
    featFlags = abs(log(featIntRat/meanIntRat)) <= log(FCLimit)
    featFlags[which(is.na(featFlags),arr.ind=T)]=FALSE
    for (b in 1:nBatch) {
      if (cvm[b, feat] == FALSE | is.na(cvm[b, feat])) {
        featFlags[, b] = featFlags[b, ] = FALSE
      }
    }
    if (any(colSums(featFlags) > 1)) {
      refBatch = min(which(colSums(featFlags) == max(colSums(featFlags))))
      refCorr = featFlags[, refBatch]
      WhichRefCorr = which(refCorr == TRUE)
      refInt = aveInt[refBatch, feat]
      refCorrIndex=WhichRefCorr[!WhichRefCorr==refBatch]
      for (b in refCorrIndex) {
        corrFact = refInt/aveInt[b, feat]
        PTcorr[batch == uniqBatch[b], feat] = PT[batch == uniqBatch[b], feat] * corrFact
      }
      featCorr=PTcorr[batch %in% uniqBatch[WhichRefCorr], feat]
      Median=median(featCorr)
      WhichPOPCorr=which(refCorr==FALSE)
      for (n in WhichPOPCorr) {
        corrFact=Median/median(PTcorr[batch == uniqBatch[n], feat])
        PTcorr[batch == uniqBatch[n], feat] = PTcorr[batch == uniqBatch[n], feat] * corrFact
      }
    }
    refCorrMat[, feat] = refCorr
  }
  return(list(PTRef = PTcorr, refCorr = refCorrMat, PTOrg = PT))
}
