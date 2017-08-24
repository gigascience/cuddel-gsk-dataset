#' DC: Grab QC samples
#'
#' Bring out 'QC' group from an XCMS object: Scaled and 100% NAs removed
#' #' **Note** that `inj` needs to be either **manually specified OR recoded** to suit your sample naming strategy.
#' @param XS An XCMS object
#' @param peakTab peaktable. Defaults to `xcms::peakTab(SX)`. Specify `PTalign` from `batchAlign` object for multibatch aligned data.
#' @param batch a batch identifier
#' @param grp an identifyer for the QC samples
#' @param inj a vector of QC injection numbers in sequence. Defaults to defined function to suit my naming strategy. Pls change accordingly.
#' @return a list containing:
#' @return batch: batch (as input)
#' @return inj: a vector of QC injection numbers in sequence (same as indata OR generated)
#' @return Feats: Scaled features (peak table)
#' @return RawFeats: Unscaled features including 100% NA features (peak table)
#' @return RawFeatsNaRm: Unscaled features excluding 100% NA features (peak table)
#' @return NAs: which features are missing from QC samples
#' @export
## Bring out 'QC' group: Scaled and NAs removed
grabQC=function(XS,peakTab,batch,grp='QC',inj) {
	incl=(XS@phenoData[,1]==batch & XS@phenoData[,2]==grp)
	if (missing(peakTab)) peakTab=peakTab(XS)
	QC=peakTab[incl,]
	QCCV=cv(QC)
	QCscale=scale(QC,center=FALSE)
	NAs=colSums(is.na(QCscale))>0
	QCRawNaRm=QC[,!NAs]
	QCFeats=QCscale[,!NAs]
	if (missing(inj)) inj=as.numeric(substr(matrix(unlist(strsplit(rownames(QC),'_')),ncol=6,byrow=TRUE)[,6],1,3))
	return(list(batch=batch,inj=inj,Feats=QCFeats,RawFeats=QC,RawFeatsNaRm=QCRawNaRm,NAs=NAs))
}

#' DC: Grab Reference samples from same batch as QC samples
#'
#' Bring out reference sample group from an XCMS object: Exclude features missing from QC samples
#' @param XS An XCMS object
#' @param peakTab peaktable. Defaults to `xcms::peakTab(SX)`. Specify `PTalign` from `batchAlign` object for multibatch aligned data.
#' @param QC a grabQC object
#' @param grp an identifyer for the reference samples
#' @param inj a vector of Ref injection numbers in sequence. Defaults to defined function to suit my naming strategy. Pls change accordingly.
#' @return a list containing:
#' @return inj: a vector of reference sample injection numbers in sequence (same as indata OR generated)
#' @return Feats: Features (peak table)
#' @export
## Bring out Reference samples from same batch as QCs
grabRef=function(XS,peakTab,QC,grp='Ref',inj) {
  incl=(XS@phenoData[,1]==QC$batch & XS@phenoData[,2]==grp)
  if (missing(peakTab)) peakTab=peakTab(XS)
  Feats=peakTab[incl,!QC$NAs]
	if (missing(inj)) inj=as.numeric(substr(matrix(unlist(strsplit(rownames(Feats),'_')),ncol=6,byrow=TRUE)[,6],1,3))
	return(list(inj=inj,Feats=Feats))
}

#' DC: Grab Batch samples
#'
#' Bring out all batch samples from an XCMS object: Exclude features missing from QC samples
#' @param XS An XCMS object
#' @param peakTab peaktable. Defaults to `xcms::peakTab(SX)`. Specify `PTalign` from `batchAlign` object for multibatch aligned data.
#' @param QC a grabQC object
#' @param inj a vector of batch injection numbers in sequence. Defaults to defined function to suit my naming strategy. Pls change accordingly.
#' @return a list containing:
#' @return inj: a vector of sample injection numbers in sequence (same as indata OR generated)
#' @return Feats: Features (peak table)
#' @export
## Bring out entire batch same as QCs
grabBatch=function(XS,peakTab,QC,inj) {
  incl=(XS@phenoData[,1]==QC$batch)
  if (missing(peakTab)) peakTab=peakTab(XS)
  Feats=peakTab[incl,!QC$NAs]
	if (missing(inj)) inj=as.numeric(substr(matrix(unlist(strsplit(rownames(Feats),'_')),ncol=6,byrow=TRUE)[,6],1,3))
	return(list(inj=inj,Feats=Feats))
}

#' Coefficient of variation (CV)
#'
#' Calculates CV per column in a matrix. Aka relative standard deviation (RSD).
#' @param mat a matrix with variables as columns
#' @return a vector of CVs
#' @export
## Simple function for calculating cv per column (ie variable)
cv=function(mat) {
  if (is.null(dim(mat))) {
    cv=sd(mat)/mean(mat)
  } else {
    mean=apply(mat,2,mean)
    sd=apply(mat,2,sd)
    cv=sd/mean
  }
	return(cv)
}

#' Root mean squared distance
#'
#'Calculate root mean squared distance from center point
#' @param mat a matrix containing observations as rows and variables as columns
#' @return a numeric rmsDist
#' @export
## Simple function for calculating root mean square distance
rmsDist=function(mat) {
	mean=colMeans(mat)
	rmsd=sqrt(sum(apply(mat,1,function(x) sum((x-mean)^2)))/nrow(mat))
	return(rmsd)
}

#' DC: Cluster features
#'
#' clust will perform clustering of features with similar drift pattern by projecting scaled QC features as coordinates in observation/injection space
#' @param QCInjs vectors with QC injections in sequence
#' @param QCFeats Feature matrix for QC injections: Scaled (but not centered) features as columns; Injections as rows
#' @param modelNames Which 'mclust' models to consider (see mclust package for details)
#' @param G Which number of clusters to consider (see mclust package for details)
#' @param report boolean whether to print a pdf report of clustering results
#' @return a clust object containing:
#' @return QCInjs: as indata
#' @return QCFeats: as indata
#' @return BIC: Bayesian Information Criteria-object for the investigated clustering models
#' @return clust: summary of the BIC-object (i.e. clustering)
#' @return BICTime: Time required for the BIC-calculation
#' @return clustTime: Time required for the clustering
#' @export
## Perform clustering
clust=function(QCInjs,QCFeats,modelNames=c('VVE'),G=seq(1,52,by=3),report=FALSE) {
	# modelNames='VVV'
	# modelNames=c('VEV','VVV')
	# modelNames=c('VEE','VEV','VVE','VVV')
	# modelNames=c('VEE','VVE')
	# modelNames=c('VII','VEI','VVI','VEE','VEV','VVE','VVV')
	# modelNames=c('EEE','EEV','EVE','EVV','VEE','VEV','VVE','VVV')
  library(mclust)
	startTime=proc.time()[3]
	mclBIC=mclustBIC(t(QCFeats),G=G,modelNames=modelNames)
	endTime=proc.time()[3]
	BICtime=endTime-startTime
	startTime=proc.time()[3]
	MC=summary(mclBIC,data=t(QCFeats))
	endTime=proc.time()[3]
	sumTime=endTime-startTime
	if (report==TRUE) {
		pdf(file=paste('cluster_BIC_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
		plot(mclBIC)
		dev.off()
	}
	return(list(QCInjs=QCInjs,QCFeats=QCFeats,BIC=mclBIC,clust=MC,BICTime=BICtime,clustTime=sumTime))
}

#' DC: Calculate intensity drift per cluster
#'
#' Clustered, scaled QC features are pooled together and intensity drift patterns are modelled per cluster
#' @param QCClust a clust object
#' @param smoothFunc choice of regression function: spline or loess (defaults to spline)
#' @param spar smoothing parameter for spline or loess
#' @param report boolean whether to print a pdf report of drift models
#' @return a driftCalc object containing:
#' @return original information from clust object (indata) and
#' @return actionInfo
#' @return ratios
#' @return corMat
#' @return deltaDist
#' @return varClust
#' @export
## Calculate drift clusters
driftCalc=function(QCClust,smoothFunc=c('spline','loess'),spar=0.2,report=FALSE) {
  library(reshape)
	if (missing(smoothFunc)) smoothFunc='spline'
	MC=QCClust$clust
	QCInjs=QCClust$QCInjs
	QCFeats=QCClust$QCFeats
	# Extract classes
		nclass=MC$G # Take out total number of identified clusters/components/groups/classes/whatever you want to call them
		classes=MC$classification # Take out the classifications for the different variables
	# Allocate variables
		cvRaw=cvCorr=deltaDist=numeric(nclass) # allocate vector for effect size of drift correction per cluster
		injs=min(QCInjs):max(QCInjs) # Make injection list
		corMat=matrix(nrow=length(injs),ncol=nclass) # Allocate matrix with correction function (rows) per cluster (column)
		cvs=varClust=list()
		ratios=matrix(nrow=nclass,ncol=4)
		rownames(ratios)=paste('cluster',1:nclass,sep='')
		colnames(ratios)=c('raw.15','corr.15','raw.2','corr.2')
	if (report==TRUE) {
		pdf(file=paste('cluster_G',nclass,'_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
		par(mfrow=c(2,1))
		par(mar=c(2,4,2,0))
	}
	### Calculate drift correction for each cluster
	# Calculate distance on scaled variables
	rmsdRaw=rmsDist(QCFeats)
		for (n in 1:nclass) {
			QCFeatsCorr=QCFeats # Allocate matrix (for drift corrected variables) for later QC distance calculation
			whichVars=which(classes==n)
			vars=QCFeats[,whichVars] # Take out cluster variables
			varClust[[n]]=colnames(QCFeats)[whichVars]
			V=as.data.frame(cbind(QCInjs,vars)) # Arrange and rearrange data
			V=melt(V,id.vars='QCInjs')
			V=V[order(V$QCInjs),]
			# Interpolations
				if (length(QCInjs)<=3) {  # 2nd degree polynomial if <4 data points
					Fit=lm(value ~ poly(QCInjs,2),data=V)
					Pred=predict(Fit,data.frame(QCInjs=injs))
					Pred=data.frame(x=injs,y=Pred)
				} else {
					if (smoothFunc=='spline') {
						splineFit=smooth.spline(V$QCInjs,V$value,spar=spar) # Cubic spline regression otherwise
						Pred=predict(splineFit,injs)  # Predict drift over all injections
					} else {
						loessFit=loess(value~QCInjs,data=V,span=spar)
						Pred=predict(loessFit,data.frame(QCInjs=injs))
						Pred=data.frame(x=injs,y=Pred)
					}
				}
			corFact=Pred$y[1]/Pred$y  # Calculate correction factors for all injections
			corMat[,n]=corFact # Store cluster correction factors in "master" matrix
			corQC=corFact[QCInjs-min(QCInjs)+1]  # Bring out correction factors for QC samples specifically
			QCFeatsCorr[,classes==n]=QCFeats[,classes==n]*corQC  # correct drift within cluster
			## Calculate rmsDist
			rmsdCorr=rmsDist(QCFeatsCorr)
			deltaDist[n]=rmsdCorr-rmsdRaw
			# deltaDist[n]=mean(dist(QCFeatsCorr))-meanDistQCFeats # Calculate change in average QC distance
			cvRaw[n]=mean(cv(QCFeats[,classes==n]))
			cvCorr[n]=mean(cv(QCFeatsCorr[,classes==n]))
			cvs[[n]]=data.frame(Raw=cv(QCFeats[,classes==n]),Corr=cv(QCFeatsCorr[,classes==n]))
			ratios[n,]=c(sum(cvs[[n]]$Raw<0.15)/nrow(cvs[[n]]),sum(cvs[[n]]$Corr<0.15)/nrow(cvs[[n]]),sum(cvs[[n]]$Raw<0.2)/nrow(cvs[[n]]),sum(cvs[[n]]$Corr<0.2)/nrow(cvs[[n]]))
			if (report==TRUE) {
				# Plot drift and drift function
				matplot(QCInjs,QCFeats[,classes==n],type='l',lty=1,col='grey',ylim=range(QCFeats[,classes==n]),main=paste('Cluster ',n,'; n=',sum(classes==n),'; Raw; Mean CV=',round(mean(cv(QCFeats[,classes==n])),3),sep=''),ylab='Scaled intensity',xlab='Injection number')
				lines(Pred,pch=2)
				matplot(QCInjs,QCFeatsCorr[,classes==n],type='l',lty=1,col='grey',ylim=range(QCFeats[,classes==n]),main=paste('Corrected; Mean CV=',round(mean(cv(QCFeatsCorr[,classes==n])),3),sep=''),ylab='Scaled intensity',xlab='Injection number')
			}
		}
	if (report==TRUE) dev.off() # Close pdf file
	clustComm=rep('None',nclass)
	actionInfo=data.frame(number=1:nclass,n=sapply(varClust,length),action=clustComm,CVRaw=cvRaw,CVCorr=cvCorr)
	QCClust$actionInfo=actionInfo
	QCClust$ratios=ratios
	QCClust$corMat=corMat
	QCClust$deltaDist=deltaDist
	QCClust$varClust=varClust
	QCDriftCalc=QCClust
	return(QCDriftCalc)
}

#' DC: Correct for intensity drift per cluster
#'
#' Perform signal intensity drift correction if resulting in increased quality of data (measured by reduced root-mean-squared distances of reference samples).
#' @param QCDriftCalc a DriftCalc object
#' @param refList a reference sample onject with 'inj' and 'Feats'
#' @param refType at present, the options "one" and "none" are supported for one or no reference samples present.
#' @param CorrObj a batch object to be corrected with 'inj' and 'Feats'. If not present defaults to the QC object
#' @param report boolean whether to print a pdf report of drift models
#' @return A Corr object consisting of:
#' @return original information from DriftCalc object (indata) and
#' @return updated actionInfo (actions taken per cluster)
#' @return QCFeatsCorr Drift-corrected QC features
#' @return refType (indata)
#' @return drift corrected reference samples
#' @return drift corrected corr/batch samples
#' @export
## Perform drift correction for clusters IF rmsdRef is improved
driftCorr=function(QCDriftCalc,refList=NA,refType=c('none','one','many'),CorrObj=NA,report=FALSE) {
	if (missing(refType)) refType='none'
	if (refType=='many') {
		cat('\nMultiple reference samples not yet implemented\n')
		break
	}
	deltaDist=QCDriftCalc$deltaDist
	varClust=QCDriftCalc$varClust
	removeFeats=QCDriftCalc$removeFeats
	if (!is.null(removeFeats)) {
	  keepClust=QCDriftCalc$keepClust
	  corrQCTemp=corrQC=QCDriftCalc$QCFeatsClean
	} else {
	  keepClust=1:(QCDriftCalc$clust$G)
	  corrQCTemp=corrQC=QCDriftCalc$QCFeats
	}
	injQC=QCDriftCalc$QCInjs
	corMat=QCDriftCalc$corMat
	clustComm=as.character(QCDriftCalc$actionInfo$action)
	ordDist=order(deltaDist)
	ordDist=ordDist[ordDist%in%keepClust]
	if (refType=='one') {
  	refClean=refList$Feats[,!colnames(refList$Feats)%in%removeFeats]
	  injRef=refList$inj
	  corrRefTemp=corrRef=refClean
	}
	if (missing(CorrObj)) {
	  injTest=injQC
	  corrTest=corrQC
	  CorrObj=list(inj=injTest,Feats=corrTest)
	} else {
	  injTest=CorrObj$inj
	  corrTest=CorrObj$Feats
	  corrTest=corrTest[,!colnames(corrTest)%in%removeFeats]
	}
	for (i in 1:length(keepClust)) {
		n=ordDist[i]
		corFact=corMat[,n] # take out cluster correction factors from "master" matrix
		corrFeats=varClust[[n]]
		if (refType=='none') { #Scheme for (suboptimal) situation without Ref samples
		  corQC=corFact[injQC-min(injQC)+1]  # Bring out correction factors for QC samples specifically
		  corrQCTemp[,colnames(corrQC)%in%corrFeats]=corrQC[,colnames(corrQC)%in%corrFeats]*corQC
		  if (rmsDist(corrQCTemp)<rmsDist(corrQC)) {
		    clustComm[n]='Corr_QC'
		    corrQC=corrQCTemp
		    corQC=corFact[injQC-min(injQC)+1]  # Bring out correction factors for QC samples specifically
		    corrQC[,colnames(corrQC)%in%corrFeats]=corrQC[,colnames(corrQC)%in%corrFeats]*corQC
		    corTest=corFact[injTest-min(injQC)+1]  # Bring out correction factors for Test samples specifically
		    corrTest[,colnames(corrTest)%in%corrFeats]=corrTest[,colnames(corrTest)%in%corrFeats]*corTest
		  } else {
		    corrQCTemp=corrQC
		  }
		}
		if (refType=='one') { # Scheme for situation with multiple injections of singular non-QC Ref samples
		  corRef=corFact[injRef-min(injQC)+1]
		  corrRefTemp[,colnames(corrRefTemp)%in%corrFeats]=corrRefTemp[,colnames(corrRefTemp)%in%corrFeats]*corRef
  		if (rmsDist(corrRefTemp)<rmsDist(corrRef)) {
  			clustComm[n]='Corr_1Ref'
  			corrRef=corrRefTemp
  			corQC=corFact[injQC-min(injQC)+1]  # Bring out correction factors for QC samples specifically
  			corrQC[,colnames(corrQC)%in%corrFeats]=corrQC[,colnames(corrQC)%in%corrFeats]*corQC
  			corTest=corFact[injTest-min(injQC)+1]  # Bring out correction factors for Test samples specifically
  			corrTest[,colnames(corrTest)%in%corrFeats]=corrTest[,colnames(corrTest)%in%corrFeats]*corTest
  		} else {
  			corrRefTemp=corrRef
  		}
		}
	}
	if (report==TRUE) {
		pdf(file=paste('Hist_Corrected_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
	  if (!is.null(removeFeats)) {
	    hist(cv(QCDriftCalc$QCFeatsClean),30,col=rgb(0,0,0,1),main='Cluster correction',xlab='CV (feature)')
	  } else {
	    hist(cv(QCDriftCalc$QCFeats),30,col=rgb(0,0,0,1),main='Cluster correction',xlab='CV (feature)')
	  }
		hist(cv(corrQC),20,col=rgb(1,1,1,.5),add=TRUE)
		legend('topright',legend=c('Clean','Corrected'),fill=c(rgb(0,0,0,1),rgb(1,1,1,0.5)))
		dev.off()
	}
	QCDriftCalc$actionInfo$action=clustComm
	QCDriftCalc$QCFeatsCorr=corrQC
	QCDriftCalc$RefType=refType
	if (refType=='one') {
  	QCDriftCalc$RefInjs=injRef
  	QCDriftCalc$RefFeats=refList$Feats
  	QCDriftCalc$RefFeatsClean=refClean
  	QCDriftCalc$RefFeatsCorr=corrRef
	}
	QCDriftCalc$TestInjs=injTest
	QCDriftCalc$TestFeats=CorrObj$Feats
	# QCDriftCalc$TestFeatsClean=
	QCDriftCalc$TestFeatsCorr=corrTest
	QCCorr=QCDriftCalc
	return(QCCorr)
}

#' DC: Remove features not passing QC test
#'
#' Remove those features with CV(QC)>limit (defined by user)
#' @param QCCorr a Corr object
#' @param CVlimit user-defined QC riterion to pass. Defaults to 0.2
#' @param report boolean whether to print a pdf report of drift models
#' @return A Clean object consisting of:
#' @return original information from Corr object (indata) and
#' @return updated actionInfo (actions taken per cluster)
#' @return QCFeatsFinal: Final peaktable (scaled) for QC samples
#' @return TestFeatsFinal: Final peaktable for batch samples
#' @return finalVars: All features kept within the final peaktable
#' @return QCcvs: CVs per features kept within the final peaktable
#' @export
## Remove individual variables with CV>limit
cleanVar=function(QCCorr,CVlimit=.2,report=FALSE) {
	QCFeats=QCCorr$QCFeats
	removeFeats=QCCorr$removeFeats
	if (!is.null(removeFeats)) {
	  QCFeatsClean=QCCorr$QCFeatsClean
	} else {
	  QCFeatsClean=QCCorr$QCFeats
	}
	QCFeatsCorr=QCCorr$QCFeatsCorr
	cvIndex=which(cv(QCFeatsCorr)>CVlimit)
	QCFeatsFinal=QCFeatsCorr[,-cvIndex]
	if (QCCorr$RefType=='one') {
	  RefFeatsFinal=QCCorr$RefFeatsCorr[,-cvIndex]
	}
	TestFeatsFinal=QCCorr$TestFeatsCorr[,-cvIndex]
	finalVars=colnames(QCFeatsFinal)
	cvFeats=mean(cv(QCFeats))      # 0.2276
	cvFeatsClean=mean(cv(QCFeatsClean))    # 0.1455
	cvFeatsCorr=mean(cv(QCFeatsCorr))   # 0.1144
	cvFeatsFinal=mean(cv(QCFeatsFinal)) # 0.1070
	QCcvs=data.frame(cvFeats=cvFeats,cvFeatsClean=cvFeatsClean,cvFeatsCorr=cvFeatsCorr,cvFeatsFinal=cvFeatsFinal)
	if (report==TRUE) {
		pdf(file=paste('Hist_Final_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
		hist(cv(QCFeatsClean),30,col=rgb(0,0,0,.5),main='Cluster cleanup',xlab='CV (feature)')
		hist(cv(QCFeatsFinal),5,col=rgb(1,1,1,.8),add=TRUE)
		legend('topright',legend=c('Clean','Final'),fill=c(rgb(0,0,0,1),rgb(1,1,1,0.5)))
		dev.off()
	}
	if (QCCorr$RefType=='one') {
	  rmsdRef=rmsDist(QCCorr$RefFeats)       # 5800887
	  rmsdRefClean=rmsDist(QCCorr$RefFeatsClean)     # 3670441
	  rmsdRefCorr=rmsDist(QCCorr$RefFeatsCorr)   # 2264153
	  rmsdRefFinal=rmsDist(RefFeatsFinal) # 2209906
	  RefRMSD=data.frame(rmsdRef=rmsdRef,rmsdRefClean=rmsdRefClean,rmsdRefCorr=rmsdRefCorr,rmsdRefFinal=rmsdRefFinal)
	  QCCorr$RefFeatsFinal=RefFeatsFinal
	  QCCorr$RefRMSD=RefRMSD
	}
	# Write algorithm to remove "removed" variables from clusters
	varClust=QCCorr$varClust
	N=length(varClust)
	CVAfter=nFeatAfter=numeric(N)
	for (n in 1:N) {
		vars=varClust[[n]]
		clustVars=as.matrix(QCFeatsFinal[,finalVars%in%vars])
		CVAfter[n]=mean(cv(clustVars))
		nFeatAfter[n]=ncol(clustVars)
	}
	actionInfo=QCCorr$actionInfo
	actionInfo=cbind(actionInfo[,1:2],nFeatAfter,actionInfo[,3:5],CVAfter)
	colnames(actionInfo)[2:3]=c('nBefore','nAfter')
	QCCorr$actionInfo=actionInfo
	QCCorr$QCFeatsFinal=QCFeatsFinal
	QCCorr$TestFeatsFinal=TestFeatsFinal
	QCCorr$finalVars=finalVars
	QCCorr$QCcvs=QCcvs
	QCFinal=QCCorr
	return(QCFinal)
}

#' DC: Grab samples for drift correction
#'
#' Wrapper function for grabbing QCs, reference and entire batch samples from XCMS-set
#' @param XS an XCMS object
#' @param batch a batch identifier
#' @param QC a QC sample identifier
#' @param Ref a reference sample identifier
#' @return An object containing:
#' @return QC: A QC object
#' @return Ref: A Ref object
#' @return Batch: A batch object
#' @export
## Wrapper function for grabbing QCs, reference and entire batch samples from XCMS-set
grabWrap=function(XS,batch,QC='QC',Ref='Ref') {
	QCObj=grabQC(XS,batch=batch,grp=QC)
	RefObj=grabRef(XS,QCObj,grp=Ref)
	BatchObj=grabBatch(XS,QCObj)
	return(list(QC=QCObj,Ref=RefObj,Batch=BatchObj))
}

#' DC: Grab samples for drift correction after batch alignment
#'
#' Wrapper function for grabbing QCs, reference and entire batch samples combining information from XCMS object and batch-aligned peaktable.
#' @param XS an XCMS object
#' @param PT a master peaktable after batch alignment
#' @param batch a batch identifier
#' @param QC a QC sample identifier
#' @param Ref a reference sample identifier
#' @return An object containing:
#' @return QC: A QC object
#' @return Ref: A Ref object
#' @return Batch: A batch object
#' @export
## Wrapper function for grabbing QCs, reference and entire batch samples from XCMS-set
grabWrapBA=function(XS,PT,batch,QC='QC',Ref='Ref') {
  QCObj=grabQC(XS,PT,batch=batch,grp=QC)
  RefObj=grabRef(XS,PT,QCObj,grp=Ref)
  BatchObj=grabBatch(XS,PT,QCObj)
  return(list(QC=QCObj,Ref=RefObj,Batch=BatchObj))
}

#' DC: Perform drift correction
#'
#' Wrapper function for all drift subfunctions (clust, driftCalc, driftCorr, cleanVar)
#' @param grabObj a grabWrap object
#' @param refType type of reference sample strategy used (at present "one" or "none" supported)
#' @param CVlimit QC feature CV limit as final feature inclusion criterion
#' @param report boolean whether to print a pdf report of drift models
#' @return A drift corrected object with QC features CV<limit containing final peak table
#' @export
## Wrapper function for all drift subfunctions (clust, driftCalc, driftCorr, cleanVar)
driftWrap=function(grabObj,refType,CVlimit=0.3,report=FALSE) {
	QCObj=grabObj$QC
	RefObj=grabObj$Ref
	CorrObj=grabObj$Batch
	if (min(RefObj$inj)<min(QCObj$inj) | min(CorrObj$inj)<min(QCObj$inj)) {
		cat('\nReference or test injection outside drift calculation region: Before first QCs injection.')
		cat('\nCalculation aborted\n')
		break
	}
	if (max(RefObj$inj)>max(QCObj$inj) | max(CorrObj$inj)>max(QCObj$inj)) {
		cat('\nReference or test injection outside drift calculation region: After last QCs injection.')
		cat('\nCalculation aborted\n')
		break
	}
	A=driftList=clust(QCObj$inj,QCObj$Feats,report=report)
	B=driftList=driftCalc(driftList,report=report)
	D=driftList=driftCorr(driftList,refList=RefObj,refType=refType,CorrObj=CorrObj,report=report)
	E=driftList=cleanVar(driftList,CVlimit=CVlimit,report=report)
}
