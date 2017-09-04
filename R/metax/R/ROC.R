
##' @title Classical univariate ROC analysis
##' @description Classical univariate ROC analysis
##' @rdname calcAUROC
##' @docType methods
##' @param x A numeric vector
##' @param y A response vector
##' @param cgroup Sample class used
##' @param plot A logical indicates whether plot
##' @param ... Additional parameter
##' @return A data.frame
##' @examples
##' \donttest{
##' x <- rnorm(50,2,2)
##' y<- rep(c("c","t"),25)
##' calcAUROC(x,y)
##' }
setGeneric("calcAUROC",function(x,y,cgroup,plot,...) 
    standardGeneric("calcAUROC"))

## classical univariate ROC analysis
##' @describeIn calcAUROC
setMethod("calcAUROC", signature(x="numeric"), function(x,y,plot=FALSE,...){
    roc.obj <- pROC::roc(y,x,percent = FALSE)
    auc.ci <- ci.auc(roc.obj, method = "bootstrap",boot.n = 500, 
                     progress = "none")
    rocdata <- data.frame(roc=roc.obj$auc,
                          lowROC=auc.ci[1],
                          upROC=auc.ci[3])
    if(plot==TRUE){
        plot.roc(roc.obj,print.auc=FALSE, legacy.axes=TRUE, 
                 col="navy", grid=TRUE,
                 xlab = "False positive rate", 
                 ylab="True positive rate")
        
    }
    return(rocdata)
    
}
)


##' @title Classical univariate ROC analysis
##' @description Classical univariate ROC analysis
##' @param para A metaXpara object
##' @param cgroup Samples used
##' @param cpu The number of CPU used
##' @param plot A logical indicates whether plot
##' @param ... Additional parameter
##' @return A metaXpara object
##' @examples
##' \dontrun{
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' para <- missingValueImpute(para)
##' addValueNorm(para) <- para
##' res <- myCalcAUROC(para,cgroup=c("S","C"))
##' }
myCalcAUROC=function(para,cgroup,cpu=0,plot=FALSE,...){
    peaksData <- dplyr::filter(para@peaksData,class %in% cgroup)
    peaksData$class <- as.character(peaksData$class)
    myROC=function(id,pData){
        dat <- dplyr::filter(pData,ID==id)
        x <- dat$valueNorm
        y <- dat$class
        roc.obj <- pROC::roc(y,x,percent = FALSE)
        auc.ci <- pROC::ci.auc(roc.obj, method = "bootstrap",boot.n = 500, 
                         progress = "none")
        rocdata <- data.frame(ID=id,roc=roc.obj$auc,lowROC=auc.ci[1],
                              upROC=auc.ci[3])
        return(rocdata)
    }
    cpu = ifelse(cpu==0,detectCores(),cpu)
    cl <- makeCluster(getOption("cl.cores", cpu))
    clusterExport(cl, c("myROC","roc","ci.auc"),envir=environment())
    clusterEvalQ(cl,library("pROC"))
    rocRes <- parLapply(cl,unique(peaksData$ID),myROC,pData=peaksData)
    stopCluster(cl)
    rocDF <- as.data.frame(rbindlist(rocRes))
}

##' @title Perform biomarker analysis and plot ROC curve
##' @description Perform biomarker analysis and plot ROC curve
##' @param para An object of \code{metaXpara}
##' @param group The sample class used
##' @param method Method for modeling
##' @param valueID The column name used
##' @param ratio (The number of sample for training)/(The number of sample for 
##' Testing)
##' @param k The number of sampling. Default is 100.
##' @param plotCICurve Whether plot confidence interval curve in the ROC figure
##' @param ... Additional parameters
##' @return The result of feature selection and modeling
##' @export
##' 
plotROC=function(para,group,fID,valueID="valueNorm",method="rf",ratio=2/3,
                 k=100,plotCICurve=TRUE,...){
    peaksData <- para@peaksData %>% filter(class %in% group, ID %in% fID) %>%
        mutate(class=as.character(class))
    
    x <- peaksData %>% 
        select_("ID","class","sample",valueID) %>% 
        spread_("ID",valueID)
    
    cl <- x$class
    
    mData <- select(x,-class,-sample)
    
    trueClassLabel <- vector(length = k, mode = "list")
    fImp <- vector(length = k, mode = "list")
    aucVec <- numeric(length = k)
    predProb <- vector(length = k, mode = "list")
    
    indList <- splitData(cl,ratio = ratio, k=k)
    accuVec <- numeric(length = k)
    
    for(i in 1:k){
        sampleTrainInd <- indList$train[i,]
        sampleTestInd <- indList$test[i,]
        
        x_train <- mData[sampleTrainInd,,drop=FALSE]
        x_test <- mData[sampleTestInd,,drop=FALSE]
        
        y_train <- cl[sampleTrainInd]
        y_test <- cl[sampleTestInd]
        
        trueClassLabel[[i]] <- y_test
        
        #save(x_train,y_train,x_test,y_test,x,mData,file="test_bug.rda")
        res <- modelEval(x_train = x_train, y_train = y_train, 
                         x_test = x_test, method = method)
        fImp[[i]] <- res$imp
        
        prob <- res$pred
        
        pred <- prediction(prob,y_test)
        aucVec[i] <- slot(performance(pred, "auc"), "y.values")[[1]]
        predProb[[i]] <- prob
        
        predNum <- as.factor(ifelse(prob>0.5, 1,0))
        atable <- table(predNum,y_test)
        accuVec[i] <- sum(diag(atable)) / sum(atable)
        
        
    }
    
    
    predProbVec <- unlist(predProb)
    trueClassLabelVec <- unlist(trueClassLabel)
    predVec <- prediction(predProbVec,trueClassLabelVec)
    meanAUC <- mean(aucVec)
    aucCI <- calculateCI(as.matrix(aucVec))
    
    ## plot ROC figure
    
    predProbList <- prediction(predProb,trueClassLabel)
    auroc <- round(meanAUC,3)
    
    perf <- performance(predProbList,"tpr","fpr")
    meanPerf <- computeAverageCurve(perf,avg.method = FALSE)
    y_all <- meanPerf@y.values[[1]]
    x_all <- meanPerf@x.values[[1]]
    
    addText <- paste("AUC = ", auroc, "\n","95% CI: ", aucCI,sep="")
    
    fig <- paste(para@outdir,"/",para@prefix,"-",paste(group,collapse = "_vs_"),"_ROC.pdf",sep="")
    message(date(),"\t","ROC curve: ",fig)
    pdf(fig,width = 4,height = 4)
    par(mar=c(3,3,1,1),mgp=c(1.6,0.6,0))
    plot(x_all,y_all, type="n", axes=F,
         xlim=c(0,1), ylim=c(0,1),cex.lab=1.1,
         xlab="1-Specificity",cex.lab=1.2,font.lab=2,cex.axis=1.1,
         ylab="Sensitivity"
    )

    box()
    axis(side=2,cex.lab=1.2,font.lab=2,cex.axis=1.1)
    lab.at <- seq(0, 1, 0.2);
    grid.at <- seq(0, 1, 0.1);
    lab.labels <- lab.at
    axis(side=1, at=lab.at,labels=as.graphicsAnnot(sprintf( "%.1f",lab.labels)),cex.lab=1.2,font.lab=2,cex.axis=1.1)
    x_all_0 <- c(0,x_all)
    y_all_0 <- c(0,y_all)
    lines(x_all_0, y_all_0, type="l", lwd=1.1, col="red")
    
    if(plotCICurve){
        res <- computeHL(perf)
        suppressWarnings(polygon(c(x_all, rev(x_all)), c(res$con.low, 
                                 rev(res$con.high)), col="#0000ff22"))
    }
    
    legend("bottomright", legend = addText,  bty="n")
    dev.off()
    
}


modelEval=function(x_train,y_train,x_test,method="rf",fImp=TRUE){
    
    # y_train <- as.factor(as.numeric(y_train)-1)
    y_train <- as.factor(as.numeric(as.factor(y_train))-1)
    res <- list()
    if(method == "rf"){
        
        model <- randomForest(x_train,y_train,ntree=300,importance = TRUE)
        res$pred <- randomForest:::predict.randomForest(model,x_test,type="prob")[,"1"]
        
        if(fImp){
            res$imp <- importance(model)[, "MeanDecreaseAccuracy"]
        }
        return(res)
    }
    
}


splitData=function(cl,ratio=0.667,k=500){
     nSample <- length(cl)
     smallClassLabel <- names(sort(table(cl)))[1]
     nSmallClass <- sum(cl == smallClassLabel)
     
     nSmallTrain <- round(ratio*nSmallClass)
     nBigTrain <- nSmallTrain
     
     nSmallTest <- nSmallClass - nSmallTrain
     nBigTest <- nSample - nSmallTrain -nBigTrain - nSmallTest
     
     
     cat("Train sample:",smallClassLabel,"=",nSmallTrain,",",
         names(sort(table(cl)))[2],"=",nBigTrain,"\n")
     
     cat("Test sample:",smallClassLabel,"=",nSmallTest,",",
         names(sort(table(cl)))[2],"=",nBigTest,"\n")
     
     ## index of each sample
     smallClassInd <- which(cl==smallClassLabel)
     bigClassInd <- c(1:nSample)[-smallClassInd]
     
     nTrain <- nSmallTrain + nBigTrain
     nTest <- nSmallTest + nBigTest
     
     indList <- list(train = matrix(0, nrow = k, ncol = nTrain),
                     test = matrix(0, nrow = k, ncol = nTest))

     for(i in 1:k){
         sampleSmallTrainInd <- sample(smallClassInd,nSmallTrain)
         sampleBigTrainInd <- sample(bigClassInd,nBigTrain)
         
         sampleTrainInd <- c(sampleSmallTrainInd,sampleBigTrainInd)
         indList$train[i,] <- sampleTrainInd
         indList$test[i,] <- c(1:nSample)[-sampleTrainInd]
         
     }
     return(indList)
}


calculateCI <- function(data, param=F){
    means <- colMeans(data, na.rm=T)
    if(param){
        sds <- apply(data, 2, sd, na.rm=T)
        nn <- nrow(data)
        std.err <- sds/sqrt(nn)
        LCL <- round(means-1.96*std.err,3)
        LCL[LCL<0] <- 0
        UCL <- round(means+1.96*std.err, 3)
        UCL[UCL>1] <- 1
        res <- paste(LCL, "-", UCL, sep="")
    }else{
        cis <- apply(data, 2, quantile, probs=c(0.025, 0.975))
        cis <- round(cis,3)
        cis[cis<0] <- 0
        cis[cis>1] <- 1
        res <- paste(cis[1,], "-", cis[2,], sep="")
    }
    return(res)
}


# now compute average curve
computeAverageCurve<-function(perf, avg.method){
    # now get the average curve
    perf.avg = perf
    if(avg.method == "vertical"){
        x.values <- seq(min(unlist(perf@x.values)), max(unlist(perf@x.values)),
                        length=max( sapply(perf@x.values, length)))
        for (i in 1:length(perf@y.values)) {
            perf.avg@y.values[[i]] <-
                approxfun(perf@x.values[[i]], perf@y.values[[i]],
                          ties=mean, rule=2)(x.values)
        }
        perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values )))
        perf.avg@x.values <- list(x.values)
    }else if(avg.method == "horizontal"){
        y.values <- seq(min(unlist(perf@y.values)), max(unlist(perf@y.values)),
                        length=max(sapply(perf@y.values, length)))
        for (i in 1:length(perf@x.values)) {
            perf.avg@x.values[[i]] <- approxfun(perf@y.values[[i]],
                                                perf@x.values[[i]],
                                                ties=mean, rule=2)(y.values)
        }
        perf.avg@x.values <- list(rowMeans( data.frame( perf.avg@x.values )));
        perf.avg@y.values <- list(y.values);
    }else{ # threshold
        all.alphas <- unlist(perf@alpha.values);
        min.alpha <- min(all.alphas);
        if(min.alpha == -Inf){
            min.alpha <- 0;
        }
        max.alpha <- max(all.alphas);
        if(max.alpha == Inf){
            max.alpha <- 1.0;
        }
        
        alpha.values <- rev(seq(min.alpha, max.alpha,length=max(sapply(perf@alpha.values, length))));
        perf.sampled <- perf;
        for (i in 1:length(perf.sampled@y.values)) {
            perf.sampled@x.values[[i]] <-
                approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                          rule=2, ties=mean)(alpha.values)
            perf.sampled@y.values[[i]] <-
                approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                          rule=2, ties=mean)(alpha.values)
        }
        ## compute average curve
        perf.avg <- perf.sampled
        perf.avg@x.values <- list(rowMeans(data.frame(perf.avg@x.values)))
        perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values)))
    }
    return(perf.avg);
}

computeHL <- function(perf){
    all.alphas <- unlist(perf@alpha.values)
    min.alpha <- min(all.alphas)
    if(min.alpha == -Inf){
        min.alpha <- 0
    }
    max.alpha <- max(all.alphas)
    if(max.alpha == Inf){
        max.alpha <- 1.0
    }
    
    alpha.values <- rev(seq(min.alpha, max.alpha,length=max(sapply(perf@alpha.values, length))))
    perf.sampled <- perf
    for (i in 1:length(perf.sampled@y.values)) {
        perf.sampled@x.values[[i]] <-
            approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                      rule=2, ties=mean)(alpha.values)
        perf.sampled@y.values[[i]] <-
            approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                      rule=2, ties=mean)(alpha.values)
    }
    y.data <- data.frame(perf.sampled@y.values)
    con.low <- apply(y.data, 1, quantile, 0.05)
    con.high <- apply(y.data, 1, quantile, 0.95)
    res <- list(con.low = con.low, con.high = con.high)
    return (res)
}






