smoothCoverageTrackNoCN0 <- function(lCT,
                                lSe,
                                lGCT,
                                lNormals=NULL,
                                method=c("loess",
                                         "lowess"))
{
    
    allRec.counts = unlist(lapply(lCT,function(x) x$records))
    allRec.counts.clusters = Mclust(allRec.counts, G = 1:10 )
    cn0.cluster = which(allRec.counts.clusters$parameters$mean < 10)

    if (length(cn0.cluster) > 0) {
            idx_excluded = allRec.counts.clusters$classification %in% cn0.cluster
            idx_clean = !idx_excluded
    } else {
            idx_excluded = rep(FALSE, length(allRec.counts))
            idx_clean = rep(TRUE, length(allRec.counts))
    }

    allRec <- unlist(lapply(lCT,function(x) log2(x$records+1)))

    if(!is.null(lNormals))
    {
        if(!"records"%in%names(lNormals[[1]]))
            allRec <- smoothNormals(allRec, lNormals)
        else
            allRec <- allRec-unlist(lapply(lNormals,function(x) log2(x$records+1)))
    }

    allGC <- unlist(lGCT)
    starts <- c(0,cumsum(sapply(lCT,nrow)[-c(length(lCT))]))+1
    ends <- cumsum(sapply(lCT,nrow))

    allGC.clean = allGC[idx_clean]
    allRec.clean = allRec[idx_clean]


    if(method[1]=="lowess")
        smoothT <- lowess(allGC.clean,allRec.clean,f=2/3,iter=10)
    if(method[1]=="loess")
        smoothT <- myloess(LL=length(allRec.clean),allRec.clean~allGC.clean)

    # 1. Get fitted/residuals for clean data
    fitted_clean <- smoothT$fitted
    resid_clean  <- smoothT$residuals
    
    # 2. Predict fitted values for excluded bins based on their GC content
    # We use the existing model (smoothT) to estimate the trend at these points
    
    fitted_excluded <- predict(smoothT, newdata = rep(median(allGC), times = sum(idx_excluded)))
    # fitted_excluded <- predict(smoothT, newdata = allGC[idx_excluded])

    # 3. Calculate residuals for excluded bins (Observed - Predicted)
    resid_excluded <- allRec[idx_excluded] - fitted_excluded

    # Initialize empty vectors for the full dataset
    final_fitted <- numeric(length(allRec))
    final_resid  <- numeric(length(allRec))

    # Fill clean data
    final_fitted[idx_clean] <- fitted_clean
    final_resid[idx_clean]  <- resid_clean
    
    # Fill reintegrated data
    final_fitted[idx_excluded] <- fitted_excluded
    final_resid[idx_excluded]  <- resid_excluded

    for(i in 1:length(lCT))
    {
        lCT[[i]] <- cbind(lCT[[i]], final_fitted[starts[i]:ends[i]],
                          final_resid[starts[i]:ends[i]])
        colnames(lCT[[i]])[(ncol(lCT[[i]])-1):ncol(lCT[[i]])] <- c("fitted","smoothed")
    }
    return(lCT)
}