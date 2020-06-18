###L/S
LS.fun<-function (dat, batch, mod = NULL,mean.only = FALSE, ref.batch = NULL) {
  if (mean.only == TRUE) {
    message("Using the 'mean only' version of ComBat")
  }
  if (length(dim(batch)) > 1) {
    stop("This version only allows one batch variable")
  }
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1 + batch)
  if (!is.null(ref.batch)) {
    if (!(ref.batch %in% levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    cat("Using batch =", ref.batch, "as a reference batch (this batch won't change)\n")
    ref <- which(levels(as.factor(batch)) == ref.batch)
    batchmod[, ref] <- 1
  }
  else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  }
  n.batches <- sapply(batches, length)
  if (any(n.batches == 1)) {
    mean.only = TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  if (!is.null(ref)) {
    check[ref] <- FALSE
  }
  design <- as.matrix(design[, !check])
  message("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n.batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if (ncol(design) > (n.batch + 1)) {
      if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[,-c(1:n.batch)]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      }
      else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  NAs <- any(is.na(dat))
  if (NAs) {
    message(c("Found", sum(is.na(dat)), "Missing Data Values"), 
            sep = " ")
  }
  cat("Standardizing Data across genes\n")
  if (!NAs) {
    B.hat <- solve(crossprod(design), tcrossprod(t(design), 
                                                 as.matrix(dat)))
  }
  else {
    B.hat <- apply(dat, 1, Beta.NA, design)
  }
  if (!is.null(ref.batch)) {
    grand.mean <- t(B.hat[ref, ])
  }
  else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  if (!NAs) {
    if (!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat - t(design[batches[[ref]], 
                                         ] %*% B.hat))^2) %*% rep(1/n.batches[ref], n.batches[ref])
    }
    else {
      var.pooled <- ((dat - t(design %*% B.hat))^2) %*% 
        rep(1/n.array, n.array)
    }
  }
  else {
    if (!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat - t(design[batches[[ref]], 
                                               ] %*% B.hat), na.rm = TRUE)
    }
    else {
      var.pooled <- rowVars(dat - t(design %*% B.hat), 
                            na.rm = TRUE)
    }
  }
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array)))
  
  message("Fitting L/S model")
  batch.design <- design[, 1:n.batch]
  if (!NAs) {
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design), 
                                                           as.matrix(s.data)))
  }
  else {
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design)
  }
  delta.hat <- NULL
  for (i in batches) {
    if (mean.only == TRUE) {
      delta.hat <- rbind(delta.hat, rep(1, nrow(s.data)))
    }
    else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[, i],
                                            na.rm = TRUE))
    }
  }
  message("Adjusting the Data\n")
  Lsdata <- s.data
  j <- 1
  for (i in batches) {
    Lsdata[, i] <- (Lsdata[, i] - t(batch.design[i,] %*% gamma.hat))/(sqrt(delta.hat[j, ]) %*% t(rep(1,n.batches[j])))
    j <- j + 1
  }
  Lsdata <- Lsdata * (sqrt(var.pooled) %*% t(rep(1,n.array))) + stand.mean
  if (!is.null(ref.batch)) {
    Lsdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  return(Lsdata)
}


###edgeR
EdgerRemoveBatchEffect<-function(data,group,batch,numBatch){
  design<-model.matrix(~factor(group)+factor(batch))
  y<-as.matrix(data)
  y <- DGEList(counts=y)
  y <- calcNormFactors(y)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit<-glmFit(y,design)
  ncol<-ncol(design)
  cutoff<-ncol-numBatch+2
  coef<-fit$coefficients[,cutoff:ncol]
  des<-design[,cutoff:ncol]
  edgeR_data<-2^(t(log2(t(data))-des%*%t(coef)))
  return(edgeR_data)
}

###limmaVoom
limmaVoomRemoveBatchEffect<-function(data,group,batch,numBatch){
  design<-model.matrix(~factor(sample)+factor(site))
  v<-voom(edata, design)
  fit <- lmFit(v,design)
  ncol<-ncol(design)
  cutoff<-ncol-numBatch+2
  coef<-fit$coefficients[,cutoff:ncol]
  des<-design[,cutoff:ncol]
  voom_data<-2^(t(log2(t(data))-des%*%t(coef)))
  return(voom_data)
}

###DESeq
DESeqRemoveBatchEffect<-function(data,group,batch,numBatch){
  modelFrame<- model.matrix(~factor(group)+factor(batch))
  cds <- newCountDataSet(data,modelFrame)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds, method="pooled")
  fit1 <- fitNbinomGLMs(cds, count ~ factor(group)+factor(batch))
  ncol<-ncol(modelFrame)
  cutoff<-ncol-numBatch+2
  coef<-fit1[,cutoff:ncol]
  des<-modelFrame[,cutoff:ncol]
  NBGLM_data<-2^(t(log2(t(data))-des%*%t(coef)))
  return(NBGLM_data)
}
