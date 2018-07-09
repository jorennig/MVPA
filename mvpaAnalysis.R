rm(list=ls(all=TRUE))
setwd("G:/MVPA")

filePath <- "G:/MVPA/betaFiles"
files <- list.files(filePath)
pids <- unique(sapply(strsplit(files, "_"), function(x) x[[1]][1]))
files <- list.files(filePath, full.names = TRUE)

i <- 1


meanAccs <- NULL
for (pid in pids){
  pidFiles <- files[grep(pid, files)]
  
  df1 <- read.table(pidFiles[1], sep="\t", header=TRUE)
  head(df1)
  
  colnames(df1) <- c(paste0("condition_",df1$Var2[1]),
                     paste0("hemisphere_",df1$Var2[1]),
                     paste0("trial_",df1$Var2[1]),
                     paste0("voxel_", df1$Var2[1], "_", 1:(ncol(df1)-3)))
  
  df2 <- NULL
  if(length(pidFiles) == 2){
    df2 <- read.table(pidFiles[2], sep="\t", header=TRUE)
    colnames(df2) <- c(paste0("condition_",df2$Var2[1]),
                       paste0("hemisphere_",df2$Var2[1]),
                       paste0("trial_",df2$Var2[1]),
                       paste0("voxel_", df2$Var2[1], "_", 1:(ncol(df2)-3)))
    df2 <- df2[order(df2[,3]),]
  }
  
  df1 <- df1[order(df1[,3]),]
  if(!is.null(df2)){
    df <- cbind(df1, df2)
  } else {
    df <- df1
  }
  
  colnames(df)
  df$trial_right
  
  y <- ifelse(df[,1] == "gestalt", 1, 0)
  y
  
  X <- df[,grep("voxel", colnames(df))]
  X
  
  # exclude NaN voxel
  X <- X[colSums(!is.na(X)) > 0]
  
  
  require(xgboost)
  
  require(caret)
  set.seed(2345)
  folds <- createFolds(as.factor(y), k = 10)
  
  acc <- NULL
  for(fold in folds){
    valid <- fold
    model <- c(1:length(y))[-valid]
    
    dmodel <- xgb.DMatrix(as.matrix(X)[model,], label = y[model], missing = NA)
    dvalid <- xgb.DMatrix(as.matrix(X)[valid,], label = y[valid], missing = NA)
    
    param <- list(objective = "binary:logistic", 
                  eta = 0.01, 
                  max_depth = 4,
                  #min_child_weight = 3,
                  subsample = 1,
                  colsample_bytree = 1
                  )
    
    xgb <- xgb.train(data = dmodel
                     , param
                     , nrounds = 100
                     , watchlist = list(val = dvalid, mod = dmodel)
                     #, eval_metric = "auc"
                     , maximize = FALSE
                     , early.stop.round = 5)
    
    
    importance_matrix <- xgb.importance(feature_names=colnames(X), model = xgb)
    xgb.plot.importance(importance_matrix)
    
    pred <- (predict(xgb, dvalid) > 0.5) * 1
    acc <- c(acc, mean(y[valid] == pred))
    
  }
  
  
  meanAccs <- c(meanAccs, mean(acc))
}

meanAccs
t.test(meanAccs, mu=0.5)
