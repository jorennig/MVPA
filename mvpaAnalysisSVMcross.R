rm(list=ls(all=TRUE))
setwd("G:/MVPA")

filePathGestalt <- "G:/MVPA/betaFilesGestalt"
filesGestalt <- list.files(filePathGestalt)
pids <- unique(sapply(strsplit(filesGestalt, "_"), function(x) x[[1]][1]))
filesGestalt <- list.files(filePathGestalt, full.names = TRUE)
filePathSubitizing <- "G:/MVPA/betaFilesSubitizing"
filesSubitizing <- list.files(filePathSubitizing, full.names = TRUE)


i <- 1


pid <- pids[2]

acc <- NULL
for (pid in pids){
  print(pid)
  pidFilesGestalt <- filesGestalt[grep(pid, filesGestalt)]
  pidFilesSubitizing <- filesSubitizing[grep(pid, filesSubitizing)]
  
  df1 <- read.table(pidFilesGestalt[1], sep="\t", header=TRUE)
  head(df1)
  
  colnames(df1) <- c(paste0("condition_",df1$Var2[1]),
                     paste0("hemisphere_",df1$Var2[1]),
                     paste0("trial_",df1$Var2[1]),
                     paste0("voxel_", df1$Var2[1], "_", 1:(ncol(df1)-3)))
  
  df2 <- NULL
  if(length(pidFilesGestalt) == 2){
    df2 <- read.table(pidFilesGestalt[2], sep="\t", header=TRUE)
    colnames(df2) <- c(paste0("condition_",df2$Var2[1]),
                       paste0("hemisphere_",df2$Var2[1]),
                       paste0("trial_",df2$Var2[1]),
                       paste0("voxel_", df2$Var2[1], "_", 1:(ncol(df2)-3)))
    df2 <- df2[order(df2[,3]),]
  }
  
  df1 <- df1[order(df1[,3]),]
  if(!is.null(df2)){
    dfGestalt <- cbind(df1, df2)
  } else {
    dfGestalt <- df1
  }
  
  df1 <- read.table(pidFilesSubitizing[1], sep="\t", header=TRUE)
  head(df1)
  
  colnames(df1) <- c(paste0("condition_",df1$Var2[1]),
                     paste0("hemisphere_",df1$Var2[1]),
                     paste0("trial_",df1$Var2[1]),
                     paste0("voxel_", df1$Var2[1], "_", 1:(ncol(df1)-3)))
  
  df2 <- NULL
  if(length(pidFilesSubitizing) == 2){
    df2 <- read.table(pidFilesSubitizing[2], sep="\t", header=TRUE)
    colnames(df2) <- c(paste0("condition_",df2$Var2[1]),
                       paste0("hemisphere_",df2$Var2[1]),
                       paste0("trial_",df2$Var2[1]),
                       paste0("voxel_", df2$Var2[1], "_", 1:(ncol(df2)-3)))
    df2 <- df2[order(df2[,3]),]
  }
  
  df1 <- df1[order(df1[,3]),]
  if(!is.null(df2)){
    dfSubitizing<- cbind(df1, df2)
  } else {
    dfSubitizing <- df1
  }
  
  colnames(dfGestalt)
  colnames(dfSubitizing)
  
  yGestalt <- ifelse(dfGestalt[,1] == "gestalt", 1, 0) # 1 = gestalt, 0 = no_gestalt
  XGestalt <- dfGestalt[,grep("voxel", colnames(dfGestalt))]

  ySubitizing <- 1 - grepl("non", dfSubitizing[,1]) * 1 # 1 = can, 0 = non_can
  XSubitizing <- dfSubitizing[,grep("voxel", colnames(dfSubitizing))]
  
  XGestalt <- XGestalt[colSums(!is.na(XGestalt)) > 0]
  XSubitizing <- XSubitizing[colSums(!is.na(XSubitizing)) > 0]
  
  colGestalt <- colnames(XGestalt)[colnames(XGestalt) %in% colnames(XSubitizing)]
  cols <- colnames(XSubitizing)[colnames(XSubitizing) %in% colGestalt]
  
  XGestalt <- XGestalt[,cols]
  XSubitizing <- XSubitizing[,cols]
  
  print(all(colnames(XGestalt) == colnames(XSubitizing)))
  
  XGestalt <-  as.matrix(XGestalt)
  XSubitizing <-  as.matrix(XSubitizing)
  
  require(xgboost)
  as.matrix(XGestalt)
  
  svm.model <- svm(yGestalt ~ ., data = XGestalt, probability=TRUE)
  pred <- (predict(svm.model, newdata=XSubitizing) > 0.5) * 1
  
  acc <- c(acc, mean(ySubitizing == pred))
  print(mean(ySubitizing == pred))
}

acc
t.test(acc, mu=0.5)
