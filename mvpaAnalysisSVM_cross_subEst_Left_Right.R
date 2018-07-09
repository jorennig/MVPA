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

accLeft <- NULL
accRight <- NULL
accGestalt <- NULL
predDf <- NULL
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
  
  yGestalt <- ifelse(dfGestalt[,1] == "gestalt", 0, 1) # 0 = gestalt, 1 = no_gestalt
  XGestaltRight <- dfGestalt[,grep("voxel_right", colnames(dfGestalt))]
  XGestaltLeft <- dfGestalt[,grep("voxel_left", colnames(dfGestalt))]

  num <- as.numeric(sapply(dfSubitizing[,1], function(c) strsplit(as.character(c), "can")[[1]][2]))
  
  ySubitizing <- ifelse(num <= 4, 0, 1) # 0 = sub, 1 = est
  XSubitizingRight <- dfSubitizing[,grep("voxel_right", colnames(dfSubitizing))]
  XSubitizingLeft <- dfSubitizing[,grep("voxel_left", colnames(dfSubitizing))]
  
  #print(all(colnames(XGestalt) == colnames(XSubitizing)))
  
  XGestaltRight <- XGestaltRight[colSums(!is.na(XGestaltRight)) > 0]
  XSubitizingRight <- XSubitizingRight[colSums(!is.na(XSubitizingRight)) > 0]
  XGestaltLeft <- XGestaltLeft[colSums(!is.na(XGestaltLeft)) > 0]
  XSubitizingLeft <- XSubitizingLeft[colSums(!is.na(XSubitizingLeft)) > 0]
  
  require(e1071)
  colGestaltRight <- colnames(XGestaltRight)[colnames(XGestaltRight) %in% colnames(XSubitizingRight)]
  colsRight <- colnames(XSubitizingRight)[colnames(XSubitizingRight) %in% colGestaltRight]
  
  colGestaltLeft <- colnames(XGestaltLeft)[colnames(XGestaltLeft) %in% colnames(XSubitizingLeft)]
  colsLeft <- colnames(XSubitizingLeft)[colnames(XSubitizingLeft) %in% colGestaltLeft]
  
  XGestaltRight <- XGestaltRight[,colsRight]
  XSubitizingRight <- XSubitizingRight[,colsRight]
  
  XGestaltLeft <- XGestaltLeft[,colsLeft]
  XSubitizingLeft <- XSubitizingLeft[,colsLeft]
  
  print(all(colnames(XGestaltRight) == colnames(XSubitizingRight)))
  print(all(colnames(XGestaltLeft) == colnames(XSubitizingLeft)))
  
  XGestaltRight <-  as.matrix(XGestaltRight)
  XSubitizingRight <-  as.matrix(XSubitizingRight)
  
  XGestaltLeft <-  as.matrix(XGestaltLeft)
  XSubitizingLeft <-  as.matrix(XSubitizingLeft)
  
  yGestalt <- factor(yGestalt)
  ySubitizing <- factor(ySubitizing)
  
  if(ncol(XGestaltRight) > 0){
    set.seed(345)
    svm_tune <- tune.svm(x=XGestaltRight, y=yGestalt, 
                     kernel="radial", cost=10^(-3:3))
    
    svm.model <- svm(yGestalt ~ ., data = XGestaltRight, probability=TRUE, kernel="radial", cross=10, cost=svm_tune$best.parameters[1])
    accGestalt <- c(accGestalt, mean(svm.model$accuracies))
    pred <- predict(svm.model, newdata=XSubitizingRight)
    
    #pred <- (predict(xgb, dtest) > 0.5) * 1
    accRight <- c(accRight, mean(ySubitizing == pred))
    print(mean(ySubitizing == pred))
    predDf <- rbind(predDf,data.frame(hemisphere=rep("right",length(pred)), 
                                      pid=rep(pid,length(pred)), 
                                      pred=(ySubitizing == pred)*1, 
                                      cor=ySubitizing, 
                                      svm=pred,
                                      condition=dfSubitizing$condition_right))
  }
  
  if(ncol(XGestaltLeft) > 0){
    set.seed(345)
    svm_tune <- tune.svm(x=XGestaltLeft, y=yGestalt, 
                         kernel="radial", cost=10^(-3:3))
    svm.model <- svm(yGestalt ~ ., data = XGestaltLeft, probability=TRUE, kernel="radial", cross=10, cost=svm_tune$best.parameters[1])
    #accGestalt <- c(accGestalt, mean(svm.model$accuracies))
    pred <- predict(svm.model, newdata=XSubitizingLeft)
    
    #pred <- (predict(xgb, dtest) > 0.5) * 1
    accLeft <- c(accLeft, mean(ySubitizing == pred))
    print(mean(ySubitizing == pred))
    predDf <- rbind(predDf,data.frame(hemisphere=rep("left",length(pred)), 
                                      pid=rep(pid,length(pred)), 
                                      pred=(ySubitizing == pred)*1, 
                                      cor=ySubitizing, 
                                      svm=pred,
                                      condition=dfSubitizing$condition_left))
  }
}


t.test(accRight, mu=0.5)
t.test(accLeft, mu=0.5)

summary(predDf)
predDf$hemisphere <- factor(predDf$hemisphere)
contrasts(predDf$hemisphere) <- contr.sum

head(predDf)
str(predDf)
predDf$condition <- as.character(predDf$condition)
predDf$numerosity <- sapply(predDf$condition, function(c) substr(c, nchar(c), nchar(c)))
predDf$arrangement <- sapply(predDf$condition, function(c) substr(c, 1, nchar(c)-1))
head(predDf)

predDf$numerosity <- factor(predDf$numerosity)
predDf$arrangement <- factor(predDf$arrangement)
contrasts(predDf$numerosity) <- contr.sum
contrasts(predDf$arrangement) <- contr.sum

require(lme4)
require(afex)
glmer1 <- mixed(cor ~ hemisphere*arrangement + (1|pid), predDf, family=binomial, method = "LRT")
glmer1$anova_table
summary(glmer1)


ls1 <- lsmeans(glmer1, c("hemisphere"))
test(ls1)

t.test(accGestalt, mu=50)
accGestalt
