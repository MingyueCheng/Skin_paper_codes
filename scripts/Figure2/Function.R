#########################
library(ranger)
library(crossRanger)
#随机森林n折交叉检验
rf.regression.cross_validation <- function(x, y, nfolds=5, ntree=500, mtry=NULL){
  folds <- balanced.folds(y, nfolds = nfolds)
  if (any(table(folds) < 5)) 
    stop("Less than 5 samples in at least one fold!\n")
  result <- list()
  result$rf.model <- list()
  if (is.factor(y)) {
    stop("the input y is factor!\n")
  }else {
    result$y <- y
    result$predicted <- NA
  }
  for (fold in sort(unique(folds))) {
    cat("The # of folds: ", fold, "\n")
    foldix <- which(folds == fold)
    y_train <- result$y[-foldix]
    data <- data.frame(y = y_train, x[-foldix, ])
    result$rf.model[[fold]] <- model <- ranger(dependent.variable.name = "y", 
                                                 data = data, keep.inbag = TRUE, importance = "permutation", 
                                                 classification = ifelse(is.factor(y_train), TRUE, FALSE), 
                                                 num.trees = ntree, mtry = mtry)
    newx <- data.frame(x[foldix, ])
    if (length(foldix) == 1) 
      newx <- data.frame(matrix(newx, nrow = 1))
    predicted_foldix <- predict(model, newx)$predictions
    
    result$predicted[foldix] <- predicted_foldix
    reg_perf <- get.reg.predict.performance(model, newx, 
                                            newy = result$y[foldix])
    result$MSE[fold] <- reg_perf$MSE
    result$RMSE[fold] <- reg_perf$RMSE
    result$nRMSE[fold] <- reg_perf$nRMSE
    result$MAE[fold] <- reg_perf$MAE
    result$MAPE[fold] <- reg_perf$MAPE
    result$MASE[fold] <- reg_perf$MASE
    result$Spearman_rho[fold] <- reg_perf$Spearman_rho
    result$R_squared[fold] <- reg_perf$R_squared
    result$Adj_R_squared[fold] <- reg_perf$Adj_R_squared
    cat("Mean squared residuals: ", reg_perf$MSE, 
        "\n")
    cat("Mean absolute error: ", reg_perf$MAE, 
        "\n")
    cat("pseudo R-squared (%explained variance): ", 
        reg_perf$R_squared, "\n")
  }
  result$importances <- matrix(0, nrow = ncol(x), ncol = nfolds+2)
  colnames(result$importances) <- c(sort(unique(folds)),"mean","rank")
  rownames(result$importances) <- names(result$rf.model[[1]]$variable.importance)
  for (fold in sort(unique(folds))) {
    result$importances[, fold] <- result$rf.model[[fold]]$variable.importance
  }
  result$importances[, "mean"] <- apply(result$importances[,1:5], 1, mean)
  result$importances[, "rank"] <- rank(-result$importances[,"mean"])
  result$params <- list(ntree = ntree, nfolds = nfolds)
  result$error.type <- "cv"
  result$nfolds <- nfolds
  cat(mean(result$MAE), "+-", sd(result$MAE), "\n")
  cat(mean(result$R_squared),"+-", sd(result$R_squared),  "\n")
  return(result)
}

#特征选择
repeated_cv_feature_select <- function(x, y, feature_importance, nfolds, repeats){
  n_features<-c(2^seq(1, log2(nrow(feature_importance))),nrow(feature_importance))
  top_n_perf<-matrix(NA, ncol=repeats+3, nrow=length(n_features))
  colnames(top_n_perf)<-c("n_features", paste(c(1:repeats),"_MAE",sep = ""), "mean", "sd")
  for (i in 1:length(n_features)) {
    select_feature <- rownames(feature_importance[which(feature_importance$rank<=n_features[i]),])
    top_n_x <- x[,select_feature]
    for (r in 1:repeats) {
      top_n_rf<-rf.regression.cross_validation(top_n_x, y, nfolds=nfolds)
      top_n_perf[i, 1] <- n_features[i]
      top_n_perf[i, r+1] <- mean(top_n_rf$MAE)
    }
  }
  top_n_perf <- as.data.frame(top_n_perf)
  top_n_perf[,"mean"] <- apply(top_n_perf[,2:(ncol(top_n_perf)-2)], 1, mean)
  top_n_perf[, "sd"] <- apply(top_n_perf[,2:(ncol(top_n_perf)-2)], 1, sd)
  return(top_n_perf)
}

#保存特征重要性表格
save_feature_importance <- function(feature_importance, outpath_dir){
  feature_importance$species <- rownames(feature_importance)
  feature_importance <- dplyr::select(feature_importance, "species", everything())
  write.table(feature_importance, outpath_dir, sep="\t", row.names = F, quote = F, col.names = T)
}

#用训练好的模型对新的输入进行年龄预测
predict_newx <- function(model, newx, true_y){
  predicted_result <- list()
  predicted_result$predicted_y <- data.frame(matrix(NA, nrow = length(true_y), ncol = length(model)),
                                             row.names = rownames(newx))
  for(i in 1:length(model)){
    predicted_result$predicted_y[,i] <- predict(model[[i]], newx)$predictions
    reg_perf <- get.reg.predict.performance(model[[i]], newx, 
                                            newy = true_y)
    predicted_result$MSE[i] <- reg_perf$MSE
    predicted_result$MAE[i] <- reg_perf$MAE
    predicted_result$R_squared[i] <- reg_perf$R_squared
  }
  cat("Predict male age MAE: ", mean(predicted_result$MAE), "+-", sd(predicted_result$MAE), "\n")
  cat("Predict male age R_squared: ", mean(predicted_result$R_squared),"+-", sd(predicted_result$R_squared),  "\n")
  predicted_result$predicted_y[,"mean"] <- apply(predicted_result$predicted_y[,1:5], 1, mean)
  return(predicted_result)
}

#皮肤年龄与实际年龄作图
plot_y_vs_predict_y <- function(y,predict_y,mae,r_squared){
  df <- data.frame(y=y,predicted_y=predict_y)
  mae_label = paste("MAE: ", as.character(round(mean(mae),2))," ± ",as.character(round(sd(mae),2))," yr",sep = "")
  r2_label = paste("R squared: ", as.character(round(mean(r_squared)*100,2))," ± ",as.character(round(sd(r_squared)*100,2))," %",sep = "")
  ggplot(df, aes(x = y, y = predicted_y)) + 
    ylab(paste("Predicted ", "age", sep = "")) + 
    xlab(paste("Observed ", "age", sep = "")) + 
    geom_point(alpha = 0.1) + 
    geom_smooth(method = "loess", span = 1) + 
    annotate(geom = "text", x = Inf, y = Inf, label = mae_label, color = "grey40", vjust = 2.4, hjust = 1.5) + 
    annotate(geom = "text", x = Inf, y = Inf, label = r2_label, color = "grey40", vjust = 4, hjust = 1.25) + 
    theme_bw()
}

#菌群预测表型与实际表型作图
plot_y_vs_predict_y_phenotype <- function(y,predict_y,mae,r_squared){
  df <- data.frame(y=y,predicted_y=predict_y)
  mae = (abs(y-predict_y))/y
  mae_label = paste("MAE: ", as.character(round(mean(mae),3))," ± ",as.character(round(sd(mae),3)),sep = "")
  r2_label = paste("R squared: ", as.character(round(mean(r_squared)*100,2))," ± ",as.character(round(sd(r_squared)*100,2))," %",sep = "")
  ggplot(df, aes(x = y, y = predicted_y)) + 
    ylab(paste("Predicted ", sep = "")) + 
    xlab(paste("Observed ", sep = "")) + 
    geom_point(alpha = 0.1) + 
    geom_smooth(method = "loess", span = 1) + 
    annotate(geom = "text", x = Inf, y = Inf, label = mae_label, color = "grey40", vjust = 2.4, hjust = 1.5) + 
    annotate(geom = "text", x = Inf, y = Inf, label = r2_label, color = "grey40", vjust = 4, hjust = 1.25) + 
    theme_bw()
}
##############################
