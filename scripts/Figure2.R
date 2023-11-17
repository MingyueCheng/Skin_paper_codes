options (stringsAsFactors = FALSE)
##1. input
species <- read.table("./data/species.txt", header = TRUE, row.names = 1, sep = "\t")
meta <- read.table("./data/metagenome_meta_246_individuals.txt", header = TRUE, row.names = 1, sep = "\t")

##2. Figure 2A_1 correlations between species beta-diversity and skin phenomes)
library(vegan)
library(ggplot2)

permanova <- adonis2(t(species[,rownames(meta)]) ~ Age+Lightening+Grease+
                        Poryphyrins+Texture+Melanin+Pore+Spot+Wrinkle+Hemoglobin+
                        UVspot,
                      data = meta, permutations=9999, method = "bray", by="margin")

permanova_results <- data.frame(row.names=rownames(permanova),"R_square" = permanova$R2, "P"=permanova$`Pr(>F)`)

#This step would cost a lot of time, we have already done it and saved the results in the file "permanova_results.txt"
permanova_results <- read.table("./intermediate_results/permanova_results.txt", header = TRUE, row.names = 1, sep = "\t")
permanova_results$var <- rownames(permanova_results)
permanova_results<-permanova_results[order(permanova_results$R_square,decreasing = T),]
permanova_results$var <- factor(permanova_results$var,levels=permanova_results$var,ordered = T)

figure2A_1 <- ggplot(permanova_results[3:nrow(permanova_results),],aes(x=var,y=R_square))+geom_bar(stat="identity",width=0.8)+
  theme(axis.text.x=element_text(angle=90,vjust = 1,hjust = 1))

sum(permanova_results[3:8,1]) #0.0748419

ggsave(figure2A_1,file="./figures/figure2A_1.pdf",width=6,height=8)
##3. Figure 2A_2: correlations between species alpha-diversity and skin phenomes
species_shannon <- as.data.frame(diversity(t(species),index='shannon'))
colnames(species_shannon) <- "shannon"
species_shannon <- cbind(species_shannon,
                                meta[rownames(species_shannon),])

get_shannon_phynotypic_correlation <- function(data){
  spearman_cor <- c()
  spearman_p <- c()
  
  for (i in (2:ncol(data))){
  spearman <-cor.test(data[,1],data[,i],method = "spearman")
  spearman_cor <- c(spearman_cor,spearman$estimate)
  spearman_p <- c(spearman_p,spearman$p.value)
  }

  results <- rbind(spearman_cor,spearman_p)
  colnames(results) <- colnames(data)[-1]
  results
}

species_shannon_phynotype_correlation <- get_shannon_phynotypic_correlation(species_shannon)

#order the phenomes according to the previous permanova results
species_shannon_phynotype_correlation<-species_shannon_phynotype_correlation[,rownames(permanova_results[3:nrow(permanova_results),])]

#plot
library(pheatmap)
figure2A_2 <- pheatmap(species_shannon_phynotype_correlation[1,],
         cluster_cols = FALSE,cluster_rows = FALSE,
         color = colorRampPalette(c("#4EA74A","white","#8E4B99"))(100),
         labels_row =colnames(species_shannon_phynotype_correlation)
         )

ggsave("./figures/figure2A_2.pdf",figure2A_2,width=2,height=4)

##4. Figure 2B: correlations between species and skin phenomes
library(Maaslin2)
library(reshape)

input_data <- as.data.frame(t(species))
input_metadata <- meta
fit_data_random_age <- Maaslin2(
  input_data, input_metadata,'./intermediate_results/Maaslin2_output_random_age_0.2prevalance',
  fixed_effects = c('Lightening','Grease','Poryphyrins',
                     'Texture','Melanin','Pore',
                     'Spot','Wrinkle','Hemoglobin',
                     'UVspot'),
  random_effects = c('Age'),
  analysis_method = "LM",
  max_significance = 0.1,
  min_prevalence = 0.2
  )

maaslin2_sig <- read.table("./intermediate_results/Maaslin2_output_random_age_0.2prevalance/all_results.tsv",header=T,row.names=NULL, dec=".", sep="\t",check.names=T,quote="")

maaslin2_sig_fdr_0.1 <- maaslin2_sig[maaslin2_sig$qval<0.1,]
maaslin2_sig_fdr_0.1_split <- split(maaslin2_sig_fdr_0.1,maaslin2_sig_fdr_0.1$metadata)

#get all the taxa in the results
maaslin2_sig_fdr_0.1_taxa_list <- maaslin2_sig_fdr_0.1$feature
maaslin2_sig_fdr_0.1_taxa_list <- maaslin2_sig_fdr_0.1_taxa_list[!duplicated(maaslin2_sig_fdr_0.1_taxa_list)]

#record q_value
maaslin2_sig_fdr_0.1_plotdata <- as.data.frame(matrix(NA,length(maaslin2_sig_fdr_0.1_taxa_list),10))
rownames(maaslin2_sig_fdr_0.1_plotdata) <- maaslin2_sig_fdr_0.1_taxa_list
colnames(maaslin2_sig_fdr_0.1_plotdata) <- rownames(permanova_results[4:nrow(permanova_results),])

for (data in maaslin2_sig_fdr_0.1_split){
  index = unique(data$metadata)
  maaslin2_sig_fdr_0.1_plotdata[data$feature,index] = data$qval
}

#record cor value
maaslin2_sig_fdr_0.1_plotdata_cor <- as.data.frame(matrix(NA,length(maaslin2_sig_fdr_0.1_taxa_list),10))
rownames(maaslin2_sig_fdr_0.1_plotdata_cor) <- maaslin2_sig_fdr_0.1_taxa_list
colnames(maaslin2_sig_fdr_0.1_plotdata_cor) <- rownames(permanova_results[4:nrow(permanova_results),])

for (data in maaslin2_sig_fdr_0.1_split){
  index = unique(data$metadata)
  maaslin2_sig_fdr_0.1_plotdata_cor[data$feature,index] = data$coef
}

scatter_plot <- function(data){
  data[data< 0.1&data>= 0.05] = 1
  data[data< 0.05&data>= 0.01] = 2
  data[data< 0.01&data>= 0.001] = 3
  data[data< 0.001] = 4
  data <- melt(t(data))
  colnames(data) <- c("skin_index","taxa","fdr")
  data$taxa2 <- sub(".+\\.s__","",data$taxa)
  data
}
maaslin2_sig_fdr_0.1_plotdata_data <- scatter_plot(maaslin2_sig_fdr_0.1_plotdata)

scatter_plot_cor <- function(data){
  data[data<0] = -1
  data[data>=0] = 1
  data <- melt(t(data))
  colnames(data) <- c("skin_index","taxa","coef")
  data$taxa2 <- sub(".+\\.s__","",data$taxa)
  data
}
maaslin2_sig_fdr_0.1_plotdata_cor_data <- scatter_plot_cor(maaslin2_sig_fdr_0.1_plotdata_cor)

#combine
maaslin2_sig_fdr_0.1_plotdata_data$coef <- maaslin2_sig_fdr_0.1_plotdata_cor_data$coef

maaslin2_sig_fdr_0.1_plotdata_data$taxa <- sub(".+\\.s__","",maaslin2_sig_fdr_0.1_plotdata_data$taxa)
figure2B <- ggplot(maaslin2_sig_fdr_0.1_plotdata_data,aes(skin_index,taxa))+
  geom_point(aes(size=as.numeric(as.factor(fdr)),color=as.factor(coef)))+
  scale_colour_manual(values=(c("#4EA74A","#8E4B99")))+
  scale_size_area(max_size = 6) +
  theme(axis.text.x=element_text(angle=90,vjust = 1,hjust = 1))

ggsave("./figures/figure2B.pdf",figure2B,width=8,height=8)

##5. Figure 2C: Prediction of shin phenomes by skin microbiome
## 5.1 functions for modeling
library(ranger)
library(crossRanger)

#RF regression function
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

#Feature selection
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

#Save feature importance
save_feature_importance <- function(feature_importance, outpath_dir){
  feature_importance$species <- rownames(feature_importance)
  feature_importance <- dplyr::select(feature_importance, "species", everything())
  write.table(feature_importance, outpath_dir, sep="\t", row.names = F, quote = F, col.names = T)
}

#Plot
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

#pipeline

mytheme <- theme(
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_blank(),
  legend.position="none")

outpath <- "./intermediate_results/species_predict_phenome"

#phenotype
phenotype <- meta[,-1]
mep <- t(species)
colnames(mep) <- gsub("\\|",".",colnames(mep))
#filter out the features with zero variance
mep <- mep[,which(apply(mep,2,var)!=0)]
cat("The number of fitlered variables (removed variables with zero variance) : ", ncol(mep) ,"\n")
#filter out the features with absent over 99% of the samples
Zero.p <- 0.99
mep <- mep[,which(colSums(mep==0)<Zero.p*nrow(mep))]
cat("The number of variables (removed variables containing over ", Zero.p," zero) in training data: ", ncol(mep) ,"\n")

#species predict phenotype (here we use Lighenting as an example)
identical(rownames(mep), rownames(phenotype))
x <- mep
y <- phenotype$Lightening

#5folds random forests
set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)

#Feature importance
feature_importance <- as.data.frame(result$importances)
feature_importance$sd <- apply(feature_importance[,1:5],1,sd)
feature_importance <- feature_importance[order(feature_importance$rank),]
#Save importance table
save_feature_importance(feature_importance, paste0(outpath, "Lightening_feature_importance.txt"))
#Plot importance
p <- ggplot(feature_importance[1:10,], aes(x = reorder(rownames(feature_importance[1:10,]),mean), y = mean))+#reorder：将x轴按mean的大小排序
  geom_bar(stat="identity", color= "#BEBDBD", fill="#BEBDBD")+#stat=identity-柱子高度代表y值大小；stat=count-柱子高度代表数据的个数
  coord_flip()+#x,y轴互换
  labs(x = 'Feature', y = 'Importance')+#设置横纵坐标轴名称
  mytheme+#设置主题:去掉灰色背景+边框+图注
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),width = .3,color="black")#添加误差线(error bar),此处用的标准差
p
#Select the most important features, re-run the RF regression
set.seed(1)
repeated_top_n_perf <- repeated_cv_feature_select(x, y, feature_importance, nfolds = 10, repeats = 5)
repeated_top_n_perf$n_features <- log2(repeated_top_n_perf$n_features)
p <- ggplot(repeated_top_n_perf)+
  geom_line(aes(x=repeated_top_n_perf[,1],y=repeated_top_n_perf[,2]), color="#BDBDBD", size=1.5)+
  geom_line(aes(x=repeated_top_n_perf[,1],y=repeated_top_n_perf[,3]), color="#BDBDBD", size=1.5)+
  geom_line(aes(x=repeated_top_n_perf[,1],y=repeated_top_n_perf[,4]), color="#BDBDBD", size=1.5)+
  geom_line(aes(x=repeated_top_n_perf[,1],y=repeated_top_n_perf[,5]), color="#BDBDBD", size=1.5)+
  geom_line(aes(x=repeated_top_n_perf[,1],y=repeated_top_n_perf[,6]), color="#BDBDBD", size=1.5)+
  geom_line(aes(x=repeated_top_n_perf[,1],y=repeated_top_n_perf[,"mean"]), size=1.5)+
  scale_x_continuous(breaks=(repeated_top_n_perf$n_features), labels = 2^(repeated_top_n_perf$n_features))+
  xlab("# of features used ") + 
  ylab("MAE") + 
  theme_bw()
p
#average MAE
p <- ggplot(repeated_top_n_perf)+
  geom_line(aes(x=repeated_top_n_perf[,1],y=repeated_top_n_perf[,"mean"]), size=1.5)+
  scale_x_continuous(breaks=(repeated_top_n_perf$n_features), labels = 2^(repeated_top_n_perf$n_features))+
  xlab("# of features used ") + 
  ylab("MAE") + 
  theme_bw()
p
#Add error bar
p <- p + geom_errorbar(aes(x = repeated_top_n_perf[,1],
                           y = repeated_top_n_perf[,"mean"],
                           ymin = mean - sd,
                           ymax = mean + sd),width = .3,color="black")#添加误差线(error bar),此处用的标准差
p
#Final select the 32 features
select_feature <- rownames(feature_importance[which(feature_importance$rank<=32),])
identical(rownames(mep),rownames(phenotype))
#选取前32个特征建模
x <- mep[, select_feature]
y <- phenotype$Lightening

set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)
#保存随机森林结果
result_file<-paste(outpath, "Lightening_select_feature.RData", sep="")
save(result, file=result_file)
#result <- get(load(result_file))
#species predict phenotype：random forests result--plot
p <- plot_y_vs_predict_y_phenotype(result$y, result$predicted, result$MAE, result$R_squared)
ggsave("./figures/figure2C.pdf",p, width=8, height=8)
