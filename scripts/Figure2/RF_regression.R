library(ggplot2)
library(ggsignif)
library(dplyr)
library(caret)

outpath <- "figures/Figure2/"
dir.create(outpath)

#figure theme
mytheme <- theme(
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_blank(),
  legend.position="none")


#998 samples
#998_samples_metadata
#data-input
metadata_998_dir <- "data/skin_imaging_phenome_998_individuals.txt"
metadata_998 <- read.table(metadata_998_dir, header=T, sep="\t", quote="", comment.char="")

#Random forest: Age prediction-based on phenotype (998 samples)
x <- metadata_998[,c(-1,-2,-13)]
y <- metadata_998$Age
set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)

#feature importance
feature_importance <- as.data.frame(result$importances)
feature_importance$sd <- apply(feature_importance[,1:5],1,sd)
feature_importance <- feature_importance[order(feature_importance$rank),]
save_feature_importance(feature_importance, paste0(outpath, "phenotype_age/998-samples/feature_importance.txt"))
#feature selection
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
#error bar (standard deviation)
p <- p + geom_errorbar(aes(x = repeated_top_n_perf[,1],
                           y = repeated_top_n_perf[,"mean"],
                           ymin = mean - sd,
                           ymax = mean + sd),width = .3,color="black")
p
#feature selection: best performance with all 10 features
select_feature <- rownames(feature_importance[which(feature_importance$rank<=10),])
#modeling with all 10 features
x <- metadata_998[, select_feature]
y <- metadata_998$Age

set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)
#save RF result
result_file<-paste(outpath, "phenotype_age/998-samples/phenotype_age_998_rf_feature_Select.RData", sep="")
save(result, file=result_file)
#result <- get(load(result_file))
#998 phenotype_age：random forests result--plot
predict_result <- data.frame(y=result$y, predict=result$predicted)
predict_result$SampleID <- rownames(metadata_998)
predict_result <- dplyr::select(predict_result, "SampleID",everything())
p <- plot_y_vs_predict_y(result$y, result$predicted, result$MAE,result$R_squared)
p
zoom=2.5
ggsave(paste0(outpath, "phenotype_age/998-samples/phenotype_predict_age_998_MAE.pdf"), 
       p, width=75*zoom, height=75*zoom, units="mm")
#Spearman correlation
cor.test(result$y, result$predicted, method="spearman")





#246 samples
#246_samples_metadata
#data-input
metadata_dir <- "./data/metagenome_meta_246_individuals.txt"
metadata <- read.table(metadata_dir, header=T, sep="\t", quote="", comment.char="",row.names = 1)

#246_species_data（annotation by metaphlan）
#data-input
metaphlan_data <- "./data/species.txt"
mep <- read.table(metaphlan_data, header=T, sep="\t", quote="", comment.char="", row.names = 1)
mep <- as.data.frame(t(mep))
colnames(mep) <- gsub("\\|",".",colnames(mep))
#Filter out features with a variance of 0
mep <- mep[,which(apply(mep,2,var)!=0)]
cat("The number of fitlered variables (removed variables with zero variance) : ", ncol(mep) ,"\n")
#Filter out features with an abundance of 0 in more than 99% of the samples
Zero.p <- 0.99
mep <- mep[,which(colSums(mep==0)<Zero.p*nrow(mep))]
cat("The number of variables (removed variables containing over ", Zero.p," zero) in training data: ", ncol(mep) ,"\n")
identical(rownames(metadata),rownames(mep))



#Random forest: Age prediction-based on species (246 samples)
x <- mep
y <- metadata$Age
set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)

#feature importance
feature_importance <- as.data.frame(result$importances)
feature_importance$sd <- apply(feature_importance[,1:5],1,sd)
feature_importance <- feature_importance[order(feature_importance$rank),]
save_feature_importance(feature_importance, paste0(outpath, "microbiota_age/feature_importance.txt"))
#feature selection
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
#error bar (standard deviation)
p <- p + geom_errorbar(aes(x = repeated_top_n_perf[,1],
                           y = repeated_top_n_perf[,"mean"],
                           ymin = mean - sd,
                           ymax = mean + sd),width = .3,color="black")
p
#feature selection: best performance with top64 features
select_feature <- rownames(feature_importance[which(feature_importance$rank<=64),])
#modeling with top64 features
x <- mep[, select_feature]
y <- metadata$Age

set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)
#save RF result
result_file<-paste(outpath, "microbiota_age/species_age_rf_feature_Select.RData", sep="")
save(result, file=result_file)
#result <- get(load(result_file))
#microbiota_age：random forests result--plot
predict_result <- data.frame(y=result$y, predict=result$predicted)
predict_result$SampleID <- rownames(mep)
predict_result <- dplyr::select(predict_result, "SampleID",everything())
p <- plot_y_vs_predict_y(result$y, result$predicted, result$MAE,result$R_squared)
p
zoom=2.5
ggsave(paste0(outpath, "microbiota_age/species_predict_age_MAE.pdf"), 
       p, width=75*zoom, height=75*zoom, units="mm")
#Spearman correlation
cor.test(result$y, result$predicted, method="spearman")




#Random forest: Age prediction-based on phenotype
x <- metadata[,c(-1)]
y <- metadata$Age
set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)

#feature importance
feature_importance <- as.data.frame(result$importances)
feature_importance$sd <- apply(feature_importance[,1:5],1,sd)
feature_importance <- feature_importance[order(feature_importance$rank),]
save_feature_importance(feature_importance, paste0(outpath, "phenotype_age/246-samples/feature_importance.txt"))
#feature selection
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
#error bar (standard deviation)
p <- p + geom_errorbar(aes(x = repeated_top_n_perf[,1],
                           y = repeated_top_n_perf[,"mean"],
                           ymin = mean - sd,
                           ymax = mean + sd),width = .3,color="black")
p
#feature selection: best performance with top8 features
select_feature <- rownames(feature_importance[which(feature_importance$rank<=8),])
#modeling with top8 features
x <- metadata[, select_feature]
y <- metadata$Age

set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)
#save RF result
result_file<-paste(outpath, "phenotype_age/246-samples/phenotype_age_rf_feature_Select.RData", sep="")
save(result, file=result_file)
#result <- get(load(result_file))
#phenotype_age：random forests result--plot
predict_result <- data.frame(y=result$y, predict=result$predicted)
predict_result$SampleID <- rownames(mep)
predict_result <- dplyr::select(predict_result, "SampleID",everything())
p <- plot_y_vs_predict_y(result$y, result$predicted, result$MAE,result$R_squared)
p
zoom=2.5
ggsave(paste0(outpath, "phenotype_age/246-samples/phenotype_predict_age_MAE.pdf"),
       p, width=75*zoom, height=75*zoom, units="mm")
#Spearman correlation
cor.test(result$y, result$predicted, method="spearman")




#Random forest: Age prediction-based on species+phenotype
mep_phenotype <- cbind(mep, metadata[,c(-1)])
x <- mep_phenotype
y <- metadata$Age
set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)
#feature importance
feature_importance <- as.data.frame(result$importances)
feature_importance$sd <- apply(feature_importance[,1:5],1,sd)
feature_importance <- feature_importance[order(feature_importance$rank),]
save_feature_importance(feature_importance, paste0(outpath, "integration_age/feature_importance.txt"))
#feature selection
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
#error bar (standard deviation)
p <- p + geom_errorbar(aes(x = repeated_top_n_perf[,1],
                           y = repeated_top_n_perf[,"mean"],
                           ymin = mean - sd,
                           ymax = mean + sd),width = .3,color="black")
p
#feature selection: best performance with top4 features
select_feature <- rownames(feature_importance[which(feature_importance$rank<=4),])
#modeling with top4 features
x <- mep_phenotype[, select_feature]
y <- metadata$Age

set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)
#save RF result
result_file<-paste(outpath, "integration_age/integration_age_rf_feature_Select.RData", sep="")
save(result, file=result_file)
#result <- get(load(result_file))
#integration_age：random forests result--plot
predict_result <- data.frame(y=result$y, predict=result$predicted)
predict_result$SampleID <- rownames(mep)
predict_result <- dplyr::select(predict_result, "SampleID",everything())
p <- plot_y_vs_predict_y(result$y, result$predicted, result$MAE,result$R_squared)
p
zoom=2.5
ggsave(paste0(outpath, "integration_age/species_phenotype_predict_age_MAE.pdf"), 
       p, width=75*zoom, height=75*zoom, units="mm")
#Spearman correlation
cor.test(result$y, result$predicted, method="spearman")






#Phenotype prediction using species
phe <- "UVspot"
x <- mep
y <- metadata[,phe]

#5folds random forests
set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)
#feature importance
feature_importance <- as.data.frame(result$importances)
feature_importance$sd <- apply(feature_importance[,1:5],1,sd)
feature_importance <- feature_importance[order(feature_importance$rank),]
save_feature_importance(feature_importance, paste0(outpath, "microbiota_phenotype/",phe,"/feature_importance.txt"))
#feature selection
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
#error bar (standard deviation)
p <- p + geom_errorbar(aes(x = repeated_top_n_perf[,1],
                           y = repeated_top_n_perf[,"mean"],
                           ymin = mean - sd,
                           ymax = mean + sd),width = .3,color="black")
p
#feature selection: best performance with topxx features
select_feature <- rownames(feature_importance[which(feature_importance$rank<=64),])
#modeling with topxx features
x <- mep[, select_feature]
y <- metadata[,phe]

set.seed(1)
result <- rf.regression.cross_validation(x, y, nfolds = 5)
#save RF result
result_file<-paste(outpath, "microbiota_phenotype/",phe,"/rf_feature_Select.RData", sep="")
save(result, file=result_file)
#result <- get(load(result_file))
#microbiota pre phenotype：random forests result--plot
predict_result <- data.frame(y=result$y, predict=result$predicted)
predict_result$SampleID <- rownames(mep)
predict_result <- dplyr::select(predict_result, "SampleID",everything())
p <- plot_y_vs_predict_y_phenotype(result$y, result$predicted, result$MAE, result$R_squared)
p
zoom=2.5
ggsave(paste0(outpath, "microbiota_phenotype/",phe,"/predict_MAE.pdf"), 
       p, width=75*zoom, height=75*zoom, units="mm")
#Spearman correlation
cor.test(result$y, result$predicted, method="spearman")

