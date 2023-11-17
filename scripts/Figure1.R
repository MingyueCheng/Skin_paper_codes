options (stringsAsFactors = FALSE)

library("factoextra")
library("FactoMineR")
library("corrplot")
library("ggplot2")
library("cowplot")


####### Figure1 B

# 0. plotting theme
age_group_color <- c("#D5221E","#3677AD","#4EA74A","#8E4B99","#EF7B19")
theme <- theme(
  plot.title=element_text(size=18, face="bold", colour="black"),
  axis.title.x=element_text(size=15, colour="black"),
  axis.title.y=element_text(size=15, angle=90, colour="black"),
  axis.text.x=element_text(size=8,angle=45,colour="black",vjust=0.5),
  axis.text.y=element_text(size=8,colour="black"),
  axis.ticks=element_line(colour="black",size=0.5),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_blank(),
  axis.line=element_line(size=0.5))

# 1. input skin phenome data
skin_1000 <- read.table("./data/skin_imaging_phenome_998_individuals.txt",header=T,row.names=1, dec=".", sep="\t",check.names=T,quote="")

# 2. separate gender, age, and age group
skin_1000_meta <- skin_1000[,1:2]
skin_1000_meta$Age_group <- skin_1000$Age_group
skin_1000_data <- skin_1000[,3:12]
skin_1000_data <- apply(skin_1000_data,1:2,function(x){x<-as.numeric(x)})

# 3. PCA analysis
PCA_data <- skin_1000_data
PCA_data <- scale(PCA_data)
rownames(PCA_data) <- rownames(skin_1000_meta)
res.pca <- PCA(PCA_data, graph = FALSE,ncp=10)
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE)
var <- get_pca_var(res.pca)

# 3.1 Plotting Figure 1B part1
p1 = fviz_pca_biplot(res.pca, 
                col.ind = skin_1000_meta[rownames(PCA_data),]$Age_group,
                palette = c("#D5221E","#3677AD","#4EA74A","#8E4B99","#EF7B19"), 
                addEllipses = FALSE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species",select.var = list(name = NULL, cos2 = 10, contrib = NULL))+
  stat_ellipse(type = "t",level=0.80,alpha=0.5,aes(x=res.pca$ind$coord[,"Dim.1"],y=res.pca$ind$coord[,"Dim.2"], 
                                                                    color=skin_1000_meta[rownames(PCA_data),]$Age_group))+
  scale_fill_manual(values = age_group_color)+
  scale_shape_manual(values=c(1,1,1,1,1))

# 3.2 Plotting Figure 1B part2
group_comparisons <- list(c("group1","group2"),c("group2","group3"),c("group3","group4"),c("group4","group5"))

library(ggpubr)
p2 = ggplot(pca_dims, aes(x=Age_group,y=Dim.1)) +
  geom_boxplot(aes(fill=Age_group))+theme+scale_fill_manual(values=age_group_color)+ stat_compare_means(comparisons = group_comparisons,methods="wilcox.test")

# store the plots of p1 p2

figure1b = plot_grid(p1,p2,ncol=1)
ggsave("./figures/figure1B.pdf",figure1b,width=8,height=12)


# 3.3 wilcox test to compare each age group against PC1

pca_dims <- as.data.frame(res.pca$ind$coord)
pca_dims$Age_group <- skin_1000_meta[rownames(pca_dims),]$Age_group
pca_dims$Gender <- skin_1000_meta[rownames(pca_dims),]$Gender
pca_dims$Product <- skin_1000_meta[rownames(pca_dims),]$Product
wilcox.test(pca_dims[pca_dims$Age_group=="group1",]$Dim.1,pca_dims[pca_dims$Age_group=="group2",]$Dim.1)
#p-value = 0.006933
wilcox.test(pca_dims[pca_dims$Age_group=="group1",]$Dim.1,pca_dims[pca_dims$Age_group=="group3",]$Dim.1)
#p-value = 2.038e-15
wilcox.test(pca_dims[pca_dims$Age_group=="group1",]$Dim.1,pca_dims[pca_dims$Age_group=="group4",]$Dim.1)
#p-value < 2.2e-16
wilcox.test(pca_dims[pca_dims$Age_group=="group1",]$Dim.1,pca_dims[pca_dims$Age_group=="group5",]$Dim.1)
#p-value = 4.476e-08


wilcox.test(pca_dims[pca_dims$Age_group=="group2",]$Dim.1,pca_dims[pca_dims$Age_group=="group3",]$Dim.1)
#p-value = 3.794e-14
wilcox.test(pca_dims[pca_dims$Age_group=="group2",]$Dim.1,pca_dims[pca_dims$Age_group=="group4",]$Dim.1)
#p-value < 2.2e-16
wilcox.test(pca_dims[pca_dims$Age_group=="group2",]$Dim.1,pca_dims[pca_dims$Age_group=="group5",]$Dim.1)
#p-value = 2.365e-08

wilcox.test(pca_dims[pca_dims$Age_group=="group3",]$Dim.1,pca_dims[pca_dims$Age_group=="group4",]$Dim.1)
#p-value = 4.26e-13
wilcox.test(pca_dims[pca_dims$Age_group=="group3",]$Dim.1,pca_dims[pca_dims$Age_group=="group5",]$Dim.1)
#p-value = 4.267e-07

wilcox.test(pca_dims[pca_dims$Age_group=="group4",]$Dim.1,pca_dims[pca_dims$Age_group=="group5",]$Dim.1)
#p-value = 0.004661

####### Figure1 C

# all individuals colored by gender
figure1c = fviz_pca_ind(res.pca, 
                col.ind = skin_1000_meta[rownames(PCA_data),]$Gender,
                palette = c("#D5221E","#3677AD"), 
                addEllipses = FALSE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species",select.var = list(name = NULL, cos2 = 10, contrib = NULL))+
  stat_ellipse(type = "t",level=0.80,alpha=0.5,aes(x=res.pca$ind$coord[,"Dim.1"],y=res.pca$ind$coord[,"Dim.2"], 
                                                                    color=skin_1000_meta[rownames(PCA_data),]$Gender))+
  scale_fill_manual(values = age_group_color)+
  scale_shape_manual(values=c(1,1))

ggsave("./figures/figure1C.pdf",figure1c,width=8,height=8)
