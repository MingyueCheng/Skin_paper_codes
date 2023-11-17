#library
library (xlsx) 
library (data.table)
library (ggplot2)
library (plyr)
library(getopt)
library(cowplot)
library(ggpubr)

options (stringsAsFactors = FALSE)

inpath = "./data/input/Figure3_input"
outpath = "./intermediate_results/Figure3_results"
Phe = "Pore"
e = 3 #least samples where MGSs/KOs occur
final.fdr.cutoffs1 = 0.01 ### FDR thresholds1
final.fdr.cutoffs2 = 0.1 ### FDR thresholds2
cor_method = "spearman"

###theme
mytheme <- theme(
  axis.title.x = element_text(size=9, colour="black"),
  axis.title.y = element_text(size=9, colour="black"),
  axis.text.x=element_text(size=9,colour="black"),
  axis.text.y=element_text(size=9,colour="black"),
  axis.ticks=element_line(colour="black",size=0.5),
  axis.ticks.length.x=unit(0.2, "cm"),
  axis.ticks.length.y=unit(0.2, "cm"),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_blank(),
  axis.line=element_line(size=0.5))

####step1_input####
### - phenotypes
phenotypes = read.table (file = paste (inpath, "phenotypes.tsv",sep="/"), sep = "\t", row.names = 1, header = T)
used_samples = rownames(phenotypes) #you could also use a subset of samples
### - MGS_abundance.tab #get
mgs_abundance = read.table (file = paste (inpath, "MGS_abundance.tsv",sep="/"), sep = "\t", row.names = 1, header = T)
mgs_abundance <- as.data.frame(t(mgs_abundance))

### - MGS_taxonomy.tab
mgs_taxonomy  = read.table (file = paste (inpath, "MGS_taxonomy.tsv",sep="/"), sep = "\t", row.names = 1, header = T)

### - KEGG_modules.tab #get
tmp = read.table (file = paste (inpath, "KEGG_modules.tsv",sep="/"), sep = "\t")
koann = strsplit (tmp[,3], split = " ")
names (koann) = tmp[,1]
module_mapping = tmp[,2] 
names (module_mapping) = tmp[,1] ; rm (tmp)

### remove KEGG references in square brackets for more clean names for plotting
module_mapping_clean = sapply(module_mapping, function(x) strsplit(x, " \\[")[[1]][1])

### - KO_abundance.tab #get
ko_abundance = read.table (file = paste (inpath, "KO_abundance.tsv",sep="/"), sep = "\t", row.names = 1, header = T)
ko_abundance <- as.data.frame(t(ko_abundance))
### - KO_to_MGS.tab #get

tmp = read.table (file = paste (inpath, "KO_to_MGS.tsv", sep="/"), sep = "\t", strip.white = T)
KO2MGS = strsplit (tmp[,2], split = " ")
names (KO2MGS) = tmp[,1] ; rm (tmp)

### - gene_to_KO.tab #get
tmp = read.table (file = paste (inpath, "KO_to_gene.tsv", sep="/"), sep = "\t", header = F)
KO2gene = strsplit (tmp[,2], split = " ")
names(KO2gene) = tmp[,1]; rm (tmp)

### - MGS_to_gene.tab #get
tmp = read.table (file = paste (inpath, "MGS_to_gene.tsv", sep="/"), sep = "\t", strip.white = T)
MGS2gene = strsplit (tmp[,2], split = " ")
names (MGS2gene) = tmp[,1] ; rm (tmp)

### - gene_abundance_sub.tab #get
gene_abundance_sub = data.frame(fread (paste (inpath, "gene_abundance_sub.tsv", sep="/"), sep = "\t", header = T), row.names = 1)
gene_abundance_sub = as.data.frame(gene_abundance_sub)

####step2_preprocess####
### MGS sparsity filter step - exclude MGSs that occur in <e of individuals
test = apply (mgs_abundance, 2, function(x) length (x [x != 0])) 
incl = names (test [test >= e])
length (incl)
MGSs = incl
mgs_abundance = mgs_abundance [,incl]
mgs_taxonomy = mgs_taxonomy [incl,]
rm (test, incl)

### KO sparsity filter step - exclude KOs that occur in <3 of the control individuals
test = apply (ko_abundance,2, function (x) length (x [x != 0]))  
incl = names (test [test >= e])
length (incl)
ko_abundance = ko_abundance [,incl]
rm (test, incl)

###make rownames consistent
mgs_abundance <- mgs_abundance[used_samples,]
ko_abundance <- ko_abundance[used_samples,]
gene_abundance_sub <- gene_abundance_sub[,used_samples]

####step3_associations_MGSs_phenotypes####
tmpMat = array (NA, c (ncol (mgs_abundance), 2))
dimnames (tmpMat) [[1]] = colnames (mgs_abundance)
dimnames (tmpMat) [[2]] = c ("estimate", "p.value")
mgs.subset = colnames(mgs_abundance)

### Associating MGSs with phenotype without de-confounding for BMI    
tmpMat [mgs.subset, c ("estimate", "p.value")] =
  t (apply (mgs_abundance [used_samples, mgs.subset], MARGIN = 2, FUN = function (x) 
    unlist (cor.test (phenotypes [used_samples, Phe], x,
                      method = cor_method, use = "pairwise.complete.obs") [c ("estimate", "p.value")])))       

cor_Phe = list ()
cor_Phe [["MGSs"]] <- tmpMat
rm (tmpMat, mgs.subset)

####step4_associations_KEGG_modules_phenotypes####
tmpMat = array (NA, c (length (koann), 2))
dimnames (tmpMat) [[1]] = names (koann)
dimnames (tmpMat) [[2]] = c ("estimate", "p.value")

ko.subset = colnames(ko_abundance)

KO_cor = apply (ko_abundance [used_samples, ko.subset], MARGIN = 2, FUN = function (x) 
  cor (phenotypes [used_samples, Phe], x,
       method = cor_method, use = "pairwise.complete.obs"))      
KO_cor_Phe = KO_cor
### Then, test for difference in correlation coefficients between KOs in the
### KEGG module and all other KOs.
for (k in names (koann)) {
  incat =    na.omit (KO_cor [   names (KO_cor) %in% koann [[k]] ]) ### select all correlations between Phe and KOs in the KEGG module
  notincat = na.omit (KO_cor [! (names (KO_cor) %in% koann [[k]])]) ### select all correlations between Phe and KOs NOT in the KEGG module
  if (length (incat) > 0 & length (notincat) > 0) {
    x = wilcox.test (incat, notincat)
    tmpMat [k,"p.value"] = x$p.value
    tmpMat [k,"estimate"] = (median (incat, na.rm = T) - median (notincat, na.rm = T))
  }
}
rm (KO_cor)
cor_Phe [["keggmodules"]] <- tmpMat
rm (tmpMat, ko.subset)

####step5_save_phenotype_association####
###store results
tmp.excel.file = paste(outpath,Phe,".xlsx",sep="")
write.xlsx (paste ("Associations with ", Phe, ""), 
            sheetName = "info", file = tmp.excel.file, row.names = F, col.names = F)

for (tmp_name in names (cor_Phe)) {
  tmpMat = cor_Phe [[paste (tmp_name)]]
  out = cbind (tmpMat [,c("estimate", "p.value")], "p.adjust" = p.adjust (tmpMat [,"p.value"], method = "BH"))
  colnames (out) = paste (rep (c (Phe), each = 3), colnames (out), sep = "_")
  write.xlsx (x = out, file = tmp.excel.file, append = T, sheetName = paste (tmp_name, sep = "_"), row.names = T, col.names = T)            
}

###extract significant results
vartotest_sig <- list ()
for (tmp_name in names (cor_Phe)) {
  tmpMat = cor_Phe [[paste (tmp_name)]]
  vartotest_sig [[paste ("fdr", final.fdr.cutoffs1, sep = "_")]] [[tmp_name]] =
    c (names (which (p.adjust (tmpMat [ , "p.value"], method = "BH") < final.fdr.cutoffs1)))
  vartotest_sig [[paste ("fdr", final.fdr.cutoffs2, sep = "_")]] [[tmp_name]] =
    c (names (which (p.adjust (tmpMat [ , "p.value"], method = "BH") < final.fdr.cutoffs2)))
}

####step6_leave_one_out####
MGS.list = colnames(mgs_abundance)
KOsets = koann [vartotest_sig[[final.fdr.cutoffs1]]$keggmodules] 

### Calculate abundance sum of ALL genes annotated to a KO
KOprofile_FUN <- function (KO) {
  KOgenes_i <- unique (unlist (KO2gene [KO])) ### get all genes that are annotated to the KO 
  return (colSums (gene_abundance_sub [KOgenes_i,])) ###求出每个样本包含这个ko里的基因总数
}

### Calculate abundance sum of genes annotated to a KO while excluding the genes
### that are part of 'LeftOutMGS'
LOOKOprofile <- function (KO, LeftOutMGS) {
  if ( any (LeftOutMGS == unlist (KO2MGS [KO])) ){	KOgenes_i <- setdiff (unique (unlist (KO2gene [KO])), MGS2gene [[LeftOutMGS]])} #得到重合部分以外的gene
  else { KOgenes_i <- unique (unlist (KO2gene [KO])) } #对于每个KO来讲，如果里面包含任何一个去除的MGS的KO包含的gene 就不要那些基因了
  return (colSums (gene_abundance_sub [KOgenes_i,])) #去除后再计算每个样本包含这个ko里的基因总数
}

### The main function calculating the effect of leaving-one-MGS-out on the 
### association between a KEGG module and a phenotype of interest
LeaveOneMGSOut_FUN <- function (KOsets = KOsets, MGS.list = MGS.list, sample.list = sample.list, y, z = NULL) {
  ### number of KOs kinds per MGS
  KO_types_per_MGS <- lapply (KOsets, function (KOset) {
    mgses <- sapply (KOset, function (KO) { intersect (unique (unlist (KO2MGS [KO])), MGS.list)})#第一步看每个KO包含哪些物种
    table (unlist (mgses))})
  
  if  (is.null (z) == T) { 
    
    print ("No 'z' is provided, calculating spearman correlation")
    
    ### Code for spearman correlation, i.e. NOT adjusting for a confounding factor
    y = y [sample.list,]
    
    ### correlating the mean KO signal (calculated with 'KOprofile_FUN') to 
    ### phenotype, and then taking the median SCC for the module
    print ("calculating spearman correlation, using all genes")
    SCCallMGS<-sapply (KOsets, function (KOset) { 
      KOprofiles<-t (sapply (KOset, KOprofile_FUN))#对KOset里的每个KO求出样本里包含这个KO基因的总数。本来样本为列，KO为行，现在KO为行，样本为列
      return (median (apply (KOprofiles [rowSums (KOprofiles) > 0, sample.list], 1, function (x) { #只计算在所有样本中至少出现1次的KO，y是表型。计算spearman，取中位数
        cor.test (x =x, y = y, method = "spearman", use = "pairwise.complete.obs", exact = FALSE)$estimate }))) 
    })
    
    ### correlating the mean KO signal while excluding contribution from 'MGS'
    ### (calculated with 'LOOKOprofile') to phenotype, and then taking the median
    ### SCC for the module.
    ### Repeating for all MGS in 'MGS.list'
    print ("calculating spearman correlation, leaving-one-MGS-out")
    SCComitingMGS <- lapply (KOsets, function (KOset) {
      sapply (intersect (unique (unlist (KO2MGS[KOset])),MGS.list), function (MGS){ 
        KOprofiles <- (sapply (KOset, function (KO) { LOOKOprofile (KO, MGS) })) 
        return (median (cor (x = KOprofiles [sample.list, colSums (KOprofiles) > 0], y = y, method = "spearman", use = "pairwise.complete.obs"))) })
    }) 
    
    ### Summarizing outout
    DeltaSCCperMGS <- lapply (names (SCComitingMGS), function (N) { 
      SCC = SCCallMGS [[N]] 
      SCC.bgadj =  (median (na.omit (KO_cor_Phe [names (KO_cor_Phe) %in% KOsets [[N]]]), na.rm = T) 
                    - median (na.omit (KO_cor_Phe[! (names (KO_cor_Phe) %in% KOsets [[N]])]), na.rm = T))
      # summary(KO_cor_Phe); adjust % SCC effect for background distribution (i.e. the fact that median(SCC for phe vs modules) is negative and not = 0)
      data.frame (SCC = SCC, SCC.bgadj = SCC.bgadj, SCC_omiting_MGS = SCComitingMGS [[N]], DeltaMGS_SCC = SCC - SCComitingMGS [[N]],
                  pctSCCeffect = 100 * (SCC - SCComitingMGS [[N]]) / SCC, pctSCCeffect.bgadj = 100 * (SCC - SCComitingMGS [[N]]) / SCC.bgadj,
                  Distinct_KOs_in_MGS = as.vector (KO_types_per_MGS [[N]][names (SCComitingMGS [[N]])]), row.names = names (SCComitingMGS [[N]])
      ) [rev (order (pctSCCeffect = 100 * (SCC - SCComitingMGS[[N]]) / SCC)),]
    })
    
    names (DeltaSCCperMGS) <- names (KOsets)
    return (DeltaSCCperMGS)
    
  }
  
  else {
    
    print ("'z' is provided, calculating partial spearman correlation, adjusting for 'z'")
    
    ### Code for partial spearman correlation, i.e. adjusting for a confounding factor
    y = y [sample.list,]
    z = z [sample.list,]
    
    ### correlating the mean KO signal (calculated with 'KOprofile_FUN') to
    ### phenotype, and then taking the median SCC for the module
    print ("calculating partial spearman correlation, using all genes")
    partialSCCallMGS <- sapply (KOsets, function (KOset) {
      KOprofiles <- t (sapply (KOset, KOprofile_FUN))
      return (median (apply (KOprofiles [rowSums (KOprofiles) > 0,], 1, function (x) { pcor.test (x = x [sample.list], y = y, z = z, method = "spearman")$estimate })))
    })
    
    ### correlating the mean KO signal while excluding contribution from 'MGS'
    ### (calculated with 'LOOKOprofile') to phenotype, and then taking the median
    ### SCC for the module.
    ### Repeating for all MGS in 'MGS.list'
    print ("calculating partial spearman correlation, leaving-one-MGS-out")
    partialSCComitingMGS <- lapply (KOsets, function (KOset) {
      sapply (intersect (unique (unlist (KO2MGS [KOset])), MGS.list), function (MGS) { # looping over all MGGs in 'MGS.list'
        KOprofiles <- t (sapply (KOset, function (KO) {LOOKOprofile (KO, MGS) }))
        return (median (apply (KOprofiles [rowSums (KOprofiles) > 0,], 1, function (x) { pcor.test (x = x [sample.list], y = y, z = z, method="spearman")$estimate })))
      })
    })
    
    ### Summarizing outout
    partialDeltaSCCperMGS<-lapply (names (partialSCComitingMGS), function (N){
      SCC = partialSCCallMGS[[N]]
      SCC.bgadj =  (median (na.omit (KO_cor_partial[names (KO_cor_partial) %in% KOsets[[N]]]), na.rm=T)
                    - median (na.omit (KO_cor_partial[! (names (KO_cor_partial) %in% KOsets[[N]])]), na.rm=T))
      ### summary(KO_cor_Phe); adjust % SCC effect for background distribution (i.e. the fact that median(SCC for phenotype vs modules) is negative and not = 0)
      data.frame (SCC =SCC, SCC.bgadj = SCC.bgadj, SCC_omiting_MGS = partialSCComitingMGS [[N]], DeltaMGS_SCC = partialSCCallMGS [[N]] - partialSCComitingMGS [[N]],
                  pctSCCeffect = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC, pctSCCeffect.bgadj = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC.bgadj,
                  Distinct_KOs_in_MGS = as.vector (KO_types_per_MGS [[N]][names (partialSCComitingMGS [[N]])]),
                  row.names = names (partialSCComitingMGS [[N]])
      ) [rev (order (pctSCCeffect = 100 * (SCC - partialSCComitingMGS [[N]]) / SCC)),]
    })
    
    names (partialDeltaSCCperMGS) <- names (KOsets)
    return (partialDeltaSCCperMGS)
    
  }
  
}
### End functions

### Next, run these functions to compute the delta SCC values for the test in
### question - change here if partial correlation to account for a confounder is
### desired.

### Performing the Leave-one-MGS-out for phenotype
DeltaSCCperMGS <- LeaveOneMGSOut_FUN (KOsets = KOsets, MGS.list = MGS.list,
                                      sample.list = used_samples, ### the subset of individuals with complete information for phenotype and metagenomic data.
                                      y = phenotypes [,Phe, drop = F]) 
save (DeltaSCCperMGS, file = paste(outpath,Phe,"_delta_SCC_per_MGS.RData",sep=""))

### Performing the Leave-one-MGS-out for phenotype.bmi.adjusted
# partialDeltaSCCperMGS = LeaveOneMGSOut_FUN (KOsets = KOsets, MGS.list = MGS.list,
#	sample.list = ctrl.no.na.4MGS, y = phenotypes [,"phenotype.IR", drop = F], z = phenotypes [,"BMI.kg.m2", drop = F])
# save (partialDeltaSCCperMGS, file = "results/partial_delta_SCC_per_MGS.RData")
#Phe = "Wrinkle"
#DeltaSCCperMGS <- load(paste(outpath,Phe,"/",Phe,"_delta_SCC_per_MGS.RData",sep=""))
#KOsets = koann[names(DeltaSCCperMGS)]
####step7_extract_top_driver####

topX = 3
output = as.data.frame (matrix (NA, nrow = length (DeltaSCCperMGS), ncol= (5 + 6 * topX)))
colnames (output) = c ("Module description","Number of genes in the KEGG module","Median SCC of KOs in the module","Median SCC of KOs out of the module", "SCC background",
                       "Top1: MGS", "Top1: number of module genes in MGS", "Top1: SCC after excluding the MGS ",  "Top1: Change of SCC after excluding the MGS", "Top1: SCC effect", "Top1: SCC effect normalized by SCC background",
                       "Top2: MGS", "Top2: number of module genes in MGS", "Top2: SCC after excluding the MGS ",  "Top2: Change of SCC after excluding the MGS", "Top2: SCC effect", "Top2: SCC effect normalized by SCC background",
                       "Top3: MGS", "Top3: number of module genes in MGS", "Top3: SCC after excluding the MGS ",  "Top3: Change of SCC after excluding the MGS", "Top3: SCC effect", "Top3: SCC effect normalized by SCC background"
                      )
rownames (output) = names (DeltaSCCperMGS)
output[,"Module description"] = module_mapping [rownames (output)]
output [,"Number of genes in the KEGG module"] = sapply (rownames (output), function (i) length (KOsets [[i]]))

for (i in names (DeltaSCCperMGS)) {
  tmp.SCC_module        = as.numeric(DeltaSCCperMGS [[i]][, "SCC"][1])
  tmp.SCC_out_of_module = as.numeric(DeltaSCCperMGS [[i]][, "SCC"][1] -  DeltaSCCperMGS [[i]][, "SCC.bgadj"][1])
  tmp.SCC_bgadj         = as.numeric(DeltaSCCperMGS [[i]][, "SCC.bgadj"][1])
  output[i,3] <- tmp.SCC_module
  output[i,4] <- tmp.SCC_out_of_module
  output[i,5] <- tmp.SCC_bgadj

  if (length (rownames (DeltaSCCperMGS[[i]])) <= topX ) { 
    
    tmp.MGSs.top = rownames (DeltaSCCperMGS [[i]]) [1:topX] 
    tmp.MGSs.top [DeltaSCCperMGS [[i]][tmp.MGSs.top, "pctSCCeffect"] <= 0] <- NA ### set top.MGSs to NA if removing them are NOT INCREASING the effect
    tmp.MGSs = c (tmp.MGSs.top) 
    
  } 
  
  else {
    
    tmp.MGSs.top = rownames (DeltaSCCperMGS [[i]]) [1:topX]
    tmp.MGSs.top [DeltaSCCperMGS [[i]][tmp.MGSs.top, "pctSCCeffect"] <= 0] <- NA ### set top.MGSs to NA if removing them are NOT INCREASING the effect
    tmp.MGSs = c (tmp.MGSs.top)
    
  }
  
  ## Add taxonomically annotation (for those MGSs that have one)
  tmp.MGSs.names = sapply (tmp.MGSs, function (i) ifelse (i %in% rownames (mgs_taxonomy),
                                                          paste (i, ": ", mgs_taxonomy [i, "Taxonomy.at.species.level"], " classified: ", mgs_taxonomy [i, "Classified.at.species.level"]),
                                                          paste (i, ": ", "Unknown")))
  
  tmp.KOsInMGS          = DeltaSCCperMGS [[i]][ tmp.MGSs, "Distinct_KOs_in_MGS"]
  tmp.SCC_omiting_MGS   = DeltaSCCperMGS [[i]][ tmp.MGSs, "SCC_omiting_MGS"]
  tmp.SCC_delta         = DeltaSCCperMGS [[i]][ tmp.MGSs, "DeltaMGS_SCC"]
  tmp.SCC_effect        = DeltaSCCperMGS [[i]][ tmp.MGSs, "pctSCCeffect"]
  tmp.SCC_effect_bgadj  = DeltaSCCperMGS [[i]][ tmp.MGSs, "pctSCCeffect.bgadj"]

  output [i, 6:ncol (output)] = c (rbind (tmp.MGSs.names, tmp.KOsInMGS, tmp.SCC_omiting_MGS,tmp.SCC_delta, tmp.SCC_effect, tmp.SCC_effect_bgadj ))
  
  rm (tmp.MGSs.names, tmp.KOsInMGS,tmp.SCC_module,
      tmp.SCC_out_of_module,tmp.SCC_bgadj,tmp.SCC_omiting_MGS,
      tmp.SCC_delta,tmp.SCC_effect,tmp.SCC_effect_bgadj)
  
}



### write 'output' to tab-delimited text file
write.table (output, file = paste(outpath,Phe,"_top_driver_species.txt",sep=""), row.names = T, col.names = NA, quote = F, sep = "\t")         

####step8_plot_leave_one_MGS#####
pdf (file = paste(outpath,Phe,"_density_plot_SCC_",Phe,".pdf",sep=""), width = 6 * 1, height = 3 * 3)
for (m in names (KOsets)) {
  
  ### Distribution of correlations for KOs in modules vs all other KOs
  incat =     na.omit (KO_cor_Phe [names (KO_cor_Phe) %in% KOsets [[m]] ]) ### select all correlations between phenotype and KOs in the KEGG module
  incat2 =    data.frame ("KO_cor_Phe" = incat, "cat" = "KOs in module")
  notincat =  na.omit (KO_cor_Phe [! (names (KO_cor_Phe) %in% KOsets [[m]])]) ### select all correlations between phenotype and KOs NOT in the KEGG module   
  notincat2 = data.frame ("KO_cor_Phe" = notincat, "cat" = "KOs not in module")
  tmp.df = rbind (incat2, notincat2)
  cdf <- ddply (tmp.df, "cat", summarise, rating.median = median (KO_cor_Phe)) ### find median for each group
  g1 <- ggplot (tmp.df, aes (x = KO_cor_Phe, fill = cat)) + 
    geom_density (alpha = 0.4) +
    geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
    ggtitle (paste0 ("KEGG module: ", m, 
                     "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                     "\n", length (KOsets [[m]]), " KOs in module vs all remaining ", (length (KO_cor_Phe) - length (KOsets[[m]])), " KOs", "\n")) +
    xlab ("SCC for KOs and phenotype") +
    theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
  
  ### Plot Distribution of correlations when leave-one-MGS-out.   
  if (nrow (DeltaSCCperMGS [[m]]) > 1) {
    g2 <- ggplot (DeltaSCCperMGS [[m]], aes (x = SCC_omiting_MGS)) + 
      geom_density (alpha = 1, fill = "grey", aes (y = ..scaled..)) + 
      geom_segment (aes (y = -0.1, yend = -0.02, x = SCC_omiting_MGS, xend = SCC_omiting_MGS)) +
      geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                       "\nSCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (DeltaSCCperMGS [[m]]))) +
      xlab ("SCC for KOs and phenotype") +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
    
  } else {
    
    ### Special circumstance when there is only one MGS (then one cannot make a 
    ### density plot, instead a black line is plotted at the value when that 
    ### one MGS is left out).
    
    g2 <- ggplot () +
      scale_x_continuous (limits = range (c (DeltaSCCperMGS [[m]][, "SCC_omiting_MGS"], cdf$rating.median))) +
      scale_y_continuous (name = "", limits = c (0, 1)) +
      geom_vline (data = cdf, aes (xintercept = rating.median, colour = cat), linetype = "dashed", size = 0.8) +
      geom_vline (data = DeltaSCCperMGS [[m]], aes (xintercept = SCC_omiting_MGS), linetype = "longdash", color = "black",size = 1) + 
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,") [[1]][1], 
                       "\nSCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (DeltaSCCperMGS [[m]]))) +
      xlab ("SCC for KOs and phenotype") +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ()) 
    
  }
  
  ### Plot distribution of correlations when leave-one-MGS-out. Bg.adjust SCC 
  ### (i.e. what we show in Figure 3c+d)
  tmp.DeltaSCCperMGS = DeltaSCCperMGS [[m]]
  ### Calculate delta-SCC in relation to the background adjusted SCC. 
  ### i.e the 'plottedvalue.s.m' shown as the second last equation in the 
  ### methods section in Pedersen et al., 2016.
  tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj = tmp.DeltaSCCperMGS$SCC.bgadj - tmp.DeltaSCCperMGS$DeltaMGS_SCC
  ### specify the range of the x-axis to include 0, i.e. [min:0] for negative 
  ### correlations and [0:max] for positive correlations
  x.range = c (min (0, tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj), max (0, tmp.DeltaSCCperMGS$DeltaMGS_SCC.bgadj))
  
  if (nrow (tmp.DeltaSCCperMGS) > 1 ) {
    
    g3 <- ggplot (tmp.DeltaSCCperMGS, aes (x = DeltaMGS_SCC.bgadj)) + 
      geom_density (alpha = 1, fill = "grey", aes (y = ..scaled..)) + 
      geom_segment (aes (y = -0.1, yend = -0.02, x = DeltaMGS_SCC.bgadj, xend = DeltaMGS_SCC.bgadj)) +
      ggtitle (paste0 ("KEGG module: ", m, 
                       "\n", strsplit (module_mapping [m], split = "\\[|,")[[1]][1], 
                       "\nbg.adj.SCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (tmp.DeltaSCCperMGS))) +
      xlab ("bg.adj.SCC for KOs and phenotype") + xlim (x.range) +
      theme_classic () + theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
    
  } else { 
    
    ### Special circumstance when there is only one MGS (then one cannot make a 
    ### density plot)
    
    g3 <- ggplot () +
      ggtitle (paste0 (m, 
                       "\n", strsplit (module_mapping [m], split="\\[|,")[[1]][1], 
                       "\nbg.adj.SCC for leave-1-MGS-out",
                       "\nNumber of MGSs = ", nrow (tmp.DeltaSCCperMGS), " - consequently not showing a density/rug plot")) 
    
  }
  
  g = plot_grid (g1, g2, g3, ncol = 1, nrow = 3, align = "v", labels = c ("a", "b", "c"), axis = "rl")
  # g = plot_grid (g1, g3, ncol = 1, nrow = 2, align = "v", labels = c ("a", "b", "c"), axis = "rl") ### use this if deleting g2
  plot (g)
  
}

dev.off ()  

####step9_plot_KEGG_module_phenotypes####
modules = KOsets
r = round(length(modules)/4) 
if (r == 0){
  r = r+1
}
c = length(modules)/r
pdf (file = paste(outpath,"plot_KEGG_module_",Phe,".pdf",sep=""), width = 8.5/4*c, height = 6.5/4*r)
modules = vartotest_sig[[final.fdr.cutoffs1]]$keggmodules

module_plot <- function(module) {
  incat = na.omit (KO_cor_Phe [ names (KO_cor_Phe) %in% koann [[module]] ]) ### select all correlations between Phe and KOs in the KEGG module
  notincat = na.omit (KO_cor_Phe [! (names (KO_cor_Phe) %in% koann [[module]])]) ### select all correlations between Phe and KOs NOT in the KEGG module
  if (length (incat) > 0 & length (notincat) > 0) {
    incat = cbind(as.numeric(incat), module)
    notincat = cbind(as.numeric(notincat), "Others")
    df = as.data.frame(rbind(incat,notincat))
    colnames(df) = c("cor","group")
    df$cor = as.numeric(df$cor)
    ggplot(df, aes(group, as.numeric(cor)))+
      geom_boxplot()+mytheme+labs(x = NULL, y = "Spearman correlation")+
      stat_compare_means(comparisons = list(c(module, "Others")),methods="wilcox.test")
  }  
}

pp <- lapply(modules, module_plot) %>% 
  plot_grid(plotlist = ., align = "h", 
            nrow = r)
plot(pp)
dev.off()
rm(r,c)


### Figure 3A
module_phe_sig_results_ann_sort <- read.table("./intermediate_results/Figure3_results/module_phe_sig_results_ann_sort_fdr_0.01.tsv",sep="\t",row.names = NULL,header=T)

figure3A <- ggplot(module_phe_sig_results_ann_sort,aes(Phenotypes,KEGG_modules))+
  geom_point(aes(size=as.numeric(as.factor(Q)),color=as.factor(Estimate2)))+
  scale_colour_manual(values=(c("#4EA74A","#8E4B99")))+
  scale_size_area(max_size = 6) + theme(legend.position = "none")

ggsave("./figures/figure3A.pdf",figure3A,width=8,height=8)
