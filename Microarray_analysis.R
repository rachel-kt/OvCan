setwd("/media/rachel/Windows/Users/All Users/Documents/project_jnu/2018_May/Ovarian Cancer/")
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("hugene10stv1cdf")
library(makecdfenv)
library(affy)
library(simpleaffy)
library(hugene10stv1cdf)
library(limma)
library(hursta2a520709cdf)
library(ComplexHeatmap)
library(circlize)
library(preprocessCore)
#Create CDF package in temporary directory 
pkgpath <- tempdir()

make.cdf.package("GPL10379_HuRSTA-2a520709_custom_MMPM.cdf", 
                 cdf.path="/media/rachel/Windows/Users/All Users/Documents/project_jnu/2018_May/Cancer/Rosetta/", 
                 compress=FALSE, 
                 species = "Homo_sapiens",
                 packagename = "hursta2a520709cdf",
                 package.path = pkgpath)
dir(pkgpath)

# USE TERMINAL to run :- sudo R CMD INSTALL pkgpath/hursta2a520709cdf
# https://www.biostars.org/p/67400/
# https://support.bioconductor.org/p/61195/

ovar.data = read.affy(covdesc = "covdesc.txt")
data.rma = rma(ovar.data)
write.exprs(data.rma, "expression_values_18-July-2018.csv")
str(ovar.data)
pheno = pData(data.rma)
#pheno$batch = c(rep(1,62), rep(2,18))
edata = exprs(data.rma)
#mod = model.matrix(~as.factor(Treatment), data=pheno)
#mod0 = model.matrix(~1,data=pheno)
#batch = pheno$batch
#batch
#modcombat = model.matrix(~1, data=pheno)
#modcombat # equal to mod0
#combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
data.rma2 = exprs(data.rma)

# http://genomicsclass.github.io/book/pages/using_limma.html
fac2 = ovar.data$Treatment
fac2
fit <- lmFit(data.rma2, design=model.matrix(~ fac2))
colnames(coef(fit))
fit <- eBayes(fit)
#tt <- topTable(fit, coef=2)
res_table = topTable(fit, coef=2, number=Inf, sort.by="none")
#write.table(res_table, "limma_results.txt")
#dim(res_table)
res_new = res_table[order(res_table$adj.P.Val),]
write.csv(res_new, "limma_results_17_July_2018_2.csv")
array_annot = read.csv("annotations.csv")
res_new$annotation = array_annot[match(rownames(res_new), array_annot$ID),10]
res_new$description = array_annot[match(rownames(res_new), array_annot$ID),2]
write.csv(res_new, "result_17_July_with annotation.csv")
res_new$threshold = as.logical((res_new$adj.P.Val < 0.01) & (abs(res_new$logFC)>2))
degs_4fc =  res_new[which(res_new$threshold),]
write.csv(degs_4fc, "four_fold_change_degs.csv")
write.csv(data.rma2[match(rownames(degs_4fc), rownames(data.rma2)),], "expression_values_of_degs.csv")
exp_data = read.csv("expression_values_of_degs.csv", row.names = 1)
exp_data.norm = normalize.quantiles(as.matrix(exp_data))

max(exp_data.norm)
min(exp_data.norm)
exp_data.norm.scaled = scale(exp_data.norm, center = T, scale = T)
max(exp_data.norm.scaled)
min(exp_data.norm.scaled)
colnames(exp_data.norm.scaled) = colnames(exp_data)
rownames(exp_data.norm.scaled) = rownames(exp_data)
ids = array_annot[match(rownames(exp_data.norm.scaled), array_annot$ID),]
#####################
# AFTER CLUST

clustered_genes = read.csv("clust_result/clust_result_genes_clusters.csv", row.names = 1)
clustered_genes$gene_symbol = array_annot[match(rownames(clustered_genes), array_annot$ID), 3]
clustered_genes$unigene = array_annot[match(rownames(clustered_genes), array_annot$ID), 7]
clustered_genes$gene_title = array_annot[match(rownames(clustered_genes), array_annot$ID), 2]
clustered_genes$gene_bank = array_annot[match(rownames(clustered_genes), array_annot$ID), 10]

write.csv(clustered_genes, "clustered_genes_annotations.csv")

#######################
#HEATMAP

heatmap_data = read.csv("clust_result/exp_for_clustering_degs_processed.csv",
                        sep = "\t",
                        row.names = 1)

f1 = colorRamp2(seq(min(heatmap_data), max(heatmap_data),length = 3),
                c("blue", "white", "red"), space = "LAB")
Heatmap(heatmap_data, 
        col = f1, 
        cluster_rows = T, 
        cluster_columns = T, 
        show_row_names = F,
        column_names_gp = gpar(fontsize = 5))
######################
#INTERACTIVE HEATMAP

f2 = colorRamp2(seq(min(exp_data.norm.scaled), max(exp_data.norm.scaled), length = 3), 
                c("green", "black", "red"),
                space = "RGB")
f3 = colorRamp2(seq(min(exp_data.norm.scaled), max(exp_data.norm.scaled), length = 3), 
                c("blue", "white", "red"),
                space = "RGB")
Heatmap(exp_data.norm.scaled,
        col = f3, 
        cluster_rows = T, 
        cluster_columns = T, 
        show_row_names = T,
        column_names_gp = gpar(fontsize = 5),
        show_column_names = T)
library(heatmaply)
dir.create("interactive_hm")
heatmaply(exp_data.norm.scaled, file = "interactive_hm/hm_21_July_2018.html")

#########################

head(heatmap_data)
col_id = colnames(heatmap_data)
# row_annot = read.csv(file = "filtered_annot.csv")
# heatmap_new_data = heatmap_data[match(row_annot$Probe_ID, rownames(heatmap_data)),]
# rownames(heatmap_new_data) = row_annot[match(rownames(heatmap_new_data),row_annot$Probe_ID),2]

############### 17 August 2018 ################
unique_genes = read.csv("unique_ids_for_network.csv", row.names = 1)
head(unique_genes)
heatmap_data_new = as.data.frame(exp_data.norm.scaled[match(rownames(unique_genes), rownames(exp_data.norm.scaled)),])
rownames(heatmap_data_new) = unique_genes$ID_annot
head(heatmap_data_new)
heatmaply(heatmap_data_new, file = "interactive_hm/hm_17_August_2018.html")
write.csv(heatmap_data_new, "data_for_network.csv")

########################################################