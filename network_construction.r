
net_dat_degs = read.csv("./new network/network_data_29-12-2019_new.csv", row.names = 1, stringsAsFactors = F)
net_dat_degs = t(scale(t(net_dat_degs), center = T, scale = T))

datExpr = t(net_dat_degs)
library(WGCNA)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.88,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

POW = 9
#write.csv(net_dat_degs, file = "data_network_pids_10_12_2019.csv")
# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr,use="p"))^POW
# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
#k=softConnectivity(datE=datExpr,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

# Turn adjacency into a measure of dissimilarity
dissADJ=1-ADJ1
dissTOM=TOMdist(ADJ1)
collectGarbage()

hierADJ=hclust(as.dist(dissADJ), method="average" )
# Plot the resulting clustering tree together with the true color assignment
sizeGrWindow(10,5);
#plotDendroAndColors(hierADJ, dendroLabels = FALSE, hang = 0.03,
#                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )
plot(hierADJ, xlab="", sub="", main = "Gene clustering on hierarchical clustering",
     labels = FALSE, hang = 0.001)

branch.number=cutreeDynamic(hierADJ,method="tree")
# This function transforms the branch numbers into colors
colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))

colorDynamicADJ=labels2colors(branch.number )

colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ,
                                                  cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))
# Plot results of all module detection methods together:
sizeGrWindow(10,5)
plotDendroAndColors(dendro = hierADJ,
                    colors=data.frame(colorStaticADJ,
                                      colorDynamicADJ, colorDynamicHybridADJ),
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Gene dendrogram and module colors")

### TOM ####
# Calculate the dendrogram

hierTOM = hclust(as.dist(dissTOM),method="average");
# The reader should vary the height cut-off parameter h1
# (related to the y-axis of dendrogram) in the following
#colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))
colorDynamicTOM = labels2colors(cutreeDynamic(hierTOM,method="tree", minClusterSize = 20, deepSplit = 4))
#colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
#                                                   deepSplit=4, pamRespectsDendro = FALSE))
# Now we plot the results
# sizeGrWindow(10,5)
# plotDendroAndColors(hierTOM,
#                     colors=data.frame(colorStaticTOM,
#                                       colorDynamicTOM, colorDynamicHybridTOM),
#                     dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
#                     main = "Gene dendrogram and module colors, TOM dissimilarity")

sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")


dynamicMods = cutreeDynamic(hierTOM,method="tree", minClusterSize = 20, deepSplit = 4)
table(dynamicMods)
mods<- table(dynamicMods)
write.csv(mods, "./new network/modules.csv")
dynamicColors = labels2colors(dynamicMods)
moduleColors = dynamicColors
table(dynamicColors)
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
ME_unmerged = MEList
ME_unmerged = MEList$eigengenes

datTraits = trait_table
moduleTraitCor = cor(ME_unmerged, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3, 3, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(ME_unmerged),
               ySymbols = names(ME_unmerged),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
"mod_trait_unmerged_big_font"
MEs = MEList$eigengenes
##############################################
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
MEDissThres = 0.25
par(mfrow=c(1,1))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
#moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(mergedColors, colorOrder)-1;
#MEs = mergedMEs;


sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,
                    colors=(cbind(moduleColors,mergedColors)),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and merged module colors, TOM dissimilarity")


nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
moduleColors = mergedColors
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

trait_table = data.frame(rownames(datExpr))
trait_table$normal = c(rep(1,12), rep(0,57))
trait_table$cancer = c(rep(0,12), rep(1,57))
rownames(trait_table) = trait_table$rownames.datExpr.
trait_table = trait_table[,-1]
datTraits = trait_table
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 9, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
cancerstatus = as.data.frame(datTraits$cancer)
names(cancerstatus) = "Cancerstatus"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, cancerstatus, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(cancerstatus), sep="");
names(GSPvalue) = paste("p.GS.", names(cancerstatus), sep="");
genes =data.frame(rownames(net_dat_degs), stringsAsFactors = F)
module = "black"
column = match(module, modNames);
moduleGenes = moduleColors==module
black_genes = genes[moduleGenes,]
column = match("blue", modNames);
moduleGenes = moduleColors=="blue"
blue_genes = genes[moduleGenes,]
column = match("lightcyan", modNames);
moduleGenes = moduleColors=="lightcyan"
lightcyan_genes = genes[moduleGenes,]
column = match("cyan", modNames);
moduleGenes = moduleColors=="cyan"
cyan_genes = genes[moduleGenes,]
column = match("grey", modNames);
moduleGenes = moduleColors=="grey"
grey_genes = genes[moduleGenes,]


write.csv(c(black_genes, blue_genes, lightcyan_genes, cyan_genes, grey_genes), "./new network/network_genes.csv")
write.csv(GSPvalue, "./Final figures/gene_Sig.csv")
write.csv(geneModuleMembership,"./Final figures/modmem.csv")
###################################3

# module membership vs gene significance p value
library(ggplot2)
scatter_data = data.frame(cbind(abs(geneModuleMembership$MMcyan),abs(geneTraitSignificance$GS.Cancerstatus)))
#scatter_data = data.frame(cbind((geneModuleMembership$MMcyan),(geneTraitSignificance$GS.Cancerstatus)))
cyan_g = read.csv("cyan genes", header = F)
x = geneModuleMembership[match(cyan_g$V1, rownames(geneModuleMembership)),]
y = geneTraitSignificance[match(cyan_g$V1, rownames(geneTraitSignificance)),]
scatter_genes = data.frame(cbind(abs(x$MMcyan),abs(y)))


ggplot(scatter_data, aes(x=X1, y=X2)) +
  geom_point(size=2, shape = 16, alpha = 1/3) +
  xlab("Module membership for Cyan module") +
  ylab("Gene significance for cancer") +
  geom_point(mapping = aes(X1, X2) ,data = scatter_genes, colour = "orange", shape = 16)

#####################################################################################
library(flashClust)
adj_matrix <- adjacency.fromSimilarity(as.matrix((cor_test)), power=9, type='unsigned')
#adj_matrix <- adjacency.fromSimilarity(as.matrix((cor_test)), power=10, type='unsigned')
gene_tree <- flashClust(as.dist(1 - adj_matrix), method="average")
plot.new()
par(mfrow=c(1,1))
plot(gene_tree, xlab="", sub="", main = "Gene clustering on hierarchical clustering",
     labels = FALSE, hang = 0.05)
w = as.dist(1 - adj_matrix)
#module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15, deepSplit=F)
module_labels <- cutreeDynamic(dendro = gene_tree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                               minClusterSize = 100, method = "tree")
#assign module colours
module_colors = labels2colors(module_labels)

#plot the dendrogram and corresponding colour bars underneath
#par(new = F)
#pdf(file = "old.pdf")
plotDendroAndColors(gene_tree, module_colors,'Module colours', dendroLabels = FALSE, hang = 0.009,
                    addGuide = TRUE, guideHang = 0.05, main='Gene clustering on hierarchical clustering')
dev.off()
#module_colors <- labels2colors(module_labels)

library(circlize)
library(ComplexHeatmap)
par(new = F)
Cor_mat_orig = cor(datExpr)
pdf(file = "./new network/COR_MATRIX_dec_2019_power9.pdf")
f1 = colorRamp2(seq(min(Cor_mat_orig), max(Cor_mat_orig),length = 3),
                c("blue", "white", "red"), space = "LAB")
Heatmap(as.data.frame(Cor_mat_orig), col = f1, cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = F, name = "correlation")
dev.off()

pdf(file = "./new network/ADJ_MATRIX_dec_2019_power9.pdf")
f1 = colorRamp2(seq(min(ADJ1), max(ADJ1),length = 3),
                c("blue", "white", "red"), space = "LAB")
Heatmap(as.data.frame(ADJ1), col = f1, cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = F, name = "correlation")
dev.off()
#################################################################

################### Exporting the network  ########################

export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.8, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
  library('igraph')
  
  # Determine filename to use
  if (is.null(filename)) {
    filename='network.graphml'
  }
  # TODO 2015/04/09
  # Add option to rescale correlations for each module before applying
  # threshold (this is simpler than the previous approach of trying to
  # determine a different threshold for each module)
  #
  # Still, modules with very low correlations should be given somewhat
  # less priority than those with very high correlations.
  
  #module_colors <- unique(nodeAttrDataFrame$color)
  #module_genes <- which(nodeAttrDataFrame$color == color)
  #module_adjmat <- adj_mat[module_genes,]
  #num_genes <- length(module_genes)
  
  # Adjust threshold if needed to limit remaining edges
  max_edges <- max_edge_ratio * nrow(adj_mat)
  
  edge_to_total_ratio <- max_edges / length(adj_mat)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
  
  # Also choose a minimum threshold to make sure that at least some edges
  # are left
  min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
  
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  # Remove edges with weights lower than the cutoff
  adj_mat[abs(adj_mat) < 0.7] <- 0
  
  # Drop any genes with no edges (TODO: Make optional)
  orphaned <- (colSums(adj_mat) == 0)
  adj_mat <- adj_mat[!orphaned, !orphaned]
  
  # Also remove annotation entries
  if (!is.null(nodeAttr)) {
    nodeAttr <- nodeAttr[!orphaned]
  }
  if (!is.null(nodeAttrDataFrame)) {
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
  }
  
  # Keep track of non-positive edges and rescale to range 0,1
  is_zero     <- adj_mat == 0
  is_negative <- adj_mat < 0
  
  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
  adj_mat[is_zero] <- 0
  adj_mat[is_negative] <- -adj_mat[is_negative]
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)))
  }
  
  # Create a new graph and add vertices
  # Weighted graph
  if (weighted) {
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  }
  
  # Add single node annotation from vector
  if (!is.null(nodeAttr)) {
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  }
  
  # Add node one or more node annotations from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
    }
  }
  
  edge_correlation_negative <- c()
  
  # neg_correlations[edge_list]
  edge_list <- get.edgelist(g)
  
  for (i in 1:nrow(edge_list)) {
    from <- edge_list[i, 1]    
    to   <- edge_list[i, 2]    
  }
  
  # Save graph to a file
  write.graph(g, filename, format='graphml')
  
  # return igraph
  return(g)
}

gene_ids <- rownames(ADJ1)
gene_info <- data.frame(cbind(gene_ids, module=moduleColors))
library(gplots)
# Include RGB versions of module colors for better assignment in Cytoscape
gene_info$color_rgb <- col2hex(gene_info$module)

# first, it's a good idea to check the distribution of edges weights in our
# correlation matrix. This will help us choose a reasonable cutoff for
# exporting the network.
g <- export_network_to_graphml(ADJ1, filename='./network-jan-2020_unsigned.graphml',
                               threshold=0.7, nodeAttrDataFrame=gene_info)

edge_list_cancer = data.frame(get.edgelist(g))
write.csv(edge_list_cancer, "./edge_list_Jan_2020.csv")
edge_list_cancer$Correlation = "NULL"
edge_list_cancer$p_n = "NULL"
i=1
Cor_mat_orig = cor(datExpr)
for (i in 1:512) {
  edge_list_cancer[i,3] = Cor_mat_orig[match(edge_list_cancer[i,1],rownames(Cor_mat_orig)), match(edge_list_cancer[i,2],colnames(Cor_mat_orig))]  
  if (edge_list_cancer[i,3] < 0) {
    edge_list_cancer[i,4] = -1  
  }   
  else{edge_list_cancer[i,4] = 1}
}
write.csv(edge_list_cancer, "./network_edges_nov.csv")
write.csv(row.names(net_dat_degs), "./network_vertices_nov.csv")

write.csv(cbind(rownames(net_dat_degs), colorDynamicTOM), "unmerged modules.csv")
