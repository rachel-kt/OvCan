#https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html
library(GOplot)
data(EC)
head(EC$david)
head(EC$genelist)
circ_2 <- circle_dat(EC$david, EC$genelist)

#david = read.csv("./goplot/david_black_2.csv", stringsAsFactors = T, sep = "\t")
david = read.csv("./../Final figures/new network/goplot/mod1david.csv", stringsAsFactors = T, sep = "\t")

genelist = read.csv("./../Final figures/new network/goplot/mod1genelist.csv", stringsAsFactors = T)
genelist$logFC = genelist$logFC*-1
circ <- circle_dat(david, genelist)

GOBar(subset(circ, category == 'BP'), title = "BP")
# Facet the barplot according to the categories of the terms 
GOBar(circ, display = 'multiple')
# Facet the barplot, add a title and change the colour scale for the z-score
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))
# Generate the bubble plot with a label threshold of 3

#GOBubble(data = circ)
#GOBubble(subset(circ, category == 'BP'))
#GOBubble(circ, ID = F, labels = 7, colour = c('orange', 'darkred', 'gold', 'blue', 'gray80','pink', 'green'))
# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
GOBubble(circ, title = 'Bubble plot', ID = T, colour = c('blue', 'green', 'red', 'yellow'), display = 'multiple', labels = 0.3)
# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)
# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)
# Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)
# Generate a circular visualization of selected terms
#IDs <- c('GO:0048407', 'GO:0061631', 'GO:0005201', 'GO:0046872', 'hsa05150', 'hsa04610', 'hsa05133', 'hsa04512', 'hsa04974', 'GO:0016328','GO:0005923','GO:0005581')
GOCircle(circ, nsub = IDs, lfc.col = "red")
# Generate a circular visualization for 10 terms
IDs <- c("GO:0019886","GO:2001199",
         'GO:0045087',
         'GO:0002504',
         'GO:0032733',
         'GO:0006911',
         'GO:0032760',
         'GO:0006954',
         'GO:0042102',
         'GO:0006955',
         'GO:0060333',
         'GO:0031295',
         'GO:0006958',
         'GO:0030666',
         'GO:0016021',
         'GO:0042613',
         'GO:0005886',
         'GO:0005765',
         'GO:0071556',
         'GO:0030669',
         'GO:0012507',
         'GO:0032588',
         'GO:0005581',
         'GO:0032395',
         'hsa05150',
         'hsa04610',
         'hsa05322',
         'hsa05152',
         'hsa05310',
         'hsa05020',
         'hsa05133')
GOCircle(circ, nsub = IDs)
