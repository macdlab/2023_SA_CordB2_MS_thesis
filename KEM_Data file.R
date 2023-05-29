#set the directory
setwd("C:/Users/hp/OneDrive/Documents/Dissertation project/KEM_Data")
#load packages
library(WGCNA)
library(tibble)
library(ggplot2)
library(tidyverse)

#do not omit this step
options(stringsAsFactors = FALSE)

#Read in the dataset
KEMData = read.csv("neonatal_expr_data_yy.csv")
#take a look in the dataset
dim(KEMData)
names(KEMData)

#make the dataframe
datExpr0 = data.frame(KEMData)
#make column to row transformation
datExpr0<-column_to_rownames(datExpr0, var = "gene_symbol")
#transpose the data
datExpr = as.data.frame(t(datExpr0))
#remove the X
rownames(datExpr) = sub("*\\X", "", rownames(datExpr))

#check the genes and samples for too many missing values
gsg = goodSamplesGenes(datExpr, verbose =3)
gsg$allOK

#cluster the samples to look for any outliers
sampleTree = hclust(dist(datExpr), method = "average")
#plot the sample tree
#open graphic output window of size 12 by 9
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
#plot a line to show the cut
abline(h = 2.0e+07, col = "red")
#determine the cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 2.0e+07, minSize = 10)
table(clust)
#we want the clust 1
keepSamples = (clust==1)
datExpr1 = datExpr[keepSamples,]
nGenes = ncol(datExpr1)
nSamples = nrow(datExpr1)
#datExpr1 now has the expression data ready for network analysis

#loading clinical data
#read in the trait data- we have two trait datasets here
#first we will complete with multimicronutrient dataset- "KEM_MMN"
#then we will work with the phenotypic dataset- "KEM_pheno"

traitData = read.csv("KEM_MMN.csv")

#traitData = read.csv("KEM_pheno.csv")

dim(traitData)
names(traitData)
#form a data frame analogous to expression data that will hold the clinical traits
KEMSamples = rownames(datExpr1)
#match the samples for which they were measured to the expression samples
traitRows = match(KEMSamples, traitData$child_no)
view(traitRows)
datTraits0 = traitData[traitRows,]
rownames(datTraits0) = traitData[traitRows, 1]
#remove the column-child_no
datTraits0<- datTraits0[,-1]

datTraits = datTraits0[,-c(8:20)] #this is for MMN dataset
datTraits = datTraits[,-c(3)]
#datTraits01<- datTraits0[, -c(1,2)]
#datTraits = datTraits01 
#the above two commands are for phenotype dataset+

collectGarbage()

#expression data is in the variable- datExpr1
#clinical traits is in the variable- datTraits
#visualize the clinical traits relating to the sample dendrogram
#recluster the samples
sampleTree2 = hclust(dist(datExpr1), method = "average")

#convert traits to color rep, white- low, red- high, grey- missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)

#plot the sample dendrogram and the colors underneath
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

#save the data
save(datExpr1, datTraits, file = "KEM-01-dataInput.RData") #this is for MMN dataset
#save(datExpr1, datTraits, file = "KEM-03-dataInput.RData") #this is for phenotype

#network construction-> 1 step

enableWGCNAThreads()
#load the data saved above
lnames = load(file = "KEM-01-dataInput.RData") #MMN
#lnames = load(file = "KEM-03-dataInput.RData") #KEM_pheno

#the variable lnames contains the names of the loaded variables
lnames

#choose the soft-thresholding powers
powers = c(c(1:20), seq(from = 12, to=20, by=2))
#call the network topology analysis function
sft = pickSoftThreshold(datExpr1, networkType = "signed", powerVector = powers, verbose = 5)

#plot the results
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.85
#scale-free topology fit index as a function of the sft power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
#the following line will correspond to using an R^2 cut-off of h
abline(h=0.85, col="red")
#mean connectivity as a function of the sft power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0, col="red")
#one-step network construction

net = blockwiseModules(datExpr1, power = 7, TOMType = "signed",
                       networkType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "KEMTraitTOM",
                       verbose = 3)

table(net$colors)

#open the graphics window
sizeGrWindow(12,9)

#convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
#plot the dendro and module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "KEM-02-networkConstruction-1step.RData")

#load the expression and trait data file
lnames = load(file = "KEM-01-dataInput.RData") #MMN
#lnames = load(file = "KEM-03-dataInput.RData") #KEM_pheno
lnames
#load the network data in the second part
lnames1 = load(file = "KEM-02-networkConstruction-1step.RData")
lnames1

#relate modules to external clinical traits

#quantify the module-trait associations
#we have eigengene for each module
#correlate the eigengenes with external traits
#define numbers of genes and samples
nGenes = ncol(datExpr1)
nSamples = nrow(datExpr1)
#recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#make the graphical representation
#color code the association with correlation value
sizeGrWindow(10,6)
#will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = " ")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 15.5, 2.5, 1))

#display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.35,
               zlim = c(-1,1))

#gene relationship to trait and important modules- GS and MM
#b2, b12, b9
#define the variable- cord blood_b2 containing the b2 levels
B2_levels = as.data.frame(datTraits$cord_b2)
names(B2_levels) = "B2_levels"
#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr1, B2_levels, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(B2_levels), sep="");
names(GSPvalue) = paste("p.GS.", names(B2_levels), sep="")

#mm-gs for modules
#all the mm-gs done by replacing module with the respective colours
module = "lightsteelblue1"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for cord blood B2 levels",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#display the gene names inside the module
#b2
#blue, brown, greenyellow
#lightsteelblue1 not for extraction
blue = colnames(datExpr1)[moduleColors=="blue"]
write.csv(blue, file = "blue-b2_signed.csv")

#lightsteelblue1
lightsteelblue1 = colnames(datExpr1)[moduleColors=="lightsteelblue1"]
write.csv(lightsteelblue1, file = "lightsteelblue1-b2_signed.csv")

#greenyellow
greenyellow = colnames(datExpr1)[moduleColors=="greenyellow"]
write.csv(greenyellow, file = "greenyellow-b2_signed.csv")

#brown
brown = colnames(datExpr1)[moduleColors=="brown"]
write.csv(brown, file = "brown-b2_signed.csv")


#identifying most important genes for one determined characteristic inside of the cluster
geneInfo0 = data.frame(KEMSamples = colnames(datExpr1),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
#write information in csv file
write.csv(geneInfo0, file = "geneInfo0.csv")

modOrder = order(-abs(cor(MEs, B2_levels, use = "p")))
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.B2_levels))
geneInfo = geneInfo0[geneOrder,]
#write the information in csv file
write.csv(geneInfo, file = "geneInfo.csv")

#visualization of the eigengenes network
#recalculate module eigengenes
MEs = moduleEigengenes(datExpr1, moduleColors)$eigengenes
#isolate cord_B2 from the clinical traits
B2levels = as.data.frame(datTraits$cord_b2)
names(B2levels) = "B2levels"
#add B2 levels to existing module eigengenes
MET = orderMEs(cbind(MEs, B2levels))
#plot the relationships among the eigengenes and the trait
#both dendrogram and heatmap will be seen together
sizeGrWindow(6,8)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(2,2,1,4), marHeatmap = c(4.3,5,0,4), cex.lab = 0.8, xLabelsAngle = 90)

#plot the dendrogram
sizeGrWindow(6,6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
#plot the heatmap matrix
#this following plot will overwrite the dendrogram plot
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5.5,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

