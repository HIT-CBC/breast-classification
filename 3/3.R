#模块
library(WGCNA)
library(stringr)
library(reshape2)

options(stringsAsFactors = FALSE)

luma_gene = read.table(file = "luma_gene.txt",header = T,sep = "\t")
datExpr <- t(luma_gene)
datExpr = as.data.frame(datExpr)
##软阈值筛选##
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##一步法网络构建：One-step network construction and module detection##
net = blockwiseModules(datExpr, power = 12, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)

##绘画结果展示##
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

##结果保存
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "AS-green-FPKM-02-networkConstruction-auto.RData")
write(colnames(datExpr),file = "gene_list.txt",sep = "\n")


#kobas,kegg通路分析
luma_kegg = read.table(file = "luma_kegg.txt",sep = "\t",header = F,quote  = "")
luma_kegg = luma_kegg[luma_kegg$V7 <= 0.05,]

luma_kegg_gene = unique(unlist(strsplit(luma_kegg[,"V8"],'\\|')))

gene_list = cbind(moduleLabels,moduleColors)

blue = rownames(gene_list)[which(gene_list[,"moduleColors"] == "blue")]
brown = rownames(gene_list)[which(gene_list[,"moduleColors"] == "brown")]
grey = rownames(gene_list)[which(gene_list[,"moduleColors"] == "grey")]
turquoise = rownames(gene_list)[which(gene_list[,"moduleColors"] == "turquoise")]
yellow = rownames(gene_list)[which(gene_list[,"moduleColors"] == "yellow")]

length(intersect(blue,luma_kegg_gene))/length(blue)
length(intersect(brown,luma_kegg_gene))/length(brown)
length(intersect(grey,luma_kegg_gene))/length(grey)
length(intersect(turquoise,luma_kegg_gene))/length(turquoise)
length(intersect(yellow,luma_kegg_gene))/length(yellow)

write(blue,file = "luma_blue_gene.txt",sep = "\n")

blue = read.table(file = "luma_blue_gene.txt",sep = "\n")

boruta_gene = read.table(file = "boruta_gene.txt",sep = "\n")
mygene = intersect(boruta_gene[,1],blue[,1])
write.table(mygene,quote = F,file = "keygene.txt",col.names = F,row.names = F)
