##########################################
#tcga数据
cancer = read.table(file = "cancer_ori.txt",header = T,sep = "\t",row.names = 1)
subtype = read.csv(file = "subtype.csv",header = T)

A_ID = gsub("-",".",subtype[which(subtype[,"Subtype"] == "BRCA_LumA"),"Patient.ID"])
test = cancer[,which(gsub(".01","",colnames(cancer)) %in% A_ID)]
library(DMwR2)
test <- knnImputation(test, k=10, scale = T, meth = "weighAvg")
write.table(test,file = "BRCA_LumA.txt",sep = "\t",quote = F)
#1000次差异位点
luma = read.table("BRCA_LumA.txt",header = T,row.names = 1,sep = "\t")
normal = read.table(file = "Normal.txt",header = T,sep = "\t")

cg = rownames(luma)

for(i in 1:1000){
  tumor = sample(luma,100)
  test = cbind(tumor,normal)
  Pvalue = apply(test,1,function(x){
    t.test(x[1:ncol(tumor)],x[(ncol(tumor)+1):(ncol(tumor)+ncol(normal))])$p.value
  })
  gap = apply(tumor,1,mean)-apply(normal,1,mean)
  fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))
  luma_cg = rownames(tumor)[which(fdr <0.05 & abs(gap) >=0.2)]
  cg = intersect(cg,luma_cg)
  
}
#位点对应基因
library(dplyr)
cgtest = read.table(file ="cgtest.txt",header = F,sep = "\t")
luma = read.table("BRCA_LumA.txt",header = T,row.names = 1,sep = "\t")

luma_cg = luma[rownames(luma)%in%cg[,1],]
luma_cg$cg = rownames(luma_cg)

luma_gene = merge(luma_cg,cgtest,by.x = "cg",by.y = "V1",all.x = F,all.y = F)
luma_gene = luma_gene[,-c(which(colnames(luma_gene) %in% c("V2","V4","V5")))]

luma_gene =  as.data.frame(luma_gene %>% group_by(V3) %>% summarise_each(funs(mean)))

luma_gene = luma_gene[,-2]

rownames(luma_gene) = luma_gene[,1]
luma_gene = luma_gene[,-1]

write.table(luma_gene,file = "luma_gene.txt",sep = "\t")


#临床数据
mmc = read.csv(file = "mmc.csv",header = T,encoding = "UTF-8")

mmc = mmc[which(mmc[,"type"] == "BRCA"),c("bcr_patient_barcode","PFS","PFS.time")]
mmc[,"bcr_patient_barcode"] = gsub("-",".",mmc[,"bcr_patient_barcode"])
write.table(mmc,file = "clinic_BRCA.txt",sep = "\t")

#geo数据

cgtest = read.table(file ="cgtest.txt",header = F,sep = "\t")
bbcg = read.table(file = "bbcg.txt",header = T,sep = "\t")
pam50 = read.table(file = "b.txt",header = T,row.names = 1)
luma = colnames(pam50)[which(pam50[1,] == "LumA")]
luma_geo = bbcg[,colnames(bbcg) %in% luma]
luma_geo = apply(luma_geo,1,function(x) {
  gsub("null",NA,x)
})

luma_geo <- knnImputation(luma_geo, k=10, scale = T, meth = "weighAvg")

luma_geo$cg = rownames(luma_geo) 
luma_geo_gene = merge(luma_geo,cgtest,by.x = "cg",by.y = "V1",all.x = T,all.y = F)
luma_geo_gene = luma_geo_gene[,-which(colnames(luma_geo_gene) %in% c("cg","V2","V4","V5"))]
luma_geo_gene[,1:(ncol(luma_geo_gene)-1)] = apply(luma_geo_gene[,1:(ncol(luma_geo_gene)-1)],2,as.numeric)
luma_geo_gene =  as.data.frame(luma_geo_gene %>% group_by(V3) %>% summarise_each(funs(mean)))
rownames(luma_geo_gene) = luma_geo_gene$V3
luma_geo_gene = luma_geo_gene[,-1]
write.table(luma_geo_gene,file = "luma_geo_gene.txt",sep = "\t")


#热图绘制，AB
load("cancer_5types.RData")
load("normal.RData")
load("BRCA_luma.RData")
load("LumB.RData")

cgtest = read.table(file ="cgtest.txt",header = F,sep = "\t")
cg = rownames(cancer)
for(i in 1:1000){
  tumor = sample(cancer,100)
  test = cbind(tumor,normal)
  Pvalue = apply(test,1,function(x){
    t.test(x[1:ncol(tumor)],x[(ncol(tumor)+1):(ncol(tumor)+ncol(normal))])$p.value
  })
  gap = apply(tumor,1,mean)-apply(normal,1,mean)
  fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))
  luma_cg = rownames(tumor)[which(fdr <0.05 & abs(gap) >=0.2)]
  cg = intersect(cg,luma_cg)
  
}
write.table(cg , file = "cancer_1000cg.txt",sep = "\t")
library(pheatmap)
library(RColorBrewer)


sam_lum = union(colnames(lumb),colnames(luma))
sam_oth = setdiff(colnames(cancer),sam_lum)
sam_normal = colnames(normal)
cg_lum = read.table(file = "lum_1000cg.txt",sep = "\t")
cg_cancer = read.table(file = "cancer_1000cg.txt",sep = "\t")
cg_oth = read.table(file = "lum_oth_1000cg.txt",sep = "\t")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

diff = cbind(cancer[which(rownames(cancer)%in% cg_cancer[,1]),sample(1:692,100)],normal[which(rownames(normal)%in% cg_cancer[,1]),])#癌症跟正常
diff = cbind(cancer[which(rownames(cancer)%in% cg_cancer[,1]),sample(sam_lum,100)],normal[which(rownames(normal)%in% cg_cancer[,1]),])#ab和正常
diff = cancer[which(rownames(cancer)%in% cg_cancer[,1]),sam_lum[1:300]]#a和b
diff = cbind(luma[which(rownames(luma) %in% cg_cancer[,1]),sample(1:349,140)],lumb[which(rownames(lumb)%in%cg_cancer[,1]),sample(1:150,60)],cancer[which(rownames(cancer)%in%cg_cancer[,1]),sam_oth])

anno_col=data.frame(sampleType=factor(rep(c("cancer","normal"),c(100,96))))
rownames(anno_col)=colnames(diff)

ann_color=list(sampleType=c(cancer='#cd0000',normal='#3a5fcd'))
xx = pheatmap(diff[1:1000,],scale="none",annotation_col=anno_col,annotation_colors=ann_color,
              clustering_method = 'ward.D2',clustering_distance_rows = 'canberra',
              show_colnames=F,show_rownames=F,fontsize_row=6,fontsize=10,
              legend_breaks=c(0.2,0.8),legend_labels=c('low','high'),
              color = colorRampPalette(c('#436eee','white','#EE0000'))(100)
)
#临床热图

library(pheatmap)

clin = read.csv("clin.csv")
clin = clin[,c("bcr_patient_barcode","A3_T","A4_N","A5_M")]
clin$bcr_patient_barcode = gsub("-",".",clin$bcr_patient_barcode)

diff = as.data.frame(t(luma_gene))
diff$bcr_patient_barcode = rownames(diff)
diff = merge(diff,clin,by = "bcr_patient_barcode",all.x = T,all.y = F)
diff$A3_T = as.factor(diff$A3_T)
diff$A4_N = as.factor(diff$A4_N)
diff$A5_M = as.factor(diff$A5_M)



anno_col=test[,c("A3_T","A4_N","A5_M")]
colnames(anno_col) = c("T_stage","N_stage","M_stage")
rownames(anno_col)=rownames(diff)
diff = as.data.frame(t(diff))
ann_colors = list(T_stage=c("white","firebrick","green","grey","black"),N_stage = c("#00447E","#F34800","red","blue","pink"),M_stage = c("white","#64A10E","#795EA2","#3370CC"))# 

pheatmap(diff[1:1000,],scale = "none",show_colnames=F,show_rownames=F,cluster_rows = F,cluster_cols = T,
         legend_breaks=c(0.2,0.4,0.6,0.8),color = colorRampPalette(c('#436eee','white','#EE0000'))(100),
         annotation_col = anno_col
)
