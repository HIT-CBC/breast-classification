library(philentropy)
library(survival)
library(survminer)
library(proxy)
library(ggplot2)

luma_gene = read.table(file = "luma_gene.txt",header = T,sep = "\t")
clinic_BRCA = read.table(file = "clinic_BRCA.txt",sep = "\t")
mygene = read.table("keygene.txt",sep = "\n")
mygene = mygene[,1]
hs6 = t(luma_gene[rownames(luma_gene) %in% mygene,])

dis_bray <- dist(hs6, method = 'chebyshev')

clust_average <- hclust(dis_bray, method = 'ward.D2')
#聚类树，默认效果
plot(clust_average,hang = -1,xlab = "sample",sub = "",labels = F)
#进行初步展示
clust_average_cut <- cutree(clust_average, k = 2)
rect.hclust(clust_average, k = 2, border = c('red', 'blue'))

hs6=cbind(hs6,clust_average_cut)
hs6=as.data.frame(hs6)

hs6$bcr_patient_barcode = rownames(hs6)
hs6 = merge(hs6,clinic_BRCA,by = "bcr_patient_barcode",all.x = T,all.y = F)

hs6$PFS=as.numeric(hs6$PFS)
hs6$PFS.time=as.numeric(hs6$PFS.time)
surv1=Surv(hs6$PFS.time,hs6$PFS)

surv_diff <- survdiff(surv1 ~ clust_average_cut,rho = 1,data = hs6)
p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)

fit=survfit(surv1~clust_average_cut,data=hs6) 
ggsurvplot( fit, 
            data = hs6, 
            size = 1, # change line size 
            palette = "aaas",
            #conf.int = TRUE, # Add confidence interval 
            pval = TRUE, # Add p-value 
            risk.table = TRUE, # Add risk table 
            risk.table.col = "clust_average_cut",# Risk table color by groups 
            legend.labs = c("cluster1","cluster2"), # Change legend labels
            risk.table.height = 0.25, # Useful to change when you have multiple groups 
            ggtheme = theme_bw() 
)
#免疫基因

immgene = read.table(file = "immGeneList.txt",sep = "\t",header = T)

immlist = intersect(rownames(luma_gene),immgene[,1])
test = luma_gene[rownames(luma_gene)%in%immlist,]

test1 = test[,colnames(test)%in%hs6$bcr_patient_barcode[hs6$clust_average_cut ==1]]
test2 = test[,colnames(test)%in%hs6$bcr_patient_barcode[hs6$clust_average_cut ==2]]

test$cluster1 = apply(test1,1,mean)
test$cluster2 = apply(test2,1,mean)
test$gene = rownames(test)

immdata = test[,c("gene","cluster1","cluster2")]
write.table(immdata,file = "immdata.txt",sep = "\t",quote = F,row.names = F)
#QDMR
immqdmr = read.table(file = "SpecificityTable.txt",sep = "\t",quote = "",header = T)
genelist = immqdmr$gene
immtable = luma_gene[rownames(luma_gene)%in%genelist,]
immtable = as.data.frame(t(immtable))
immtable$Labels = ifelse(rownames(immtable)%in%hs6$bcr_patient_barcode[hs6$clust_average_cut ==1],"cluster1","cluster2")

my_comparisons <- list(c("cluster1", "cluster2"))


e<-immtable %>% 
  dplyr::filter(Labels %in% c("cluster1","cluster2")) %>% 
  ggviolin(x = "Labels", y = c(colnames(immtable)[1:19]),
           combine=T, select=c("cluster1","cluster2"),order=c("cluster1","cluster2"),
           palette ="aaas", color="Labels",shape="Labels",
           ylab="beta_value",xlab=FALSE,
           panel.labs = list(Labels=c("cluster1","cluster2")),
           font.label = list(size = 14, face = "bold", color ="red"),
           add = "jitter", add.params = list(fill = "white"))
e+stat_compare_means(method = "t.test",
                     label = "p.signif",##星号设置
                     comparisons = my_comparisons)+theme_gray(base_size = 14)
##tcga和geo箱型图
library(ggplot2)
library(ggpubr)

tcgavio = hs6
rownames(tcgavio) = tcgavio[,1]
tcgavio = tcgavio[,-1]

tcgavio$Labels = ifelse(tcgavio[,"clust_average_cut"] == "1","cluster1","cluster2")
#tcgavio$Labels = ifelse(tcgavio[,"clust_average_cut"] == "2","cluster1","cluster2")

tcgavio = tcgavio[,c(1:22,ncol(tcgavio))]#注意修改

my_comparisons <- list(c("cluster1", "cluster2"))


e<-tcgavio %>% 
  dplyr::filter(Labels %in% c("cluster1","cluster2")) %>% 
  ggboxplot(x = "Labels", y = c(colnames(tcgavio)[1:22]),
            combine=T, select=c("cluster1","cluster2"),order=c("cluster1","cluster2"),
            palette ="aaas", color="Labels",shape="Labels",
            ylab="beta_value",xlab=FALSE,
            panel.labs = list(Labels=c("cluster1","cluster2")),
            font.label = list(size = 14, face = "bold", color ="red"),
            add = "jitter", add.params = list(fill = "white"))
e+stat_compare_means(method = "t.test",
                     label = "p.signif",##星号设置
                     comparisons = my_comparisons)+theme_gray(base_size = 14)





hs6$clust_average_cut = ifelse(hs6[,"clust_average_cut"] == "2",1,2)
##
#tsne
library(Rtsne)

write.table(clust_average_cut,"clust_average_cut.txt",sep='\t')
group=read.table("clust_average_cut.txt",sep='\t',header=TRUE)

cluster = group$x

tsne_obj <- Rtsne(dis_bray, is_distance = TRUE,perplexity = 4)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(group$x),name = rownames(group))
ggplot(tsne_data,aes(x = X, y = Y,color = cluster)) +
  geom_point()+
  scale_color_manual(values = c( "#EE0000FF","#3B4992FF"))

#免疫浸润
count = read.csv(file = "TCGA-BRCA.htseq_counts.tsv",header=TRUE,sep='\t')

gene=count[,1]
gene=gsub("\\..*","",gene)

library(clusterProfiler)
library(org.Hs.eg.db)
library(proxy)
library(philentropy)
library(vegan)  
library(magrittr)
library(dplyr)
##文件为一列gene entrez ID, 文件名设为 enterz.txt
eg <- bitr(gene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")

count$Ensembl_ID = gene
test = merge(x = eg,y = count,by.x = "ENSEMBL",by.y = "Ensembl_ID",all.x = T,all.y = F)
write.table(test,file = "ensembl_symbol.txt",sep = "\t",quote = F)
cluster1 = names(clust_average_cut)[which(clust_average_cut==1)]
cluster2 = names(clust_average_cut)[which(clust_average_cut==2)]

tumor = test[,c(2,which(as.numeric(substr(colnames(test),14,15)) < 10))]

cluster1fpkm = tumor[,c(1,which(substr(colnames(tumor),1,12) %in% cluster1))]
cluster2fpkm = tumor[,c(1,which(substr(colnames(tumor),1,12) %in% cluster2))]

cluster1fpkm =  cluster1fpkm %>% group_by(SYMBOL) %>% summarise_each(funs(mean))
cluster2fpkm =  cluster2fpkm %>% group_by(SYMBOL) %>% summarise_each(funs(mean))
#t = cluster2fpkm%>%group_by(SYMBOL)%>%summarise(y = )
cluster1fpkm = as.data.frame(cluster1fpkm)
cluster2fpkm = as.data.frame(cluster2fpkm)


rownames(cluster1fpkm) = cluster1fpkm$SYMBOL
rownames(cluster2fpkm) = cluster2fpkm$SYMBOL
cluster1fpkm = cluster1fpkm[,-1]
cluster2fpkm = cluster2fpkm[,-1]


cluster1=2^cluster1fpkm-1
cluster2=2^cluster2fpkm-1
tpm1 = t(t(cluster1)/colSums(cluster1))*10^6
tpm2 = t(t(cluster2)/colSums(cluster2))*10^6
write.table(tpm1,"tpm1.txt",sep='\t')
write.table(tpm2,"tpm2.txt",sep='\t')
CIBERSORT("LM22.txt","tpm1.txt",perm = 100,QN = F,"cluster1_CIBERSORT-Results.txt")
CIBERSORT("LM22.txt","tpm2.txt",perm = 100,QN = F,"cluster2_CIBERSORT-Results.txt")

######画图
ciber1=read.table("cluster1_CIBERSORT-Results.txt",sep='\t',header=TRUE)
ciber2=read.table("cluster2_CIBERSORT-Results.txt",sep='\t',header=TRUE)

ciber1=ciber1[,1:23]
ciber2=ciber2[,1:23]

ciber=rbind(ciber1,ciber2)

ciber=cbind(rep(c("cluster1","cluster2"),c(259,92)),ciber)
colnames(ciber)[1]=c("cluster")
ciber=as.matrix(ciber)
ciberplot=ciber[,3:24]
dim(ciberplot) <- c(351*22,1)
ciberplot=cbind(ciberplot,rep(ciber[,1],22))
ciberplot=cbind(ciberplot,rep(colnames(ciber)[3:24],each=351))




colnames(ciberplot)=c("cell_value","cluster","cell_type")
ciberplot=as.data.frame(ciberplot)
ciberplot$cell_value=as.numeric(ciberplot$cell_value)

library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)

P<-ciberplot %>% 
  ## 确定x,y
  ggplot(aes(x = cell_type, y = cell_value, fill = cluster)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "cell_value")+
  scale_x_discrete(name = "cell_type") +
  ggtitle("Comparison of immune cell infiltration") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1,size = 11)) 
P
p2<-P+scale_fill_lancet()
p2
p2+stat_compare_means(method = "wilcox.test",
                      label = "p.signif",##星号设置
                      hide.ns=TRUE       )+
  theme_gray(base_size = 14)+theme(plot.title = element_text(hjust = 0.5,size = 14, face =  "bold"),
                                   text = element_text(size = 12),
                                   axis.title = element_text(face="bold"),
                                   axis.text.x=element_text(angle=45, hjust=1, vjust=1,size = 11)) 
#######训练集测试集的相关性

r1=read.table("cluster2_CIBERSORT-Results_7test.txt",sep='\t',header=TRUE)
r2=read.table("cluster2_CIBERSORT-Results_7test.txt",sep='\t',header=TRUE)

library(vegan)
library(dplyr)
library(corrplot)
corrplot(cor(r1[,2:8],r2[,2:8]), method = "circle", type = 'upper')

#TNM堆积图
library(ggplot2)
library(reshape2)
clin = read.csv("clin.csv")
clin = clin[,c("bcr_patient_barcode","A3_T","A4_N","A5_M")]
clin$bcr_patient_barcode = gsub("-",".",clin$bcr_patient_barcode)
test = hs6[,c("bcr_patient_barcode","clust_average_cut")]
test = merge(test,clin,by = "bcr_patient_barcode",all.x = T,all.y = F)

table(test[test$clust_average_cut == 1,"A3_T"])
table(test[test$clust_average_cut == 1,"A4_N"])
table(test[test$clust_average_cut == 1,"A5_M"])
table(test[test$clust_average_cut == 2,"A3_T"])
table(test[test$clust_average_cut == 2,"A4_N"])
table(test[test$clust_average_cut == 2,"A5_M"])

mydata<-read.csv("TNMdata.csv",sep=",",na.strings="NA",stringsAsFactors=FALSE)

sum<-sort((rowSums(mydata[,2:ncol(mydata)]))/3,index.return=TRUE)

mydata$cluster <- factor(mydata$cluster, levels = mydata$cluster[order(sum$ix)])
mydata<-melt(mydata,id.vars='cluster')
colnames(mydata) = c("cluster","Stage","value")
library(RColorBrewer)
ggplot(data=mydata,aes(cluster,value,fill=Stage))+
  geom_bar(stat="identity",position="fill", color="black", width=0.65,size=0.25)+
  scale_fill_hue(h=c(2,200))+
  coord_flip()+scale_fill_brewer(palette = 'Set3')+
  theme(axis.title=element_text(size=15,face="plain",color="black"),
        axis.text = element_text(size=12,face="plain",color="black"),
        legend.title=element_text(size=13,face="plain",color="black"),
        legend.position = "right")

#雷达图
library(fmsb) 
library(ggsci)
set.seed(99)
tcgavio=hs6[,colnames(hs6)!=c("PFS","PFS.time")]
rownames(tcgavio) = tcgavio$bcr_patient_barcode
tcgavio = tcgavio[,-1]
colnames(tcgavio)[ncol(tcgavio)] = "Labels"
clustermean=aggregate(tcgavio[,1:20],by=list(gene=tcgavio$Labels),FUN=mean)

rownames(clustermean)=clustermean[,1]
clustermean=clustermean[,-1]
data=rbind(rep(0.8,20) , rep(0,20) ,clustermean)
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=pal_aaas("default", alpha = 0.5)(7)
colors_border=pal_aaas("default")(7)
radarchart( data  , axistype=1 , 
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
            
            vlcex=0.8 
)
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)

data1=data[1:3,]
colors_in=pal_aaas("default", alpha = 0.5)(1)
colors_border=pal_aaas("default")(1)
radarchart( data1  ,  
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            
            cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
            
            vlcex=0.8 
)
legend(x=0.7, y=1, legend = rownames(data1[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3) 

data2=data[c(1,2,4),]
colors_in=c("#EE000099")

colors_border=c("#EE0000FF")
radarchart( data2 ,  
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            
            cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
            
            vlcex=0.8 
)
legend(x=0.7, y=1, legend = rownames(data2[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3) 


data3=data[c(1,2,5),]
colors_in=c("#008B4599")

colors_border=c("#008B45FF")
radarchart( data3,  
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            
            cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
            
            vlcex=0.8 
)
legend(x=0.7, y=1, legend = rownames(data3[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)


#annotation生存
annotation = read.table(file = "annotation-tcga.tsv",header = T,sep = "\t")
ann_brca = annotation[annotation$TCGA_project == "BRCA",]

ann_brca = ann_brca[,c("X","MFP","Immune.Subtype","TCGA.Subtype","PFS","PFS_Flag")]
ann_brca = ann_brca[!is.na(ann_brca$PFS),]
ann_brca$X = gsub("-",".",ann_brca$X)
load("BRCA_luma.RData")
load("LumB.RData")

ann_luma = ann_brca[ann_brca$TCGA.Subtype == "BRCA.LumA",]
ann_lumb = ann_brca[ann_brca$TCGA.Subtype == "BRCA.LumB",]

ann_lumb = rbind(ann_lumb,ann_luma[ann_luma$X %in% colnames(lumb),])
ann_luma = ann_luma[-c(which(ann_luma$X %in% colnames(lumb))),]

library(survival)
library(survminer)

hs6 = ann_luma
hs6 = ann_lumb

surv1=Surv(hs6$PFS,hs6$PFS_Flag)
fit=survfit(surv1~Immune.Subtype,data=hs6) 
ggsurvplot( fit, 
            data = hs6, 
            size = 1, # change line size 
            palette = "aaas",
            #conf.int = TRUE, # Add confidence interval 
            pval = TRUE, # Add p-value 
            risk.table = TRUE, # Add risk table 
            risk.table.col = "Immune.Subtype",# Risk table color by groups 
            #legend.labs = c("D","F","IE","IE/F"), # Change legend labels
            legend.labs = c("C1","C2","C3","C4","C6"), # Change legend labels
            risk.table.height = 0.25, # Useful to change when you have multiple groups 
            ggtheme = theme_bw() ,
            pval.method = F,
            test.for.trend = F
)

#CYT评分GZMA和PRF1
test = read.table(file = "ensembl_symbol.txt",sep = "\t")
test = test[,c(2,which(as.numeric(substr(colnames(test),14,15)) <= 9))]

cluster_a = read.table(file = "cluster.txt",sep = "\t",header = T)
cyt = test[test$SYMBOL %in% c("GZMA","PRF1"),which(substr(colnames(test),1,12) %in% cluster_a$bcr_patient_barcode)]
rownames(cyt) = c("GZMA","PRF1")
cyt = as.data.frame(t(cyt))
cyt$bcr_patient_barcode = substr(rownames(cyt),1,12)

out = merge(cyt,cluster_a,by = "bcr_patient_barcode",all.x = T,all.y = F)
out$cyt = rowMeans(out[,2:3])
out$Labels = as.character(out$Labels)
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("1", "2"))
ggplot(out,aes(Labels,cyt,fill = Labels))+
  geom_boxplot()+
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     comparisons = my_comparisons)

