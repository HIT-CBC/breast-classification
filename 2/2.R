##Boruta
library(Boruta)
load("BRCA_luma.RData")
normal = read.table(file = "Normal.txt",header = T,sep = "\t")
train = sample(luma,244,replace = F)
test = luma[,which(colnames(luma) %in% setdiff(colnames(luma),colnames(train)))]

diff = cbind(train,normal)

Pvalue = apply(diff,1,function(x){
        t.test(x[1:ncol(train)],x[(ncol(train)+1):ncol(diff)])$p.value
})
gap = apply(train,1,mean)-apply(normal,1,mean)

fdr = p.adjust(Pvalue, method = "BH",length(Pvalue))

train_cg = train[which(fdr <= 0.05 & abs(gap) >=0.2),]

type=rep(c("luma","normal"),c(ncol(train),ncol(normal)))

normal = cbind(rownames(normal),normal)
colnames(normal)[1] = "cg"

train_cg = cbind(rownames(train_cg),train_cg)
colnames(train_cg)[1] = "cg"

tdifftrain1 = merge(train_cg,normal,by = "cg",all.x = F,all.y = F)
normal = normal[,-1]

rownames(tdifftrain1) = tdifftrain1[,"cg"]
tdifftrain1 = as.data.frame(t(tdifftrain1[,-1]))
tdifftrain1=cbind(tdifftrain1,type)
tdifftrain1[,ncol(tdifftrain1)]=as.factor(tdifftrain1[,ncol(tdifftrain1)])

Boruta(type~.,data=tdifftrain1,doTrace=2)->Boruta.train.extended

boruta_signif <- names(Boruta.train.extended$finalDecision[Boruta.train.extended$finalDecision %in% c("Confirmed", "Tentative")])

#分类器
load(file = "boruta.RData")
load(file = "cancer_5types.RData")
unluma = setdiff(colnames(cancer),union(colnames(test),colnames(train)))
unluma_train = sample(unluma,244,replace = F)
unluma_test = setdiff(unluma,unluma_train)

#非a样本不够，也按7：3分类，均分到a样本中

traindata = cancer[which(rownames(cancer) %in% boruta_signif),which(colnames(cancer) %in% union(colnames(train),unluma_train))]
testdata = cancer[which(rownames(cancer) %in% boruta_signif),which(colnames(cancer) %in% union(colnames(test),unluma_test))]

traindata = as.data.frame(t(rbind(ifelse(colnames(traindata) %in% unluma ,"unluma","luma"),traindata)))
testdata = as.data.frame(t(rbind(ifelse(colnames(testdata) %in% unluma ,"unluma","luma"),testdata)))

traindata[,2:ncol(traindata)] = apply(traindata[,2:ncol(traindata)],2,as.numeric)
testdata[,2:ncol(testdata)] = apply(testdata[,2:ncol(testdata)],2,as.numeric)

colnames(traindata)[1] = "type"
colnames(testdata)[1] = "type"

traindata$type=as.factor(traindata$type)     #样本集
testdata$type=as.factor(testdata$type)       #验证集

library(caret)

inTrain<- createDataPartition(y=traindata$type,p=0.75,list=FALSE)

training<- traindata[inTrain, ]              #训练集
testing<- traindata[-inTrain, ]              #测试集

fitControl <- trainControl(## 10-fold CV
        method = "repeatedcv",
        number = 10,
        ## repeated ten times
        repeats = 3,
        classProbs = TRUE,
        summaryFunction = twoClassSummary)

modelFit=train(type~.,data=training,method="svmRadial",metric = "ROC",trControl = fitControl) 

atesting = confusionMatrix(predict(modelFit,newdata=testing),testing$type)$byClass['Balanced Accuracy']
atraining = confusionMatrix(predict(modelFit,newdata=training),training$type)$byClass['Balanced Accuracy']
atest = confusionMatrix(predict(modelFit,newdata=testdata),testdata$type)$byClass['Balanced Accuracy']
#AUC值ROC曲线

library(pROC)

prediction1 <- predict(modelFit,newdata=training)

roctest=roc(training$type,factor(prediction1,ordered=T))
plot.roc(roctest,col="red",print.auc =TRUE,print.auc.col = "black",lwd = 3,identity.lwd = 1.5)

#keygene
cgtest = read.table(file ="cgtest.txt",header = F,sep = "\t")
boruta_gene = unique(cgtest$V3[which(cgtest$V1 %in% boruta_signif)])

write.table(boruta_gene,file = "boruta_gene.txt",sep = "\t",row.names = F,col.names = F,quote = F)






