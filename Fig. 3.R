########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
##SBS compare
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- SBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
resistant_matrix <- merge_matrix[which(merge_matrix$y==1),]
sensitive_matrix <- merge_matrix[which(merge_matrix$y==0),]

###############################################################################################################################################
CVgroup <- function(k,datasize,seed){
    set.seed(seed)
	n <- rep(1:k,ceiling(datasize/k))[1:datasize]
	temp <- sample(n,datasize)
	x <- 1:k
	dataseq <- 1:datasize
	cvlist <- lapply(x,function(x) dataseq[temp==x])
	return(cvlist)
}
k=10
datasize <- nrow(merge_matrix)
library(e1071)
datasize <- nrow(merge_matrix)
result <- c()
time <- 1
seed=123

###############################################################################################################################################

###############################################################################################################################################
##RF
library(randomForest)
k=10
result <- c()
time <- 1
seed=123
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			test <- as.data.frame(test)
			train <- rbind(resistant_train,sensitive_train)
			train <- as.data.frame(train)
			train$y <- factor(train$y)
			otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
            pred <- predict(otu_train.forest,test,type="prob")
		    each_result <- data.frame(prediction=pred[,2],real_data=test$y,time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

##Van Allen
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- SBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_SBS_assign_all/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- SBS_Exposure
###############################################################################################################################################
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- patient_1[which(patient_1[,8]!="long-survival"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,8]=="nonresponse")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

###############################################################################################################################################
##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


######################################################################################################################
##Hellmann
######################################################################################################################
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- SBS_Exposure
###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- SBS_Exposure
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
y <- rep(0,nrow(patient))
y[which(patient[,30]=="noresponse")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy
	   
######################################################################################################################
##Rizvi
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- SBS_Exposure
###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- SBS_Exposure
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
y <- rep(0,nrow(patient))
y[which(patient[,15]=="NDB")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy
	   
######################################################################################################################################################
##DBS
######################################################################################################################################################
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- DBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
resistant_matrix <- merge_matrix[which(merge_matrix$y==1),]
sensitive_matrix <- merge_matrix[which(merge_matrix$y==0),]

###############################################################################################################################################
CVgroup <- function(k,datasize,seed){
    set.seed(seed)
	n <- rep(1:k,ceiling(datasize/k))[1:datasize]
	temp <- sample(n,datasize)
	x <- 1:k
	dataseq <- 1:datasize
	cvlist <- lapply(x,function(x) dataseq[temp==x])
	return(cvlist)
}
k=10
datasize <- nrow(merge_matrix)
library(e1071)
datasize <- nrow(merge_matrix)
result <- c()
time <- 1
seed=123

###############################################################################################################################################

##RF
library(randomForest)
k=10
result <- c()
time <- 1
seed=123
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			test <- as.data.frame(test)
			train <- rbind(resistant_train,sensitive_train)
			train <- as.data.frame(train)
			train$y <- factor(train$y)
			otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
            pred <- predict(otu_train.forest,test,type="prob")
		    each_result <- data.frame(prediction=pred[,2],real_data=test$y,time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

##Van Allen
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- DBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_DBS_assign_all/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- DBS_Exposure
###############################################################################################################################################
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- patient_1[which(patient_1[,8]!="long-survival"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,8]=="nonresponse")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

###############################################################################################################################################
##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

######################################################################################################################
##Hellmann
######################################################################################################################
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- DBS_Exposure
###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix

DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_DBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- DBS_Exposure
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
y <- rep(0,nrow(patient))
y[which(patient[,30]=="noresponse")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy
######################################################################################################################
##Rizvi
######################################################################################################################
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- DBS_Exposure
###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_DBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
merge_matrix <- DBS_Exposure
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
y <- rep(0,nrow(patient))
y[which(patient[,15]=="NDB")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

##Fig.3A
SBScompare <- read.table("G:/all_work/mutation signature/immune/new_outputdata/SBS_compare_auc.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
SBScompare[,1] <- round(SBScompare[,1],2)
SBScompare$type <- factor(SBScompare$type,levels=c("SBS+DBS","SBS","DBS"))
SBScompare$data <- factor(SBScompare$data,levels=c("Miao","Allen","Hellmann","Rizvi"))
ggplot(SBScompare, aes(fill=type, y=AUROC, x=data)) + 
    geom_bar(position="dodge", stat="identity")+theme_classic()+
	scale_fill_manual(values=c("#F24F4F","#FCAEAD","#FDD4D0"))+
	theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,color = "black", size = 10))+ylab("AUROC")
##Fig.3B	
SBScompare <- read.table("G:/all_work/mutation signature/immune/new_outputdata/SBS_compare_accuracy.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
SBScompare[,1] <- round(SBScompare[,1],2)
SBScompare$type <- factor(SBScompare$type,levels=c("SBS+DBS","SBS","DBS"))
SBScompare$data <- factor(SBScompare$data,levels=c("Miao","Allen","Hellmann","Rizvi"))
ggplot(SBScompare, aes(fill=type, y=Accuracy, x=data)) + 
    geom_bar(position="dodge", stat="identity")+theme_classic()+
	scale_fill_manual(values=c("#F24F4F","#FCAEAD","#FDD4D0"))+
	theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,color = "black", size = 10))+ylab("Accuracy")
	


########################################################################################################################################
##all mutation
##signature
library(sigminer)
laml <- read_maf(maf = "G:/all_work/mutation signature/immune/Allen_mutation.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
write.table(t(mt_tally_all$all_matrices$SBS_96),"G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_SBS_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
write.table(t(mt_tally_all$all_matrices$DBS_78),"G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_DBS_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "ID")
write.table(t(mt_tally_all$all_matrices$ID_83),"G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_ID_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")
#SBS
from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix","/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_SBS_Roh",
"/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_SBS_Roh.txt",
opportunity_genome="GRCh37",
context_type="96")

#DBS
#from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix","/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_DBS_Roh",
"/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_DBS_Roh.txt",
opportunity_genome="GRCh37",
context_type="DINUC")

#ID
#from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix","/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_ID_Roh",
"/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_ID_Roh.txt",
opportunity_genome="GRCh37",
context_type="ID")


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_SBS_Roh/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_DBS_Roh/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
resistant_matrix <- merge_matrix[which(merge_matrix$y==1),]
sensitive_matrix <- merge_matrix[which(merge_matrix$y==0),]

###############################################################################################################################################
CVgroup <- function(k,datasize,seed){
    set.seed(seed)
	n <- rep(1:k,ceiling(datasize/k))[1:datasize]
	temp <- sample(n,datasize)
	x <- 1:k
	dataseq <- 1:datasize
	cvlist <- lapply(x,function(x) dataseq[temp==x])
	return(cvlist)
}
k=10
datasize <- nrow(merge_matrix)
library(e1071)
datasize <- nrow(merge_matrix)
result <- c()
time <- 1
seed=123

###############################################################################################################################################
##RF
library(randomForest)
k=10
result <- c()
time <- 1
seed=123
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			test <- as.data.frame(test)
			train <- rbind(resistant_train,sensitive_train)
			train <- as.data.frame(train)
			train$y <- factor(train$y)
			otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
            pred <- predict(otu_train.forest,test,type="prob")
		    each_result <- data.frame(prediction=pred[,2],real_data=test$y,time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

######################################################################################################################
##Van Allen
######################################################################################################################
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,8]!="long-survival"),]
mutation <- read.table("G:/all_work/mutation signature/immune/Levi_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
each_mutation <- mutation[which(!is.na(match(mutation[,1],patient[,1]))),]
colnames(each_mutation)[1] <- "Tumor_Sample_Barcode"
write.table(each_mutation,"G:/all_work/mutation signature/immune/Levi_mutation_all/fisher_mutation_all_Roh.txt",col.names=T,row.names=F,quote=F,sep="\t")
library(sigminer)
laml <- read_maf(maf = "G:/all_work/mutation signature/immune/Levi_mutation_all/fisher_mutation_all_Roh.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
matrix_1 <- t(mt_tally_all$all_matrices$SBS_96)[,which(colSums(t(mt_tally_all$all_matrices$SBS_96))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/Levi_mutation_all/fisher_SBS_all_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
matrix_1 <- matrix(t(mt_tally_all$all_matrices$DBS_78)[,which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)],nrow=78)
rownames(matrix_1) <- rownames(t(mt_tally_all$all_matrices$DBS_78))
colnames(matrix_1) <- colnames(t(mt_tally_all$all_matrices$DBS_78))[which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/Levi_mutation_all/fisher_DBS_all_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")

from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/Levi_mutation_all/fisher_SBS_all_Roh.txt",
                    output="/data/home/huangyue/mutation_signature/immune/Levi_mutation_all/fisher_SBS_assign_all_Roh", 
					context_type="96",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_SBS_Roh/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/Levi_mutation_all/fisher_DBS_all_Roh.txt",
                    output="/data/home/huangyue/mutation_signature/immune/Levi_mutation_all/fisher_DBS_assign_all_Roh", 
					context_type="DINUC",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_DBS_Roh/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)

SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_SBS_Roh/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_DBS_Roh/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_all/fisher_SBS_assign_all_Roh/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_all/fisher_DBS_assign_all_Roh/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

###############################################################################################################################################
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- patient_1[which(patient_1[,8]!="long-survival"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,8]=="nonresponse")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

###############################################################################################################################################
##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy
	      
######################################################################################################################
######################################################################################################################
######################################################################################################################
##Hellmann
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
mutation <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
each_mutation <- mutation[which(!is.na(match(mutation[,12],patient[,1]))),]
#each_mutation$End_position <- each_mutation$End_position+1
write.table(each_mutation,"G:/all_work/mutation signature/immune/independent_data/lung_1_all/fisher_mutation_all_Roh.txt",col.names=T,row.names=F,quote=F,sep="\t")
library(sigminer)
laml <- read_maf(maf = "G:/all_work/mutation signature/immune/independent_data/lung_1_all/fisher_mutation_all_Roh.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
matrix_1 <- t(mt_tally_all$all_matrices$SBS_96)[,which(colSums(t(mt_tally_all$all_matrices$SBS_96))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/independent_data/lung_1_all/fisher_SBS_all_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
matrix_1 <- matrix(t(mt_tally_all$all_matrices$DBS_78)[,which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)],nrow=78)
rownames(matrix_1) <- rownames(t(mt_tally_all$all_matrices$DBS_78))
colnames(matrix_1) <- colnames(t(mt_tally_all$all_matrices$DBS_78))[which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/independent_data/lung_1_all/fisher_DBS_all_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")

from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/independent_data/lung_1_all/fisher_SBS_all_Roh.txt",
                    output="/data/home/huangyue/mutation_signature/immune/independent_data/lung_1_all/fisher_SBS_assign_all_Allen_Roh", 
					context_type="96",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_SBS_Roh/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/independent_data/lung_1_all/fisher_DBS_all_Roh.txt",
                    output="/data/home/huangyue/mutation_signature/immune/independent_data/lung_1_all/fisher_DBS_assign_all_Allen_Roh", 
					context_type="DINUC",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_DBS_Roh/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)

######################################################################################################################
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_SBS_Roh/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_DBS_Roh/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1
###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_all/fisher_SBS_assign_all_Allen_Roh/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_all/fisher_DBS_assign_all_Allen_Roh/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1

patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
y <- rep(0,nrow(patient))
y[which(patient[,30]=="noresponse")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

###############################################################################################################################################
##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy	   
######################################################################################################################
##Rizvi
patient <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
mutation <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
each_mutation <- mutation[which(!is.na(match(mutation[,12],patient[,1]))),]
write.table(each_mutation,"G:/all_work/mutation signature/immune/independent_data/lung_2_all/fisher_mutation_all_Roh.txt",col.names=T,row.names=F,quote=F,sep="\t")
library(sigminer)
laml <- read_maf(maf = "G:/all_work/mutation signature/immune/independent_data/lung_2_all/fisher_mutation_all_Roh.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
matrix_1 <- t(mt_tally_all$all_matrices$SBS_96)[,which(colSums(t(mt_tally_all$all_matrices$SBS_96))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/independent_data/lung_2_all/fisher_SBS_all_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
matrix_1 <- matrix(t(mt_tally_all$all_matrices$DBS_78)[,which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)],nrow=78)
rownames(matrix_1) <- rownames(t(mt_tally_all$all_matrices$DBS_78))
colnames(matrix_1) <- colnames(t(mt_tally_all$all_matrices$DBS_78))[which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/independent_data/lung_2_all/fisher_DBS_all_Roh.txt",col.names=T,row.names=T,quote=F,sep="\t")


from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/independent_data/lung_2_all/fisher_SBS_all_Roh.txt",
                    output="/data/home/huangyue/mutation_signature/immune/independent_data/lung_2_all/fisher_SBS_assign_all_Allen_Roh", 
					context_type="96",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_SBS_Roh/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/independent_data/lung_2_all/fisher_DBS_all_Roh.txt",
                    output="/data/home/huangyue/mutation_signature/immune/independent_data/lung_2_all/fisher_DBS_assign_all_Allen_Roh", 
					context_type="DINUC",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_all/fisher_DBS_Roh/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)

######################################################################################################################
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_SBS_Roh/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_all/fisher_DBS_Roh/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_all/fisher_SBS_assign_all_Allen_Roh/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_all/fisher_DBS_assign_all_Allen_Roh/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1

patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
y <- rep(0,nrow(patient))
y[which(patient[,15]=="NDB")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

###############################################################################################################################################
##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy
################################################
##Fig.3C
library(ggradar)
library(dplyr)
library(scales)
library(tibble)
data <- as.data.frame(matrix(c(0.8288,0.6036,0.8175,0.6613,0.7245,0.721,0.8699,0.6719), ncol=4))
data <- cbind(c("MSmodel","allmutation"),data)
colnames(data) <- c("group","Miao" , "Allen" , "Hellmann" , "Rizvi")
ggradar(data,values.radar=c("0","0.5","1"),background.circle.colour="white",group.colours=c("#3D7671","#F24F4F"))
################################################
##Fig.3D
library(ggradar)
library(dplyr)
library(scales)
library(tibble)
data <- as.data.frame(matrix(c(0.7657613,0.6275502,0.80,0.70,0.7123288,0.7066667,0.8214286,0.6333333), ncol=4))
data <- cbind(c("MSmodel","allmutation"),data)
colnames(data) <- c("group","Miao" , "Allen" , "Hellmann" , "Rizvi")
ggradar(data,values.radar=c("0","0.5","1"),background.circle.colour="white",group.colours=c("#3D7671","#F24F4F"))


##################################################################
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

##no normalization
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
#Melanoma_1 <- Melanoma[which(Melanoma[,18]!="stable disease"),]
#Melanoma_1 <- Melanoma[which(Melanoma$cancer_type=="Lung"|Melanoma$cancer_type=="Melanoma"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
resistant_matrix <- merge_matrix[which(merge_matrix$y==1),]
sensitive_matrix <- merge_matrix[which(merge_matrix$y==0),]

###############################################################################################################################################
CVgroup <- function(k,datasize,seed){
    set.seed(seed)
	n <- rep(1:k,ceiling(datasize/k))[1:datasize]
	temp <- sample(n,datasize)
	x <- 1:k
	dataseq <- 1:datasize
	cvlist <- lapply(x,function(x) dataseq[temp==x])
	return(cvlist)
}
k=10
datasize <- nrow(merge_matrix)
library(e1071)
datasize <- nrow(merge_matrix)
result <- c()
time <- 1
seed=123

###############################################################################################################################################
##SVM
library(e1071)
seed=123
result <- c()
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			test$y <- factor(test$y)
			train <- rbind(resistant_train,sensitive_train)
			train$y <- factor(train$y)
			svr.model <- svm(y ~. , train,probability = TRUE)
			summary(svr.model)
            pred1 <- predict(svr.model, test,probability =TRUE)
			pred <- attr(pred1,"probabilities")[,1]
		    each_result <- data.frame(prediction=pred,real_data=test$y,time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)
SVM_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

###############################################################################################################################################
##xgboost
library(xgboost)
seed=123
result <- c()
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			#test$y <- factor(test$y)
			test <- as.matrix(test)
			test_1 <- list(data=test[,-7],label=test[,7])
            dtest <- xgb.DMatrix(data=test_1$data,label=test_1$label)
			train <- rbind(resistant_train,sensitive_train)
			#train$y <- factor(train$y)
			train <- as.matrix(train)
			train_1 <- list(data=train[,-7],label=train[,7])
		    dtrain <- xgb.DMatrix(data=train_1$data,label=train_1$label)
			xgb <- xgboost(data=dtrain,nround=100,objective="binary:logistic",learning_rate=0.1)
			pred <- predict(xgb, newdata=dtest,type="response")
		    each_result <- data.frame(prediction=pred,real_data=test[,7],time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)
xgboost_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)

###############################################################################################################################################
##RF
library(randomForest)
k=10
result <- c()
time <- 1
seed=123
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			test <- as.data.frame(test)
			train <- rbind(resistant_train,sensitive_train)
			train <- as.data.frame(train)
			train$y <- factor(train$y)
			otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
            pred <- predict(otu_train.forest,test,type="prob")
		    each_result <- data.frame(prediction=pred[,2],real_data=test$y,time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy



##Enet
library(glmnet)
library(caret)
k=10
seed=123
result <- c()
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			test <- as.data.frame(test)
			train <- rbind(resistant_train,sensitive_train)
			train <- as.data.frame(train)
			train$y <- factor(train$y)
			model <- train(
                     y~.,data = train, method="glmnet", 
		             trControl = trainControl("cv",number=10),
                     tuneLength = 10
		    )
            pred <- predict(model,test,type="prob")
		    each_result <- data.frame(prediction=pred[,2],real_data=test$y,time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Ridge
library(magrittr)
library(glmnet)
library(caret)
k=10
datasize <- nrow(merge_matrix)
result <- c()
time <- 1
seed=123
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			test <- as.matrix(test)
			train <- rbind(resistant_train,sensitive_train)
			train <- as.matrix(train)
			#train[,18] <- as.character(train[,18])
			set.seed(123456)
            glmnet1 <- cv.glmnet(x=train[,-7],y=as.character(train[,7]),alpha=0,family="binomial")
		    model <- glmnet(x=train[,-7],y=as.character(train[,7]),alpha=0,lambda = glmnet1$lambda.min,family="binomial")
            pred <- model %>% predict(test[,-7],type="response")%>% as.vector()
		    each_result <- data.frame(prediction=pred,real_data=test[,7],time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)


library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

##Logistic
k=10
datasize <- nrow(merge_matrix)
result <- c()
time <- 1
seed=123
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			train <- rbind(resistant_train,sensitive_train)
			set.seed(123456)
            model <- glm(y~., data = train, family = "binomial")
            pred <- model %>% predict(test[,-7],type="response")%>% as.vector()
		    each_result <- data.frame(prediction=pred,real_data=test[,7],time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

##KNN
library(glmnet)
library(caret)
k=10
seed=123
result <- c()
for(i in 1:10){
	cvlist_1 <- CVgroup(k=k,datasize=nrow(resistant_matrix),seed=seed)
	cvlist_2 <- CVgroup(k=k,datasize=nrow(sensitive_matrix),seed=seed)
	#print(cvlist)
	for(j in 1:k){
	    for(m in 1:k){
		    resistant_test <- resistant_matrix[cvlist_1[[j]],]
			resistant_train <- resistant_matrix[-cvlist_1[[j]],]
			sensitive_test <- sensitive_matrix[cvlist_2[[m]],]
			sensitive_train <- sensitive_matrix[-cvlist_2[[m]],]
			test <- rbind(resistant_test,sensitive_test)
			test_1 <- apply(test[,-7],2,scale)
			test_1[,which(colSums(test_1)=="NaN")] <- 0
			test_1 <- as.data.frame(test_1)
			test_1$y <- factor(test$y)
			rownames(test_1) <- rownames(test)
			test <- test_1
			train <- rbind(resistant_train,sensitive_train)
			train_1 <- apply(train[,-7],2,scale)
			train_1[,which(colSums(train_1)=="NaN")] <- 0
			train_1 <- as.data.frame(train_1)
			train_1$y <- factor(train$y)
			rownames(train_1) <- rownames(train)
			train <- train_1
			grid1 <- expand.grid(.k = seq(2, 20, by = 1))
			model <- train(
                     y~.,data = train, method="knn", 
		             trControl = trainControl("cv",number=10),
                     tuneGrid = grid1
		    )
            pred <- predict(model,test,type="prob")
		    each_result <- data.frame(prediction=pred[,2],real_data=test$y,time=time)
		    result <- rbind(result,each_result)
		    time <- time+1
		}
	}
	seed <- seed+1
	print(seed)
}
result <- as.data.frame(result)


library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


######################################################################################################################
##Levi
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
#Melanoma_1 <- patient[which(patient[,18]!="stable disease"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_SBS_assign_all/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_DBS_assign_all/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

##no normalization
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

###############################################################################################################################################
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- patient_1[which(patient_1[,8]!="long-survival"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,8]=="nonresponse")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

library(e1071)
datasize <- nrow(merge_matrix)
result <- c()
##SVM
library(e1071)
set.seed(123)
train$y <- factor(train$y)
svr_model <- svm(y~., data = train,probability = TRUE)
pred <- predict(svr_model, test,probability =TRUE)
pred1 <- attr(pred,"probabilities")[,2]
result <- data.frame(prediction=pred1,real_data=test$y)
result$sample <- rownames(test)
svm_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


###############################################################################################################################################
##xgboost
library(xgboost)
train <- as.matrix(train)
test <- as.matrix(test)
library(xgboost)
train_1 <- list(data=train[,-7],label=train[,7])
dtrain <- xgb.DMatrix(data=train_1$data,label=train_1$label)
test_1 <- list(data=test[,-7],label=test[,7])
dtest <- xgb.DMatrix(data=test_1$data,label=test_1$label)
xgb <- xgboost(data=dtrain,nround=100,objective="binary:logistic",learning_rate=0.1)
pred <- predict(xgb, newdata=dtest,type="class")
result <- data.frame(prediction=pred,real_data=test[,7])
result <- as.data.frame(result)
result$sample <- rownames(test)
xgboost_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy



##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Enet
set.seed(123)
library(glmnet)
library(caret)
train <- as.data.frame(train)
train$y <- factor(train$y)
test <- as.data.frame(test)
model <- train(
         y~.,data = train, method="glmnet",family="binomial",
		 trControl = trainControl("cv",number=10),
         tuneLength = 10
		 )
prediction <- predict(model,test[,-7],type="prob")
result <- data.frame(prediction=prediction[,2],real_data=test$y)
result <- as.data.frame(result)
result$sample <- rownames(test)
Enet_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Ridge
set.seed(123)
train <- as.matrix(train)
test <- as.matrix(test)
glmnet1 <- cv.glmnet(x=train[,-7],y=train[,7],alpha=0,family="binomial")
model <- glmnet(x=train[,-7],y=train[,7],alpha=0,lambda = glmnet1$lambda.min,family="binomial")
pred1 <- predict(model, test[,-7],type="response")
result <- data.frame(prediction=pred1,real_data=test[,7])
result$sample <- rownames(test)
linear_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##logistic
train <- as.data.frame(train)
test <- as.data.frame(test)
model <- glm(y~., data = train, family = "binomial")
pred1 <- predict(model, test[,-7],type="response")
result <- data.frame(prediction=pred1,real_data=test[,7])
result$sample <- rownames(test)
logist_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy



SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
#Melanoma_1 <- patient[which(patient[,18]!="stable disease"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_SBS_assign_all/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_DBS_assign_all/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1
###############################################################################################################################################
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- patient_1[which(patient_1[,8]!="long-survival"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,8]=="nonresponse")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

##KNN
set.seed(123)
library(glmnet)
library(caret)
train <- as.data.frame(train)
train$y <- factor(train$y)
test <- as.data.frame(test)
grid1 <- expand.grid(.k = seq(2, 20, by = 1))
model <- train(
         y~.,data = train, method="knn",
		 trControl = trainControl("cv",number=10),
         tuneGrid = grid1
		 )
prediction <- predict(model,test[,-7],type="prob")
result <- data.frame(prediction=prediction[,2],real_data=test$y)
result <- as.data.frame(result)
result$sample <- rownames(test)

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


######################################################################################################################
##Hellmann
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
#Melanoma_1 <- patient[which(patient[,18]!="stable disease"),]
#Melanoma_1 <- Melanoma[which(Melanoma$cancer_type=="Lung"|Melanoma$cancer_type=="Melanoma"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_DBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
y <- rep(0,nrow(patient))
y[which(patient[,30]=="noresponse")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

###############################################################################################################################################
##SVM
library(e1071)
set.seed(123)
train$y <- factor(train$y)
svr_model <- svm(y~., data = train,probability = TRUE)
pred <- predict(svr_model, test,probability =TRUE)
pred1 <- attr(pred,"probabilities")[,2]
result <- data.frame(prediction=pred1,real_data=test$y)
result$sample <- rownames(test)
svm_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


###############################################################################################################################################
##xgboost
library(xgboost)
train <- as.matrix(train)
test <- as.matrix(test)
library(xgboost)
train_1 <- list(data=train[,-7],label=train[,7])
dtrain <- xgb.DMatrix(data=train_1$data,label=train_1$label)
test_1 <- list(data=test[,-7],label=test[,7])
dtest <- xgb.DMatrix(data=test_1$data,label=test_1$label)
xgb <- xgboost(data=dtrain,nround=100,objective="binary:logistic",learning_rate=0.1)
pred <- predict(xgb, newdata=dtest,type="class")
result <- data.frame(prediction=pred,real_data=test[,7])
result <- as.data.frame(result)
result$sample <- rownames(test)
xgboost_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy



##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Enet
set.seed(123)
library(glmnet)
library(caret)
train <- as.data.frame(train)
train$y <- factor(train$y)
test <- as.data.frame(test)
model <- train(
         y~.,data = train, method="glmnet",family="binomial",
		 trControl = trainControl("cv",number=10),
         tuneLength = 10
		 )
prediction <- predict(model,test[,-7],type="prob")
result <- data.frame(prediction=prediction[,2],real_data=test$y)
result <- as.data.frame(result)
result$sample <- rownames(test)
Enet_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Ridge
set.seed(123)
train <- as.matrix(train)
test <- as.matrix(test)
glmnet1 <- cv.glmnet(x=train[,-7],y=train[,7],alpha=0,family="binomial")
model <- glmnet(x=train[,-7],y=train[,7],alpha=0,lambda = glmnet1$lambda.min,family="binomial")
pred1 <- predict(model, test[,-7],type="response")
result <- data.frame(prediction=pred1,real_data=test[,7])
result$sample <- rownames(test)
linear_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##logistic
train <- as.data.frame(train)
test <- as.data.frame(test)
model <- glm(y~., data = train, family = "binomial")
pred1 <- predict(model, test[,-7],type="response")
result <- data.frame(prediction=pred1,real_data=test[,7])
result$sample <- rownames(test)
logist_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy




##KNN
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1
###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
#Melanoma_1 <- patient[which(patient[,18]!="stable disease"),]
#Melanoma_1 <- Melanoma[which(Melanoma$cancer_type=="Lung"|Melanoma$cancer_type=="Melanoma"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_DBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
y <- rep(0,nrow(patient))
y[which(patient[,30]=="noresponse")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix


set.seed(123)
library(glmnet)
library(caret)
train <- as.data.frame(train)
train$y <- factor(train$y)
test <- as.data.frame(test)
grid1 <- expand.grid(.k = seq(2, 20, by = 1))
model <- train(
         y~.,data = train, method="knn",
		 trControl = trainControl("cv",number=10),
         tuneGrid = grid1
		 )
prediction <- predict(model,test[,-7],type="prob")
result <- data.frame(prediction=prediction[,2],real_data=test$y)
result <- as.data.frame(result)
result$sample <- rownames(test)

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy
	   
######################################################################################################################
##Rizvi
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

##no normalization
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
#Melanoma_1 <- patient[which(patient[,18]!="stable disease"),]
#Melanoma_1 <- Melanoma[which(Melanoma$cancer_type=="Lung"|Melanoma$cancer_type=="Melanoma"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_DBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

##no normalization
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
y <- rep(0,nrow(patient))
y[which(patient[,15]=="NDB")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix

###############################################################################################################################################
##SVM
library(e1071)
set.seed(123)
train$y <- factor(train$y)
svr_model <- svm(y~., data = train,probability = TRUE)
pred <- predict(svr_model, test,probability =TRUE)
pred1 <- attr(pred,"probabilities")[,2]
result <- data.frame(prediction=pred1,real_data=test$y)
result$sample <- rownames(test)
svm_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


###############################################################################################################################################
##xgboost
library(xgboost)
train <- as.matrix(train)
test <- as.matrix(test)
library(xgboost)
train_1 <- list(data=train[,-7],label=train[,7])
dtrain <- xgb.DMatrix(data=train_1$data,label=train_1$label)
test_1 <- list(data=test[,-7],label=test[,7])
dtest <- xgb.DMatrix(data=test_1$data,label=test_1$label)
xgb <- xgboost(data=dtrain,nround=100,objective="binary:logistic",learning_rate=0.1)
pred <- predict(xgb, newdata=dtest,type="class")
result <- data.frame(prediction=pred,real_data=test[,7])
result <- as.data.frame(result)
result$sample <- rownames(test)
xgboost_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy



##RF
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],real_data=test$y)
result$sample <- rownames(test)
RF_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Enet
set.seed(123)
library(glmnet)
library(caret)
train <- as.data.frame(train)
train$y <- factor(train$y)
test <- as.data.frame(test)
model <- train(
         y~.,data = train, method="glmnet",family="binomial",
		 trControl = trainControl("cv",number=10),
         tuneLength = 10
		 )
prediction <- predict(model,test[,-7],type="prob")
result <- data.frame(prediction=prediction[,2],real_data=test$y)
result <- as.data.frame(result)
result$sample <- rownames(test)
Enet_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Ridge
set.seed(123)
train <- as.matrix(train)
test <- as.matrix(test)
glmnet1 <- cv.glmnet(x=train[,-7],y=train[,7],alpha=0,family="binomial")
model <- glmnet(x=train[,-7],y=train[,7],alpha=0,lambda = glmnet1$lambda.min,family="binomial")
pred1 <- predict(model, test[,-7],type="response")
result <- data.frame(prediction=pred1,real_data=test[,7])
result$sample <- rownames(test)
linear_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##logistic
train <- as.data.frame(train)
test <- as.data.frame(test)
model <- glm(y~., data = train, family = "binomial")
pred1 <- predict(model, test[,-7],type="response")
result <- data.frame(prediction=pred1,real_data=test[,7])
result$sample <- rownames(test)
logist_result <- result

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##KNN
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

##no normalization
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure

merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1
###############################################################################################################################################
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
#Melanoma_1 <- patient[which(patient[,18]!="stable disease"),]
#Melanoma_1 <- Melanoma[which(Melanoma$cancer_type=="Lung"|Melanoma$cancer_type=="Melanoma"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_DBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)

##no normalization
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
merge_matrix_1 <- apply(merge_matrix,2,scale)
rownames(merge_matrix_1) <- rownames(merge_matrix)
merge_matrix <- merge_matrix_1
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
y <- rep(0,nrow(patient))
y[which(patient[,15]=="NDB")] <- 1
patient$y <- y
each_drug_1 <- patient[match(rownames(merge_matrix),patient[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
test <- merge_matrix
test[,5] <- 0

set.seed(123)
library(glmnet)
library(caret)
train <- as.data.frame(train)
train$y <- factor(train$y)
test <- as.data.frame(test)
grid1 <- expand.grid(.k = seq(2, 20, by = 1))
model <- train(
         y~.,data = train, method="knn",
		 trControl = trainControl("cv",number=10),
         tuneGrid = grid1
		 )
prediction <- predict(model,test[,-7],type="prob")
result <- data.frame(prediction=prediction[,2],real_data=test$y)
result <- as.data.frame(result)
result$sample <- rownames(test)

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1],direction="<")
auc <- round(auc(rocobj),4)
auc
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy
		
###################################################
##Fig.3E
mydata <- read.table("G:/all_work/mutation signature/immune/new_outputdata/model_auc.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
AUC_matrix <- matrix(rep(0,length(unique(mydata[,2]))*length(unique(mydata[,3]))),ncol=length(unique(mydata[,3])))
rownames(AUC_matrix) <- unique(mydata[,2])
colnames(AUC_matrix) <- unique(mydata[,3])
AUC_matrix <- AUC_matrix[,c(3,2,1,4,5,6,7)]
for(i in 1:nrow(AUC_matrix)){
    each_data <- mydata[which(mydata[,2]==rownames(AUC_matrix)[i]),]
	AUC_matrix[i,] <- each_data[match(colnames(AUC_matrix),each_data[,3]),1]
}

library(pheatmap)
pheatmap(AUC_matrix, display_numbers = TRUE,number_color = "black",number_format="%.4f",border_color="white",
         cluster_rows=FALSE, cluster_cols=FALSE,color = colorRampPalette(c("#ABD9C3", "white", "#D6A67F"))(50))
###################################################
##Fig.3F		 
mydata <- read.table("G:/all_work/mutation signature/immune/new_outputdata/model_acc.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
mydata[,1] <- round(mydata[,1],4)
AUC_matrix <- matrix(rep(0,length(unique(mydata[,2]))*length(unique(mydata[,3]))),ncol=length(unique(mydata[,3])))
rownames(AUC_matrix) <- unique(mydata[,2])
colnames(AUC_matrix) <- unique(mydata[,3])
AUC_matrix <- AUC_matrix[,c(3,2,1,4,5,6,7)]
for(i in 1:nrow(AUC_matrix)){
    each_data <- mydata[which(mydata[,2]==rownames(AUC_matrix)[i]),]
	AUC_matrix[i,] <- each_data[match(colnames(AUC_matrix),each_data[,3]),1]
}

library(pheatmap)
pheatmap(AUC_matrix, display_numbers = TRUE,number_color = "black",number_format="%.4f",border_color="white",
         cluster_rows=FALSE, cluster_cols=FALSE,color = colorRampPalette(c("#DEC4C3", "white", "#D6A67F"))(50))
