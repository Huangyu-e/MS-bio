##Fisher
Allen_mutation_1 <- read.table("G:/all_work/mutation signature/immune/Allen_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
resistant_patient <- Melanoma_1[which(Melanoma_1[,19]=="no clinical benefit"),1]
sensitive_patient <- Melanoma_1[which(Melanoma_1[,19]=="clinical benefit"),1]
Allen_mutation <- Allen_mutation_1[which(!is.na(match(Allen_mutation_1[,24],Melanoma_1[,1]))),]
sel_gene <- c()
for(gene in unique(Allen_mutation[,1])){
    each_gene <- Allen_mutation[which(Allen_mutation[,1]==gene),]
	mut_patient <- unique(each_gene[,24])
	nonmutat_patient <- setdiff(Allen_mutation[,24],mut_patient)
	a <- length(intersect(resistant_patient,mut_patient))
	b <- length(intersect(sensitive_patient,mut_patient))
	c <- length(intersect(resistant_patient,nonmutat_patient))
	d <- length(intersect(sensitive_patient,nonmutat_patient))
	p <- fisher.test(matrix(c(a,b,c,d),ncol=2))$p.value
	if(p < 0.05){
	    each_gene <- c(gene,p)
		sel_gene <- rbind(sel_gene,each_gene)
	}
}
sel_gene <- as.data.frame(sel_gene)##719
sel_gene[,2] <- as.numeric(sel_gene[,2])
#sel_gene[,2] <- round(sel_gene[,2],3)
write.table(sel_gene,"G:/all_work/mutation signature/immune/Allen_mutation_Roh/sel_gene.txt",col.names=F,row.names=F,quote=F,sep="\t")
sel_mutation <- Allen_mutation[which(!is.na(match(Allen_mutation[,1],sel_gene[,1]))),]
write.table(sel_mutation,"G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_mutation_all.txt",col.names=T,row.names=F,quote=F,sep="\t")
##signature
library(sigminer)
laml <- read_maf(maf = "G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_mutation_all.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
write.table(t(mt_tally_all$all_matrices$SBS_96),"G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all.txt",col.names=T,row.names=T,quote=F,sep="\t")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
write.table(t(mt_tally_all$all_matrices$DBS_78),"G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_DBS_all.txt",col.names=T,row.names=T,quote=F,sep="\t")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "ID")
write.table(t(mt_tally_all$all_matrices$ID_83),"G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_ID_all.txt",col.names=T,row.names=T,quote=F,sep="\t")
#SBS
from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix","/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_SBS_all_100",
"/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_SBS_all.txt",
opportunity_genome="GRCh37",
context_type="96")

#DBS
#from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix","/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_DBS_all_100",
"/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_DBS_all.txt",
opportunity_genome="GRCh37",
context_type="DINUC")

#ID
#from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor("matrix","/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_ID_all_100",
"/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_ID_all.txt",
opportunity_genome="GRCh37",
context_type="ID")


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
			otu_train.forest <- randomForest(y~., data = train, importance = TRUE,)
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
rocobj <- roc(result[,2],result[,1])
######################################
##Fig.2A
ggroc(rocobj,linetype=1,size=2,legacy.axes = T,color="#F24F4F")+
geom_abline(intercept = 0,#截距
              slope = 1,#斜率
              color = "black",
              size = 1,
              linetype = "dashed")+#类型
labs(x = "1 - Specificity",
y = "Sensivity")+
ggtitle("ROC Curve")+
				  #scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+
				  #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+#x,y轴范围
  theme_classic()+
  #scale_colour_manual(values=c("#F24E4E","#2C64A0","#5D4F8E"))+#换金钱色
  
  theme_bw()+
  theme(panel.grid=element_blank())+
  annotate("text",x=0.75,y=0.125,label=paste("AUROC = ", round(rocobj$auc,4)))


##Van Allen
######################################################################################################################
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,8]!="long-survival"),]
mutation <- read.table("G:/all_work/mutation signature/immune/Levi_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
each_mutation <- mutation[which(!is.na(match(mutation[,1],patient[,1]))),]
sel_mutation <- each_mutation[which(!is.na(match(each_mutation[,2],sel_gene[,1]))),]
colnames(sel_mutation)[1] <- "Tumor_Sample_Barcode"
write.table(sel_mutation,"G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_mutation_all.txt",col.names=T,row.names=F,quote=F,sep="\t")
library(sigminer)
laml <- read_maf(maf = "G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_mutation_all.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
matrix_1 <- t(mt_tally_all$all_matrices$SBS_96)[,which(colSums(t(mt_tally_all$all_matrices$SBS_96))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_SBS_all.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
matrix_1 <- matrix(t(mt_tally_all$all_matrices$DBS_78)[,which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)],nrow=78)
rownames(matrix_1) <- rownames(t(mt_tally_all$all_matrices$DBS_78))
colnames(matrix_1) <- colnames(t(mt_tally_all$all_matrices$DBS_78))[which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_DBS_all.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "ID")
matrix_1 <- t(mt_tally_all$all_matrices$ID_83)[,which(colSums(t(mt_tally_all$all_matrices$ID_83))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/Levi_mutation_Roh/fisher_ID_all.txt",col.names=T,row.names=T,quote=F,sep="\t")

from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/Levi_mutation_Roh/fisher_SBS_all.txt",
                    output="/data/home/huangyue/mutation_signature/immune/Levi_mutation_Roh/fisher_SBS_assign_all", 
					context_type="96",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/Levi_mutation_Roh/fisher_DBS_all.txt",
                    output="/data/home/huangyue/mutation_signature/immune/Levi_mutation_Roh/fisher_DBS_assign_all", 
					context_type="DINUC",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/Levi_mutation_Roh/fisher_ID_all.txt",
                    output="/data/home/huangyue/mutation_signature/immune/Levi_mutation_Roh/fisher_ID_assign_all", 
					context_type="ID",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_ID_all_100/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
					

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

###############################################################################################################################################
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

library(pROC)
library(ggplot2)
result[,1] <- as.numeric(result[,1])
rocobj <- roc(result[,2],result[,1])
######################################
##Fig.2B
ggroc(rocobj,linetype=1,size=2,legacy.axes = T,color="#F24F4F")+
geom_abline(intercept = 0,#截距
              slope = 1,#斜率
              color = "black",
              size = 1,
              linetype = "dashed")+#类型
labs(x = "1 - Specificity",
y = "Sensivity")+
ggtitle("ROC Curve")+
				  #scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+
				  #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+#x,y轴范围
  theme_classic()+
  #scale_colour_manual(values=c("#F24E4E","#2C64A0","#5D4F8E"))+#换金钱色
  
  theme_bw()+
  theme(panel.grid=element_blank())+
  annotate("text",x=0.75,y=0.125,label=paste("AUROC = ", round(rocobj$auc,4)))

######################################
##Fig.2E
result$type <- "nonresponse"
result[which(result[,1]<0.5),5] <- "response"
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- patient_1[which(patient_1[,8]!="long-survival"),]
follow_up <- merge(result,Melanoma_1,by.x="sample",by.y="patient")
library(survival)
library(survminer)
fit <- survfit(Surv(overall_survival, dead) ~ type, data = follow_up)
ggsurvplot(fit,
       pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", # Specify median survival
       ggtheme = theme_classic(), # Change ggplot2 theme
       palette = c("#F7B882","#7FC17F")
       )
######################################
##Fig.2F
fit <- survfit(Surv(progression_free, progression) ~ type, data = follow_up)
ggsurvplot(fit,
       pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", # Specify median survival
       ggtheme = theme_classic(), # Change ggplot2 theme
       palette = c("#F7B882","#7FC17F")
       )

######################################################################################################################
##Hellmann
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
mutation <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
each_mutation <- mutation[which(!is.na(match(mutation[,12],patient[,1]))),]
sel_mutation <- each_mutation[which(!is.na(match(each_mutation[,1],sel_gene[,1]))),]
write.table(sel_mutation,"G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_mutation_all.txt",col.names=T,row.names=F,quote=F,sep="\t")
library(sigminer)
laml <- read_maf(maf = "G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_mutation_all.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
matrix_1 <- t(mt_tally_all$all_matrices$SBS_96)[,which(colSums(t(mt_tally_all$all_matrices$SBS_96))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_SBS_all.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
matrix_1 <- matrix(t(mt_tally_all$all_matrices$DBS_78)[,which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)],nrow=78)
rownames(matrix_1) <- rownames(t(mt_tally_all$all_matrices$DBS_78))
colnames(matrix_1) <- colnames(t(mt_tally_all$all_matrices$DBS_78))[which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_DBS_all.txt",col.names=T,row.names=T,quote=F,sep="\t")


from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/independent_data/lung_1_Roh/fisher_SBS_all.txt",
                    output="/data/home/huangyue/mutation_signature/immune/independent_data/lung_1_Roh/fisher_SBS_assign_all_Allen", 
					context_type="96",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/independent_data/lung_1_Roh/fisher_DBS_all.txt",
                    output="/data/home/huangyue/mutation_signature/immune/independent_data/lung_1_Roh/fisher_DBS_assign_all_Allen", 
					context_type="DINUC",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)

######################################################################################################################
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
######################################
##Fig.2C
ggroc(rocobj,linetype=1,size=2,legacy.axes = T,color="#F24F4F")+
geom_abline(intercept = 0,#截距
              slope = 1,#斜率
              color = "black",
              size = 1,
              linetype = "dashed")+#类型
labs(x = "1 - Specificity",
y = "Sensivity")+
ggtitle("ROC Curve")+
				  #scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+
				  #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+#x,y轴范围
  theme_classic()+
  #scale_colour_manual(values=c("#F24E4E","#2C64A0","#5D4F8E"))+#换金钱色
  
  theme_bw()+
  theme(panel.grid=element_blank())+
  annotate("text",x=0.75,y=0.125,label=paste("AUROC = ", round(rocobj$auc,4)))

######################################
##Fig.2G
result$type <- "nonresponse"
result[which(result[,1]<0.5),5] <- "response"
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
follow_up <- merge(result,patient,by.x="sample",by.y="Primary_Sample_Barcode")
library(survival)
library(survminer)
fit <- survfit(Surv(PFS.months., Progression_status) ~ type, data = follow_up)
ggsurvplot(fit,
       pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", # Specify median survival
       ggtheme = theme_classic(), # Change ggplot2 theme
       palette = c("#F7B882","#7FC17F")
       )
	   
	   
######################################################################################################################
##Rizvi
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
mutation <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
each_mutation <- mutation[which(!is.na(match(mutation[,12],patient[,1]))),]
sel_mutation <- each_mutation[which(!is.na(match(each_mutation[,1],sel_gene[,1]))),]
write.table(sel_mutation,"G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_mutation_all.txt",col.names=T,row.names=F,quote=F,sep="\t")
library(sigminer)
laml <- read_maf(maf = "G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_mutation_all.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
matrix_1 <- t(mt_tally_all$all_matrices$SBS_96)[,which(colSums(t(mt_tally_all$all_matrices$SBS_96))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_SBS_all.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
matrix_1 <- matrix(t(mt_tally_all$all_matrices$DBS_78)[,which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)],nrow=78)
rownames(matrix_1) <- rownames(t(mt_tally_all$all_matrices$DBS_78))
colnames(matrix_1) <- colnames(t(mt_tally_all$all_matrices$DBS_78))[which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)]
write.table(matrix_1,"G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_DBS_all.txt",col.names=T,row.names=T,quote=F,sep="\t")


from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/independent_data/lung_2_Roh/fisher_SBS_all.txt",
                    output="/data/home/huangyue/mutation_signature/immune/independent_data/lung_2_Roh/fisher_SBS_assign_all_Allen", 
					context_type="96",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/independent_data/lung_2_Roh/fisher_DBS_all.txt",
                    output="/data/home/huangyue/mutation_signature/immune/independent_data/lung_2_Roh/fisher_DBS_assign_all_Allen", 
					context_type="DINUC",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)

######################################################################################################################
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
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,19]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix


SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_DBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
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
######################################
##Fig.2D
ggroc(rocobj,linetype=1,size=2,legacy.axes = T,color="#F24F4F")+
geom_abline(intercept = 0,#截距
              slope = 1,#斜率
              color = "black",
              size = 1,
              linetype = "dashed")+#类型
labs(x = "1 - Specificity",
y = "Sensivity")+
ggtitle("ROC Curve")+
				  #scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+
				  #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+#x,y轴范围
  theme_classic()+
  #scale_colour_manual(values=c("#F24E4E","#2C64A0","#5D4F8E"))+#换金钱色
  
  theme_bw()+
  theme(panel.grid=element_blank())+
  annotate("text",x=0.75,y=0.125,label=paste("AUROC = ", round(rocobj$auc,4)))

######################################
##Fig.2H
result$type <- "nonresponse"
result[which(result[,1]<0.5),5] <- "response"
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
follow_up <- merge(result,patient,by.x="sample",by.y="Study.ID")
library(survival)
library(survminer)
fit <- survfit(Surv(PFS..mos., Event....) ~ type, data = follow_up)
ggsurvplot(fit,
       pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", # Specify median survival
       ggtheme = theme_classic(), # Change ggplot2 theme
       palette = c("#F7B882","#7FC17F")
       )