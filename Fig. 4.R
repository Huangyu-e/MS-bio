
##########################################################################################################################################################
##TMB
##########################################################################################################################################################
##Miao
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
library(maftools)
laml <- read.maf("G:/all_work/mutation signature/immune/Allen_mutation.txt")
TMB <- tmb(laml)
TMB <- as.data.frame(TMB)
patient <- merge(Melanoma_1,TMB,by.x="pair_id",by.y="Tumor_Sample_Barcode")
patient$TMB <- patient$total/38
library(pROC)
library(ggplot2)
patient$y <- rep(0,nrow(patient))
patient[which(patient[,19]=="no clinical benefit"),27] <- 1
patient[,25] <- as.numeric(patient[,25])
rocobj <- roc(patient[,27],patient[,25],direction=">")
auc <- round(auc(rocobj),4)
auc
##Van Allen
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,8]!="long-survival"),]
patient$TMB <- patient$nonsynonymous/38
	   
library(pROC)
library(ggplot2)
patient$y <- rep(0,nrow(patient))
patient[which(patient[,8]=="nonresponse"),43] <- 1
patient[,41] <- as.numeric(patient[,41])
rocobj <- roc(patient[,43],patient[,41],direction=">")
auc <- round(auc(rocobj),4)
auc

##Fig.4B
new_patient <- merge(patient,RF_result,by.x="patient",by.y="sample")
new_patient$log10TMB <- log10(new_patient$TMB)
new_patient[which(new_patient$log10TMB>1&new_patient$prediction>0.5),47] <- "HT-NR"
new_patient[which(new_patient$log10TMB>1&new_patient$prediction<=0.5),47] <- "HT-R"
new_patient[which(new_patient$log10TMB<=1&new_patient$prediction<=0.5),47] <- "LT-R"
new_patient[which(new_patient$log10TMB<=1&new_patient$prediction>0.5),47] <- "LT-NR"
new_patient$V47 <- factor(new_patient$V47,levels=c("HT-R","HT-NR","LT-R","LT-NR"))

ggplot(data=new_patient, aes(log10TMB, prediction,color=V47)) +
  geom_point(alpha = 0.6,size=3)+
  geom_vline(xintercept = 1, color="grey40", 
             linetype="longdash", lwd = 0.5) +
  geom_hline(yintercept = 0.5, color="grey40", 
             linetype="longdash", lwd = 0.5)+
  theme_bw(base_size = 12#, base_family = "Times" #修改字体
           ) +
  theme(panel.grid=element_blank())+
  scale_color_manual(values=c("#438963","#AD223E","#E6AA3C","#2C65A0"))
##Fig.4E  
new_patient_1 <- new_patient[which(new_patient$V47=="HT-NR"|new_patient$V47=="LT-R"),]
library(survival)
library(survminer)
fit <- survfit(Surv(overall_survival, dead) ~ V47, data = new_patient_1)
ggsurvplot(fit,
       pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #surv.median.line = "hv", # Specify median survival
       ggtheme = theme_classic(), # Change ggplot2 theme
       palette = c("#AD223E", "#E6AA3C")
       )
##Fig.4F	   
library(survival)
library(survminer)
fit <- survfit(Surv(progression_free, progression) ~ V47, data = new_patient_1)
ggsurvplot(fit,
       pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       #surv.median.line = "hv", # Specify median survival
       ggtheme = theme_classic(), # Change ggplot2 theme
       palette = c("#AD223E", "#E6AA3C")
       )


##lHellmann
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
patient$TMB1 <- patient$TMB/38
	   
library(pROC)
library(ggplot2)
patient$y <- rep(0,nrow(patient))
patient[which(patient[,30]=="noresponse"),35] <- 1
patient[,33] <- as.numeric(patient[,33])
rocobj <- roc(patient[,35],patient[,33],direction=">")
auc <- round(auc(rocobj),4)
auc


##Fig.4C
new_patient <- merge(patient,RF_result,by.x="Primary_Sample_Barcode",by.y="sample")
new_patient$log10TMB <- log10(new_patient$TMB1)
new_patient[which(new_patient$log10TMB>1&new_patient$prediction>0.5),39] <- "HT-NR"
new_patient[which(new_patient$log10TMB>1&new_patient$prediction<=0.5),39] <- "HT-R"
new_patient[which(new_patient$log10TMB<=1&new_patient$prediction<=0.5),39] <- "LT-R"
new_patient[which(new_patient$log10TMB<=1&new_patient$prediction>0.5),39] <- "LT-NR"
new_patient$V39 <- factor(new_patient$V39,levels=c("HT-R","HT-NR","LT-R","LT-NR"))

ggplot(data=new_patient, aes(log10TMB, prediction,color=V39)) +
  geom_point(alpha = 0.6,size=3)+
  geom_vline(xintercept = 1, color="grey40", 
             linetype="longdash", lwd = 0.5) +
  geom_hline(yintercept = 0.5, color="grey40", 
             linetype="longdash", lwd = 0.5)+
  theme_bw(base_size = 12#, base_family = "Times" #修改字体
           ) +
  theme(panel.grid=element_blank())+
  scale_color_manual(values=c("#438963","#E6AA3C","#2C65A0"))


##Rizvi
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
patient$TMB <- patient$Nonsyn./38
	   
library(pROC)
library(ggplot2)
patient$y <- rep(0,nrow(patient))
patient[which(patient[,15]=="NDB"),22] <- 1
patient[,20] <- as.numeric(patient[,20])
rocobj <- roc(patient[,22],patient[,20],direction=">")
auc <- round(auc(rocobj),4)
auc


##Fig.4D
new_patient <- merge(patient,RF_result,by.x="Study.ID",by.y="sample")
new_patient$log10TMB <- log10(new_patient$TMB)
new_patient[which(new_patient$log10TMB>1&new_patient$prediction>0.5),26] <- "HT-NR"
new_patient[which(new_patient$log10TMB>1&new_patient$prediction<=0.5),26] <- "HT-R"
new_patient[which(new_patient$log10TMB<=1&new_patient$prediction<=0.5),26] <- "LT-R"
new_patient[which(new_patient$log10TMB<=1&new_patient$prediction>0.5),26] <- "LT-NR"
new_patient$V26 <- factor(new_patient$V26,levels=c("HT-R","HT-NR","LT-R","LT-NR"))
ggplot(data=new_patient, aes(log10TMB, prediction,color=V26)) +
  geom_point(alpha = 0.6,size=3)+
  geom_vline(xintercept = 1, color="grey40", 
             linetype="longdash", lwd = 0.5) +
  geom_hline(yintercept = 0.5, color="grey40", 
             linetype="longdash", lwd = 0.5)+
  theme_bw(base_size = 12#, base_family = "Times" #修改字体
           ) +
  theme(panel.grid=element_blank())+
  scale_color_manual(values=c("#438963","#E6AA3C","#2C65A0"))
  
##Fig.4A
library(ggradar)
library(dplyr)
library(scales)
library(tibble)
data <- as.data.frame(matrix(c(0.8288,0.6399,0.8175,0.6748,0.7245,0.7484,0.8699,0.8403), ncol=4))
data <- cbind(c("MSmodel","TMB"),data)
colnames(data) <- c("group","Miao" , "Allen" , "Hellmann" , "Rizvi")
ggradar(data,values.radar=c("0","0.5","1"),background.circle.colour="white",group.colours=c("#F24F4F","#2C649F"))
##########################################################################################################################################################
##mutation
##########################################################################################################################################################
##Miao
Allen_mutation_1 <- read.table("G:/all_work/mutation signature/immune/Allen_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient$response <- 0
patient[which(patient[,19]=="no clinical benefit"),22] <- 1

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="B2M"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 1
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="JAK1"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 1
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="JAK2"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 1
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="PTEN"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 1
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="STK11"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 1
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="MUC16"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="NOTCH4"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="POLE"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="POLD1"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="TP53"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="KRAS"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="EGFR"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 1
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="COL3A1"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- Allen_mutation_1[which(Allen_mutation_1[,1]=="NRAS"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 1
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

DDR_gene <- c("MSH2","MSH6","PMS2","POLE","BRCA2")
mut_patient <- Allen_mutation_1[which(!is.na(match(Allen_mutation_1[,1],DDR_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

NER_gene <- c("ERCC3","ERCC2","ERCC5","ERCC4")
mut_patient <- Allen_mutation_1[which(!is.na(match(Allen_mutation_1[,1],NER_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

HR_gene <- c("RAD51","RAD51B","RAD51D","XRCC2","RAD52","RAD54L","RAD50","MRE11","NBN","ARID1A","BLM","BRCA2")
mut_patient <- Allen_mutation_1[which(!is.na(match(Allen_mutation_1[,1],HR_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

gene52 <- read.table("G:/all_work/mutation signature/immune/gene52.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
mut_patient <- Allen_mutation_1[which(!is.na(match(Allen_mutation_1[,1],gene52[,1]))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,24]))),23] <- 0
TP <- length(which(patient[,22]==1&patient[,23]==1))
FP <- length(which(patient[,22]==0&patient[,23]==1))
FN <- length(which(patient[,22]==1&patient[,23]==0))
TN <- length(which(patient[,22]==0&patient[,23]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Van Allen
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,8]!="long-survival"),]
mutation <- read.table("G:/all_work/mutation signature/immune/Levi_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)

patient$response <- 0
patient[which(patient[,8]=="nonresponse"),41] <- 1

mut_patient <- mutation[which(mutation[,2]=="B2M"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 1
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="JAK1"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 1
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="JAK2"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 1
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="PTEN"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 1
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="STK11"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 1
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="MUC16"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="NOTCH4"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="POLE"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="POLD1"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="TP53"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="KRAS"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 1
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="EGFR"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 1
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="COL3A1"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,2]=="NRAS"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 1
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

DDR_gene <- c("MSH2","MSH6","PMS2","POLE","BRCA2")
mut_patient <- mutation[which(!is.na(match(mutation[,2],DDR_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

NER_gene <- c("ERCC3","ERCC2","ERCC5","ERCC4")
mut_patient <- mutation[which(!is.na(match(mutation[,2],NER_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

HR_gene <- c("RAD51","RAD51B","RAD51D","XRCC2","RAD52","RAD54L","RAD50","MRE11","NBN","ARID1A","BLM","BRCA2")
mut_patient <- mutation[which(!is.na(match(mutation[,2],HR_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

gene52 <- read.table("G:/all_work/mutation signature/immune/gene52.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
mut_patient <- mutation[which(!is.na(match(mutation[,2],gene52[,1]))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,1]))),42] <- 0
TP <- length(which(patient[,41]==1&patient[,42]==1))
FP <- length(which(patient[,41]==0&patient[,42]==1))
FN <- length(which(patient[,41]==1&patient[,42]==0))
TN <- length(which(patient[,41]==0&patient[,42]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy



##Hellmann 
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
mutation_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
mutation <- mutation_1[which(!is.na(match(mutation_1[,12],patient[,1]))),]

patient$response <- 0
patient[which(patient[,30]=="noresponse"),33] <- 1

mut_patient <- mutation[which(mutation[,1]=="B2M"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 1
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="JAK1"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 1
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="JAK2"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 1
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="PTEN"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 1
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="STK11"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 1
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="MUC16"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="NOTCH4"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="POLE"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="POLD1"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="TP53"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="KRAS"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="EGFR"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 1
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="COL3A1"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="NRAS"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 1
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

DDR_gene <- c("MSH2","MSH6","PMS2","POLE","BRCA2")
mut_patient <- mutation[which(!is.na(match(mutation[,1],DDR_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

NER_gene <- c("ERCC3","ERCC2","ERCC5","ERCC4")
mut_patient <- mutation[which(!is.na(match(mutation[,1],NER_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

HR_gene <- c("RAD51","RAD51B","RAD51D","XRCC2","RAD52","RAD54L","RAD50","MRE11","NBN","ARID1A","BLM","BRCA2")
mut_patient <- mutation[which(!is.na(match(mutation[,1],HR_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

gene52 <- read.table("G:/all_work/mutation signature/immune/gene52.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
mut_patient <- mutation[which(!is.na(match(mutation[,1],gene52[,1]))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),34] <- 0
TP <- length(which(patient[,33]==1&patient[,34]==1))
FP <- length(which(patient[,33]==0&patient[,34]==1))
FN <- length(which(patient[,33]==1&patient[,34]==0))
TN <- length(which(patient[,33]==0&patient[,34]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Rizvi
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
mutation_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
mutation <- mutation_1[which(!is.na(match(mutation_1[,12],patient[,1]))),]

patient$response <- 0
patient[which(patient[,15]=="NDB"),20] <- 1

mut_patient <- mutation[which(mutation[,1]=="B2M"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 1
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="JAK1"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 1
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="JAK2"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 1
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="PTEN"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 1
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="STK11"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 1
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="MUC16"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="NOTCH4"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="POLE"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="POLD1"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="TP53"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="KRAS"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="EGFR"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 1
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="COL3A1"),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

mut_patient <- mutation[which(mutation[,1]=="NRAS"),]
patient$mutation_type <- 0
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 1
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

DDR_gene <- c("MSH2","MSH6","PMS2","POLE","BRCA2")
mut_patient <- mutation[which(!is.na(match(mutation[,1],DDR_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

NER_gene <- c("ERCC3","ERCC2","ERCC5","ERCC4")
mut_patient <- mutation[which(!is.na(match(mutation[,1],NER_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

HR_gene <- c("RAD51","RAD51B","RAD51D","XRCC2","RAD52","RAD54L","RAD50","MRE11","NBN","ARID1A","BLM","BRCA2")
mut_patient <- mutation[which(!is.na(match(mutation[,2],HR_gene))),]
patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy

gene52 <- read.table("G:/all_work/mutation signature/immune/gene52.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
mut_patient <- mutation[which(!is.na(match(mutation[,1],gene52[,1]))),]

patient$mutation_type <- 1
patient[which(!is.na(match(patient[,1],mut_patient[,12]))),21] <- 0
TP <- length(which(patient[,20]==1&patient[,21]==1))
FP <- length(which(patient[,20]==0&patient[,21]==1))
FN <- length(which(patient[,20]==1&patient[,21]==0))
TN <- length(which(patient[,20]==0&patient[,21]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy


##Fig.4G
mutation_compare <- read.table("G:/all_work/mutation signature/immune/new_outputdata/mutation_compare.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
data <- mutation_compare
data$group <- factor(data$group,levels=c("Miao","Allen","Hellmann","Rizvi"))
library(tidyverse)
data = data %>% arrange(group, value)
# library
library(tidyverse)

empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))
 
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
 
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
 
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  geom_segment(data=grid_data, aes(x = end, y = 0.8, xend = start, yend = 0.8), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.6, xend = start, yend = 0.6), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.4, xend = start, yend = 0.4), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = end, y = -0.02, xend = start, yend = -0.02), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(data$id),4), y = c(0.2, 0.4, 0.6, 0.8), label = c("0.2", "0.4", "0.6", "0.8") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=1) +
  ylim(-0.5,1) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  scale_fill_manual(values=c("#C7DAB3","#FEB281","#88C8C9","#A0C3B9","#C9DC80"))+
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.1, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -0.1, label=group), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
 
p


##expression
aa <- readRDS("ICB_Van_Allen.rds")
bb <- aa@ExperimentList@listData$expr@assays@data@listData$expr
write.table(bb,"G:/all_work/mutation signature/immune/Levi_exp.txt",col=T,row=T,sep="\t",quote=F)


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

Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")

merge_matrix_1 <- merge_matrix[which(!is.na(match(rownames(merge_matrix),Melanoma_1[,1]))),]

###############################################################################################################################################
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,8]=="nonresponse")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix_1),Melanoma_1[,1]),]
merge_matrix_1 <- as.data.frame(merge_matrix_1)
merge_matrix_1$y <- each_drug_1$y
test <- merge_matrix_1

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
result$predcted_type <- 0
result[which(result[,1]>0.5),4] <- 1
TP <- length(which(result[,2]==1&result[,4]==1))
FP <- length(which(result[,2]==0&result[,4]==1))
FN <- length(which(result[,2]==1&result[,4]==0))
TN <- length(which(result[,2]==0&result[,4]==0))
accuracy <- (TP+TN)/(TP+TN+FP+FN)
accuracy
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1
Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]
Melanoma_1$express <- as.numeric(Levi_exp_1[which(rownames(Levi_exp_1)=="ENSG00000120217"),])
library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc

Melanoma_1$express <- as.numeric(Levi_exp_1[which(rownames(Levi_exp_1)=="ENSG00000163599"),])
library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc

Melanoma_1$express <- as.numeric(Levi_exp_1[which(rownames(Levi_exp_1)=="ENSG00000188389"),])
library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc


##CD8 T cell
CD8T <- read.table("G:/all_work/mutation signature/immune/signature/CD8_T_cell.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1

Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]
cd8texp <- Levi_exp_1[which(!is.na(match(rownames(Levi_exp_1),CD8T[,2]))),]
cd8texp_1 <- colMeans(cd8texp)
Melanoma_1$express <- as.numeric(cd8texp_1)

library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc

##CYT
CYT <- read.table("G:/all_work/mutation signature/immune/signature/CYT.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1

Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]
CYTtexp <- Levi_exp_1[which(!is.na(match(rownames(Levi_exp_1),CYT[,2]))),]
library(psych)
CYT_score <- unlist(lapply(CYTtexp,geometric.mean,2))
Melanoma_1$express <- as.numeric(CYT_score)
Melanoma_2 <- Melanoma_1[which(Melanoma_1$express!="NaN"),]

library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_2[,72],Melanoma_2[,73])
auc <- round(auc(rocobj),4)
auc

##CTL
CTL <- read.table("G:/all_work/mutation signature/immune/signature/CTL.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1

Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]
CTLtexp <- Levi_exp_1[which(!is.na(match(rownames(Levi_exp_1),CTL[,2]))),]
CTLtexp_1 <- colMeans(CTLtexp)
Melanoma_1$express <- as.numeric(CTLtexp_1)

library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc

##IFNR
IFN <- read.table("G:/all_work/mutation signature/immune/signature/IFN.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1

Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]
IFNtexp <- Levi_exp_1[which(!is.na(match(rownames(Levi_exp_1),IFN[,2]))),]
IFNtexp_1 <- colMeans(IFNtexp)
Melanoma_1$express <- as.numeric(IFNtexp_1)

library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc

##T cell inflamed
IFN <- read.table("G:/all_work/mutation signature/immune/signature/T_cell_inflamed.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1

Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]
IFNtexp <- Levi_exp_1[which(!is.na(match(rownames(Levi_exp_1),IFN[,2]))),]
IFNtexp_1 <- colMeans(IFNtexp)
Melanoma_1$express <- as.numeric(IFNtexp_1)

library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc



##LRRC15 CAF
IFN <- read.table("G:/all_work/mutation signature/immune/signature/LRRC15.CAF",header=F,sep="\t",quote="",stringsAsFactors=F)
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1

library(multiGSEA)
Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]
IFNtexp <- Levi_exp_1[which(!is.na(match(rownames(Levi_exp_1),IFN[,2]))),]
IFNtexp_1 <- eigenWeightedMean(IFNtexp)
Melanoma_1$express <- as.numeric(IFNtexp_1$score)

library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc

##NLRP3
IFN <- read.table("G:/all_work/mutation signature/immune/signature/NLRP3.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1
Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]

library(GSVA)
signature <- list(IFN[,2])
GSVA_res <-gsva(as.matrix(Levi_exp_1),
                   signature,
                  method ="ssgsea")

Melanoma_1$express <- as.numeric(GSVA_res)

library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc

##CD8
IFN <- read.table("G:/all_work/mutation signature/immune/signature/cd8.txt",header=F,sep="\t",quote="",stringsAsFactors=F)
Levi_exp <- read.table("G:/all_work/mutation signature/immune/Levi_exp_1.csv",header=T,sep=",",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/Levi_mutation_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_11 <- patient_1[which(patient_1[,8]!="long-survival"),]
patient_2 <- read.table("G:/all_work/mutation signature/immune/Levi_exp_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_21 <- patient_2[which(patient_2[,17]!="long-survival"),]
Melanoma_1 <- merge(patient_11,patient_21,by="patient")
Melanoma_1$y <- 0
Melanoma_1[which(Melanoma_1[,8]=="nonresponse"),72] <- 1
Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]

Levi_exp_1 <- Levi_exp[,match(Melanoma_1[,1],colnames(Levi_exp))]
IFNtexp <- Levi_exp_1[which(!is.na(match(rownames(Levi_exp_1),IFN[,2]))),]
IFNtexp_1 <- colSums(IFNtexp)
Melanoma_1$express <- as.numeric(IFNtexp_1)

library(pROC)
library(ggplot2)
rocobj <- roc(Melanoma_1[,72],Melanoma_1[,73])
auc <- round(auc(rocobj),4)
auc



##Fig.4H
mydata <- read.table("G:/all_work/mutation signature/immune/new_outputdata/expression_compare.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
rx<-10.5 * 35.43307
ry<-7.5 * 35.43307
rr<-5.00 * 35.43307
circos <- data.frame(paste("<circle cx=\"", rx, "\" cy=\"", ry, "\" r=\"", .2 * rr, "\" style=\"stroke:#b7b7b7; stroke-width:1; stroke-dasharray: 10 5; fill:none\"/>", sep = ""),
                     paste("<circle cx=\"", rx, "\" cy=\"", ry, "\" r=\"", .4 * rr, "\" style=\"stroke:#6c6b6b; stroke-width:2; fill:none\"/>", sep = ""),
                     paste("<circle cx=\"", rx, "\" cy=\"", ry, "\" r=\"", .6 * rr, "\" style=\"stroke:#b7b7b7; stroke-width:1; stroke-dasharray: 10 5; fill:none\"/>", sep = ""),
                     paste("<circle cx=\"", rx, "\" cy=\"", ry, "\" r=\"", .8*rr, "\" style=\"stroke:#6c6b6b; stroke-width:2; fill:none\"/>", sep = "")
)
angle_base <- 360/nrow(mydata)
for (i in 1: nrow(mydata)) {
  mydata[i, 6] <- angle_base * (i-1)
} 
names(mydata)[6] <- "angle"

cx1<-rx - 35/(nrow(mydata)/6) - 20 * mydata$AUC

cy1<-ry - mydata$AUC * 59 / 0.25
cx2<-rx + 35/(nrow(mydata)/6) + 20 * mydata$AUC
cy2<-cy1
mydata$bezier <- paste("<path d=\"M",rx,",",ry," C", cx1, ",", cy1, ",", cx2, ",", cy2, ",",rx,",",ry,"\" transform=\"rotate(", mydata$angle, " ",rx," ",ry,")\"", " fill=\"#", mydata$fill, "\" stroke = \"black\" stroke-width = \"0.5px\" />", sep = "")
tx<-rx - nchar(mydata$Name) * 3
ty<-ry - rr - 6
mydata$text <- paste("<text x=\"", tx, "\" y=\"", ty, "\" font-size=\"", mydata$text_size, "\" transform=\"rotate(", mydata$angle, " ",rx," ",ry,")\"", " fill=\"#", mydata$text_color, "\" >", mydata$Name, "</text>", sep = "")
first_line <- data.frame("<?xml version=\"1.0\" standalone=\"no\"?>",
                         "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"",
                         "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">",
                         "",
                         paste("<svg id=\"svg\" width=\"744.0945\" height=\"1052.362\">", "\t")
)
write.table(first_line[1, 1], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(first_line[1, 2], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 3], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 4], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(first_line[1, 5], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(circos[1, 1], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(circos[1, 2], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(circos[1, 3], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(circos[1, 4], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(mydata$bezier, "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
write.table(mydata$text, "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
last_line <- data.frame(paste("</svg>"))
write.table(last_line[1, 1], "G:/all_work/mutation signature/immune/new_graph/radar.svg", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
