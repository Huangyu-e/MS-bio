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
sel_gene <- as.data.frame(sel_gene)##880
write.table(sel_gene,"G:/all_work/mutation signature/immune/sel_gene.txt",col=F,row=F,sep="\t",quote=F)
##TCGA
sel_gene <- read.table("/data/home/huangyue/mutation_signature/immune/MC3/sel_gene.txt",header=FALSE,sep="\t",quote="",stringsAsFactors=F)
library(maftools)
mutation <- read.maf("/data/home/huangyue/mutation_signature/immune/MC3/mc3.v0.2.8.PUBLIC.maf.gz")
mutation_1 <- as.data.frame(mutation@data)
sample_type <- function(x){
    split_name <- unlist(strsplit(x,split="-"))
	type <- split_name[4]
	return(type)
}
mutation_1$type <- unlist(lapply(as.character(mutation_1[,16]),sample_type))
mutation_2 <- mutation_1[which(mutation_1$type=="01A"),]
sel_mutation <- mutation_2[which(!is.na(match(mutation_2[,1],sel_gene[,1]))),]

mutation_transform <- function(x){
    split_name <- unlist(strsplit(x,split="-"))
	merge_1 <- paste(split_name[1],split_name[2],sep="-")
	merge_2 <- paste(merge_1,split_name[3],sep="-")
	return(merge_2)
}
sel_mutation[,16] <- as.character(sel_mutation[,16])
sel_mutation$sample_1 <- unlist(lapply(sel_mutation[,16],mutation_transform))

Thorsson <- read.table("/data/home/huangyue/mutation_signature/immune/MC3/Thorsson.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
delete_sample <- Thorsson[which(Thorsson[,2]=="DLBC"|Thorsson[,2]=="LAML"|Thorsson[,2]=="THYM"),]##删掉这三个主要包含免疫细胞的癌型
sel_mutation_1 <- sel_mutation[which(is.na(match(sel_mutation$sample_1,delete_sample[,1]))),]
write.table(sel_mutation_1,"/data/home/huangyue/mutation_signature/immune/MC3/TCGA_sel_mutation.txt",col.names=T,row.names=F,quote=F,sep="\t")

library(sigminer)
library(maftools)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
laml <- read.maf(maf = "/data/home/huangyue/mutation_signature/immune/MC3/TCGA_sel_mutation.txt")
mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "SBS")
matrix_1 <- t(mt_tally_all$all_matrices$SBS_96)[,which(colSums(t(mt_tally_all$all_matrices$SBS_96))>0)]
write.table(matrix_1,"/data/home/huangyue/mutation_signature/immune/MC3/fisher_SBS_TCGA.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "DBS")
matrix_1 <- t(mt_tally_all$all_matrices$DBS_78)[,which(colSums(t(mt_tally_all$all_matrices$DBS_78))>0)]
write.table(matrix_1,"/data/home/huangyue/mutation_signature/immune/MC3/fisher_DBS_TCGA.txt",col.names=T,row.names=T,quote=F,sep="\t")

mt_tally_all <- sig_tally(laml,
                ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                mode = "ID")
matrix_1 <- t(mt_tally_all$all_matrices$ID_83)[,which(colSums(t(mt_tally_all$all_matrices$ID_83))>0)]
write.table(matrix_1,"/data/home/huangyue/mutation_signature/immune/MC3/fisher_ID_TCGA.txt",col.names=T,row.names=T,quote=F,sep="\t")

from SigProfilerAssignment import Analyzer as Analyze
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/MC3/fisher_SBS_TCGA.txt",
                    output="/data/home/huangyue/mutation_signature/immune/MC3/fisher_SBS_assign", 
					context_type="96",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/MC3/fisher_DBS_TCGA.txt",
                    output="/data/home/huangyue/mutation_signature/immune/MC3/fisher_DBS_assign", 
					context_type="DINUC",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)
Analyze.denovo_fit( samples="/data/home/huangyue/mutation_signature/immune/MC3/fisher_ID_TCGA.txt",
                    output="/data/home/huangyue/mutation_signature/immune/MC3/fisher_ID_assign", 
					context_type="ID",
                    signatures="/data/home/huangyue/mutation_signature/immune/Allen_mutation_Roh/fisher_ID_all_100/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt",
                    genome_build="GRCh37", 
                    verbose=True)


###############################################################################################################################################					
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation/fisher_DBS_all_100/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
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
patient <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- patient[which(patient[,17]!="long-term survival"),]
y <- rep(0,nrow(Melanoma_1))
y[which(Melanoma_1[,17]=="no clinical benefit")] <- 1
Melanoma_1$y <- y
each_drug_1 <- Melanoma_1[match(rownames(merge_matrix),Melanoma_1[,1]),]
merge_matrix <- as.data.frame(merge_matrix)
merge_matrix$y <- each_drug_1$y
train <- merge_matrix

SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/fisher_SBS_assign/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
DBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/fisher_DBS_assign/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
all_signature <- union(colnames(SBS_Exposure),colnames(DBS_Exposure))
all_sample <- union(rownames(SBS_Exposure),rownames(DBS_Exposure))
merge_matrix <- matrix(rep(0,length(all_signature)*length(all_sample)),ncol=length(all_signature))
colnames(merge_matrix) <- all_signature
rownames(merge_matrix) <- all_sample
SBS_Exposure <- as.matrix(SBS_Exposure)
DBS_Exposure <- as.matrix(DBS_Exposure)
merge_matrix[match(rownames(SBS_Exposure),rownames(merge_matrix)),match(colnames(SBS_Exposure),colnames(merge_matrix))] <- SBS_Exposure
merge_matrix[match(rownames(DBS_Exposure),rownames(merge_matrix)),match(colnames(DBS_Exposure),colnames(merge_matrix))] <- DBS_Exposure
test <- merge_matrix
###############################################################################################################################################

###############################################################################################################################################
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
pred1 <- predict(otu_train.forest, test,type="prob")
result <- data.frame(prediction=pred1[,2],sample=rownames(test))
write.table(result,"G:/all_work/mutation signature/immune/independent_data/MC3/result_RF.txt",col.names=T,row.names=F,quote=F,sep="\t")

##expression
merge_exp <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/TCGA_exp.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
exp_type <- function(x){
    split_name <- unlist(strsplit(x,split="\\."))
	return(split_name[5])
}
sample_type <- unlist(lapply(colnames(merge_exp),exp_type))
merge_exp_1 <- merge_exp[,which(sample_type=="01A")]
exp_transform <- function(x){
    split_name <- unlist(strsplit(x,split="\\."))
	merge_1 <- paste(split_name[2],split_name[3],sep=".")
	merge_2 <- paste(merge_1,split_name[4],sep=".")
	merge_3 <- paste(merge_2,split_name[5],sep=".")
	return(merge_3)
}
aa <- unlist(lapply(colnames(merge_exp_1),exp_transform))
colnames(merge_exp_1) <- aa


rownames(merge_exp_1)[which(grepl("SLC35E2",rownames(merge_exp_1)))] <- c("\"SLC35E2B|728661\"","\"SLC35E2A|9906\"")

gene_split <- function(x){
    split_name <- unlist(strsplit(x,split="\\|"))
	return(split_name[1])
}
bb <- unlist(lapply(rownames(merge_exp_1),gene_split))

gene_split_1 <- function(x){
    split_name <- unlist(strsplit(x,split="\""))
	return(split_name[2])
}
cc <- unlist(lapply(bb,gene_split_1))

rownames(merge_exp_1) <- cc
write.table(merge_exp_1,"G:/all_work/mutation signature/immune/independent_data/MC3/TCGA_exp_deal.txt",col.names=T,row.names=T,quote=F,sep="\t")


#######################################################################################################################################
result <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/result_RF.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
mutation_transform <- function(x){
    split_name <- unlist(strsplit(x,split="-"))
	merge_1 <- paste(split_name[1],split_name[2],sep="-")
	merge_2 <- paste(merge_1,split_name[3],sep="-")
	return(merge_2)
}
result$sample_1 <- unlist(lapply(result[,2],mutation_transform))
result$type <- "nonresponse"
result[which(result[,1]<=0.5),4] <- "response"

##Fig.5A
cbioportal <- read.table("G:/all_work/mutation signature/immune/cbioportal.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
result_merge <- merge(result,score_table,by.x="sample_1",by.y="TCGA.Participant.Barcode")
result_merge$TMB <- result_merge$Mutation.Count/38
Nonsilent <- result_merge[,c(4,21)]
Nonsilent_1 <- Nonsilent[which(!is.na(Nonsilent[,2])),]
Nonsilent <- Nonsilent_1
Nonsilent[,2] <- log10(Nonsilent[,2]+1)
Nonsilent$type <- factor(Nonsilent$type,levels=c("response","nonresponse"))
ggplot(Nonsilent, aes(x=type, y=TMB,fill=type)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+scale_fill_manual(values=c("#7FC17F","#F7B882"))+theme_classic()+ylab("log10 Nonsilent.Mutation.Rate")
compare_means(TMB ~ type,  data = Nonsilent)##****
##Fig.5B
MSI <- result_merge[,c(4,22)]
MSI <- MSI[which(!is.na(MSI[,2])),]
MSI[,2] <- log10(MSI[,2]+1)
MSI$type <- factor(MSI$type,levels=c("response","nonresponse"))
ggplot(Nonsilent, aes(x=type, y=MSI,fill=type)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+scale_fill_manual(values=c("#7FC17F","#F7B882"))+theme_classic()+ylab("log10 Nonsilent.Mutation.Rate")
compare_means(MSI ~ type,  data = Nonsilent)##****

##Fig.5C
score_table <- read.table("G:/all_work/mutation signature/immune/Thorsson.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
result_merge <- merge(result,score_table,by.x="sample_1",by.y="TCGA.Participant.Barcode")
SNV.Neoantigens_map <- result_merge[,c(4,18)]
SNV.Neoantigens_map_1 <- SNV.Neoantigens_map[which(!is.na(SNV.Neoantigens_map[,2])),]
SNV.Neoantigens_map <- SNV.Neoantigens_map_1
SNV.Neoantigens_map[,2] <- log10(SNV.Neoantigens_map[,2]+1)
SNV.Neoantigens_map$type <- factor(SNV.Neoantigens_map$type,levels=c("response","nonresponse"))
ggplot(SNV.Neoantigens_map, aes(x=type, y=SNV.Neoantigens,fill=type)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+scale_fill_manual(values=c("#7FC17F","#F7B882"))+theme_classic()+ylab("log10 SNV.Neoantigens")
compare_means(SNV.Neoantigens ~ type,  data = SNV.Neoantigens_map)##***

##Fig.5D
Indel.Neoantigens_map <- result_merge[,c(4,19)]
Indel.Neoantigens_map_1 <- Indel.Neoantigens_map[which(!is.na(Indel.Neoantigens_map[,2])),]
Indel.Neoantigens_map <- Indel.Neoantigens_map_1
Indel.Neoantigens_map[,2] <- log10(Indel.Neoantigens_map[,2]+1)
Indel.Neoantigens_map$type <- factor(Indel.Neoantigens_map$type,levels=c("response","nonresponse"))
ggplot(Indel.Neoantigens_map, aes(x=type, y=Indel.Neoantigens,fill=type)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+scale_fill_manual(values=c("#7FC17F","#F7B882"))+theme_classic()+ylab("log10 Indel.Neoantigens")	
compare_means(Indel.Neoantigens ~ type,  data = Indel.Neoantigens_map)##****

##Fig.5E
TCR_map <- result_merge[,c(4,30)]
TCR_map_1 <- TCR_map[which(!is.na(TCR_map[,2])),]
TCR_map <- TCR_map_1
TCR_map[,2] <- log10(TCR_map[,2]+1)
TCR_map$type <- factor(TCR_map$type,levels=c("response","nonresponse"))
ggplot(TCR_map, aes(x=type, y=TCR.Richness,fill=type)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+scale_fill_manual(values=c("#7FC17F","#F7B882"))+theme_classic()+ylab("log10 TCR")
compare_means(TCR.Richness ~ type,  data = TCR_map)##****

##Fig.5F
BCR_map <- result_merge[,c(4,28)]
BCR_map_1 <- BCR_map[which(!is.na(BCR_map[,2])),]
BCR_map <- BCR_map_1
BCR_map[,2] <- log10(BCR_map[,2]+1)
BCR_map$type <- factor(BCR_map$type,levels=c("response","nonresponse"))
ggplot(BCR_map, aes(x=type, y=BCR.Richness,fill=type)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+scale_fill_manual(values=c("#7FC17F","#F7B882"))+theme_classic()+ylab("log10 BCR")
compare_means(BCR.Richness ~ type,  data = BCR_map)##****

##Fig.5G
merge_exp_2 <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/TCGA_exp_deal.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
result <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/result_RF.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
##Estimate
library(IOBR)
cibersort <- deconvo_tme(eset = merge_exp_2, method = "estimate")
cibersort <- as.data.frame(cibersort)
write.table(cibersort,"G:/all_work/mutation signature/immune/independent_data/MC3/Estimate.txt",col.names=T,row.names=T,quote=F,sep="\t")

result <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/result_RF.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Estimate <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/Estimate.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
mutation_transform <- function(x){
    split_name <- unlist(strsplit(x,split="-"))
	merge_1 <- paste(split_name[1],split_name[2],sep="-")
	merge_2 <- paste(merge_1,split_name[3],sep="-")
	merge_3 <- paste(merge_2,split_name[4],sep="-")
	return(merge_3)
}
result$sample_1 <- unlist(lapply(result[,2],mutation_transform))
result$type <- "nonresponse"
result[which(result[,1]<=0.5),4] <- "response"
Estimate_1 <- merge(Estimate,result,by.x="ID",by.y="sample_1")
Estimate_1$type <- factor(Estimate_1$type,levels=c("response","nonresponse"))
library(ggpubr)
ggboxplot(
  Estimate_1, x = "type", y = "ESTIMATEScore_estimate",fill="type" ,palette = c("#7FC17F", "#F7B882"),width=0.5
  )+stat_compare_means()+ylab("ImmuneScore_estimate")
compare_means(ImmuneScore_estimate ~ type,data = Estimate_1)

##Fig.5H
Lymphocytes_map <- result_merge[,c(4,62)]
Lymphocytes_map_1 <- Lymphocytes_map[which(!is.na(Lymphocytes_map[,2])),]
Lymphocytes_map <- Lymphocytes_map_1
Lymphocytes_map$type <- factor(Lymphocytes_map$type,levels=c("response","nonresponse"))
ggboxplot(
  Lymphocytes_map, x = "type", y = "Lymphocytes",fill="type" ,palette = c("#7FC17F", "#F7B882"),width=0.5
  )+stat_compare_means()+ylab("Lymphocytes.Fraction")
compare_means(Lymphocytes ~ type,data = Lymphocytes_map)

##Fig.5I
Leukocyte_map <- result_merge[,c(4,8)]
Leukocyte_map_1 <- Leukocyte_map[which(!is.na(Leukocyte_map[,2])),]
Leukocyte_map <- Leukocyte_map_1
Leukocyte_map$type <- factor(Leukocyte_map$type,levels=c("response","nonresponse"))
library(ggpubr)
ggboxplot(
  Leukocyte_map, x = "type", y = "Leukocyte.Fraction",fill="type" ,palette = c("#7FC17F", "#F7B882"),width=0.5
  )+stat_compare_means()+ylab("Leukocyte.Fraction")
compare_means(Leukocyte.Fraction ~ type,data = Leukocyte_map)

##Fig.5J
TIL_map <- result_merge[,c(4,11)]
TIL_map_1 <- TIL_map[which(!is.na(TIL_map[,2])),]
TIL_map <- TIL_map_1
TIL_map$type <- factor(TIL_map$type,levels=c("response","nonresponse"))
ggboxplot(
  TIL_map, x = "type", y = "TIL.Regional.Fraction",fill="type" ,palette = c("#7FC17F", "#F7B882"),width=0.5
  )+stat_compare_means()+ylab("TILs")
compare_means(TIL.Regional.Fraction ~ type,data = TIL_map)

##Fig.5K
coinhibitor <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/coinhibitor.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
inhibitor <- coinhibitor[which(coinhibitor[,2]=="Coinhibitors"),]
stimulators <- coinhibitor[which(coinhibitor[,2]=="Costimulators"),]

merge_exp_2 <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/TCGA_exp_deal.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
result <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/result_RF.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
merge_exp_1 <- log2(merge_exp_2+1)

MHC <- inhibitor
MHC_exp <- c()
for(i in 1:nrow(MHC)){
    each_exp <- data.frame(expre=as.numeric(merge_exp_1[which(rownames(merge_exp_1)==MHC[i,1]),][1,]),sample_1=colnames(merge_exp_1),gene=MHC[i,1])
	MHC_exp <- rbind(MHC_exp,each_exp)
}
MHC_exp <- as.data.frame(MHC_exp)
MHC_exp[,1] <- as.numeric(MHC_exp[,1])
MHC_exp_1 <- MHC_exp[which(!is.na(MHC_exp[,1])),]

mutation_transform <- function(x){
    split_name <- unlist(strsplit(x,split="-"))
	merge_1 <- paste(split_name[1],split_name[2],sep=".")
	merge_2 <- paste(merge_1,split_name[3],sep=".")
	merge_3 <- paste(merge_2,split_name[4],sep=".")
	return(merge_3)
}
result$sample_1 <- unlist(lapply(result[,2],mutation_transform))
result$type <- "nonresponse"
result[which(result[,1]<=0.5),4] <- "response"
MHC_exp_2 <- merge(MHC_exp_1,result,by="sample_1")
MHC_exp_3 <- merge(MHC_exp_2,MHC,by.x="gene",by.y="symbol")
MHC_exp_3$type <- factor(MHC_exp_3$type,levels=c("response","nonresponse"))
MHC_exp_3 <- unique(MHC_exp_3)

MHC_exp_3$expre <- as.numeric(MHC_exp_3$expre)

library(ggunchained)
ggplot(MHC_exp_3, aes(x = gene, y = expre, fill = type)) + 
geom_split_violin(trim = F,colour=NA,scale="width")+ 
geom_point(stat = 'summary',fun=mean,
             position = position_dodge(width = 0.9))+
					xlab("") +theme_classic()+ylab("expression")+
  scale_fill_manual(values = c("#7FC17F","#F7B882"))+
  theme(axis.text.x = element_text(angle = 45, hjust = .5, vjust = 0.5)) + 
  
  guides(fill = guide_legend(title = NULL)) + 
  theme(legend.background = element_blank(),
        legend.position=c(0,0),legend.justification = c(0,0))+
		stat_compare_means(aes(group=type))
stat <- c()
for(i in unique(MHC_exp_3$gene)){
    each_map <- MHC_exp_3[which(MHC_exp_3[,1]==i),]
	each_stat <- c(i,compare_means(expre ~ type, data = each_map)$p)
	stat <- rbind(stat,each_stat)
}
stat <- as.data.frame(stat)
stat[,2] <- as.numeric(stat[,2])
stat[which(stat[,2]>0.05),]

##Fig.5L
merge_exp <- read.table("/data/home/huangyue/mutation_signature/immune/MC3/TCGA_exp_deal.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
result <- read.table("/data/home/huangyue/mutation_signature/immune/MC3/result_RF.txt",header=T,sep="\t",quote="",stringsAsFactors=F)

library(GSVA)
GSVA_res <-gsva(as.matrix(merge_exp),
                   signature,
                  method ="gsva")
write.table(GSVA_res,"/data/home/huangyue/mutation_signature/immune/MC3/Heetal_gsva_1.txt",col.names=T,row.names=T,quote=F,sep="\t")

import <- read.table("/data/home/huangyue/mutation_signature/immune/MC3/immport_pathway.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
immport_pathway <- list()
for(j in unique(import[,6])){
    immport_pathway[[j]] <- import[which(import[,6]==j),1]
}
GSVA_res <-gsva(as.matrix(merge_exp),
                   immport_pathway,
                  method ="gsva")
write.table(GSVA_res,"/data/home/huangyue/mutation_signature/immune/MC3/import_gsva_1.txt",col.names=T,row.names=T,quote=F,sep="\t")

##Fig.5M
library(GSVA)
pathway <- list.files("/data/home/huangyue/mutation_signature/immune/MC3/MsigDB")
signature <- list()
for(i in 1:length(pathway)){
    each_pathway_path <- paste("/data/home/huangyue/mutation_signature/immune/MC3/MsigDB",pathway[[i]],sep="/")
    each_pathway <- read.table(each_pathway_path,header=F,sep="\t",quote="",stringsAsFactors=F)
    signature[[pathway[[i]]]]= as.vector(as.character(each_pathway[1,-c(1,2)]))

}
GSVA_res <-gsva(as.matrix(merge_exp),
                   signature,
                  method ="gsva")
write.table(GSVA_res,"/data/home/huangyue/mutation_signature/immune/MC3/MsigDB_gsva_1.txt",col.names=T,row.names=T,quote=F,sep="\t")

GSVA_res <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/import_gsva_1.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
result <- read.table("G:/all_work/mutation signature/immune/independent_data/MC3/result_RF.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
mutation_transform <- function(x){
    split_name <- unlist(strsplit(x,split="-"))
	merge_1 <- paste(split_name[1],split_name[2],sep=".")
	merge_2 <- paste(merge_1,split_name[3],sep=".")
	merge_3 <- paste(merge_2,split_name[4],sep=".")
	return(merge_3)
}
result$sample_1 <- unlist(lapply(result$sample,mutation_transform))
result$type <- "nonresponse"
result[which(result[,1]<=0.5),4] <- "response"
response_sample <- result[which(result[,4]=="response"),3]
nonresponse_sample <- result[which(result[,4]=="nonresponse"),3]
mean_all <- c()
for(i in 1:nrow(GSVA_res)){
    each_pathway <- as.matrix(GSVA_res[i,])
	rownames(each_pathway) <- rownames(GSVA_res)[i]
	colnames(each_pathway) <- colnames(GSVA_res)
	mean_response <- mean(as.numeric(each_pathway[1,which(!is.na(match(colnames(each_pathway),response_sample)))]))
	mean_nonresponse <- mean(as.numeric(each_pathway[1,which(!is.na(match(colnames(each_pathway),nonresponse_sample)))]))
	each_mean <- c(rownames(each_pathway),mean_response,mean_nonresponse)
	mean_all <- rbind(mean_all,each_mean)
}
mean_all <- as.data.frame(mean_all)
colnames(mean_all) <- c("pathway","response","nonresponse")
mean_all[,2] <- as.numeric(mean_all[,2])
mean_all[,3] <- as.numeric(mean_all[,3])
ggplot(mean_all,aes(x=nonresponse,y=response))+geom_point(size=5,color="#447EAE")+geom_abline(intercept=0,slope=1,color="#DC2020",lty=5,lwd=1)+
geom_text(aes(nonresponse,response,label=pathway, vjust = 1))+
  theme_classic()+ylim(c(-0.2,0.3))+xlim(c(-0.2,0.3))#Immport
  
score_all <- c()
for(i in 1:ncol(GSVA_res)){
    each_score <- data.frame(signature = rownames(GSVA_res),score = as.numeric(GSVA_res[,i]),sample_1 = colnames(GSVA_res)[i])
	score_all <- rbind(score_all,each_score)
}

score_all <- as.data.frame(score_all)
score_all_1 <- merge(score_all,result,by="sample_1")
score_all_1$type <- "nonresponse"
score_all_1[which(score_all_1[,4]<=0.5),6] <- "response"
score_all_1$type <- factor(score_all_1$type,levels=c("response","nonresponse"))
score_all_2 <- score_all_1[which(!is.na(score_all_1[,3])),]
stat <- c()
for(i in unique(score_all_1$signature)){
    each_map <- score_all_1[which(score_all_1[,2]==i),]
	each_stat <- c(i,compare_means(score ~ type, data = each_map)$p)
	stat <- rbind(stat,each_stat)
}
stat <- as.data.frame(stat)
stat[,2] <- as.numeric(stat[,2])
stat[which(stat[,2]>0.05),]