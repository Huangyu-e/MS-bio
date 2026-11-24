##Fig.6A
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
library(randomForest)
set.seed(123)
train$y <- factor(train$y)
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
varImpPlot(otu_train.forest)

##Fig.6B
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_2 <- Melanoma_1[which(!is.na(match(Melanoma_1[,1],rownames(SBS_Exposure)))),]
resistant_patient <- Melanoma_2[which(Melanoma_2[,19]=="no clinical benefit"),1]
sensitive_patient <- Melanoma_2[which(Melanoma_2[,19]=="clinical benefit"),1]
SBS_Exposure$type <- "nonresponse"
SBS_Exposure[which(!is.na(match(rownames(SBS_Exposure),sensitive_patient))),5] <- "response"
SBS_Exposure$SBS96C <- log10(SBS_Exposure$SBS96C+1)
library(ggpubr)
ggboxplot(
  SBS_Exposure, x = "type", y = "SBS96C",color="type" ,palette = c("#274FA2","#DD6D6D"),add = "jitter",
  width=0.4
  )+stat_compare_means()
  
##Fig.6C
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_1_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,3]=="10.1016/j.ccell.2018.03.018"),]
Melanoma_2 <- patient[which(!is.na(match(patient[,1],rownames(SBS_Exposure)))),]
Melanoma_3 <- Melanoma_2[match(rownames(SBS_Exposure),Melanoma_2[,1]),]
cor.test(SBS_Exposure[,3],Melanoma_3$neoantigen_load,method="spearman")

data <- data.frame(SBS96C=log10(SBS_Exposure$SBS96C+1),neoantigen_load=log10(Melanoma_3$neoantigen_load+1))
library(ggplot2)
library(hrbrthemes)
ggplot(data, aes(x=SBS96C, y=neoantigen_load)) +
  geom_point() +
  geom_smooth(method=lm , color="black", fill="#69b3a2", se=TRUE) +
  theme_classic()
cor.test(data[,1],data[,2],method="spearman")

##Fig.6D
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_Roh/fisher_SBS_assign_all_Allen/De_Novo_Solution/Activities/De_Novo_Activities.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
patient_1 <- read.table("G:/all_work/mutation signature/immune/independent_data/lung_2_noSD/patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient <- patient_1[which(patient_1[,15]!="NR"),]
Melanoma_2 <- patient[which(!is.na(match(patient[,1],rownames(SBS_Exposure)))),]
Melanoma_3 <- Melanoma_2[match(rownames(SBS_Exposure),Melanoma_2[,1]),]
cor.test(SBS_Exposure[,3],Melanoma_3$Neoant.,method="spearman")
data <- data.frame(SBS96C=log10(SBS_Exposure$SBS96C+1),neoantigen_load=log10(Melanoma_3$Neoant.+1))
library(ggplot2)
library(hrbrthemes)
ggplot(data, aes(x=SBS96C, y=neoantigen_load)) +
  geom_point() +
  geom_smooth(method=lm , color="black", fill="#69b3a2", se=TRUE) +
  theme_classic()
cor.test(data[,1],data[,2],method="spearman")

  
##Fig.6E
SBS_Exposure <- read.table("G:/all_work/mutation signature/immune/Allen_mutation_Roh/fisher_SBS_all_100/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
#library(maftools)
#laml <- read.maf("G:/all_work/mutation signature/immune/Allen_mutation.txt")
#Allen_mutation_1 <- as.data.frame(laml@data)
Allen_mutation_1 <- read.table("G:/all_work/mutation signature/immune/Allen_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
Melanoma_1 <- read.table("G:/all_work/mutation signature/immune/Allen_patient.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
#Melanoma_1 <- Melanoma[which(Melanoma[,18]!="stable disease"),]
resistant_patient <- Melanoma_1[which(Melanoma_1[,19]=="no clinical benefit"),1]
sensitive_patient <- Melanoma_1[which(Melanoma_1[,19]=="clinical benefit"),1]
Allen_mutation <- Allen_mutation_1[which(!is.na(match(Allen_mutation_1[,24],Melanoma_1[,1]))),]

fisher_gene <- c()
for(gene in unique(Allen_mutation[,1])){
    each_gene <- Allen_mutation[which(Allen_mutation[,1]==gene),]
	mut_sample <- unique(each_gene[,24])
	nonmutat_sample <- setdiff(Allen_mutation[,24],mut_sample)
	mut_96C <- SBS_Exposure[which(!is.na(match(rownames(SBS_Exposure),mut_sample))),3]
	if(length(mut_96C)>2){
	nonmut_96C <- SBS_Exposure[which(!is.na(match(rownames(SBS_Exposure),nonmutat_sample))),3]
    p <- wilcox.test(mut_96C,nonmut_96C)$p.value
	if(p < 0.01){
	    each_gene <- c(gene,p)
		fisher_gene <- rbind(fisher_gene,each_gene)
	}
	}
}
fisher_gene <- as.data.frame(fisher_gene)
colnames(fisher_gene) <- c("gene","wilcox")
write.table(fisher_gene,"G:/all_work/mutation signature/immune/Allen_mutation_Roh/wilcox_gene_614.txt",col.names=T,row.names=F,quote=F,sep="\t")

fisher_gene <- read.table("/data/home/huangyue/mutation_signature/wilcox_gene_614.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
library(clusterProfiler)
eg <- bitr(fisher_gene[,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
KEGG_pathway <- enrichKEGG(eg[,2], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,pAdjustMethod = 'none', 
                                minGSSize = 1,maxGSSize = 500,qvalueCutoff = 1,use_internal_data = FALSE)
pdf("/data/home/huangyue/mutation_signature/kegg.pdf")
p <- barplot(KEGG_pathway,showCategory = 10,color = "pvalue")
print(p)
dev.off()
#hsa_kegg <- clusterProfiler::download_KEGG("hsa")##353
#hsa_GO <- select(org.Hs.eg.db,keys(org.Hs.eg.db),columns = c('GO'))#18349
#GO_BP <- hsa_GO[which(hsa_GO[,4]=="BP"),]#12342

ego_BP <- enrichGO(
  gene          = eg[,2],
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff = 1,
  minGSSize = 1,
  readable=TRUE)
ego_BP <- as.data.frame(ego_BP)
ego_BP$log10p <- (-log10(ego_BP$pvalue))
ego_BP_1 <- ego_BP[1:10,]
ego_BP_1$Description <- factor(ego_BP_1$Description,levels=c("peptide cross-linking","positive regulation of peptide secretion",
"microtubule-based movement","action potential","multicellular organismal signaling","regulation of receptor binding",
"regulation of T cell differentiation in thymus","T cell differentiation in thymus","positive regulation of synaptic vesicle exocytosis",
"membrane depolarization during action potential"))
pdf("/data/home/huangyue/mutation_signature/GOBP_1.pdf")
p <- ggplot(data=ego_BP_1, aes(x=Description, y=log10p)) +
geom_bar(stat="identity", width=0.8,fill="#2EAADF") + coord_flip() + 
theme_bw() + 
xlab("GO term") + 
theme(axis.text=element_text(face = "bold", color="gray50")) +
labs(title = "The Most Enriched GO Terms")
print(p)
dev.off()

  
  
  
  
  
  