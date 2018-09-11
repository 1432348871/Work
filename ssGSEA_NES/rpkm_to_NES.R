#input: bmsbms038_rpkm_normalized.txt  gene_set_GMX.txt
#output: Melanoma_TILs_tpm_01.txt
#-----data processing------#
rpkm<-as.matrix(read.table("bms038_rpkm_normalized.txt",sep="\t",header=T))
rpkm1=apply(rpkm[,-1],2,as.numeric)
rpkm_to_tpm<-function(rpkm) TPM=rpkm*10^6/sum(rpkm)#rpkm convert tpm according to 生信菜鸟团biotrainee.com
logTPM=log(rpkm_to_tpm(rpkm1)+1,2)
mean_TPM=apply(logTPM,1,mean)
loa=which(mean_TPM<0.1) 
rownames(logTPM)=rpkm[,1]
data<-logTPM[-loa,]
#--------ssGSEA NES--------#
#packages
source("https://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("GSVA")
biocLite("scales")
biocLite("XML")
library("scales")
library("XML")
library("GSVA")
#GeneSet
gene_set=as.matrix(read.table("gene_set_GMX.txt",sep="\t",header=TRUE))
GeneSet<-list(Activated_B_cell=gene_set[ ,1],
Activated_CD4_T_cell=gene_set[ ,2],
Activated_CD8_T_cell=gene_set[ ,3],
Activated_dendritic_cell=gene_set[ ,4],
CD56bright_natural_killer_cell=gene_set[ ,5],
CD56dim_natural_killer_cell=gene_set[ ,6],
Central_memory_CD4_T_cell=gene_set[ ,7],
Central_memory_CD8_T_cell=gene_set[ ,8],
Effector_memeory_CD4_T_cell=gene_set[ ,9],
Effector_memeory_CD8_T_cell=gene_set[ ,10],
Eosinophil=gene_set[ ,11],
Gamma_delta_T_cell=gene_set[ ,12],
Immature_B_cell=gene_set[ ,13],
Immature_dendritic_cell=gene_set[ ,14],
Macrophage=gene_set[ ,15],
Mast_cell=gene_set[ ,16],
MDSC=gene_set[ ,17],
Memory_B_cell=gene_set[ ,18],
Monocyte=gene_set[ ,19],
Natural_killer_cell=gene_set[ ,20],
Natural_killer_T_cell=gene_set[ ,21],
Neutrophil=gene_set[ ,22],
Plasmacytoid_dendritic_cell=gene_set[ ,23],
Regulatory_T_cell=gene_set[ ,24],
T_follicular_helper_cell=gene_set[ ,25],
Type_1_T_helper_cell=gene_set[ ,26],
Type_17_T_helper_cell=gene_set[ ,27],
Type_2_T_helper_cell=gene_set[ ,28]
)
#GSVA ssgsea
NES=gsva(data,GeneSet, method="ssgsea",mx.diff=TRUE,kcdf="Gaussian",min.sz=1,max.sz=Inf,tau=0.25,ssgsea.norm=TRUE)#
write.table(NES,"Melanoma_TILs_tpm_01.txt",sep="\t",col.names=T,row.names=T,quote=F)













