# GEO data preprocess
#input  TAMs_GPL6947_GSE35449_series_matrix GPL6947-13512.txt
#output  CAF_exp.txt MDSC_exp.txt TAM_exp.txt
setwd("/mnt/md0/Work_xulei/data/GEO_series_matrix")
data<-read.table("TAMs_GPL6947_GSE35449_series_matrix.txt",header=T,sep="\t");#GEO series-matix data/ have been remove the annotation information
data<-na.omit(data)
id<-read.table("GPL6947-13512.txt",header=T,sep="\t");#GPL data/ incloding probe_ID(col 1),Gene_name(col 2),Gene_id(col 3)
ind<-match(data[,1],id[,1])
length(which(is.na(ind)))
exp<-cbind(Name=id[ind,2],data[,-1])
exp1<-na.omit(exp)
drop<-grep(pattern='///',exp1[,1])
exp2<-exp1[-drop,]
exp3<-as.numeric(as.character(unlist(exp2)))
exp3<-matrix(exp3,ncol=ncol(exp2))
geneidfactor<-factor(exp3[,1])
gene_exp_matrix<-apply(exp3,2,function(x) tapply(x,geneidfactor,mean))#对相同因子取jun值;
rownames(gene_exp_matrix)<-unique(exp2[,1])
#-------CAF #dim 19922 24
CAF_exp<-gene_exp_matrix[,-1]
colname<-read.csv("GSE39396_colnames.csv",sep="",header=F)
colnames(CAF_exp)=paste(colname[,1],colname[,2],colname[,3])
#-----MDSC  #dim 19044 9 
MDSC_exp<-gene_exp_matrix[,-1]
colnames(MDSC_exp)=c(paste("A375-MDSC",1:3,sep=""),paste("MDSC+4-IPP",1:3,sep=""),paste("Cultured Monocytes",1:3,sep=""))
#-----TAMs  #dim 25158 21
TAM_exp<-gene_exp_matrix[,-1]
colnames(TAM_exp)<-c("M0_1","M2_1" ,"M1_1" ,"M1_2", "M2_2" ,"M0_3" ,"M2_3" ,"M1_3" ,"M0_4" ,"M2_4" ,"M1_4" ,"M0_5" ,"M2_5"  ,"M1_5" ,"M1_6" ,"M0_6" ,"M2_7" ,"M0_7","M0_2","M2_6","M1_7")
#------intersect 15791 
A<-intersect(rownames(CAF_exp),rownames(MDSC_exp))#17758
B<-intersect(A,rownames(TAM_exp))#15791
CAF<-CAF_exp[match(B,rownames(CAF_exp)),]
MDSC<-MDSC_exp[match(B,rownames(MDSC_exp)),]
TAM<-TAM_exp[match(B,rownames(TAM_exp)),]
#------cbind 3 cell_exp 
CAF1<-CAF[order(rownames(CAF)),]	
MDSC1<-MDSC[order(rownames(MDSC)),]	
TAM1<-TAM[order(rownames(TAM)),]	
Cell_exp<-cbind(CAF1,MDSC1,TAM1)	
write.table(TAM1,"TAM_exp.txt",col.names=T,row.names=T,sep="\t")
write.table(CAF1,"CAF_exp.txt",col.names=T,row.names=T,sep="\t")
write.table(MDSC1,"MDSC_exp.txt",col.names=T,row.names=T,sep="\t")
write.table(Cell_exp,"Cell_exp.txt",col.names=T,row.names=T,sep="\t")
