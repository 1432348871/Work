# GEO data preprocess
#data prepare
data<-read.table("/mnt/md0/Work_xulei/data/GEO_series_matrix/TAMs_GPL6947_GSE35449_series_matrix.txt",header=T,sep="\t");
data<-na.omit(data);
id<-read.table("/mnt/md0/Work_xulei/data/GEO_series_matrix/GPL6947-13512.txt",header=T,sep="\t");
ind<-match(data[,1],id[,1]);
length(which(is.na(ind)));
exp<-cbind(Name=id[ind,2],data[,-1]);
exp1<-na.omit(exp)
drop<-grep(pattern='///',exp1[,1])
exp2<-exp1[-drop,]
exp3=as.numeric(as.character(unlist(exp2)));
exp3=matrix(exp3,ncol=ncol(exp2));
geneidfactor=factor(exp3[,1]);
gene_exp_matrix=apply(exp3,2,function(x) tapply(x,geneidfactor,mean))#对相同因子取jun值;
rownames(gene_exp_matrix)<-unique(exp2[,1])
