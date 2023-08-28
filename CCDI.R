#################pyroptosis#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/pyroptosis.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=pyroptosis
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_pyroptosis.Rdata")
###
rm(list = ls())
load("1sur/surv_pyroptosis.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_pyroptosis.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_pyroptosis.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_pyroptosis.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_pyroptosis_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_pyroptosis_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_pyroptosis.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_pyroptosis.txt",
            sep="\t",
            quote=F,
            row.names=F)



rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_pyroptosis.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:6) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/pyroptosis_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/pyroptosis_level1.pdf",width = 3,height = 3)
p7
dev.off()  
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_pyroptosis.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/multiCox_pyroptosis.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()
#################ferroptosis#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/ferroptosis.txt',sep="\t",header=T,check.names=F)#283
names(rt)="symbol"
ferroptosis =merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=ferroptosis 
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_ferroptosis.Rdata")
###
rm(list = ls())
load("1sur/surv_ferroptosis.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_ferroptosis.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_ferroptosis.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_ferroptosis.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_ferroptosis_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_ferroptosis_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_ferroptosis.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_ferroptosis.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_ferroptosis.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:6) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/ferroptosis_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/ferroptosis_level1.pdf",width = 3,height = 3)
p7
dev.off()  
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_ferroptosis.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/multiCox_ferroptosis.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()
#################autophagy#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/autophagy.txt',sep="\t",header=T,check.names=F)#222
names(rt)="symbol"
autophagy=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=autophagy
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_autophagy.Rdata")
###
rm(list = ls())
load("1sur/surv_autophagy.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_autophagy.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_autophagy.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_autophagy.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_autophagy_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_autophagy_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_autophagy.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_autophagy.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_autophagy.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:3) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/autophagy_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/autophagy_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_autophagy.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_autophagy.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()
#################cuproptosis#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/cuproptosis.txt',sep="\t",header=T,check.names=F)#17
names(rt)="symbol"
cuproptosis=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=cuproptosis
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_cuproptosis.Rdata")
###
rm(list = ls())
load("1sur/surv_cuproptosis.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_cuproptosis.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_cuproptosis.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.5#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#1
write.table(outTab,file="2unicox/unicox_cuproptosis.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_cuproptosis_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_cuproptosis_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
#multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
#outTab=data.frame(outTab)
write.table(outTab,file="3multicox/multiCox_cuproptosis.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_cuproptosis.txt",
            sep="\t",
            quote=F,
            row.names=F)

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_cuproptosis.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
library(tidyverse)
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:2) %>% 
  ##
  dplyr::select(gene,value,everything()) 
library(ggridges)
pdf("3multicox/cuproptosis_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/cuproptosis_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_cuproptosis.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))



pdf(file="3multicox/mmultiCox_cuproptosis.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()
#################parthanatos#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/parthanatos.txt',sep="\t",header=T,check.names=F)#23
names(rt)="symbol"
parthanatos=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=parthanatos
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_parthanatos.Rdata")
###
rm(list = ls())
load("1sur/surv_parthanatos.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_parthanatos.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_parthanatos.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_parthanatos.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_parthanatos_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_parthanatos_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=colnames(rt[3]),outTab)
# outTab=gsub("`","",outTab)
# dir.create("3multicox")

write.table(outTab,file="3multicox/multiCox_parthanatos.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_parthanatos.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_parthanatos.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/parthanatos_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/parthanatos_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_parthanatos.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_parthanatos.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################Entotic_cell_death#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/Entotic cell death .txt',sep="\t",header=T,check.names=F)#23
names(rt)="symbol"
Entotic_cell_death =merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=Entotic_cell_death 
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_Entotic_cell_death .Rdata")
###
rm(list = ls())
load("1sur/surv_Entotic_cell_death .Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_Entotic_cell_death .txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_Entotic_cell_death .txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_Entotic_cell_death .txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_Entotic_cell_death _exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_Entotic_cell_death _exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_Entotic_cell_death .txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_Entotic_cell_death .txt",
            sep="\t",
            quote=F,
            row.names=F)




rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_Entotic_cell_death .txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:2) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/Entotic_cell_death_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/Entotic_cell_death_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_Entotic_cell_death .txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_Entotic_cell_death.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################necrosis#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/necrosis.txt',sep="\t",header=T,check.names=F)#498
names(rt)="symbol"
necrosis=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=necrosis
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_necrosis.Rdata")
###
rm(list = ls())
load("1sur/surv_necrosis.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_necrosis.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_necrosis.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_necrosis.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_necrosis_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_necrosis_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_necrosis.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_necrosis.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_necrosis.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:6) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/necrosis_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/necrosis_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_necrosis.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_necrosis.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################necroptosis#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/necroptosis.txt',sep="\t",header=T,check.names=F)#500
names(rt)="symbol"
necroptosis=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=necroptosis
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_necroptosis.Rdata")
###
rm(list = ls())
load("1sur/surv_necroptosis.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_necroptosis.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_necroptosis.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_necroptosis.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_necroptosis_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_necroptosis_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_necroptosis.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_necroptosis.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_necroptosis.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:9) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/necroptosis_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/necroptosis_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_necroptosis.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_necroptosis.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################lysosome_dependent_cell_death#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/lysosome dependent cell death.txt',sep="\t",header=T,check.names=F)#194
names(rt)="symbol"
lysosome_dependent_cell_death=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=lysosome_dependent_cell_death
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_lysosome_dependent_cell_death.Rdata")
###
rm(list = ls())
load("1sur/surv_lysosome_dependent_cell_death.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_lysosome_dependent_cell_death.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_lysosome_dependent_cell_death.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_lysosome_dependent_cell_death.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_lysosome_dependent_cell_death_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_lysosome_dependent_cell_death_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_lysosome_dependent_cell_death.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_lysosome_dependent_cell_death.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_lysosome_dependent_cell_death.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:4) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/lysosome_dependent_cell_death_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/lysosome_dependent_cell_death_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_lysosome_dependent_cell_death.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_lysosome_dependent_cell_death.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################intrinsic_apoptosis#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/intrinsic apoptosis.txt',sep="\t",header=T,check.names=F)#500
names(rt)="symbol"
intrinsic_apoptosis=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=intrinsic_apoptosis
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_intrinsic_apoptosis.Rdata")
###
rm(list = ls())
load("1sur/surv_intrinsic_apoptosis.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_intrinsic_apoptosis.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_intrinsic_apoptosis.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_intrinsic_apoptosis.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_intrinsic_apoptosis_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_intrinsic_apoptosis_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_intrinsic_apoptosis.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_intrinsic_apoptosis.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_intrinsic_apoptosis.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:6) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/intrinsic_apoptosis_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/intrinsic_apoptosis_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_intrinsic_apoptosis.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_intrinsic_apoptosis.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################immunogenic_cell_death#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/immunogenic cell death.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
immunogenic_cell_death=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=immunogenic_cell_death
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_immunogenic_cell_death.Rdata")
###
rm(list = ls())
load("1sur/surv_immunogenic_cell_death.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_immunogenic_cell_death.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_immunogenic_cell_death.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_immunogenic_cell_death.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_immunogenic_cell_death_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_immunogenic_cell_death_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)


multiCox=step(multiCox,direction = "both")
####评估cox基本假设###
library("survival")
library("survminer")
test.ph <- cox.zph(multiCox)
pdf("2unicox/cox.pdf",15,10)
ggcoxzph(test.ph)
dev.off()
#####
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_immunogenic_cell_death.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_immunogenic_cell_death.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_immunogenic_cell_death.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:8) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/immunogenic_cell_death_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/immunogenic_cell_death_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_immunogenic_cell_death.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_immunogenic_cell_death.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################extrinsic_apoptosis#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/extrinsic apoptosis.txt',sep="\t",header=T,check.names=F)#500
names(rt)="symbol"
extrinsic_apoptosis=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=extrinsic_apoptosis
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_extrinsic_apoptosis.Rdata")
###
rm(list = ls())
load("1sur/surv_extrinsic_apoptosis.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_extrinsic_apoptosis.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_extrinsic_apoptosis.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_extrinsic_apoptosis.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_extrinsic_apoptosis_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_extrinsic_apoptosis_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_extrinsic_apoptosis.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_extrinsic_apoptosis.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_extrinsic_apoptosis.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:9) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/extrinsic_apoptosis_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/extrinsic_apoptosis_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_extrinsic_apoptosis.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_extrinsic_apoptosis.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################Apoptosis_like_morphology#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/Apoptosis like morphology.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
Apoptosis_like_morphology=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=Apoptosis_like_morphology
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_Apoptosis_like_morphology.Rdata")
###
rm(list = ls())
load("1sur/surv_Apoptosis_like_morphology.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_Apoptosis_like_morphology.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_Apoptosis_like_morphology.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_Apoptosis_like_morphology.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_Apoptosis_like_morphology_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_Apoptosis_like_morphology_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_Apoptosis_like_morphology.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_Apoptosis_like_morphology.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_Apoptosis_like_morphology.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:4) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/Apoptosis_like_morphology_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/Apoptosis_like_morphology_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_Apoptosis_like_morphology.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_Apoptosis_like_morphology.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

#################necrosis_like_morphology#################
rm(list = ls())
load("data/GSE53625_symbol.Rdata")
rt=read.table('phenotype/necrosis like morphology.txt',sep="\t",header=T,check.names=F)#73
names(rt)="symbol"
necrosis_like_morphology=merge(rt,data1,by="symbol")[,-c(2,3)]
###
library(tidyverse)
clin_data <- read.table('data/time.txt',sep="\t",header=T,check.names=F) %>%
  dplyr::select(id,
                fustat,
                futime
  ) %>%
  mutate(id = id,
         futime =futime/12,
         fustat = fustat
  ) %>%
  dplyr::select(1,3,2)
matrix_miRNA=necrosis_like_morphology
matrix_sig_miRNA <- matrix_miRNA %>%  
  column_to_rownames("symbol") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-"),##变“.”为-
         id = str_sub(sample_id, 1, 12)) %>%
  relocate(id) %>%
  dplyr::select(-2) %>%
  as.data.frame()

surv_immuneRNA <- inner_join(clin_data, matrix_sig_miRNA, by = "id") %>%
  as.data.frame()

dir.create("1sur")
table(duplicated(surv_immuneRNA$id))
#FALSE  
#179
surv_immuneRNA=surv_immuneRNA[!duplicated(surv_immuneRNA$id),]
save(surv_immuneRNA,file = "1sur/surv_necrosis_like_morphology.Rdata")
###
rm(list = ls())
load("1sur/surv_necrosis_like_morphology.Rdata")
rt=surv_immuneRNA
library(survival)
library(survminer)
write.table(rt,file="1sur/surv_necrosis_like_morphology.txt",sep="\t",row.names=F,quote=F)
rt=read.table('1sur/surv_necrosis_like_morphology.txt',sep="\t",header=T,check.names=F,row.names = 1)
pFilter=0.05#过滤条件

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    # diff=survdiff(Surv(futime, fustat) ~rt[,gene],data = rt)
    # pValue=1-pchisq(diff$chisq,df=1)
    # if(pValue<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       # KMPvalue=pValue,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       coxPvalue=coxP) )
    # }
  }
}

dir.create("2unicox")
#13
write.table(outTab,file="2unicox/unicox_necrosis_like_morphology.txt",sep="\t",row.names=F,quote=F)  
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="2unicox/unicox_necrosis_like_morphology_exp.txt",sep="\t",row.names=F,quote=F)
###
rm(list = ls())
outTab=read.table("2unicox/unicox_necrosis_like_morphology_exp.txt",header=T,sep="\t",row.names=1,check.names=F)      #读取文件
rt=outTab
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
dir.create("3multicox")
write.table(outTab,file="3multicox/multiCox_necrosis_like_morphology.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
id=rownames(rt)
write.table(cbind(id,rt[,outCol],riskScore),
            file="3multicox/risk_necrosis_like_morphology.txt",
            sep="\t",
            quote=F,
            row.names=F)
rm(list = ls())

load("data/GSE53625_symbol.Rdata")
rt=read.table('3multicox/multiCox_necrosis_like_morphology.txt',sep="\t",header=T,check.names=F)#388
names(rt)="symbol"
pyroptosis=merge(rt,data1,by="symbol")[-c(2:8)]
data=pyroptosis%>%
  column_to_rownames("symbol")%>%
  t()%>%
  data.frame()
data$group=rep(c("tumor","normal"),179)

mydata<-data %>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="value",1:2) %>% 
  ##
  dplyr::select(gene,value,everything()) 

pdf("3multicox/necrosis_like_morphology_level2.pdf",width = 3,height = 3)

ggplot(mydata,aes(y=gene,x=value,fill=group))+
  geom_density_ridges() +
  scale_fill_manual(values = c("#0099B499","#ED000099"))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off() 
P5 <- mydata %>%
  ## 确定x,y
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7,outlier.size = 0) +
  scale_alpha_manual(values = c( "#00468B99","#ED000099"))+
  scale_y_continuous(name = " Exprssion")+
  scale_x_discrete(name = "gene") +
  # ggtitle("Boxplot of MMR") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 90, hjust = 1))

library(ggpubr)
p7 = P5 + stat_compare_means(method = "wilcox.test", label = "p.signif" , label.y = 13.5, size = 4)
pdf("3multicox/necrosis_like_morphology_level1.pdf",width = 3,height = 3)
p7
dev.off() 
##############unicox图#############
rm(list = ls())
rt <- read.table("3multicox/multiCox_necrosis_like_morphology.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))



pdf(file="3multicox/mmultiCox_necrosis_like_morphology.pdf", width = 5,height =3)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))


##绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(3.5,n+1,cex=text.cex,font=2,adj=1)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00468B99",lwd=0.5)
abline(v=1,col="#1B191999",lty=2,lwd=0.5)
boxcolor = ifelse(as.numeric(hr) > 1, "#ED000099", "#00468B99")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.5)
axis(1)
dev.off()

