setwd("~")
population=read.table("~/population1.txt", header=TRUE)

args <- commandArgs(trailingOnly = TRUE)
name<- args[]

data=read.table(paste0(name,".vcf.ped"), colClasses = "character")

Maternal=data[ , !c(TRUE,FALSE) ]
Paternal=data[ , !c(FALSE,TRUE) ]

Maternalp=cbind(population,Maternal)
Paternalp=cbind(population,Paternal)
vcf=read.table(paste0(name,".vcf"), colClasses = "character")
vcf.ref=t(vcf[,2:4])
vcf.ref=apply(format(vcf.ref), 2, paste, collapse="_")
colnames(Maternalp) <- c("pop", "superpop","sampleId",vcf.ref)
colnames(Paternalp) <- c("pop", "superpop","sampleId",vcf.ref)

Maternalp<- Maternalp[ , !grepl( "esv" , names( Maternalp ) ) ]
Maternalp.dd=Maternalp[,-c(1:3)]
dd <- apply( Maternalp.dd ,1, paste , collapse = "" )
dd<-as.data.frame(dd)
matern=cbind(Maternalp.dd=Maternalp,dd)

Paternalp<- Paternalp[ , !grepl( "esv" , names( Paternalp ) ) ]
Paternalp.dd=Paternalp[,-c(1:3)]
dd <- apply( Paternalp.dd ,1, paste , collapse = "" )
dd<-as.data.frame(dd)
patern=cbind(Paternalp.dd=Paternalp,dd)
vcf.ref=t(vcf[,2:4])
vcf.ref<- vcf.ref[ , !grepl( "esv" , vcf.ref[2,]) ]
vcf.ref=apply(format(vcf.ref), 2, paste, collapse="_")
colnames(matern) <- c("pop", "superpop","sampleId",vcf.ref,"dd")
colnames(patern) <- c("pop", "superpop","sampleId",vcf.ref,"dd")

combined=rbind(matern,patern)

freq.all <- data.frame(table(combined$dd)) ##Retreive haplotype numbers for 2504 individuals all together
freq.pop <- as.data.frame.matrix(table(combined$dd,combined$pop)) ##Retreive haplotype numbers for populations
freq.suppop <- as.data.frame.matrix(table(combined$dd,combined$superpop)) ##Retreive haplotype numbers for superpopulations


######Representing the data
freq.pop1=freq.pop[apply(freq.pop,1,function(x) !all(x<4)),] ###removing rows with all populations <0.005 AF

op1=merge(freq.pop1,freq.suppop,by="row.names",all.x=TRUE)
op2=merge(op1,freq.all, by.x="Row.names", by.y="Var1", all.x = TRUE)
comb = combined[!duplicated(combined$dd),]
op3=merge(op2,comb[,-c(1:3)],by.x="Row.names", by.y="dd",all.x=TRUE)

aa=op3[vapply(op3, function(x) length(unique(x)) > 1, logical(1L))]
aa=aa[order(aa$Freq, decreasing = TRUE), ]
write.table(format(aa,digits=3),paste0(name,"-Hap.xls"), quote=FALSE, sep="\t", row.names = FALSE)
