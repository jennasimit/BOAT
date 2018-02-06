args=(commandArgs(TRUE))
trait1.file		=	args[1] 
trait2.file		=	args[2] 
outputT1.file	=	args[3]
outputT2.file	=	args[4]
output1.file	=	args[5]
R			=	as.numeric(args[6])
pi0			=	as.numeric(args[7])
alpha		=	as.numeric(args[8])

#R=20
#pi0=0.99
#alpha=0.001

W=0.04
ABF=TRUE
PV=TRUE

source("pre-overlap-analysis.R")

t1 <- read.table(trait1.file,header=TRUE)
t2 <- read.table(trait2.file,header=TRUE)

abf.data <- data.abf.fn(t1=t1,t2=t2,W=W)

write.table(abf.data$trait1,file=outputT1.file,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(abf.data$trait2,file=outputT2.file,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

abf.measure <- clump.measure.fn(t1=abf.data$trait1,t2=abf.data$trait2,R=R,pi0=pi0,alpha=alpha,ABF=ABF,PV=PV)

write.table(abf.measure,file=output1.file,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

