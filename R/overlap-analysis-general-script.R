args=(commandArgs(TRUE))
trait1.file			=	args[1] 			#trait1 output from summary-stat-overlap.R
trait2.file			=	args[2] 			# trait2 output from summary-stat-overlap.R
R					=	as.numeric(args[3]) # ratio of costs = type II error cost/ type I error cost
pi0					=	as.numeric(args[4]) # Pr(H_0); e.g. 0.99, 0.999
alpha				=	as.numeric(args[5]) # p-value sig threshold
ABForPV				=	as.numeric(args[6])	# takes value 1 to use ABFs and value 0 to use p-values
clump.file			=	args[7] 			#pruned SNPs using ABF or PV measure measure; plink output
main.output.file1	=	args[8]				# output list of significant overlap variants and their marginal statistics
main.output.file2	=	args[9]				# output file with parameter information and p-value for test of more overlap variants than expected by chance

ABF <- 1*(ABForPV==1)
PV <- 1*(ABForPV==0)

# output files
fpv.out=args[8]
flist.abf=args[9]
flist.pv=args[10]

source("overlap-analysis.R")

prep.out <- analysis.prep.fn(trait1.file,trait2.file,clump.file,pi0,R,alpha,ABF=ABF,PV=PV)

out <- analysis.fn(prep.out,alpha,R,pi0,ABF=ABF,PV=PV)

write.table(out$overlap.pv,main.output.file2,quote=FALSE,row.names=FALSE)
if(ABF) write.table(out$sigABF,main.output.file1,quote=FALSE,row.names=FALSE,col.names=TRUE)
if(PV) write.table(out$sigpv,main.output.file1,quote=FALSE,row.names=FALSE,col.names=TRUE)							

q(save="no");
