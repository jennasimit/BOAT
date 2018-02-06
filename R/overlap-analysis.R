ftrait1=args[1] #trait1 output from summary-stat-overlap.R
ftrait2=args[2] # trait2 output from summary-stat-overlap.R
R=as.numeric(args[3]) # ratio of costs = type II error cost/ type I error cost
pi0=as.numeric(args[4]) # Pr(H_0); enter as a constant initially, but can also form it as a function 
alpha=as.numeric(args[5]) # p-value sig threshold
fclump=args[6] #pruned SNPs using ABF or pv measure
fpv.out=args[7]
flist.abf=args[8]
flist.pv=args[9]


analysis.prep.fn <- function(ftrait1,ftrait2,fclump,pi0,R,alpha,ABF=TRUE,PV=TRUE) 	{

t1 <- read.table(ftrait1,header=TRUE) # output from summary-stat-overlap.R
t2 <- read.table(ftrait2,header=TRUE) # output from summary-stat-overlap.R

if(ABF) 	{
CLabf <- read.table(fclump,header=TRUE) # plink clump output -> abf

# extract CLabf SNPs from each trait for abf analysis 
t1abf <- t1[match(CLabf$SNP,t1$SNP),]
t2abf <- t2[match(CLabf$SNP,t2$SNP),]

## abf analysis	prep								
PO <- pi0/(1-pi0) # prior odds of no association
theta <- PO/R
# reject H0 of no assocn at SNP if ABF > PO/R
sigabf1 <- 1*(t1abf$abf > theta)
sigabf2 <- 1*(t2abf$abf > theta)
return(list(t1abf=t1abf,t2abf=t2abf,sigabf1=sigabf1,sigabf2=sigabf2))

			}
if(!ABF) {
sigabf1 <- NA
sigabf2 <- NA
		}
		
if(PV)		{	
CLpv <- read.table(fclump,header=TRUE) # plink clump output -> pv

# extract CLpv SNPs from each trait for pv analysis
t1pv <- t1[match(CLpv$SNP,t1$SNP),]
t2pv <- t2[match(CLpv$SNP,t2$SNP),]

## pv analysis prep
sigpv1 <- 1*(t1pv$P.value < alpha)
sigpv2 <- 1*(t2pv$P.value < alpha)
return(list(t1pv=t1pv,t2pv=t2pv,sigpv1=sigpv1,sigpv2=sigpv2))
			}

if(!PV) {
sigpv1 <- NA
sigpv2 <- NA
		}

#return(list(t1abf=t1abf,t2abf=t2abf,t1pv=t1pv,t2pv=t2pv,sigabf1=sigabf1,sigabf2=sigabf2,sigpv1=sigpv1,sigpv2=sigpv2))
		}

				
# McNemar mid p-value				
McNemarMidp.fn <- function(contable) {

 n <- contable[1,2] + contable[2,1]

 # exact one-sided McNemar p-value
 l <- min(c(contable[1,2], contable[2,1]))
 print(c(n,l))
 onepv <- pbinom(l,size=n,prob=.5)
 pv <- 2*onepv - dbinom(l,size=n,prob=.5)

  return(pv)
									}


analysis.fn <- function(prep.out,alpha,R,pi0,ABF=TRUE,PV=TRUE)		{

attach(prep.out)

if(ABF) 	{
contable.abf <- table(sigabf1,sigabf2)
pv.abf <- McNemarMidp.fn(contable.abf)									

# overlap snp lists
ind <- which(sigabf1*sigabf2 == 1)
sigABF <- cbind(t1abf[ind,],t2abf[ind,])

			}

if(PV)		{			
contable.pv <- table(sigpv1,sigpv2) 
pv.pv <- McNemarMidp.fn(contable.pv)									

ind <- which(sigpv1*sigpv2 == 1)
sigpv <- cbind(t1pv[ind,],t2pv[ind,])

			}
			
PO <- pi0/(1-pi0) # prior odds of no association
Ltheta <- log10(PO/R)

if(!ABF) {
sigABF <- NA
pv.abf <- NA
		}
		
if(!PV) {
sigpv <- NA
pv.pv <- NA
		}		


pv.out <- list(ABFpv=pv.abf,PVpv=pv.pv,LABFthreshold=Ltheta,PVthreshold=alpha,R=R,pi0=pi0)

detach(prep.out)


return(list(overlap.pv=pv.out,sigABF=sigABF,sigpv=sigpv))

																}
