

# abf calculation
abf.fn <- function(i,beta,se,W)	{
 V <- se[i]^2
 z <- beta[i]/sqrt(V)
 t1 <- V/(V+W)
 t2 <- W/(V+W)
 abf <- sqrt(t1)*exp(t2*(z^2)/2)
 return(abf)
								}

data.abf.fn <- function(t1,t2,W=0.04)	{

# only consider SNPs present in both data sets by rsid and on chr
intsnp <- intersect(t1$SNP,t2$SNP)
trait1 <- t1[match(intsnp,t1$SNP),]
trait2 <- t2[match(intsnp,t2$SNP),]
m <- dim(trait1)[1]

attach(trait1)
nmat <- matrix(1:m,nrow=m)
abf <- apply(nmat,1,abf.fn,beta=Effect,se=StdErr,W) # 2 by m output
t1 <- data.frame(trait1,abf=abf) 
detach(trait1)							
trait1 <- t1
rm(t1)

attach(trait2)
nmat <- matrix(1:m,nrow=m)
abf <- apply(nmat,1,abf.fn,beta=Effect,se=StdErr,W) # 2 by m output
t2 <- data.frame(trait2,abf=abf) 
detach(trait2)				
trait2 <- t2
rm(t2)	

return(list(trait1=trait1,trait2=trait2))

								}
								

clump.measure.fn <- function(t1,t2,R,pi0,alpha,ABF=ABF,PV=PV)	{

sigt1t2 <- t1
								
if(ABF)	{
PO <- pi0/(1-pi0) # prior odds of no association
theta <- PO/R

M <- max(t1$abf,t2$abf,na.rm=TRUE)
sigabfBoth <- 1*(t1$abf > theta)*(t2$abf > theta)
tabf <- pmax(t1$abf,t2$abf) + sigabfBoth*M  # use most sig value of two traits at the SNP and if both are sig then boost value to increase chance of not pruning out;
# this accounts for the scenarios where both traits are marginally sig at one SNP, and at a SNP in LD with it, only one trait meets sig and has a very large sig value -->
# taking the most sig value at each SNP would mean we keep the SNP with only one trait assoc and prune out the one that is sig for both traits
sigt1t2$abf <- tabf
sigt1t2$neg.abf <- -tabf
		} 

if(PV)	{
sigpvboth <- 1*(t1$P.value < alpha)*(t2$P.value < alpha)
tpv <- pmin(t1$P.value,t2$P.value) - sigpvboth # may get a -ve value here but that is only for pruning purposes
sigt1t2$P.value <- tpv
		}


return(sigt1t2)	
#write.table(sigt1t2,file=fabfpv,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

	}							
