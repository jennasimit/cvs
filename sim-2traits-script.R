#OR1a=as.numeric(args[1]) # OR for cv1 of trait 1
#OR2a=as.numeric(args[2]) # OR for cv2 of trait 1
#OR1b=as.numeric(args[3]) # OR for cv1 of trait 2
#OR2b=as.numeric(args[4]) # OR for cv2 of trait 2
#N=as.numeric(args[5])  # number of controls=cases1=cases2
#rep=as.numeric(args[6])  # rep number 
#DIRout=args[7]
#DIRtmp=args[8]

mydir <- "tmp"


OR1a=1.2; OR2a=1.2; OR1b=1.2; OR2b=1.2

trait1="AD"
trait2="AC"

r2=0.9
mppthr=0.001

mT1=1
MT1=3
mT2=1
MT2=3

#fABFout <- paste(DIRout,"/BF_",rep,".txt",sep="")

prev=0.1

N0=N
N1=N
N2=N

kappa <- c(1,10,50,500,1000)
nk <- length(kappa)

library(MTFM)
source("sim-2traits-functions.R")

snpG <-	read.table(gzfile("/scratch/wallace/IL2RA/null_100k.controls_1.tags.Gmat.txt.gz"),header=TRUE) # simulated null data in hapgen2 

A <- c("rs12722563","rs12722522","rs12722508","rs7909519","rs61839660","rs12722496","rs12722495","rs41295049")
B <- "rs2104286"
C <- c("rs11597367","rs35285258","rs11594656")
D <- c("rs3118475","rs62626317","rs41295055","rs62626325","rs41295105","rs56382813","rs41295079")
E <- "rs6602437"
F <- c("rs41295121","rs41295159")
sAITD <- "rs706779" # AITD
sAA <- "rs3118470" # AA
sUC <- "rs4147359" # UC
sRA <- "rs10795791" # RA



c12 <- grep(sample(A,1),rownames(snpG))
c1 <- grep(sample(D,1),rownames(snpG))
c2 <- grep(sample(C,1),rownames(snpG))

causals1.ind <- c(c12,c1)
causals2.ind <- c(c12,c2)


sim <- phen.gen.fn(beta1=c(log(prev),log(OR1a),log(OR2a)),beta2=c(log(prev),log(OR1b),log(OR2b)),snpG=snpG,N0=N0,N1=N1,N2=N2,causals1.ind,causals2.ind)
Gm <- new("SnpMatrix",(sim$G+1)) # snp cols, indivs rows
data1 <- data.frame(Y=sim$y,sim$G)


c0 <- grep("control.",rownames(Gm))
c1 <- grep("case1.",rownames(Gm))
c2 <- grep("case2.",rownames(Gm))

G1 <- Gm[c(c0,c1),]
G2 <- Gm[c(c0,c2),]
t1pheno <- c(rep(0,N0),rep(1,N1))
t2pheno <- c(rep(0,N0),rep(1,N2))


tags <- make.tags.fn(r2=r2,G0=Gm[c0,],mydir=NA) # set mydir=NA if do not want to save tags to an Rdata file

pp <- sharedmPP.fn(kappa,G1,t1pheno,G2,t2pheno,tags,mppthr,mT1,MT1,mT2,MT2,trait1,trait2,mydir)

ppc <- vector("list",2)
for(j in 1:2) { for(k in 1:nk) { ppc[[j]] <- rbind(ppc[[j]],PP2collapse.fn(pp[[j]][,k]))}}

# save marginal model PP
#write.table(pp$pp1,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
#write.table(pp$pp2,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)

# save marginal model PP in terms of snp group
#write.table(ppc[[1]],fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
#write.table(ppc[[2]],fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
