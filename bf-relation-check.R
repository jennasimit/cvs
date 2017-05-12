args=(commandArgs(TRUE))
#OR1=as.numeric(args[1])
#OR2=as.numeric(args[2])
N0=as.numeric(args[1])
N1=as.numeric(args[2])
N2=as.numeric(args[3])  
#rep=as.numeric(args[6])  # rep number 
#DIRout=args[7]
#DIRtmp=args[8]

#DIRout="/home/ja628/causal_sim/results/IL2RA/model-AD-AC-ORs-T1-1.2-1.2-T2-1.2-1.2/N1a_1000-N0_1000-N1b_1000"
#DIRtmp <- "/home/ja628/scratch/causal_sim/IL2RA/model-AD-AC-ORs-T1-1.2-1.2-T2-1.2-1.2/N1a_1000-N0_1000-N1b_1000"
# assume causals A+D, A+C

#fABFout <- paste(DIRout,"/BF_",rep,".txt",sep="")

OR1=1.2
OR2=1.2
prev=0.1



# max no. of causal snps for each trait
mT1=3 
mT2=3


library(GUESSFM)
library(snpStats)
library(mlogitBMA)
library(BMA)
library(Rcpp)
library(RcppArmadillo)

source("/home/ja628/scratch/scripts/IL2RA_general_scripts/myglib.R")
options(scipen=999)

convert.fn <- function(gcalls) {
 N <- dim(gcalls)[2]-5 # number of individuals*3; first 5 cols are snp info
 m <- dim(gcalls)[1] # number of snps
 snames <- gcalls[,2]
 gcalls <- gcalls[,-(1:5)]
 G <- matrix(0,nrow=m,ncol=N/3)
 rownames(G) <- snames
 colnames(G) <- paste("indiv",1:(N/3),sep="")
 
 gen.fn <- function(i,gcalls) {
  N <- dim(gcalls)[2]
  G <- numeric(N/3)
  c1 <- seq(1,N,by=3) #AA
  c2 <- seq(2,N,by=3) #AB
  c3 <- seq(3,N,by=3) #BB
  G[which(gcalls[i,c1]==1)] <- 0
  G[which(gcalls[i,c2]==1)] <- 1
  G[which(gcalls[i,c3]==1)] <- 2
  return(G)
  }
 
 Gmat <- apply(matrix(1:m,ncol=1),1,gen.fn,gcalls)
 G <- t(Gmat)
 rownames(G) <- snames
 colnames(G) <- paste("indiv",1:(N/3),sep="")
 
  return(G)
				}
   

 
##
inv.logit.fn <-function(x) return(exp(x)/(1+exp(x)))
##

phen.gen.fn <-function(beta1=c(-2.3,.2,.2),beta2=c(-2.3,.2,.2),snpG,N0=100,N1=100,N2=100,causals1.ind,causals2.ind) {
# beta[1] = log(prev)
 N <-dim(snpG)[2] # number of indivs
 #causal1 <- snpG[causals1.ind,]
 #causal2 <- snpG[causals2.ind,]
 n0=0
 n1=0
 n2=0
 G0 <-NULL
 G1 <- NULL
 G2 <- NULL
 j=1
 while( (n0<=N0 | n1 <=N1 | n2 <=N2) & (j <= N) ) {
  p1 <-inv.logit.fn(beta1[1]+sum(beta1[-1]*snpG[causals1.ind,j]))
  p2 <-inv.logit.fn(beta2[1]+sum(beta2[-1]*snpG[causals2.ind,j]))
  u<-runif(2)
   if(u[1]>p1 & u[1]>p2) {
     n0<-n0 + 1
     G0 <-rbind(G0,snpG[,j])
   } else if(u[1]<=p1 & u[1]>p2) {
    n1<-n1 + 1
     G1 <-rbind(G1,snpG[,j])
   } else if(u[1]<=p2 & u[1]>p1) {
     n2<-n2 + 1
     G2 <-rbind(G2,snpG[,j])
    } else {
     if(u[2]<0.5 & n1<N1 ) { n1<-n1 + 1; G1 <-rbind(G1,snpG[,j]) 
     } else { n2<-n2 + 1; G2 <-rbind(G2,snpG[,j]) }
    }
 #print(c(u[1],p1,p2,n0,n1,n2,u[2]))
j <- j+1
	}
#print(dim(G0))
#print(dim(G1))
#print(dim(G2))
G<-rbind(G0[1:N0,],G1[1:N1,],G2[1:N2,])
rownames(G)<-c(paste("control.",1:N0,sep=""),paste("case1.",1:N1,sep=""),paste("case2.",1:N2,sep=""))
colnames(G) <- rownames(snpG)
#colnames(G)<-paste("rs",1:dim(G)[2],sep="")

#return(list(G0=G0,G1=G1,n0=n0,n1=n1,causal.ind=causal.ind))		
 return(list(G=G,y=c(rep(0,N0),rep(1,N1),rep(2,N2)) ))
# output: G=genotype matrix (rows=indiv, cols=snps), y=case(1);control(0), causal.ind =indices of casual variants
  	  	     }

###

T1mods.fn <- function(mt1,msnps) {
 myLetters <- letters[1:26]
match("a", myLetters)
s=length(msnps)
modlist.fn <- function(mT1,s) {
 modT1 <- vector("list",mT1)
  for(i in 1:mT1) {
   mci <- combn(letters[1:s],i,simplify=TRUE)  # all combinations of i msnps
   nci <- dim(mci)[2]
   modT1[[i]] <- matrix(0,nrow=nci,ncol=s)
   for(j in 1:nci)  modT1[[i]][j,match(mci[,j],myLetters)] <- 1 
  }
 modsT1 <- rep(0,s)
 for(i in 1:mT1) modsT1 <- rbind(modsT1,modT1[[i]])
 return(modsT1)
}

T1mods <- modlist.fn(mT1=mT1,s=s)

 
return(T1mods)
}


##

abfT1.fn <- function(data1,mT1=3,msnps) {
# data1 <- data.frame(Y=t1pheno,G1.mat)
# approximates bfs for single trait, considering all "best" tag models from each of T1 and T2
# outputs to file fabf

mods <- T1mods.fn(mT1,msnps)
colnames(mods) <- msnps

mod1 <- glib.1(x=data1[,-1],y=data1$Y,error="binomial", link="logit",models=mods)

logABF <- mod1$bf$twologB10[,1]*0.5

ind <- order(logABF,decreasing=TRUE)
out <- data.frame(logABF=logABF[ind],M=mod1$models[ind,],row.names=NULL)
cnames <- c("logABF",names(data1)[-1])
names(out) <- cnames
#write.table(out,fabf,row.names=FALSE,col.names=cnames,quote=FALSE)
return(out)
}

##

abf.fn <- function(s,mT1,mT2,m1) {
#s=length(m1snps) # number of snps considered in models for the 2 traits
#mTi=max no. of causal snps for Ti; i=1,2
#traits="AD-AC"

#fabf <- paste("/home/ja628/causal_sim/results/IL2RA/",traits,"-abf.txt",sep="")

myLetters <- letters[1:26]
match("a", myLetters)

modlist.fn <- function(mT1,s) {
 modT1 <- vector("list",mT1)
  for(i in 1:mT1) {
   mci <- combn(letters[1:s],i,simplify=TRUE)  # all combinations of i msnps
   nci <- dim(mci)[2]
   modT1[[i]] <- matrix(0,nrow=nci,ncol=s)
   for(j in 1:nci)  modT1[[i]][j,match(mci[,j],myLetters)] <- 1 
  }
 modsT1 <- rep(0,s)
 for(i in 1:mT1) modsT1 <- rbind(modsT1,modT1[[i]])
 return(modsT1)
}

T1mods <- modlist.fn(mT1=mT1,s=s)
if(mT2==mT1) {T2mods <- T1mods} else {T2mods <- modlist.fn(mT1=mT2,s=s)}

nT1 <- dim(T1mods)[1]
nT2 <- dim(T2mods)[1]

T1modsrep <- matrix(rep(t(T1mods),nT2),ncol=ncol(T1mods),byrow=TRUE)
T2modsrep <- matrix(rep(T2mods,each=nT1),ncol=ncol(T2mods),byrow=FALSE)
T1T2mods <- cbind(T1modsrep,T2modsrep)

# add column of 1's to the models for the trait2*effect variable
T1T2mods1 <- cbind(T1T2mods,1)
#T1T2mods0 <- cbind(T1T2mods,0)

#T1T2mods <- rbind(T1T2mods1,T1T2mods0)
T1T2mods <- T1T2mods1

mod1 <- glib.1(x=m1$data[,(4+s):(4+3*s)],y=m1$data$Y.star,error="binomial", link="logit",models=T1T2mods)

logABF <- mod1$bf$twologB10[,1]*0.5

ind <- order(logABF,decreasing=TRUE)

out <- data.frame(logABF=logABF[ind],M=mod1$models[ind,])


cnames <- c("logABF",names(m1$data[,(4+s):(4+3*s)]))
#write.table(out,fabf,row.names=FALSE,col.names=cnames,quote=FALSE)
names(out) <- cnames

out <- out[,-dim(out)[2]] # rm last column, which is 1 for z_2



return(out) # logABF,mod
}

###


#/home/ja628/scratch/software/impute2-1KG/null_100k_1.controls.tags.gen 
######

#g0 <- read.table(paste("/home/ja628/scratch/software/impute2-1KG/N1a_",N1,"-N0_",N0,"-N1b_",N2,"/null_100k_",rep,".controls.tags.gen",sep=""),header=FALSE,as.is=TRUE)
#Nn <- (dim(g0)[2]-5)/3
#snpG <- convert.fn(g0) # genotype matrix (snp rows, indiv cols)
#write.table(snpG,"/home/ja628/scratch/software/impute2-1KG/null_100k.controls.tags.Gmat.txt",quote=FALSE,row.names=TRUE,col.names=TRUE)
snpG <- read.table(gzfile("/home/ja628/scratch/software/impute2-1KG/null_100k.controls.tags.Gmat.txt.gz"),header=TRUE)
Nn <- dim(snpG)[2]
#gcalls <- snpG[,sample(1:Nn,replace=FALSE)] # randomly re-order indivs


A <- c("rs12722563","rs12722522","rs12722508","rs7909519","rs61839660","rs12722496","rs12722495","rs41295049")
B <- "rs2104286"
C <- c("rs11597367","rs35285258","rs11594656")
D <- c("rs3118475","rs62626317","rs41295055","rs62626325","rs41295105","rs56382813")
E <- "rs6602437"
F <- c("rs41295121","rs41295159")
sAITD <- "rs706779" # AITD
sAA <- "rs3118470" # AA
sUC <- "rs4147359" # UC
sRA <- "rs10795791" # RA


msnps <- c(A[1],B[1],C[1],D[1],E[1])
snpG <- snpG[msnps,]

#c12 <- grep(sample(A,1),rownames(snpG))
#c1 <- grep(sample(D,1),rownames(snpG))
#c2 <- grep(sample(C,1),rownames(snpG))

#causals1.ind <- c(c12,c1)
#causals2.ind <- c(c12,c2)

causals1.ind <- 1
causals2.ind <- 1

sim <- phen.gen.fn(beta1=c(log(prev),log(OR1)),beta2=c(log(prev),log(OR2)),snpG=snpG,N0=N0,N1=N1,N2=N2,causals1.ind,causals2.ind)
Gm <- new("SnpMatrix",(sim$G+1)) # snp cols, indivs rows
data1 <- data.frame(Y=sim$y,sim$G)


c0 <- grep("control.",rownames(Gm))
c1 <- grep("case1.",rownames(Gm))
c2 <- grep("case2.",rownames(Gm))

G1 <- Gm[c(c0,c1),]
G2 <- Gm[c(c0,c2),]

data1.1 <- data.frame(Y=sim$y[c(c0,c1)],sim$G[c(c0,c1),])
data1.2 <- data.frame(Y=sim$y[c(c0,c2)],sim$G[c(c0,c2),])
data1.2$Y[data1.2$Y==2] <- 1

m1snps <- msnps
print(m1snps)
s=length(m1snps) # number of snps considered in models for the 2 traits
print(data1[1,])
data.m1 <- data1[,c("Y",m1snps)]

m1 <- mlogit2logit(Y ~ 1|. -Y,data.m1,choices=0:2,base.choice=1)

bft1t2 <- abf.fn(s=s,mT1=mT1,mT2=mT2,m1=m1)

#write.table(logBF,fABFout,quote=FALSE,col.names=TRUE,row.names=FALSE)

bft1 <- abfT1.fn(data1.1,mT1=3,m1snps)
bft2 <- abfT1.fn(data1.2,mT1=3,m1snps)

t1mods <- as.matrix(bft1[,-1])
t2mods <- as.matrix(bft2[,-1])
s <- dim(t1mods)[2]
t1t2mods <- as.matrix(bft1t2[,-1])
colnames(t1mods) <- NULL
colnames(t2mods) <- NULL
colnames(t1t2mods) <- NULL

bfall <- c()

for(k in 1:dim(t1mods)[1]) {
 for(j in 1:dim(t2mods)[1]) {
ind <- which(apply(t1t2mods, 1, identical, c(t1mods[k,],t2mods[j,])))
tmp <- cbind(bft1t2[ind,],logBF1=bft1[k,1],logBF2=bft2[j,1])
bfall <- rbind(bfall, tmp)
}
}
fdir <- "/home/ja628/causal_sim/results/IL2RA/"
write.table(bfall,paste(fdir,"bf-compare-bin-multinom-N0-",N0,"-N1-",N1,"-N2-",N2,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
 
 
model1 <- lm(bfall[,1]~bfall[,"logBF1"]+bfall[,"logBF2"])
 s12 <- bfall[,"logBF1"]+bfall[,"logBF2"]
model2 <- lm(bfall[,1]~s12) 
 
fname <- paste(fdir,"results-bf-compare-bin-multinom-N0-",N0,"-N1-",N1,"-N2-",N2,".RData",sep="")
save(model1,model2,file=fname)

#PP1 <- PP.fn(out,shared=1)
#PP5 <- PP.fn(out,shared=5)
#PP10 <- PP.fn(out,shared=10)
#PP20 <- PP.fn(out,shared=20)


#####

#unlink(g0)

