args=(commandArgs(TRUE))
#r2=as.numeric(args[1])
#mppthr=as.numeric(args[2])
t1=as.numeric(args[1])
t2=as.numeric(args[2])

source("/home/ja628/scratch/scripts/IL2RA_general_scripts/multi-trait/multi-trait-fine-map.R")

r2=0.9
mppthr=0.01

traits <- c("T1D","GRAVES","MS","JIA")
mtraits <- c(4,1,1,0)
Mtraits <- c(6,3,3,3)


trait1=traits[t1]
trait2=traits[t2]
mT1=mtraits[t1]
MT1=Mtraits[t1]
mT2=mtraits[t2]
MT2=Mtraits[t2]


# directory to save files 
fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"

# do only one time: create tag snps object from common controls
# make.tags.fn(r2=0.9)
load(file.path(fdir,paste("tags-r2-",r2,".RData",sep="") ) )


# do once: save "best" tags SNPs for each trait
#sel.snps.fn(trait1,mppthr=0.01,r2=0.9,fdir)
#sel.snps.fn(trait2,mppthr=0.01,r2=0.9,fdir)

fname <- file.path(fdir,paste(trait1,"-gfm-on-tags-r2-",r2,"-best-snps-mpp-",mppthr,".RData",sep=""))
load(fname)
t1snps <- bestsnps

fname <- file.path(fdir,paste(trait2,"-gfm-on-tags-r2-",r2,"-best-snps-mpp-",mppthr,".RData",sep=""))
load(fname)
t2snps <- bestsnps

msnps <- union(t1snps[,"var"],t2snps[,"var"])
s <- length(msnps)

T1mods <- T1mods.fn(mT1,MT1,msnps)
T2mods <- T1mods.fn(mT2,MT2,msnps)


f1 <- paste("/scratch/wallace/IL2RA/full-data/",trait1,sep="")
load(file.path(f1,"data.RData"))
t1snp.data <- snp.data
t1pheno <- cc.ph
t1strat <- strat

f1 <- paste("/scratch/wallace/IL2RA/full-data/",trait2,sep="")
load(file.path(f1,"data.RData"))
t2snp.data <- snp.data
t2pheno <- cc.ph
t2strat <- strat


G1.mat <- as(t1snp.data[,msnps],"numeric")
data1 <- data.frame(Y=t1pheno,G1.mat)

G2.mat <- as(t2snp.data[,msnps],"numeric")
data2 <- data.frame(Y=t2pheno,G2.mat)


traits <- paste(trait1,"-",trait2,sep="")

fabf <- file.path(fdir,paste(trait1,"-bin-tags-abf-",traits,".txt",sep=""))
bf1 <- abfT1.fn(data1,trait1,T1mods,fabf)

fabf <- file.path(fdir,paste(trait2,"-bin-tags-abf",traits,".txt",sep=""))
bf2 <- abfT1.fn(data2,trait2,T2mods,fabf)


#nullmod <- matrix(c(rep(0,(2*s)),1),nrow=1) # no snp effects and only trait effect
#t2m1snp <- t2snp.data[t2pheno==1,]
#allsnps <- rbind(t1snp.data,t2m1snp)
#allpheno <- c(t1pheno,rep(2,sum(t2pheno)))
#y <- character(length(allpheno))
#y[allpheno==0] <- "CONTROL"
#y[allpheno==1] <- trait1
#y[allpheno==2] <- trait2
#G1.mat <- as(allsnps,"numeric")
#data1 <- data.frame(Y=y,G1.mat)
#m1 <- mlogit2logit(Y ~ 1|. -Y,data1,choices=c("CONTROL",trait1,trait2),base.choice=1)
#Nullmod <- glib.1(x=m1$data[,(4+s):dim(m1$data)[2]],y=m1$data$Y.star,error="binomial", link="logit",models=nullmod,post.bymodel = FALSE)
#logBF0 <- Nullmod$bf$twologB10[,1]*0.5

T1T2abf <- T1T2abf.fn(T1abf=bf1,T2abf=bf2)
fname <- file.path(fdir,paste(traits,"-tags-bf.txt",sep=""))
write.table(T1T2abf$out12,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)


indbest <- which(T1T2abf$out12[,1]>log(3))
best1 <- unique(T1T2abf$out1[indbest,])
best2 <- unique(T1T2abf$out2[indbest,])
best12 <- T1T2abf$out12[indbest,]

e1 <- expand.tags.bf(best1,tags)
e2 <- expand.tags.bf(best2,tags)

e1modmat <- mod2matrix.fn(e1)
e2modmat <- mod2matrix.fn(e2)

m2snps <- colnames(e1modmat)

G1.mat <- as(t1snp.data[,m2snps],"numeric")
data1.2 <- data.frame(Y=t1pheno,G1.mat)

G2.mat <- as(t2snp.data[,m2snps],"numeric")
data2.2 <- data.frame(Y=t2pheno,G2.mat)


traits <- paste(trait1,"-",trait2,sep="")

fabf <- file.path(fdir,paste(trait1,"-bin-abf-",traits,".txt",sep=""))
bf1 <- abfT1.fn(data1.2,trait1,e1modmat,fabf)

fabf <- file.path(fdir,paste(trait2,"-bin-abf",traits,".txt",sep=""))
bf2 <- abfT1.fn(data2.2,trait2,e2modmat,fabf)

T1bfexp <- e1
T1bfexp$logbf <- bf1$logABF

T2bfexp <- e2
T2bfexp$logbf <- bf2$logABF

out <- T1T2bfexp.fn(T1T2abf,T1bfexp,T2bfexp,indbest) 
fname <- file.path(fdir,paste(traits,"-bf-final.txt",sep=""))
write.table(out,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)

