args=(commandArgs(TRUE))
#r2=as.numeric(args[1])
#mppthr=as.numeric(args[2])
t1=as.numeric(args[1])
t2=as.numeric(args[2])
t3=as.numeric(args[3])

source("/home/ja628/scratch/scripts/IL2RA_general_scripts/multi-trait/multi-trait-fine-map-v5.R")

r2=0.9
mppthr=0.001

traits <- c("T1D","GRAVES","MS","JIA")
mtraits <- c(4,1,1,0)
Mtraits <- c(5,3,2,2)


trait1=traits[t1]
trait2=traits[t2]
trait3=traits[t3]

mT1=mtraits[t1]
MT1=Mtraits[t1]
mT2=mtraits[t2]
MT2=Mtraits[t2]
mT3=mtraits[t3]
MT3=Mtraits[t3]


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

fname <- file.path(fdir,paste(trait3,"-gfm-on-tags-r2-",r2,"-best-snps-mpp-",mppthr,".RData",sep=""))
load(fname)
t3snps <- bestsnps


msnps <- union(union(t1snps[,"var"],t2snps[,"var"]),t3snps[,"var"])
s <- length(msnps)

T1mods <- T1mods.fn(mT1,MT1,msnps)
T2mods <- T1mods.fn(mT2,MT2,msnps)
T3mods <- T1mods.fn(mT3,MT3,msnps)

j.mat <- matrix(1:dim(T1mods)[1],ncol=1) 
t1mod <- apply(j.mat,1,format.mod.fn,T1mods)
t1size <- apply(T1mods,1,sum)
T1mod <- data.frame(mod=t1mod,size=t1size)

j.mat <- matrix(1:dim(T2mods)[1],ncol=1) 
t2mod <- apply(j.mat,1,format.mod.fn,T2mods)
t2size <- apply(T2mods,1,sum)
T2mod <- data.frame(mod=t2mod,size=t2size)

j.mat <- matrix(1:dim(T3mods)[1],ncol=1) 
t3mod <- apply(j.mat,1,format.mod.fn,T3mods)
t3size <- apply(T3mods,1,sum)
T3mod <- data.frame(mod=t3mod,size=t3size)


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

f1 <- paste("/scratch/wallace/IL2RA/full-data/",trait3,sep="")
load(file.path(f1,"data.RData"))
t3snp.data <- snp.data
t3pheno <- cc.ph
t3strat <- strat


G1.mat <- as(t1snp.data[,msnps],"numeric")
data1 <- data.frame(Y=t1pheno,G1.mat)

G2.mat <- as(t2snp.data[,msnps],"numeric")
data2 <- data.frame(Y=t2pheno,G2.mat)

G3.mat <- as(t3snp.data[,msnps],"numeric")
data3 <- data.frame(Y=t3pheno,G3.mat)


traits <- paste(trait1,"-",trait2,"-",trait3,sep="")

bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=T1mod$mod,family="binomial")[[1]] 
bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=T2mod$mod,family="binomial")[[1]] 
bf3 <- abf.calc(y=data3[,1],x=data3[,-1],models=T3mod$mod,family="binomial")[[1]] 


bf1$size <- t1size
bf2$size <- t2size
bf3$size <- t3size


nullmod <- matrix(c(rep(0,(2*s)),1,rep(0,s),1),nrow=1) # no snp effects and only trait effects
#t2m1snp <- t2snp.data[t2pheno==1,]
allsnps <- rbind(G1.mat[,msnps],G2.mat[t2pheno==1,msnps],G3.mat[t3pheno==1,msnps])
allpheno <- c(t1pheno,rep(2,sum(t2pheno)),rep(3,sum(t3pheno)))
y <- character(length(allpheno))
y[allpheno==0] <- "CONTROL"
y[allpheno==1] <- trait1
y[allpheno==2] <- trait2
y[allpheno==3] <- trait3
#G1.mat <- as(allsnps,"numeric")
data1 <- data.frame(Y=y,allsnps)
m1 <- mlogit2logit(Y ~ 1|. -Y,data1,choices=c("CONTROL",trait1,trait2,trait3),base.choice=1)
Nullmod <- glib.1(x=m1$data[,(4+s):dim(m1$data)[2]],y=m1$data$Y.star,error="binomial", link="logit",models=nullmod,post.bymodel = FALSE)
logBF0 <- Nullmod$bf$twologB10[,1]*0.5


T1T2T3abf <- T1T2T3abfexp.fn(bf1,bf2,bf3,logBF0)


# calculate PP and only expand joint models with PP>0.001
pp<-PP123.fn(T1T2T3abf,shared=1,s=length(unique(tags(tags))),mT1=3,mT2=3,mT3=3)

fname <- file.path(fdir,paste(traits,"-PP-tags.RData",sep=""))
save(pp,T1T2T3abf,file=fname)

indbest <- which(pp$PP>0.0001)

BFkeep <- pp[-indbest,]

BFexp.notbest <- T1T2T3bfexp.low.fn(BFkeep,tags,logBF0)


# expand tags at all marginal models that are part of a "best" joint model
tmp <- data.frame(logbf=pp$t1.logbf,mod=pp$t1.mod,size=pp$t1.size)
best1 <- unique(tmp[indbest,])

tmp <- data.frame(logbf=pp$t2.logbf,mod=pp$t2.mod,size=pp$t2.size)
best2 <- unique(tmp[indbest,])

tmp <- data.frame(logbf=pp$t3.logbf,mod=pp$t3.mod,size=pp$t3.size)
best3 <- unique(tmp[indbest,])

e1 <- expand.tags.bf(best1,tags)
e2 <- expand.tags.bf(best2,tags)
e3 <- expand.tags.bf(best3,tags)

G1.mat <- as(t1snp.data,"numeric")
data1 <- data.frame(Y=t1pheno,G1.mat)

G2.mat <- as(t2snp.data,"numeric")
data2 <- data.frame(Y=t2pheno,G2.mat)

G3.mat <- as(t3snp.data,"numeric")
data3 <- data.frame(Y=t3pheno,G3.mat)


bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=e1$str,family="binomial") 
bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=e2$str,family="binomial") 
bf3 <- abf.calc(y=data3[,1],x=data3[,-1],models=e3$str,family="binomial") 

T1bfexp <- e1
T1bfexp$logbf <- bf1[[1]][,"lBF"]

T2bfexp <- e2
T2bfexp$logbf <- bf2[[1]][,"lBF"]

T3bfexp <- e3
T3bfexp$logbf <- bf3[[1]][,"lBF"]



# at expanded tags, re-fit joint models 

T1T2T3bfexp <- T1T2T3bfexp.best.fn(T1T2T3abf,T1bfexp,T2bfexp,T3bfexp,indbest,logBF0)

# merge expanded re-fit joint with non-expanded non-re-fit joint (PP<0.001) 
tmp<-BFexp.notbest[,c("t1.logbf","t1.str","t1.size","t2.logbf","t2.str","t2.size","t3.logbf","t3.str","t3.size","t1t2t3.logbf","t1t2t3.Nexpmod")]
bf <- rbind(T1T2T3bfexp,tmp)

fname <- file.path(fdir,paste(traits,"-bf-final.txt",sep=""))
write.table(bf,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)

##### PP calc

## binomial calc
# pp1 <- PP.bin.fn(T1bfexp,s=length(snps(tags)),mT1=3)
# pp2 <- PP.bin.fn(T2bfexp,s=length(snps(tags)),mT1=3)

s=length(snps(tags))

shared <- c(1,5,10,20,100,1000)
ns <- length(shared)

PP <- NULL
for(k in 1:ns) PP <- rbind(PP,PP123.fn(bf,shared=shared[k],s=s,mT1=3,mT2=3,mT3=3,details=FALSE))
row.names(PP) <- shared
fname <- file.path(fdir,paste(traits,"-PP-shared.txt",sep=""))
PP <- t(PP)
write.table(PP,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)

#PP <- read.table(fname,header=TRUE,row.names=1)

PPall <- data.frame(PP,bf$t1t2t3.Nexpmod)

alltraits <- c(trait1,trait2,trait3)
K <- length(alltraits)

PPmarg <- PP.marg.fn(PP,K)
for(k in 1:K) {
fname <- file.path(fdir,paste(alltraits[k],"-PP-shared-",traits,".txt",sep=""))
write.table(PPmarg[[k]],fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
mpp <- MPP.fn(PPmarg[[k]])
fname <- file.path(fdir,paste(alltraits[k],"-MPP-shared-",traits,".txt",sep=""))
write.table(mpp,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)

}





