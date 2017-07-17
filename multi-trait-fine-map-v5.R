library(annotSnpStats)
library(snpStats)
library(GUESSFM)
library(mlogitBMA)
library(BMA)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(data.table)
library(gtools)
library(fields)
library(dplyr)


source("/home/ja628/scratch/scripts/IL2RA_general_scripts/myglib.R")
options(scipen=999)

make.tags.fn <- function(r2){
load("/scratch/wallace/IL2RA/full-data/MS/data.RData")
load("/scratch/wallace/IL2RA/full-data/MS/tags.RData")
DATA <- new("SnpMatrix",snp.data[cc.ph==0,snps(tags)]) 
## DATA <- as(snp.data[cc.ph==0,snps(tags)],"SnpMatrix")
tfile <- file.path(fdir,paste("tags-r2-",r2,".RData",sep="") )
tags <- tag(DATA, tag.threshold = r2)
message("saving tags object to ", tfile)
save(tags, file = tfile)
}

##

gfm.on.tags.fn <- function(t1,r2,tags,fdir) {
 #fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
  
 tsnp <- snp.data[,unique(tags@tags)] 
 
 mydir <- file.path(fdir,paste(t1,"-gfm-on-tags-r2-",r2,sep=""))
 run.bvs(X=tsnp,Y=cc.ph,tag.r2=NA,nexp=3,nsave=1000,gdir=mydir,wait=TRUE) 

}

gfm.on.tags.final.fn <- function(t1,r2,tags,fdir) {
 #fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 
 mydir <- file.path(fdir,paste(t1,"-gfm-on-tags-r2-",r2,sep=""))
 d <- read.snpmod(mydir)
 dx <- expand.tags(d,tags)

 fname <- file.path(fdir,paste(t1,"-dx-gfm-r2-",r2,"-expanded.Rdata",sep=""))
 save(dx,file=fname)
 
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
  
 X <- snp.data[,snps(tags)] 
 
 best <- best.models(dx,pp.thr=0.0001)
 abf <- abf.calc(y=cc.ph,x=X,models=best$str,family="binomial") 
 sm <- abf2snpmod(abf,expected=3,nsnps=abf$nsnp)
 best2 <- best.models(sm,pp.thr=0.0001) # best models after re-fit
 bestsnps <- best.snps(sm,pp.thr=0)
fname <- file.path(fdir,paste(t1,"-gfm-r2-",r2,"-final-best-mods.RData",sep=""))
save(best2,file=fname)
fname <- file.path(fdir,paste(t1,"-gfm-r2-",r2,"-final-best-snps.RData",sep=""))
save(bestsnps,file=fname)
}

##

pp.nsnp.mods.fn <- function(t1,r2,tags,fdir) {
# fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
fname <- file.path(fdir,paste(t1,"-dx-gfm-r2-",r2,"-expanded.Rdata",sep=""))
 load(fname)
 
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
  
 X <- snp.data[,snps(tags)] 
 
 best <- best.models(dx,pp.thr=0.0001)
 abf <- abf.calc(y=cc.ph,x=X,models=best$str,family="binomial") 
 sm <- abf2snpmod(abf,expected=3,nsnps=abf$nsnp)
out <- pp.nsnp(sm,expected=3)
print(t1)
return(out)
}

###


gfm.adj.fn <- function(t1,msnps,r2,tags,traits,fdir) {
 #fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 # msnps is union of bestsnps from traits 1 and 2
 # traits <- paste(t1,"-",t2,sep="")
  
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
  
 X <- snp.data[,msnps] 
 X.all <- snp.data[,snps(tags)] 
 
 mydir <- file.path(fdir,paste(t1,"-gfm-on-",traits,"-bestsnps-r2-",r2,sep=""))
 run.bvs(X=X,Y=cc.ph,tag.r2=NA,nexp=3,nsave=1000,gdir=mydir,wait=TRUE) 
 d <- read.snpmod(mydir)
 dx <- expand.tags(d,tags)

 fname <- file.path(fdir,paste(t1,"-dx-gfm-on-",traits,"-bestsnps-r2-",r2,"-expanded.Rdata",sep=""))
 save(dx,file=fname)

 
 best <- best.models(dx,pp.thr=0.0001)
 abf <- abf.calc(y=cc.ph,x=X.all,models=best$str,family="binomial") 
 sm <- abf2snpmod(abf,expected=3,nsnps=abf$nsnp)
 best2 <- best.models(sm,pp.thr=0.0001) # best models after re-fit
 bestsnps <- best.snps(sm,pp.thr=0)
fname <- file.path(fdir,paste(t1,"-gfm-on-",traits,"-bestsnps-r2-",r2,"-final-best-mods.RData",sep=""))
save(best2,file=fname)
fname <- file.path(fdir,paste(t1,"-gfm-on-",traits,"-bestsnps-r2-",r2,"-final-best-snps.RData",sep=""))
save(bestsnps,file=fname)
}
###


###

# find tag SNPs with MPP > mppthr after running GUESSFM on tag SNPs and save to file
# print table of posterior distribution for number of causal variants in GFM models
sel.snps.fn <- function(t1,mppthr=0.01,r2=0.9,fdir) {
 #fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 
 mydir <- file.path(fdir,paste(t1,"-gfm-on-tags-r2-",r2,sep=""))
 #mydir <- file.path(fdir,paste(t1,"-gfm-on-tags",sep="")) #r2=0.99

 d <- read.snpmod(mydir)  
 bestsnps <- best.snps(d,pp.thr=mppthr)

 print(t1)
 print(pp.nsnp(d,expected=3))
 fname <- file.path(fdir,paste(t1,"-gfm-on-tags-r2-",r2,"-best-snps-mpp-",mppthr,".RData",sep=""))
 save(bestsnps,file=fname)
}

###

T1mods.fn <- function(mT1,MT1,msnps) {
 # mT1 and MT1 are the min and max number of causals in a model for trait t1
 # msnps is the vector of snps to consider in the models
 s <- length(msnps)
 if(mT1>0) mc <- mT1:MT1
 if(mT1==0) mc <- (mT1+1):MT1 # start at 1 instead of 0 and do 0 snp model later
 modT1 <- vector("list",length(mc))
  for(i in 1:length(mc)) {
   mci <- combn(msnps,mc[i],simplify=TRUE)  # all combinations of i msnps
   nci <- dim(mci)[2]
   modT1[[i]] <- matrix(0,nrow=nci,ncol=s,dimnames=list(NULL,msnps))
   for(j in 1:nci)  modT1[[i]][j,match(mci[,j],msnps)] <- 1 
  }
 if(mT1 > 0) modsT1 <- NULL
 if(mT1 == 0) modsT1 <- matrix(rep(0,s),nrow=1)
 for(i in 1:length(mc)) modsT1 <- rbind(modsT1,modT1[[i]])
 return(modsT1)
}


###

abfT1.fn <- function(data1,trait1,mods,fabf) {
# data1 <- data.frame(Y=t1pheno,G1.mat)
# approximates bfs for single trait, considering models based on all "best" tag snps from both T1 and T2
# outputs to file fabf


mod1 <- glib.1(x=data1[,-1],y=data1$Y,error="binomial", link="logit",models=mods)

logABF <- mod1$bf$twologB10[,1]*0.5

#ind <- order(logABF,decreasing=TRUE)
#out <- data.frame(logABF=logABF[ind],M=mod1$models[ind,],row.names=NULL)
out <- data.frame(logABF=logABF,M=mod1$models,row.names=NULL)

cnames <- c("logABF",names(data1)[-1])
names(out) <- cnames
out$size <- apply(out[,-1],1,sum)


if(!is.na(fabf)) write.table(out,fabf,row.names=FALSE,col.names=TRUE,quote=FALSE)
return(out)
}

##

format.mod.fn <- function(k,out) {
 ind <- which(out[k,]==1)
 if(length(ind)>0) {
  mod <- paste(names(out[k,][ind]),sep="",collapse="%")
  } else {mod <- ""}
  return(mod)
	}



##
mod.fn <- function(k,out) {
 ind <- which(out[k,-c(1,dim(out)[2])]==1)
 if(length(ind)>0) {
  mod <- paste(names(out[k,-1][ind]),sep="",collapse="%")
  } else {mod <- 0}
  return(mod)
	}

##
options(stringsAsFactors = FALSE)
format.bf <- function(bf,s) {
 j.mat <- matrix(1:dim(bf)[1],ncol=1) 
 mod <- apply(j.mat,1,mod.fn,bf)
 out <- data.frame(logbf=bf[,1],mod,size=bf[,dim(bf)[2]])
 return(out)
 }

##

T1T2abf.fn <- function(T1abf,T2abf) {

nT1 <- dim(T1abf)[1]
nT2 <- dim(T2abf)[1]
s <- dim(T1abf)[2]-2 # number of snps (2 columns have bf and size)

T1modsrep <- matrix(rep(t(T1abf),nT2),ncol=ncol(T1abf),byrow=TRUE)
T2modsrep <- matrix(rep(as.matrix(T2abf),each=nT1),ncol=ncol(T2abf),byrow=FALSE)
T1T2mods <- cbind(T1modsrep,T2modsrep)

tmp <- T1T2mods[,c(1:(s+2))]
colnames(tmp) <- names(T1abf)
out1 <- format.bf(tmp,s)
tmp <- T1T2mods[,c((s+3):(2*s+4))]
colnames(tmp) <- names(T2abf)
out2 <- format.bf(tmp,s)


colnames(T1T2mods) <- c(paste("t1",names(T1abf),sep="."),paste("t2",names(T2abf),sep="."))

t1t2 <- T1T2mods[,"t1.logABF"]+T1T2mods[,"t2.logABF"] 
out12 <- cbind(t1t2.logbf=t1t2,T1T2mods)



return(list(out1=out1,out2=out2,out12=out12))
}

##

T1T2abfexp.fn <- function(T1bfexp,T2bfexp,logBF0) {
# use to combine expanded model bfs from traits; all models

nT1 <- dim(T1bfexp)[1]
nT2 <- dim(T2bfexp)[1]
T1abf <- T1bfexp
T2abf <- T2bfexp

T1modsrep <- matrix(rep(t(T1abf),nT2),ncol=ncol(T1abf),byrow=TRUE)
T2modsrep <- matrix(rep(as.matrix(T2abf),each=nT1),ncol=ncol(T2abf),byrow=FALSE)
T1T2mods <- cbind(T1modsrep,T2modsrep)


colnames(T1T2mods) <- c(paste("t1",names(T1abf),sep="."),paste("t2",names(T2abf),sep="."))
T1T2mods <- as.data.frame(T1T2mods,stringsAsFactors =FALSE)

T1T2mods$t1.lBF <- as.numeric(T1T2mods$t1.lBF)
T1T2mods$t2.lBF <- as.numeric(T1T2mods$t2.lBF)


t1t2 <- T1T2mods$t1.lBF+T1T2mods$t2.lBF+logBF0
out12 <- data.frame(t1t2.logbf=t1t2,T1T2mods)
names(out12) <- c("t1t2.logbf", "t1.str"  , "t1.tag"  ,   "t1.logbf" ,"t1.size"  ,  "t2.str" ,"t2.tag"   ,  "t2.logbf","t2.size")


out12$t1t2.Nexpmod <- 1
return(out12)
}


###

T1T2bfexp.best.fn <- function(T1T2abf,T1bfexp,T2bfexp,indbest,logBF0) {
# on expanded marginal models, re-calculate joint bfs at "best" joint models so have joint bfs on expanded joint models

outkeep <- T1T2abf[indbest,]
names(T1bfexp) <- c("et1.logbf","t1.str","et1.size","et1.str","et1.rank")
b1 <- inner_join(outkeep,T1bfexp) # merge on t1.str to give all t1  expanded models 
names(T2bfexp) <- c("et2.logbf","t2.str","et2.size","et2.str","et2.rank")
b12 <- inner_join(b1,T2bfexp) # merge on t2.str

bf.all <- b12[,c("et1.logbf","et1.str", "et1.size","et2.logbf","et2.str", "et2.size")]
names(bf.all) <- c("t1.logbf","t1.str", "t1.size","t2.logbf","t2.str", "t2.size")

    
bf.all$t1t2.logbf <- bf.all$t1.logbf +bf.all$t2.logbf+logBF0

bf.all$t1t2.Nexpmod <- 1

out <- bf.all[order(bf.all$t1t2.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0

return(out)
}
##

T1T2bfexp.best.sim.v1.fn <- function(T1T2abf,T1bfexp,T2bfexp,indbest) {
# on expanded marginal models, re-calculate joint bfs at "best" joint models so have joint bfs on expanded joint models

 out1keep <- T1T2abf$out1[indbest,] # re-fit these
 out2keep <- T1T2abf$out2[indbest,]


 bfall <- c()

k.mat <- matrix(1:dim(T1T2abf$out12[indbest,])[1],ncol=1)

exp.fn <- function(k,T1bfexp,T2bfexp,out1keep,out2keep) {
#  for(k in 1:dim(T1T2abf$out12)[1]) {   
    ind1 <- which(apply(as.matrix(T1bfexp[,"mod"]), 1, identical, out1keep[k,"mod"]))
    tmp1 <- T1bfexp[ind1,c("logbf","str","size")]
    
    ind2 <- which(apply(as.matrix(T2bfexp[,"mod"]), 1, identical, out2keep[k,"mod"]))
    tmp2 <- T2bfexp[ind2,c("logbf","str","size")]
    
    T1modsrep <- matrix(rep(t(tmp1),dim(tmp2)[1]),ncol=ncol(tmp1),byrow=TRUE)
	T2modsrep <- matrix(rep(as.matrix(tmp2),each=dim(tmp1)[1]),ncol=ncol(tmp2),byrow=FALSE)
	tmp <- cbind(T1modsrep,T2modsrep)
	T1T2mods <- data.frame(tmp)
	names(T1T2mods) <- c("t1.logbf","t1.str", "t1.size","t2.logbf","t2.str", "t2.size")
   return(T1T2mods)
#    bfall <- rbind(bfall, T1T2mods)
    }

bfall <- apply(k.mat,1,exp.fn,T1bfexp,T2bfexp,out1keep,out2keep)
bf.all <- do.call(rbind,bfall)

bf.all$t1.logbf <- as.numeric(bf.all$t1.logbf)
bf.all$t2.logbf <- as.numeric(bf.all$t2.logbf)
    
t1t2bf <- bf.all$t1.logbf +bf.all$t2.logbf

bf.all$t1t2.logbf <- t1t2bf

out <- bf.all[order(bf.all$t1t2.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0

return(out)
}

####

T1T2bfexp.low.fn <- function(BFkeep,tags,logBF0) {
# at expanded marginal models with low PP (keep same bf as tag model) use counts of no. of expanded models to get bfs (and counts of these models) at  expanded joint models

tmp <- data.frame(logbf=BFkeep$t1.logbf,mod=BFkeep$t1.mod,size=BFkeep$t1.size)
nbest1 <- unique(tmp)
en1 <- expand.tags.bf(nbest1,tags)
T1lowexp <- en1[match(unique(en1$mod),en1$mod),c("logbf","mod","size")]
test<-table(en1$mod)
T1lowexp$Nexpmod <- test[match(T1lowexp$mod,names(test))]

tmp <- data.frame(logbf=BFkeep$t2.logbf,mod=BFkeep$t2.mod,size=BFkeep$t2.size)
nbest2 <- unique(tmp)
en2 <- expand.tags.bf(nbest2,tags)
T2lowexp <- en2[match(unique(en2$mod),en2$mod),c("logbf","mod","size")]
test<-table(en2$mod)
T2lowexp$Nexpmod <- test[match(T2lowexp$mod,names(test))]

b1 <- T1lowexp[match(BFkeep[,"t1.mod"],as.matrix(T1lowexp[,"mod"])),c("logbf","mod","size","Nexpmod")]
b2 <- T2lowexp[match(BFkeep[,"t2.mod"],as.matrix(T2lowexp[,"mod"])),c("logbf","mod","size","Nexpmod")]
bf.all  <- data.frame(b1,b2)
names(bf.all) <- c("t1.logbf","t1.str", "t1.size","t1.Nexpmod","t2.logbf","t2.str", "t2.size","t2.Nexpmod")
    
t1t2bf <- bf.all$t1.logbf +bf.all$t2.logbf+logBF0
Nt1t2 <- bf.all$t1.Nexpmod *bf.all$t2.Nexpmod

bf.all$t1t2.logbf <- t1t2bf
bf.all$t1t2.Nexpmod <- Nt1t2

out <- bf.all[order(bf.all$t1t2.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0

return(out)
}
##

T1T2T3bfexp.low.fn <- function(BFkeep,tags,logBF0) {
# at expanded marginal models with low PP (keep same bf as tag model) use counts of no. of expanded models to get bfs (and counts of these models) at  expanded joint models

tmp <- data.frame(logbf=BFkeep$t1.logbf,mod=BFkeep$t1.mod,size=BFkeep$t1.size)
nbest1 <- unique(tmp)
en1 <- expand.tags.bf(nbest1,tags)
T1lowexp <- en1[match(unique(en1$mod),en1$mod),c("logbf","mod","size")]
test<-table(en1$mod)
T1lowexp$Nexpmod <- test[match(T1lowexp$mod,names(test))]

tmp <- data.frame(logbf=BFkeep$t2.logbf,mod=BFkeep$t2.mod,size=BFkeep$t2.size)
nbest2 <- unique(tmp)
en2 <- expand.tags.bf(nbest2,tags)
T2lowexp <- en2[match(unique(en2$mod),en2$mod),c("logbf","mod","size")]
test<-table(en2$mod)
T2lowexp$Nexpmod <- test[match(T2lowexp$mod,names(test))]

tmp <- data.frame(logbf=BFkeep$t3.logbf,mod=BFkeep$t3.mod,size=BFkeep$t3.size)
nbest3 <- unique(tmp)
en3 <- expand.tags.bf(nbest3,tags)
T3lowexp <- en3[match(unique(en3$mod),en3$mod),c("logbf","mod","size")]
test<-table(en3$mod)
T3lowexp$Nexpmod <- test[match(T3lowexp$mod,names(test))]

b1 <- T1lowexp[match(BFkeep[,"t1.mod"],as.matrix(T1lowexp[,"mod"])),c("logbf","mod","size","Nexpmod")]
b2 <- T2lowexp[match(BFkeep[,"t2.mod"],as.matrix(T2lowexp[,"mod"])),c("logbf","mod","size","Nexpmod")]
b3 <- T3lowexp[match(BFkeep[,"t3.mod"],as.matrix(T3lowexp[,"mod"])),c("logbf","mod","size","Nexpmod")]

bf.all  <- data.frame(b1,b2,b3)
names(bf.all) <- c("t1.logbf","t1.str", "t1.size","t1.Nexpmod","t2.logbf","t2.str", "t2.size","t2.Nexpmod","t3.logbf","t3.str", "t3.size","t3.Nexpmod")
     
bf.all$t1t2t3.logbf  <- bf.all$t1.logbf +bf.all$t2.logbf +bf.all$t3.logbf +logBF0
bf.all$t1t2t3.Nexpmod <- bf.all$t1.Nexpmod *bf.all$t2.Nexpmod*bf.all$t3.Nexpmod


out <- bf.all[order(bf.all$t1t2t3.logbf, decreasing = TRUE), ]
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t3.size==0)
if(length(check0)>0) out[check0,"t3.str"] <- 0

return(out)
}



###

T1T2T3abf.fn <- function(T1abf,T2abf,T3abf) {

nT1 <- dim(T1abf)[1]
nT2 <- dim(T2abf)[1]
nT3 <- dim(T3abf)[1]
s <- dim(T1abf)[2]-2 # number of snps (2 columns have bf and size)

T1modsrep <- matrix(rep(t(T1abf),nT2),ncol=ncol(T1abf),byrow=TRUE)
T2modsrep <- matrix(rep(as.matrix(T2abf),each=nT1),ncol=ncol(T2abf),byrow=FALSE)
T1T2mods <- cbind(T1modsrep,T2modsrep)

nT1T2 <- dim(T1T2mods)[1]

T1T2modsrep <- matrix(rep(t(T1T2mods),nT3),ncol=ncol(T1T2mods),byrow=TRUE)
T3modsrep <- matrix(rep(as.matrix(T3abf),each=nT1T2),ncol=ncol(T3abf),byrow=FALSE)
T1T2T3mods <- cbind(T1T2modsrep,T3modsrep) 

tmp <- T1T2T3mods[,c(1:(s+2))]
colnames(tmp) <- names(T1abf)
out1 <- format.bf(tmp,s)
tmp <- T1T2T3mods[,c((s+3):(2*s+4))]
colnames(tmp) <- names(T2abf)
out2 <- format.bf(tmp,s)
tmp <- T1T2T3mods[,c((2*s+5):(3*s+6))]
colnames(tmp) <- names(T3abf)
out3 <- format.bf(tmp,s)


colnames(T1T2T3mods) <- c(paste("t1",names(T1abf),sep="."),paste("t2",names(T2abf),sep="."),paste("t3",names(T3abf),sep="."))

t1t2t3 <- T1T2T3mods[,"t1.logABF"]+T1T2T3mods[,"t2.logABF"] ++T1T2T3mods[,"t3.logABF"] 
out123 <- cbind(t1t2t3.logbf=t1t2t3,T1T2T3mods)


return(list(out1=out1,out2=out2,out3=out3,out123=out123))
}
###

T1T2T3bfexp.best.fn <- function(T1T2T3abf,T1bfexp,T2bfexp,T3bfexp,indbest,logBF0) {
# on expanded marginal models, re-calculate joint bfs at "best" joint models so have joint bfs on expanded joint models

outkeep <- T1T2T3abf[indbest,]
names(T1bfexp) <- c("et1.logbf","t1.str","et1.size","et1.str","et1.rank")
b1 <- inner_join(outkeep,T1bfexp) # merge on t1.str to give all t1  expanded models 
names(T2bfexp) <- c("et2.logbf","t2.str","et2.size","et2.str","et2.rank")
b12 <- inner_join(b1,T2bfexp) # merge on t2.str
names(T3bfexp) <- c("et3.logbf","t3.str","et3.size","et3.str","et3.rank")
b123 <- inner_join(b12,T3bfexp) # merge on t2.str


bf.all <- b123[,c("et1.logbf","et1.str", "et1.size","et2.logbf","et2.str", "et2.size","et3.logbf","et3.str", "et3.size")]
names(bf.all) <- c("t1.logbf","t1.str", "t1.size","t2.logbf","t2.str", "t2.size","t3.logbf","t3.str", "t3.size")
    
bf.all$t1t2t3.logbf <- bf.all$t1.logbf +bf.all$t2.logbf|+bf.all$t3.logbf+logBF0
bf.all$t1t2t3.Nexpmod <- 1

out <- bf.all[order(bf.all$t1t2t3.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0
check0 <- which(out$t3.size==0)
if(length(check0)>0) out[check0,"t3.str"] <- 0

return(out)
}




###

T1T2T3abfexp.fn <- function(T1abf,T2abf,T3abf,logBF0) {

nT1 <- dim(T1abf)[1]
nT2 <- dim(T2abf)[1]
nT3 <- dim(T3abf)[1]

T1modsrep <- matrix(rep(t(T1abf),nT2),ncol=ncol(T1abf),byrow=TRUE)
T2modsrep <- matrix(rep(as.matrix(T2abf),each=nT1),ncol=ncol(T2abf),byrow=FALSE)
T1T2mods <- cbind(T1modsrep,T2modsrep)

colnames(T1T2mods) <- c(paste("t1",names(T1abf),sep="."),paste("t2",names(T2abf),sep="."))
T1T2mods <- as.data.frame(T1T2mods,stringsAsFactors =FALSE)

nT1T2 <- dim(T1T2mods)[1]

#T1T2modsrep <- big.matrix(nrow=nT1T2*nT3,ncol=ncol(T1T2mods))
#a <- seq(1,(nT1T2*nT3),by=nT1T2)
#b <- seq((nT1T2+1),(nT1T2*nT3),by=nT1T2)

T1T2modsrep <- matrix(rep(t(T1T2mods),nT3),ncol=ncol(T1T2mods),byrow=TRUE)

T3modsrep <- matrix(rep(as.matrix(T3abf),each=nT1T2),ncol=ncol(T3abf),byrow=FALSE)
T1T2T3mods <- cbind(T1T2modsrep,T3modsrep) 

colnames(T1T2T3mods) <- c(paste("t1",names(T1abf),sep="."),paste("t2",names(T2abf),sep="."),paste("t3",names(T3abf),sep="."))

T1T2T3mods <- as.data.frame(T1T2T3mods,stringsAsFactors =FALSE)

T1T2T3mods$t1.lBF <- as.numeric(T1T2T3mods$t1.lBF)
T1T2T3mods$t2.lBF <- as.numeric(T1T2T3mods$t2.lBF)
T1T2T3mods$t3.lBF <- as.numeric(T1T2T3mods$t3.lBF)


t1t2t3 <- T1T2T3mods$t1.lBF + T1T2T3mods$t2.lBF + T1T2T3mods$t3.lBF + logBF0
out123 <- data.frame(t1t2T3.logbf=t1t2t3,T1T2T3mods)
names(out123) <- c("t1t2t3.logbf", "t1.str"  , "t1.tag"  ,   "t1.logbf" ,"t1.size"  ,  "t2.str" ,"t2.tag"   ,  "t2.logbf","t2.size", "t3.str"  , "t3.tag"  ,   "t3.logbf" ,"t3.size")


out123$t1t2t3.Nexpmod <- 1
return(out123)
}






###

prior.fn <- function(k,bf,shared=1,s,mT1=3,mT2=3) {
#  s <- bf[k,"s"]
  m1snp <- unlist(strsplit(as.character(bf$t1.str[k]),"%"))
  m2snp <- unlist(strsplit(as.character(bf$t2.str[k]),"%"))
  m12snp <- intersect(m1snp,m2snp)
  nT1 <- as.numeric(bf$t1.size[k])
  nT2 <- as.numeric(bf$t2.size[k])
  nT1T2 <- length(m12snp)
  #print(c(nT1,nT2,nT1T2))
  prior12 <- shared*(nT1T2 > 0) + 1*(nT1T2 == 0)   # >1 if at least one overlap, 1 otherwise
  #prior1 <- dbinom(nT1,size=s,prob=mT1/s)
  #prior2 <- dbinom(nT2,size=s,prob=mT2/s)
  prior1 <- dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
  prior2 <- dbinom(nT2,size=s,prob=mT2/s)/choose(s,nT2)
  p <- prior1*prior2*prior12
  return(p)
  }

adjprior.fn <- function(prior,out) {
 adjprior <- prior
 ind10 <- which(out$t2.size==0 & out$t1.size>0)
 ind01 <- which(out$t1.size==0 & out$t2.size>0)
 ind00 <- which(out$t2.size==0 & out$t1.size==0)
 ind11 <- which(out$t2.size>0 & out$t1.size>0)
 den <- sum(prior[ind11])
 den0 <- 0
 if(length(ind10) >0) den0 <- den0 + sum(prior[ind10])
 if(length(ind01) >0) den0 <- den0 + sum(prior[ind01])
 if(length(ind00) >0) den0 <- den0 + prior[ind00]
 adjprior[ind11] <- (1-den0)*prior[ind11]/den
 return(adjprior) 
}

##
PP.fn <- function(out,shared,s=length(snps(tags)),mT1=3,mT2=3,details=TRUE) {
  k.mat <- matrix(1:dim(out)[1],ncol=1)
  prior <- apply(k.mat,1,prior.fn,bf=out,shared=shared,s=s,mT1=mT1,mT2=mT2)
  if(min(out$t1.size)>0 & min(out$t2.size)>0 ) {
  adjprior <- prior/sum(prior*out$t1t2.Nexpmod)  
  } else { adjprior <- adjprior.fn(prior,out)
  }
  post <- exp(log(adjprior)+out[,"t1t2.logbf"]-max(out[,"t1t2.logbf"]))
  ppost <- post/sum(post*out$t1t2.Nexpmod)
  outPP <- matrix(ppost,nrow=1,dimnames=list(NULL,paste(out$t1.str,out$t2.str,sep="NEXT")))
 # indC <- order(colnames(tmp))
  #outPP <- tmp[,indC]
if(details) outPP <- data.frame(PP=ppost,prior=adjprior,logprior=log(adjprior),logbf=out[,"t1t2.logbf"],t1.mod=out$t1.str, t2.mod=out$t2.str,t1.logbf=out$t1.logbf,t2.logbf=out$t2.logbf,t1.size=out$t1.size,t2.size=out$t2.size)
  

  #ind <- order(ppost,decreasing=TRUE)
  #outpp <- data.frame(ppost[ind],logABF=log10(exp(out[ind,1])),out[ind,-1])  
  #outpp <- data.frame(ppost[ind],logABF=log10(exp(out[ind,1])),out[ind,-1]) 
  return(outPP)
}
###

prior12.fn <- function(k,bf,shared=1,s,mT1=3,mT2=3) {
  m1snp <- unlist(strsplit(as.character(bf$out1$mod[k]),"%"))
  m2snp <- unlist(strsplit(as.character(bf$out2$mod[k]),"%"))
  m12snp <- intersect(m1snp,m2snp)
  nT1 <- length(m1snp)
  nT2 <- length(m2snp)
  nT1T2 <- length(m12snp)
  #print(c(nT1,nT2,nT1T2))
  prior12 <- shared*(nT1T2 > 0) + 1*(nT1T2 == 0)   # >1 if at least one overlap, 1 otherwise
  prior1 <- dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
  prior2 <- dbinom(nT2,size=s,prob=mT2/s)/choose(s,nT2)
  #prior1 <- dbinom(nT1,size=s,prob=mT1/s)
  #prior2 <- dbinom(nT2,size=s,prob=mT2/s)
  p <- prior1*prior2*prior12
  return(p)
  }

##
PP12.fn <- function(out,shared,s=length(snps(tags)),mT1=3,mT2=3) {
  k.mat <- matrix(1:dim(out$out12)[1],ncol=1)
  prior <- apply(k.mat,1,prior12.fn,bf=out,shared=shared,s=s,mT1=mT1,mT2=mT2)
  adjprior <- prior/sum(prior)
  post <- exp(log(adjprior)+out$out12[,"t1t2.logbf"]-max(out$out12[,"t1t2.logbf"]))
  ppost <- post/sum(post)
  outPP <- data.frame(PP=ppost,prior=adjprior,logprior=log(adjprior),logbf=out$out12[,"t1t2.logbf"],t1.mod=out$out1$mod, t2.mod=out$out2$mod,t1.logbf=out$out1$logbf,t2.logbf=out$out2$logbf,t1.size=out$out1$size,t2.size=out$out2$size)
  return(outPP)
}
####

####

prior.bin.fn <- function(k,bf,s,mT1=3) {
#  s <- bf[k,"s"]
  m1snp <- unlist(strsplit(as.character(bf$str[k]),"%"))
   nT1 <- length(m1snp)
   prior1 <- dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
   return(prior1)
  }

##
PP.bin.fn <- function(out,s=length(snps(tags)),mT1=3) {
  k.mat <- matrix(1:dim(out)[1],ncol=1)
  prior <- apply(k.mat,1,prior.bin.fn,bf=out,s=s,mT1=mT1)
  adjprior <- prior/sum(prior)
  post <- exp(log(adjprior)+out[,1]-max(out[,1]))
  ppost <- post/sum(post)
  tmp <- matrix(ppost,nrow=1,dimnames=list(NULL,out$str))
  indC <- order(tmp[1,],decreasing=TRUE)
  outPP <- tmp[,indC]

  #ind <- order(ppost,decreasing=TRUE)
  #outpp <- data.frame(ppost[ind],logABF=log10(exp(out[ind,1])),out[ind,-1])  
  #outpp <- data.frame(ppost[ind],logABF=log10(exp(out[ind,1])),out[ind,-1]) 
  return(outPP)
}

###


PP123.fn <- function(out,shared,s=length(snps(tags)),mT1=3,mT2=3,mT3=3,details=TRUE) {
  k.mat <- matrix(1:dim(out)[1],ncol=1)
  prior <- apply(k.mat,1,prior3T.fn,bf=out,shared=shared,s=s,mT1=mT1,mT2=mT2,mT3=mT3)
  adjprior <- prior/sum(prior*out$t1t2t3.Nexpmod)
  post <- exp(log(adjprior)+out[,"t1t2t3.logbf"]-max(out[,"t1t2t3.logbf"]))
  ppost <- post/sum(post*out$t1t2t3.Nexpmod)  
  if(!details) outPP <- matrix(ppost,nrow=1,dimnames=list(NULL,paste(out$t1.str,out$t2.str,out$t3.str,sep="NEXT")))
  if(details) outPP <- data.frame(PP=ppost,prior=adjprior,logprior=log(adjprior),logbf=out[,"t1t2t3.logbf"],t1.mod=out$t1.str, t2.mod=out$t2.str,t3.mod=out$t3.str,t1.logbf=out$t1.logbf,t2.logbf=out$t2.logbf,t3.logbf=out$t3.logbf,t1.size=out$t1.size,t2.size=out$t2.size,t3.size=out$t3.size)
  return(outPP)
}
#

###
prior3T.fn <- function(k,bf,shared=1,s,mT1=3,mT2=3,mT3=3) {
#  s <-length(snps(tags))
# mTi = Expected no. causals
  m1snp <- unlist(strsplit(as.character(bf$t1.str[k]),"%"))
  m2snp <- unlist(strsplit(as.character(bf$t2.str[k]),"%"))
  m3snp <- unlist(strsplit(as.character(bf$t3.str[k]),"%"))
  
  m12snp <- intersect(m1snp,m2snp)
  m13snp <- intersect(m1snp,m3snp)
  m23snp <- intersect(m2snp,m3snp)
   
  n1 <- length(m1snp)
  n2 <- length(m2snp)
  n3 <- length(m3snp)

  K=3
  m <- choose(K,2)
  n12 <- 1*(length(m12snp)>0)
  n13 <- 1*(length(m13snp)>0)
  n23 <- 1*(length(m23snp)>0)
  
  C <- (n12+n13+n23)/m
 
  prior1 <- dbinom(n1,size=s,prob=mT1/s)/choose(s,n1)
  prior2 <- dbinom(n2,size=s,prob=mT2/s)/choose(s,n2)
  prior3 <- dbinom(n3,size=s,prob=mT3/s)/choose(s,n3)
  p <- prior1*prior2*prior3*shared^C
  return(p)
  }
###





mod.split.fn <- function(k,PP,tnum) {
# when two traits, tnum is 1 or 2 , i.e. 1st or 2nd trait
# called by PP.marg.fn
tmp <- unlist(strsplit(row.names(PP)[k],"NEXT"))
msnp <- tmp[tnum]
if(is.na(msnp)) msnp <- "0" # null model for trait tnum
return(msnp)
}


## after joint models are split into models for each trait, combine common model PP to get PP for that model; called by PP.marg.fn
 mergePP.fn <- function(pp) {
     mnames <- unique(pp$mod)
     m<- dim(pp)[2]
     g <- length(mnames)
     ns <- m-2 # number of values for sharing scale (1st col is model, last col is number of models)
     p1 <- matrix(0,nrow=g,ncol=ns,dimnames=list(mnames,names(pp)[-c(1,m)])) 
     
     for(j in 1:g) { 
     ind1 <- which(pp$mod == mnames[j])
     tmp <- pp[ind1,-c(1,m)]*pp[ind1,m] 
     p1[j,] <- apply(tmp,2,sum)      
     }
    
     return(p1)
    	}


# PP for each trait,splitting joint models
#PP.marg.fn <-function(PP) {
 
# mod <- apply(matrix(1:dim(PP)[1],ncol=1),1,mod.split.fn,PP,1)
#  pp1 <- data.frame(mod,PP,row.names=NULL)
#  mod <- apply(matrix(1:dim(PP)[1],ncol=1),1,mod.split.fn,PP,2)
#  pp2 <- data.frame(mod,PP,row.names=NULL)
    	
#   PP1 <- mergePP.fn(pp1)  
#   PP2 <- mergePP.fn(pp2)  

#spp1 <- PP1[order(PP1[,1],decreasing=TRUE),]
#spp2 <- PP2[order(PP2[,1],decreasing=TRUE),]
               
#return(list(PP1=spp1,PP2=spp2))
#}

#

PP.marg.fn <-function(PPall,K) {
 m<- dim(PPall)[2]
 PP <- PPall[,-m] # last column is no. of joint models with same PP 
PPm <- vector("list",K) 
 for(k in 1:K) {
  mod <- apply(matrix(1:dim(PP)[1],ncol=1),1,mod.split.fn,PP,k)
  tmp <- data.frame(mod,PP,Nmod=PPall[,m],row.names=NULL)
  pp <-  mergePP.fn(tmp)
  PPm[[k]] <- pp[order(pp[,1],decreasing=TRUE),]
    	}
          
return(PPm)
}

###

MPP.fn<-function(PP1) {
 mnames <- rownames(PP1)
 sep.fn <- function(k,mnames) {
   msep <- unlist(strsplit(as.character(mnames[k]),"%"))
   return(msep)
   }
  msep <- apply(matrix(1:length(mnames),ncol=1),1,sep.fn,mnames)
 
 #out <- matrix(PP,nrow=1,dimnames=list(NULL,mod.names))
 #indC <- order(colnames(out))
 #outPP <- matrix(out[indC],nrow=1,dimnames=list(NULL,colnames(out)[indC]))
 
 #out <- outPP
 #out<-PP
 # mnames <- colnames(out)
 # gnames <- c(LETTERS[1:6],"sAITD","sAA","sUC","sRA","G8","G11","G12")
  gnames <- unique(unlist(msep)) # snps 
  
  
 
  check.fn <- function(k,msep,out,gnames) {
     g <- length(gnames)
     p1 <- numeric(g) 
     for(j in 1:g) { 
     ind1 <- gnames[j] %in% msep[[k]]
     if(ind1) p1[j] <- out[k] 
     }
     return(p1)
    	}
   
   mpp1 <- NULL
   for(k in 1:dim(PP1)[2]) {
    tmp1 <- apply(matrix(1:length(mnames),ncol=1),1,check.fn,msep,PP1[,k],gnames)  
    mpp1 <- rbind(mpp1,apply(tmp1,1,sum) )
    }
    
   mpp1 <- data.frame(mpp1,row.names=colnames(PP1))  
   names(mpp1)<-gnames    
return(t(mpp1))
}



###
expand.tags.bf <- function(best, tags) {
   
    bsnps <- unique(unlist(strsplit(as.character(best$mod),"%")))    
    check <- which(bsnps=="0" | bsnps=="1") 
    if(length(check)>0) bsnps <- bsnps[-check]
    B <- dim(best)[1]
    wh <- which(make.names(tags(tags)) %in% bsnps)
           
    
    if (!length(wh)) 
        stop("none of the supplied tags are amongst the best SNPs in d")
    proxies <- split(make.names(snps(tags)[wh]), make.names(tags(tags)[wh]))
   
    if (!all(bsnps %in% names(proxies))) 
        stop("not all model SNPs found in tags object")
    message("expanding tags for ", B, " models over ", length(proxies), 
        " tag SNPs, tagging a total of ", length(unlist(proxies)), 
        " SNPs.")
    pm <- mclapply(as.list(1:B), function(i) {
        if (best[i, "size"] == 0) {
            pm.str <- ""
        }
        else {
            pm.str <- apply(do.call(expand.grid, proxies[unlist(strsplit(as.character(best$mod[i]),"%"))]), 
                1, makestr)
        }
        gc()
        
        return(pm.str)
    })
    npm <- sapply(pm, length)
    neighb <- as.data.table(best)[rep(1:length(pm), times = npm), 
        ]
    setnames(neighb, sub("str", "index.str", names(neighb)))
    neighb[, `:=`(str, unlist(pm))]
    #neighb[, `:=`(index, str == index.str)]
   # neighb[, `:=`(logPP, logPP - logsum(logPP))]
# neighb[, `:=`(PP, exp(neighb$logPP))]
#    neighb <- neighb[order(neighb$PP, decreasing = TRUE), ]
    neighb <- neighb[order(neighb$logbf, decreasing = TRUE), ]
    neighb[, `:=`(rank, 1:nrow(neighb))]
    best <- as.data.frame(neighb)
   # d@model.snps <- strsplit(neighb$str, "%")
   # return(marg.snps(d))
   return(best)
}

###
modmat.fn  <- function(k,mods,msnps) {
  ind <- match(mods[[k]],msnps)
  out <- numeric(length(msnps))
  out[ind] <- 1
  return(out)
  }
  
####

mod2matrix.fn <- function(best1) {

 msnps <- unique(unlist(strsplit(best1$str,"%"))) 
 mods1 <- best1$str
 
 Nmod <- length(mods1)
 Nsnp <- length(msnps)
 
 mods <- strsplit(mods1,"%")
  
T1mods <- t(apply(matrix(1:Nmod,ncol=1),1,modmat.fn,mods,msnps))
colnames(T1mods) <- msnps
return(T1mods)
}

####

T1T2bfexp.fn <- function(T1T2abf,T1bfexp,T2bfexp) {
# use expanded marginal models, to get expanded joint models



out1keep <- T1T2abf$out1 
 out2keep <- T1T2abf$out2


 bfall <- c()

k.mat <- matrix(1:dim(T1T2abf$out12)[1],ncol=1)


exp.fn <- function(k,T1bfexp,T2bfexp,out1keep,out2keep) {
#  for(k in 1:dim(T1T2abf$out12)[1]) {   
    ind1 <- which(apply(as.matrix(T1bfexp[,"mod"]), 1, identical, out1keep[k,"mod"]))
    tmp1 <- T1bfexp[ind1,c("logbf","str","size")]
    
    ind2 <- which(apply(as.matrix(T2bfexp[,"mod"]), 1, identical, out2keep[k,"mod"]))
    tmp2 <- T2bfexp[ind2,c("logbf","str","size")]
    
    T1modsrep <- matrix(rep(t(tmp1),dim(tmp2)[1]),ncol=ncol(tmp1),byrow=TRUE)
	T2modsrep <- matrix(rep(as.matrix(tmp2),each=dim(tmp1)[1]),ncol=ncol(tmp2),byrow=FALSE)
	tmp <- cbind(T1modsrep,T2modsrep)
	T1T2mods <- data.frame(tmp)
	names(T1T2mods) <- c("t1.logbf","t1.str", "t1.size","t2.logbf","t2.str", "t2.size")
   return(T1T2mods)
#    bfall <- rbind(bfall, T1T2mods)
    }

bfall <- apply(k.mat,1,exp.fn,T1bfexp,T2bfexp,out1keep,out2keep)
bf.all <- do.call(rbind,bfall)

bf.all$t1.logbf <- as.numeric(bf.all$t1.logbf)
bf.all$t2.logbf <- as.numeric(bf.all$t2.logbf)
    
t1t2bf <- bf.all$t1.logbf +bf.all$t2.logbf

bf.all$t1t2.logbf <- t1t2bf

out <- bf.all[order(bf.all$t1t2.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0

return(out)
}

####

T1T2T3bfexp.fn <- function(T1T2T3abf,T1bfexp,T2bfexp,T3bfexp) {
# on expanded marginal models, re-calculate joint bfs at "best" joint models so have joint bfs on expanded joint models; on "not best"  models
# keep bf of model for tag snps

 #out1keep <- T1T2abf$out1[indbest,] # re-fit these
 #out2keep <- T1T2abf$out2[indbest,]

out1keep <- T1T2T3abf$out1 # re-fit these
 out2keep <- T1T2T3abf$out2
out3keep <- T1T2T3abf$out3

 bfall <- c()

k.mat <- matrix(1:dim(T1T2abf$out12)[1],ncol=1)

exp.fn <- function(k,T1bfexp,T2bfexp,out1keep,out2keep) {
#  for(k in 1:dim(T1T2abf$out12)[1]) {   
    ind1 <- which(apply(as.matrix(T1bfexp[,"mod"]), 1, identical, out1keep[k,"mod"]))
    tmp1 <- T1bfexp[ind1,c("logbf","str","size")]
    
    ind2 <- which(apply(as.matrix(T2bfexp[,"mod"]), 1, identical, out2keep[k,"mod"]))
    tmp2 <- T2bfexp[ind2,c("logbf","str","size")]
    
    T1modsrep <- matrix(rep(t(tmp1),dim(tmp2)[1]),ncol=ncol(tmp1),byrow=TRUE)
	T2modsrep <- matrix(rep(as.matrix(tmp2),each=dim(tmp1)[1]),ncol=ncol(tmp2),byrow=FALSE)
	tmp <- cbind(T1modsrep,T2modsrep)
	T1T2mods <- data.frame(tmp)
	names(T1T2mods) <- c("t1.logbf","t1.str", "t1.size","t2.logbf","t2.str", "t2.size")
   return(T1T2mods)
#    bfall <- rbind(bfall, T1T2mods)
    }

bfall <- apply(k.mat,1,exp.fn,T1bfexp,T2bfexp,out1keep,out2keep)
bf.all <- do.call(rbind,bfall)

bf.all$t1.logbf <- as.numeric(bf.all$t1.logbf)
bf.all$t2.logbf <- as.numeric(bf.all$t2.logbf)
    
t1t2bf <- bf.all$t1.logbf +bf.all$t2.logbf

bf.all$t1t2.logbf <- t1t2bf

out <- bf.all[order(bf.all$t1t2.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0

return(out)
}


###########


make.snp.groups.fn <- function() {
#load("/scratch/wallace/IL2RA/full-data/snp-picker-grouped.RData")
#allgroups <- as(gSP2,"groups") 
#load("/home/ja628/causal_sim/results/IL2RA/real-data/new/snpgroups.RData")
load("/home/ja628/causal_sim/results/IL2RA/real-data/new/newgroups.RData")

#allgroups <- groupsx
#Gtmp <- allgroups@.Data
Gtmp <- newg
snpGroups <- list(A=Gtmp[[1]],B=Gtmp[[5]],C=Gtmp[[2]],D=Gtmp[[6]],E=Gtmp[[8]],F=Gtmp[[7]],S=Gtmp[[3]],G11=Gtmp[[4]])

#snpGroups <- list(A=Gtmp[[4]],B=Gtmp[[14]],C=c(Gtmp[[6]],Gtmp[[13]]),D=Gtmp[[15]],E=Gtmp[[1]],F=Gtmp[[2]],S=Gtmp[[5]],G11=Gtmp[[3]],sAITD=c(Gtmp[[11]],"rs706779.6098824.T.C"),S1=Gtmp[[7]],S3=Gtmp[[8]],S15=Gtmp[[9]],S2=Gtmp[[10]],S1.1=Gtmp[[12]],S2.1=Gtmp[[16]],S12=Gtmp[[17]],S5=Gtmp[[18]])
return(snpGroups)
}

## plots

pp.filter.fn <- function(pp,pthr) {
 # rows=kappa; cols=snps
 m <- dim(pp)[1]
 ind <- NULL
 for(k in 1:m) {ind <- union(ind,which(pp[k,]>=pthr)); print(ind)}
 out <- pp[,ind]
 return(out)
 }

MPP.plot.fn <- function(mpp1,shared,pthr) { 
 if(!is.na(pthr)) mpp1 <- pp.filter.fn(mpp1,pthr=pthr)
 m <- dim(mpp1)[1]
 p1 <- dim(mpp1)[2]  	
 par(las=1) 
 #par(cex.axis=0.8) 
 par(cex.axis=0.5) 
 image(1:m,1:p1,mpp1,xlab="",ylab="",col=rev(heat.colors(128)),zlim=c(0,1.01),axes=FALSE)
 axis(2,at=1:p1,labels=colnames(mpp1))
 box()
 axis(1,at=1:m,shared)
 image.plot(1:m,1:p1,mpp1,xlab="",ylab="",col=rev(heat.colors(128)),zlim=c(0,1.01),legend.only=TRUE )
 return(mpp1)
 		      	       		 }


##

t1t2plot.fn <- function(alltraits,shared,pthr) {
#  ti <- c(t1,t2)
 ti=alltraits
  traits <- paste(alltraits[1],"-",alltraits[2],sep="")
 k=alltraits[1]
  fname <- file.path(fdir,paste(k,"-MPP-shared-",traits,".txt",sep=""))
  MPP1 <- t(read.table(fname,header=TRUE,row.names=1))
  mpp1 <- pp.filter.fn(MPP1,pthr=pthr)
 k=alltraits[2] 
  fname <- file.path(fdir,paste(k,"-MPP-shared-",traits,".txt",sep=""))
  MPP2 <- t(read.table(fname,header=TRUE,row.names=1))
  mpp2 <- pp.filter.fn(MPP2,pthr=pthr)
  MPP <- smartbind(mpp1,mpp2,fill=0) 
 
 if(sum(ti == "GRAVES")>0) ti[which(ti=="GRAVES")]<-"AITD"
rownames(MPP) <- c(paste(ti[1],rownames(mpp1),sep="."),paste(ti[2],rownames(mpp2),sep="."))
colnames(MPP) <- unlist(strsplit(colnames(MPP),"[.]"))[c(TRUE,FALSE,FALSE,FALSE)]
pname <- file.path(fdir,paste("MPP-",traits,".pdf",sep=""))
pdf(pname)
par(mfrow=c(1,2))
mpp1<-MPP.plot.fn(as.matrix(MPP[grep(ti[1],rownames(MPP)),]),shared,pthr=NA)
title(ti[1])
mpp2<-MPP.plot.fn(as.matrix(MPP[grep(ti[2],rownames(MPP)),]),shared,pthr=NA)
title(ti[2])
dev.off()
return(list(mpp1=mpp1,mpp2=mpp2))
}