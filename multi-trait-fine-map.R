library(annotSnpStats)
library(snpStats)
library(GUESSFM)
library(mlogitBMA)
library(BMA)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(data.table)




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
 if(mT1==0) mc <- (mT1+1):MT1
 modT1 <- vector("list",length(mc))
  for(i in mc) {
   mci <- combn(msnps,i,simplify=TRUE)  # all combinations of i msnps
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

ind <- order(logABF,decreasing=TRUE)
out <- data.frame(logABF=logABF[ind],M=mod1$models[ind,],row.names=NULL)
cnames <- c("logABF",names(data1)[-1])
names(out) <- cnames
out$size <- apply(out[,-1],1,sum)


write.table(out,fabf,row.names=FALSE,col.names=TRUE,quote=FALSE)
return(out)
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





###

prior.fn <- function(k,bf,shared=1,s,mT1=3,mT2=3) {
#  s <- bf[k,"s"]
  m1snp <- unlist(strsplit(as.character(bf$t1.str[k]),"%"))
  m2snp <- unlist(strsplit(as.character(bf$t2.str[k]),"%"))
  m12snp <- intersect(m1snp,m2snp)
  nT1 <- length(m1snp)
  nT2 <- length(m2snp)
  nT1T2 <- length(m12snp)
  #print(c(nT1,nT2,nT1T2))
  prior12 <- shared*(nT1T2 > 0) + 1*(nT1T2 == 0)   # >1 if at least one overlap, 1 otherwise
  prior1 <- dbinom(nT1,size=s,prob=mT1/s)
  prior2 <- dbinom(nT2,size=s,prob=mT2/s)
  p <- prior1*prior2*prior12
  return(p)
  }

##
PP.fn <- function(out,shared,s=length(snps(tags)),mT1=3,mT2=3) {
  k.mat <- matrix(1:dim(out)[1],ncol=1)
  prior <- apply(k.mat,1,prior.fn,bf=out,shared=shared,s=s,mT1=mT1,mT2=mT2)
  post <- exp(log(prior)+out[,"t1t2.logbf"]-max(out[,1]))
  ppost <- post/sum(post)
  tmp <- matrix(ppost,nrow=1,dimnames=list(NULL,paste(out$t1.str,out$t2.str,sep="NEXT")))
  indC <- order(colnames(tmp))
  outPP <- tmp[,indC]

  #ind <- order(ppost,decreasing=TRUE)
  #outpp <- data.frame(ppost[ind],logABF=log10(exp(out[ind,1])),out[ind,-1])  
  #outpp <- data.frame(ppost[ind],logABF=log10(exp(out[ind,1])),out[ind,-1]) 
  return(outPP)
}


###

mod.split.fn <- function(k,PP,tnum) {
# when two traits, tnum is 1 or 2 , i.e. 1st or 2nd trait
# called by PP.marg.fn
tmp <- unlist(strsplit(row.names(PP)[k],"NEXT"))
msnp <- tmp[tnum]
return(msnp)
}


## after joint models are split into models for each trait, combine common model PP to get PP for that model; called by PP.marg.fn
 mergePP.fn <- function(pp) {
     mnames <- unique(pp$mod)
     g <- length(mnames)
     ns <- dim(pp)[2]-1 # number of values for sharing scale
     p1 <- matrix(0,nrow=g,ncol=ns,dimnames=list(mnames,names(pp)[-1])) 
     for(j in 1:g) { 
     ind1 <- which(pp$mod == mnames[j])
     p1[j,] <- apply(pp[ind1,-1],2,sum) 
     }
     return(p1)
    	}


# PP for each trait,splitting joint models
PP.marg.fn <-function(PP) {
 
  mod <- apply(matrix(1:dim(PP)[1],ncol=1),1,mod.split.fn,PP,1)
  pp1 <- data.frame(mod,PP,row.names=NULL)
  mod <- apply(matrix(1:dim(PP)[1],ncol=1),1,mod.split.fn,PP,2)
  pp2 <- data.frame(mod,PP,row.names=NULL)
    	
   PP1 <- mergePP.fn(pp1)  
   PP2 <- mergePP.fn(pp2)  

spp1 <- PP1[order(PP1[,1],decreasing=TRUE),]
spp2 <- PP2[order(PP2[,1],decreasing=TRUE),]
               
return(list(PP1=spp1,PP2=spp2))
}

###

MPP.fn<-function(PP) {
 mnames <- names(PP)
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
  gnames <- c(LETTERS[1:6],"sAITD","sAA","sUC","sRA","G8","G11","G12")
  g <- length(gnames)
  gnames1 <- paste("z_1.",gnames,sep="")
  gnames2 <- paste("z_2.",gnames,sep="")
 
  check.fn <- function(k,msep,out,gnames1) {
     p1 <- numeric(g) 
     for(j in 1:g) { 
     ind1 <- gnames1[j] %in% msep[[k]]
     if(ind1) p1[j] <- out[k] 
     }
     return(p1)
    	}
    tmp1 <- apply(matrix(1:length(mnames),ncol=1),1,check.fn,msep,PP,gnames1)  
    mpp1 <- apply(tmp1,1,sum) 
    tmp2 <- apply(matrix(1:length(mnames),ncol=1),1,check.fn,msep,PP,gnames2)  
    mpp2 <- apply(tmp2,1,sum) 
    
   mpp <- data.frame(p1=mpp1,p2=mpp2,row.names=gnames)      
return(t(mpp))
}




###
expand.tags.bf <- function(best, tags) {
   
    bsnps <- unique(unlist(strsplit(as.character(best$mod),"%")))    
    check <- which(bsnps=="0") 
    if(length(check)>0) bsnps <- bsnps[-check]
    B <- dim(best)[1]
    wh <- which(make.names(tags(tags)) %in% bsnps)
         
   # t2n1 <- setdiff(make.names(tags(tags2)),make.names(tags(tags1)))
    #tags2n1 <- taggedby(tags2,t2n1)
    
    
    if (!length(wh)) 
        stop("none of the supplied tags are amongst the best SNPs in d")
    proxies <- split(make.names(snps(tags)[wh]), make.names(tags(tags)[wh]))
   # proxies2 <- split(make.names(snps(tags2)[wh]), make.names(tags(tags2)[wh]))
 
 # if (!all(bsnps %in% names(proxies)))    {
 #   ind <- which(!(bsnps %in% names(proxies)))
 #   proxnew <- as.character(tagsof(tags,bsnps[ind])[,1])
    #bsnps[ind] <- proxnew
  #  wh <- which(make.names(tags(tags)) %in% bsnps)
  #  proxies <- split(make.names(snps(tags)[wh]), make.names(tags(tags)[wh]))
 #}
 
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

T1T2bfexp.fn <- function(T1T2abf,T1bfexp,T2bfexp,indbest) {
# on expanded marginal models, re-calculate joint bfs so have joint bfs on expanded joint models
 out1keep <- T1T2abf$out1[indbest,]
 out2keep <- T1T2abf$out2[indbest,]

 bfall <- c()

  for(k in 1:length(indbest)) {   
    ind1 <- which(apply(as.matrix(T1bfexp[,"mod"]), 1, identical, out1keep[k,"mod"]))
    tmp1 <- T1bfexp[ind1,c("logbf","str","size")]
    
    ind2 <- which(apply(as.matrix(T2bfexp[,"mod"]), 1, identical, out2keep[k,"mod"]))
    tmp2 <- T2bfexp[ind2,c("logbf","str","size")]
    
    T1modsrep <- matrix(rep(t(tmp1),dim(tmp2)[1]),ncol=ncol(tmp1),byrow=TRUE)
	T2modsrep <- matrix(rep(as.matrix(tmp2),each=dim(tmp1)[1]),ncol=ncol(tmp2),byrow=FALSE)
	tmp <- cbind(T1modsrep,T2modsrep)
	T1T2mods <- data.frame(tmp)
	names(T1T2mods) <- c("t1.logbf","t1.str", "t1.size","t2.logbf","t2.str", "t2.size")
   
    bfall <- rbind(bfall, T1T2mods)
    }

bfall$t1.logbf <- as.numeric(bfall$t1.logbf)
bfall$t2.logbf <- as.numeric(bfall$t2.logbf)
    
t1t2bf <- bfall$t1.logbf +bfall$t2.logbf

bfall$t1t2.logbf <- t1t2bf
out <- bfall[order(bfall$t1t2.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0

return(out)
}

###########
