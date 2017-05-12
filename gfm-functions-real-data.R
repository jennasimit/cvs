library(annotSnpStats)
library(GUESSFM)

# snp groups 

make.snp.groups.fn <- function() {
#load("/scratch/wallace/IL2RA/snp-picker-grouped.RData")
load("/scratch/wallace/IL2RA/full-data/snp-picker-grouped.RData")
allgroups <- as(gSP2,"groups") 
Gtmp <- allgroups@.Data
d="/scratch/wallace/IL2RA/full-data/RA" 
dfile <- file.path(d,"data.RData") 
tfile <- file.path(d,"tags.RData")
(load(dfile)) 
(load(tfile)) 

Gtmp <- allgroups@.Data
Gd <-  c("rs3118475","rs62626317","rs41295055","rs41295079","rs62626325","rs41295105","rs56382813")

gd <- c()
for(k in 1:length(Gd)) {
 snps <- colnames(snp.data)
 ind <- grep(Gd[k],snps)
 if(length(ind)>0) gd <- c(gd,snps[ind])
 		      }

gb <- c(snps[grep("rs2104286",snps)],snps[grep("rs12722488",snps)],snps[grep("rs12722504",snps)])

snpGroups <- list(A=Gtmp[[1]],B=gb,C=c(Gtmp[[3]],Gtmp[[6]][1:4]),D=gd,E=Gtmp[[5]],F=Gtmp[[4]],sAITD=Gtmp[[7]],sAA=Gtmp[[11]],sUC=Gtmp[[9]],sRA=Gtmp[[2]],G8=Gtmp[[6]][-(1:4)],G11=Gtmp[[10]],G12=Gtmp[[8]])


#snpGroups <- list(A=Gtmp[[1]],B=Gtmp[[7]],C=c(Gtmp[[3]],Gtmp[[8]][1:4]),D=gd,E=Gtmp[[6]],F=Gtmp[[5]],sAITD=Gtmp[[10]],sAA=Gtmp[[16]],sUC=Gtmp[[14]],sRA=Gtmp[[2]],G4=Gtmp[[4]],G8=Gtmp[[8]][-(1:4)],G9=Gtmp[[9]],G11=Gtmp[[11]],G12=Gtmp[[12]],G13=Gtmp[[13]],G15=Gtmp[[15]],G17=Gtmp[[17]])
#snpGroups <- list(A=Gtmp[[1]],B=Gtmp[[7]],C=Gtmp[[3]],D=gd,E=Gtmp[[6]],F=Gtmp[[5]],sAITD=Gtmp[[10]],sAA=Gtmp[[16]],sUC=Gtmp[[14]],sRA=Gtmp[[2]],G4=Gtmp[[4]],G8=Gtmp[[8]],G9=Gtmp[[9]],G11=Gtmp[[11]],G12=Gtmp[[12]],G13=Gtmp[[13]],G15=Gtmp[[15]],G17=Gtmp[[17]])


return(snpGroups)
}



## group MPP fn and best tag 
Ginfo.fn <- function(G=snpGroups,gfmbestsnp) {
 Ng <- length(G)
 tags <-  numeric(Ng)
 gmpp <-  numeric(Ng)
 for(k in 1:Ng) {
  snps <- G[[k]]
  ind <- match(snps,gfmbestsnp[,3],nomatch=0)
  tmp <- gfmbestsnp[ind,]
  if(dim(tmp)[1]>0) {
  tags[k] <- tmp[which.max(tmp[,"Marg_Prob_Incl"]),3]
  gmpp[k] <- sum(tmp[,"Marg_Prob_Incl"]) 
  					}
  }   
  out <- data.frame(tags=tags,gmpp=gmpp)
  row.names(out) <- names(snpGroups)
  return(out)
}

## group PP fn

summarise.fn <- function(bestmod,snpGroups) {

m <- dim(bestmod)[1]
B <- snpGroups
b <- length(snpGroups)
gnames <- names(snpGroups)

check <- matrix(0,nrow=m,ncol=b)
for(i in 1:m) {
  bsnp <- unlist(strsplit(as.character(bestmod$snps[i]),"%"))
  checkb <- c()
  for(k in 1:b) checkb <- c(checkb, 1*(length(intersect(unlist(B[k]),bsnp))>0))
  check[i,] <- checkb
	}


c2 <- combn(gnames,2,simplify=TRUE)
c2names <- apply(c2,2,name.fn <- function(cc) return(paste(cc,collapse=".")))
c3 <- combn(gnames,3,simplify=TRUE)
c3names <- apply(c3,2,name.fn <- function(cc) return(paste(cc,collapse=".")))
c4 <- combn(gnames,4,simplify=TRUE)
c4names <- apply(c4,2,name.fn <- function(cc) return(paste(cc,collapse=".")))
c5 <- combn(gnames,5,simplify=TRUE)
c5names <- apply(c5,2,name.fn <- function(cc) return(paste(cc,collapse=".")))
c6 <- combn(gnames,6,simplify=TRUE)
c6names <- apply(c6,2,name.fn <- function(cc) return(paste(cc,collapse=".")))



m2 <- matrix(0,nrow=length(c2names),ncol=b); for(i in 1:length(c2names)) { m2[i,match(c2[,i],gnames)] <- 1 }
m3 <- matrix(0,nrow=length(c3names),ncol=b); for(i in 1:length(c3names)) { m3[i,match(c3[,i],gnames)] <- 1 }
m4 <- matrix(0,nrow=length(c4names),ncol=b); for(i in 1:length(c4names)) { m4[i,match(c4[,i],gnames)] <- 1 }
m5 <- matrix(0,nrow=length(c5names),ncol=b); for(i in 1:length(c5names)) { m5[i,match(c5[,i],gnames)] <- 1 }
m6 <- matrix(0,nrow=length(c6names),ncol=b); for(i in 1:length(c6names)) { m6[i,match(c6[,i],gnames)] <- 1 }

m1 <- diag(b)
M <- rbind(m1,m2,m3,m4,m5,m6,rep(0,b))
c1names <- gnames
rownames(M) <- c(c1names,c2names,c3names,c4names,c5names,c6names,"M0")

out <- data.frame(PP=bestmod$PP,check)
names(out) <- c("PP",paste("B",1:b,sep=""))

outPP <- as.data.frame(matrix(0,nrow = 1, ncol = dim(M)[1],dimnames = list(NULL, rownames(M))))

indr <- NULL
ind <- 1:m
r <- 0
indr <- NULL
while(length(ind)>0 & r<dim(outPP)[2])  {
  while(length(indr)==0) {
  r <- r+1; 
  indr <- which(apply(out[,-1], 1, function(x) all(x == M[r,])) )
  rM <- rownames(M)[r] 
			}
 outPP[,rM] <- sum(bestmod$PP[indr])
 ind <- setdiff(ind,indr)
 indr <- NULL
 }

return(outPP[which(outPP>0)])
}


##



# full guessfm analysis
gfm.fn <- function(t1) {
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
 fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 tfile <- file.path(fdir,"tags.RData") 
# tfile <- file.path(f1,"tags.RData") 
 load(tfile)
 X <- snp.data[,snps(tags)]
 mydir <- file.path(fdir,paste(t1,"-gfm",sep=""))
 run.bvs(X=X,Y=cc.ph,tag.r2=NA,nexp=3,nsave=1000,gdir=mydir,wait=TRUE) 
}

gfm.expanded.fn <- function(t1) {
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 tfile <- file.path(f1,"tags.RData") 
 #tfile <- "/home/ja628/causal_sim/results/IL2RA/real-data/new/tags.RData"
 load(dfile)
 load(tfile)
 
 d <- read.snpmod(f1)
 dx <- expand.tags(d,tags)

 fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-dx-gfm-expanded.Rdata",sep="")
 save(dx,file=fname)
}


gfm.by.groups.fn <- function(t1,snpGroups) {
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
 tfile <- file.path(f1,"tags.RData") 
 load(tfile)
 X <- snp.data[,snps(tags)]
 
 fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-dx-gfm-expanded.Rdata",sep="")
 (load(fname))
 best <- best.models(dx,pp.thr=0.0001)
 abf <- abf.calc(y=cc.ph,x=X,models=best$str,family="binomial") 
 sm <- abf2snpmod(abf,expected=3,nsnps=abf$nsnp)
 best2 <- best.models(sm,pp.thr=0.0001) # best models after re-fit
 bestsnps <- best.snps(sm,pp.thr=0)
 
 PPgfm <- summarise.fn(best2,snpGroups)
 fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-gfm-groupPP.txt",sep="")
 write.table(PPgfm,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)
 MPPgfm <- Ginfo.fn(snpGroups,bestsnps)
 fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-gfm-groupMPP.txt",sep="")
 write.table(MPPgfm,fname,quote=FALSE,row.names=TRUE,col.names=TRUE)
}



sel.models.fn <- function(t1) {
 fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 tfile <- file.path(fdir,"tags.RData") 
 load(dfile)
 load(tfile)
 
 tsnp <- snp.data[,unique(tags@tags)] 
 
 mydir <- file.path(fdir,paste(t1,"-gfm-on-tags",sep=""))
 #run.bvs(X=tsnp,Y=cc.ph,tag.r2=NA,nexp=3,nsave=1000,gdir=mydir,wait=TRUE) 

 d <- read.snpmod(mydir) 
 #bestmods <- best.models(d,pp.thr=0.0001) 
 #fname <- file.path(fdir,paste(t1,"-gfm-on-tags-best-models0.0001.RData",sep=""))
 #save(bestmods,file=fname)
 bestsnps <- best.snps(d,pp.thr=0)
 fname <- file.path(fdir,paste(t1,"-gfm-on-tags-best-snps.RData",sep=""))
 save(bestsnps,file=fname)
}

sel.models.by.group.fn <- function(t1,snpGroups) {
 fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 fname <- file.path(fdir,paste(t1,"-gfm-on-tags-best-models0.0001.RData",sep=""))
 load(fname)
 fname <- file.path(fdir,paste(t1,"-gfm-on-tags-best-snps.RData",sep=""))
 load(fname)

 PPgfm <- summarise.fn(bestmods,snpGroups)
 fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-gfm-on-tags-groupPP.txt",sep="")
 write.table(PPgfm,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)
 
 MPPgfm <- Ginfo.fn(snpGroups,bestsnps)
 fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-gfm-on-tags-groupMPP.txt",sep="")
 write.table(MPPgfm,fname,quote=FALSE,row.names=TRUE,col.names=TRUE)
 }


