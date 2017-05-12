args=(commandArgs(TRUE))
r2=as.numeric(args[1])
mppthr=as.numeric(args[2])
tk=as.numeric(args[3])

library(annotSnpStats)
library(GUESSFM)


traits <- c("T1D","GRAVES","MS","JIA")
t1=traits[tk]

load("/scratch/wallace/IL2RA/full-data/MS/data.RData")
load("/scratch/wallace/IL2RA/full-data/MS/tags.RData")
DATA <- new("SnpMatrix",snp.data[cc.ph==0,snps(tags)]) 

tags <- tag(DATA, tag.threshold = r2)

fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"


gfm.on.tags.fn <- function(t1,r2,tags,fdir) {
 #fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
  
 tsnp <- snp.data[,unique(tags@tags)] 
 
 mydir <- file.path(fdir,paste(t1,"-gfm-on-tags-r2-",r2,sep=""))
 run.bvs(X=tsnp,Y=cc.ph,tag.r2=NA,nexp=3,nsave=1000,gdir=mydir,wait=TRUE) 

}

sel.snps.fn <- function(t1,mppthr=0.01,r2,fdir) {
 #fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 
 #mydir <- file.path(fdir,paste(t1,"-gfm-on-tags-r2-",r2,sep=""))
 mydir <- file.path(fdir,paste(t1,"-gfm-on-tags",sep=""))

 d <- read.snpmod(mydir)  
 bestsnps <- best.snps(d,pp.thr=mppthr)
 fname <- file.path(fdir,paste(t1,"-gfm-on-tags-r2-",r2,"-best-snps-mpp-",mppthr,".RData",sep=""))
 save(bestsnps,file=fname)
}


gfm.on.tags.fn(t1=t1,r2=r2,tags=tags,fdir=fdir)
sel.snps.fn(t1=t1,mppthr=mppthr,r2=r2,fdir=fdir)


load("/home/ja628/causal_sim/results/IL2RA/real-data/new/tags.RData") # r2=0.99