args=(commandArgs(TRUE))
t1=args[1]

library(annotSnpStats)
library(GUESSFM)

# full guessfm analysis
gfm.fn <- function(t1) {
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
 #load("/scratch/wallace/IL2RA/full-data/MS/tags.RData")
 fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 tfile <- file.path(fdir,"tags.RData") 
 load(tfile)
 X <- snp.data[,snps(tags)]
 mydir <- file.path(fdir,paste(t1,"-gfm",sep=""))
 run.bvs(X=X,Y=cc.ph,tag.r2=NA,nexp=3,nsave=1000,gdir=mydir,wait=TRUE) 

 d <- read.snpmod(mydir)
 dx <- expand.tags(d,tags)
 fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-dx-gfm-expanded.Rdata",sep="")
 save(dx,file=fname)
 #(load(fname))
 #best <- best.models(dx,pp.thr=0.0001)
 #abf <- abf.calc(y=cc.ph,x=X,models=best$str,family="binomial") 
 #sm <- abf2snpmod(abf,expected=3,nsnps=abf$nsnp)
 #best2 <- best.models(sm,pp.thr=0.0001) # best models after re-fit
 #bestsnps <- best.snps(sm,pp.thr=0)
 
# PPgfm <- summarise.fn(best2,snpGroups)
# fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-gfm-PP.txt",sep="")
# write.table(PPgfm,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)
# MPPgfm <- Ginfo.fn(snpGroups,bestsnps)
# fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-gfm-MPP.txt",sep="")
 #write.table(MPPgfm,fname,quote=FALSE,row.names=TRUE,col.names=TRUE)
}


#load("/scratch/wallace/IL2RA/full-data/MS/data.RData")
#load("/scratch/wallace/IL2RA/full-data/MS/tags.RData")
#DATA <- new("SnpMatrix",snp.data[cc.ph==0,snps(tags)]) 
#DATA <- as(snp.data[cc.ph==0,snps(tags)],"SnpMatrix")
#tfile <- "/home/ja628/causal_sim/results/IL2RA/real-data/new/tags.RData"
#tags0 <- tag(DATA, tag.threshold = 0.99)
#message("saving tags object to ", tfile)
#save(tags, file = tfile)



sel.models.fn <- function(t1) {
# fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
# tfile <- file.path(f1,"tags.RData") 
# tfile <- file.path(fdir,"tags.RData") 
load("/scratch/wallace/IL2RA/full-data/MS/tags.RData")
 load(dfile)
# load(tfile)
 
 tsnp <- snp.data[,unique(tags@tags)] 
 
 mydir <- file.path(fdir,paste(t1,"-gfm-on-tags",sep=""))
 run.bvs(X=tsnp,Y=cc.ph,tag.r2=NA,nexp=3,nsave=1000,gdir=mydir,wait=TRUE) 

 d <- read.snpmod(mydir) 
 bestmods <- best.models(d,pp.thr=0.0001)
 
 fname <- file.path(fdir,paste(t1,"-gfm-on-tags-best-models0.0001.txt",sep=""))
 save(bestmods,file=fname)
}

#fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"
#sel.models.fn("MS",fdir)
#sel.models.fn("JIA",fdir)
#sel.models.fn("GRAVES",fdir)
#sel.models.fn("T1D",fdir)


gfm.adj.fn <- function(t1,fdir) {
 
 tfile <- file.path(fdir,"tags.RData") 
 load(tfile)
 
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData")  
 load(dfile)
 t1data <- snp.data
 fname <- file.path(fdir,paste(t1,"-gfm-on-tags-best-models0.0001.txt",sep=""))
 (load(fname))
 best1 <- bestmods
 

 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t2,sep="")
 dfile <- file.path(f1,"data.RData") 
 load(dfile)
 t2data <- snp.data
 fname <- file.path(fdir,paste(t2,"-gfm-on-tags-best-models0.0001.txt",sep=""))
 (load(fname))
 best2 <- bestmods
 m2snps <- unique(unlist(strsplit(best2$str,"%")))
 
 msnps <- union(m1snps,m2snps)

fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-dx-gfm-expanded.Rdata",sep="")
 (load(fname))
 best <- best.models(dx,pp.thr=0.0001)
 abf <- abf.calc(y=cc.ph,x=snp.data,models=best1$str,family="binomial") 
 sm <- abf2snpmod(abf,expected=3,nsnps=abf$nsnp)
 best2 <- best.models(sm,pp.thr=0.0001) # best models after re-fit
 bestsnps <- best.snps(sm,pp.thr=0)
 
 



}

gfm.fn(t1)
#sel.models.fn(t1)