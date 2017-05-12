## load DATA 
## dx is a list of expanded GUESSFM outputs (one per disease) 
## dirs is a list of directories used for the different GUESSFM runs (one per disease) 
## load the control data - you will just have one set 

library(GUESSFM)
library(annotSnpStats)

############

gfm.expanded.fn <- function(t1) {
 f1 <- paste("/scratch/wallace/IL2RA/full-data/",t1,sep="")
 dfile <- file.path(f1,"data.RData") 
 tfile <- "/home/ja628/causal_sim/results/IL2RA/real-data/new/tags.RData"
 load(dfile)
 load(tfile)
 
 d <- read.snpmod(f1)
 dx <- expand.tags(d,tags)

 fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",t1,"-dx-gfm-expanded.Rdata",sep="")
 save(dx,file=fname)
}

#gfm.fn(t1,snpGroups)
gfm.expanded.fn("JIA")
gfm.expanded.fn("MS")
gfm.expanded.fn("GRAVES")
gfm.expanded.fn("T1D")
gfm.expanded.fn("RA")

###

traits <- c("JIA","GRAVES","MS","T1D","RA")

#dirs <- as.list(paste(rep("/scratch/wallace/IL2RA/full-data/",length(traits)),traits,sep=""))
#names(dirs) <- traits
#dfiles <-unlist(sapply(dirs,list.files,pattern="data.RData",full=TRUE)) 
#DATA <- lapply(dfiles, function(f) { (load(f)); return(as(snp.data[cc.ph==0,],"SnpMatrix"))}) 

load("/scratch/wallace/IL2RA/full-data/MS/data.RData")
#DATA <- as(snp.data[cc.ph==0,],"SnpMatrix")
DATA <- new("SnpMatrix",snp.data[cc.ph==0,])

m <- length(traits)
dx.all <- vector("list",m)
names(dx.all) <- traits
for(k in 1:m) {
fname <- paste("/home/ja628/causal_sim/results/IL2RA/real-data/new/",traits[k],"-dx-gfm-expanded.Rdata",sep="")
(load(fname))
dx.all[[k]] <- dx
}

dx <- dx.all

## you can use lapply here 
SPx <- lapply(dx, snp.picker, DATA,start.thr=0.01,r2.gap=0.05,nochange.thr = 0.001,nochange.run=2) 
fname <- "/home/ja628/causal_sim/results/IL2RA/real-data/new/snp-picker-groups.RData"
save(SPx,file=fname)

groupsx <- as(SPx[[1]],"groups") # first group
if(length(SPx)>1)        
for(i in 2:length(SPx)) groupsx <- union(groupsx,as(SPx[[i]],"groups"))
fname <- "/home/ja628/causal_sim/results/IL2RA/real-data/new/snp-picker-groups-merged.RData"
save(groupsx,file=fname)


