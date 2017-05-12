args=(commandArgs(TRUE))
#r2=as.numeric(args[1])
#mppthr=as.numeric(args[2])
t1=as.numeric(args[1])
t2=as.numeric(args[2])

source("/home/ja628/scratch/scripts/IL2RA_general_scripts/multi-trait/multi-trait-fine-map.R")

r2=0.9
mppthr=0.01

traits.all <- c("T1D","GRAVES","MS","JIA")

trait1=traits.all[t1]
trait2=traits.all[t2]
traits <- paste(trait1,"-",trait2,sep="")

# directory to save files 
fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"

load(file.path(fdir,paste("tags-r2-",r2,".RData",sep="") ) )
s=length(snps(tags))

fname <- file.path(fdir,paste(traits,"-bf-final.txt",sep=""))
bf <- read.table(fname,header=TRUE,as.is=TRUE)

shared <- c(1,5,10,20,100)
ns <- length(shared)

PP <- NULL

for(k in 1:ns) PP <- rbind(PP,PP.fn(bf,shared=shared[k],s=s,mT1=3,mT2=3))
row.names(PP) <- shared
fname <- file.path(fdir,paste(traits,"-PP-shared.txt",sep=""))
write.table(t(PP),fname,row.names=TRUE,col.names=TRUE,quote=FALSE)

#PP <- read.table(fname,header=TRUE,row.names=1)

PPmarg <- PP.marg.fn(PP)
fname <- file.path(fdir,paste(trait1,"-PP-shared-",traits,".txt",sep=""))
write.table(PPmarg$PP1,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
fname <- file.path(fdir,paste(trait2,"-PP-shared-",traits,".txt",sep=""))
write.table(PPmarg$PP2,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)