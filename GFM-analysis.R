args=(commandArgs(TRUE))
#r2=as.numeric(args[1])
#mppthr=as.numeric(args[2])
t1=as.numeric(args[1])

source("/home/ja628/scratch/scripts/IL2RA_general_scripts/multi-trait/multi-trait-fine-map.R")

r2=0.9
mppthr=0.01

traits <- c("T1D","GRAVES","MS","JIA")


trait1=traits[t1]

# directory to save files 
fdir <- "/home/ja628/causal_sim/results/IL2RA/real-data/new"

# do only one time: create tag snps object from common controls
# make.tags.fn(r2=0.9)
load(file.path(fdir,paste("tags-r2-",r2,".RData",sep="") ) )


# do once: save "best" tags SNPs for each trait
#sel.snps.fn(trait1,mppthr=0.01,r2=0.9,fdir)
#sel.snps.fn(trait2,mppthr=0.01,r2=0.9,fdir)

gfm.on.tags.final.fn(t1=trait1,r2=r2,tags=tags,fdir=fdir)

