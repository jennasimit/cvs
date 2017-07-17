args=(commandArgs(TRUE))
#r2=as.numeric(args[1])
#mppthr=as.numeric(args[2])
t1=as.numeric(args[1])
t2=as.numeric(args[2])

source("/home/ja628/scratch/scripts/IL2RA_general_scripts/multi-trait/multi-trait-fine-map-v5.R")

r2=0.9
mppthr=0.001

traits <- c("T1D","GRAVES","MS","JIA")
mtraits <- c(4,1,1,0)
Mtraits <- c(5,3,2,2)


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

j.mat <- matrix(1:dim(T1mods)[1],ncol=1) 
t1mod <- apply(j.mat,1,format.mod.fn,T1mods)
t1size <- apply(T1mods,1,sum)
T1mod <- data.frame(mod=t1mod,size=t1size)

j.mat <- matrix(1:dim(T2mods)[1],ncol=1) 
t2mod <- apply(j.mat,1,format.mod.fn,T2mods)
t2size <- apply(T2mods,1,sum)
T2mod <- data.frame(mod=t2mod,size=t2size)


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


G1.mat <- as(t1snp.data,"numeric")
data1 <- data.frame(Y=t1pheno,G1.mat)

G2.mat <- as(t2snp.data,"numeric")
data2 <- data.frame(Y=t2pheno,G2.mat)


traits <- paste(trait1,"-",trait2,sep="")


bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=T1mod$mod,family="binomial")[[1]] 
bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=T2mod$mod,family="binomial")[[1]] 

bf1$size <- t1size
bf2$size <- t2size


nullmod <- matrix(c(rep(0,(2*s)),1),nrow=1) # no snp effects and only trait effect
#t2m1snp <- t2snp.data[t2pheno==1,]
allsnps <- rbind(G1.mat[,msnps],G2.mat[t2pheno==1,msnps])
allpheno <- c(t1pheno,rep(2,sum(t2pheno)))
y <- character(length(allpheno))
y[allpheno==0] <- "CONTROL"
y[allpheno==1] <- trait1
y[allpheno==2] <- trait2
#G1.mat <- as(allsnps,"numeric")
data1 <- data.frame(Y=y,allsnps)
m1 <- mlogit2logit(Y ~ 1|. -Y,data1,choices=c("CONTROL",trait1,trait2),base.choice=1)
Nullmod <- glib.1(x=m1$data[,(4+s):dim(m1$data)[2]],y=m1$data$Y.star,error="binomial", link="logit",models=nullmod,post.bymodel = FALSE)
logBF0 <- Nullmod$bf$twologB10[,1]*0.5

T1T2abf <- T1T2abfexp.fn(bf1,bf2,logBF0)
#fname <- file.path(fdir,paste(traits,"-tags-bf.txt",sep=""))
#write.table(T1T2abf$out12,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)

# calculate PP and only expand joint models with PP>0.001
#pp<-PP12.fn(T1T2abf,shared=1,s=length(msnps),mT1=3,mT2=3)
pp<-PP.fn(T1T2abf,shared=1,s=length(unique(tags(tags))),mT1=3,mT2=3)



indbest <- which(pp$PP>0.00001)

BFkeep <- pp[-indbest,] # expand low pp models to see how many exist for each model and each of these will have the same bf as the tag model


BFexp.notbest <- T1T2bfexp.low.fn(BFkeep,tags,logBF0)




# expand tags at all marginal models that are part of a "best" joint model
tmp <- data.frame(logbf=pp$t1.logbf,mod=pp$t1.mod,size=pp$t1.size)
best1 <- unique(tmp[indbest,])

tmp <- data.frame(logbf=pp$t2.logbf,mod=pp$t2.mod,size=pp$t2.size)
best2 <- unique(tmp[indbest,])

e1 <- expand.tags.bf(best1,tags)
e2 <- expand.tags.bf(best2,tags)


G1.mat <- as(t1snp.data,"numeric")
data1 <- data.frame(Y=t1pheno,G1.mat)

G2.mat <- as(t2snp.data,"numeric")
data2 <- data.frame(Y=t2pheno,G2.mat)

bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=e1$str,family="binomial") 
bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=e2$str,family="binomial") 


T1bfexp <- e1
T1bfexp$logbf <- bf1[[1]][,"lBF"]

T2bfexp <- e2
T2bfexp$logbf <- bf2[[1]][,"lBF"]


#fname <- file.path(fdir,paste(trait1,"-bin-bf-final-",traits,".txt",sep=""))
#write.table(T1bfexp,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)

#fname <- file.path(fdir,paste(trait2,"-bin-bf-final-",traits,".txt",sep=""))
#write.table(T2bfexp,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)

# at expanded tags, re-fit joint models 

T1T2bfexp <- T1T2bfexp.best.fn(T1T2abf,T1bfexp,T2bfexp,indbest,logBF0)


# merge expanded re-fit joint with expanded non-re-fit joint 
tmp<-BFexp.notbest[,c("t1.logbf","t1.str","t1.size","t2.logbf","t2.str","t2.size","t1t2.logbf","t1t2.Nexpmod")]

BF <- rbind(T1T2bfexp,tmp)

fname <- file.path(fdir,paste(traits,"-bf-final.txt",sep=""))
write.table(BF,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)

##### PP calc

## binomial calc
# pp1 <- PP.bin.fn(T1bfexp,s=length(snps(tags)),mT1=3)
# pp2 <- PP.bin.fn(T2bfexp,s=length(snps(tags)),mT1=3)

s=length(snps(tags))

#fname <- file.path(fdir,paste(traits,"-bf-final.txt",sep=""))
#bf <- read.table(fname,header=TRUE,as.is=TRUE)
bf <- BF

shared <- c(1,5,10,20,100,1000)
ns <- length(shared)

PP <- NULL

for(k in 1:ns) PP <- rbind(PP,PP.fn(bf,shared=shared[k],s=s,mT1=3,mT2=3,details=FALSE))
row.names(PP) <- shared
fname <- file.path(fdir,paste(traits,"-PP-shared.txt",sep=""))
PP <- t(PP)
write.table(PP,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)

#PP <- read.table(fname,header=TRUE,row.names=1)

#PP <- vector("list",ns)
#for(k in 1:ns) PP[[k]] <- PP.fn(bf,shared=shared[k],s=s,mT1=3,mT2=3)

PPall <- data.frame(PP,BF$t1t2.Nexpmod)

alltraits <- c(trait1,trait2)
K <- length(alltraits)
traits <- paste(trait1,"-",trait2,sep="")
#fname <- file.path(fdir,paste(traits,"-PP-shared.txt",sep=""))
#PP <- read.table(fname,header=TRUE,row.names=1)

mpp <- vector("list",K)
PPmarg <- PP.marg.fn(PPall,K)
for(k in 1:K) {
fname <- file.path(fdir,paste(alltraits[k],"-PP-shared-",traits,".txt",sep=""))
write.table(PPmarg[[k]],fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
mpp[[k]] <- MPP.fn(PPmarg[[k]])
fname <- file.path(fdir,paste(alltraits[k],"-MPP-shared-",traits,".txt",sep=""))
write.table(mpp[[k]],fname,row.names=TRUE,col.names=TRUE,quote=FALSE)

}


#calc.fn <- function(t1,t2) {
#trait1=traits[t1];trait2=traits[t2];traitsb<-paste(trait1,"-",trait2,sep="")
#PP1 <- read.table(file.path(fdir,paste(trait1,"-PP-shared-",traitsb,".txt",sep="")),header=TRUE,row.names=1)
#PP2 <- read.table(file.path(fdir,paste(trait2,"-PP-shared-",traitsb,".txt",sep="")),header=TRUE,row.names=1)
#mpp1 <- MPP.fn(PP1)
#fname <- file.path(fdir,paste(trait1,"-MPP-shared-",traitsb,".txt",sep=""))
#write.table(mpp1,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
#mpp2 <- MPP.fn(PP2)
#fname <- file.path(fdir,paste(trait2,"-MPP-shared-",traitsb,".txt",sep=""))
#write.table(mpp2,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
#}


snpGroups <- make.snp.groups.fn()

MPP <- smartbind(t(mpp[[1]]),t(mpp[[2]]),fill=0) 
colnames(MPP)[which(colnames(MPP)=="0")] <- "m0" 

t1t2plot.fn(alltraits,shared=shared,pthr=0.01)

G <- character(length(colnames(MPP)))
for(k in 1:length(colnames(MPP))) {
check <- names(snpGroups)[grep(colnames(MPP)[k],snpGroups)];print(c(colnames(MPP)[k],check))
if(length(check)!=0) {
G[k] <- check} else {G[k] <- colnames(MPP)[k]}
}

snpG <- unique(G)
ng <- length(snpGroups)
mppG <- matrix(0,ncol=ng,nrow=dim(MPP)[1],dimnames=list(rownames(MPP),names(snpGroups)))
for(k in 1:ng) mppG[,k] <- apply(as.matrix(MPP[,which(G==names(snpGroups)[k])]),1,sum)

mppGS <- mppG
notinG <- setdiff(snpG,names(snpGroups))
if(length(notinG)>0) {
mppS <- matrix(0,ncol=length(notinG),nrow=dim(MPP)[1],dimnames=list(rownames(MPP),notinG))
for(k in 1:length(notinG)) mppS[,k] <- apply(as.matrix(MPP[,which(G==notinG[k])]),1,sum)
colnames(mppS) <- unlist(strsplit(colnames(mppS),"[.]"))[c(TRUE,FALSE,FALSE,FALSE)]
mppGS <- cbind(mppG,mppS)
}
rownames(mppGS) <- c(paste(alltraits[1],shared,sep="."),paste(alltraits[2],shared,sep="."))

fname <- file.path(fdir,paste(traits,"-group-MPP-shared.txt",sep=""))
write.table(mppGS,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)

pname <- file.path(fdir,paste("MPP-by-group-",traits,".pdf",sep=""))
pdf(pname)
par(mfrow=c(1,2))
MPP.plot.fn(as.matrix(mppGS[grep(alltraits[1],rownames(mppGS)),]),shared=shared,pthr=0.001)
if(alltraits[1] != "GRAVES") title(alltraits[1])
if(alltraits[1] == "GRAVES") title("AITD")

MPP.plot.fn(as.matrix(mppGS[grep(alltraits[2],rownames(mppGS)),]),shared=shared,pthr=0.001)
if(alltraits[2] != "GRAVES") title(alltraits[2])
if(alltraits[2] == "GRAVES") title("AITD")

dev.off()


colnames(MPP) <- unlist(strsplit(colnames(MPP),"[.]"))[c(TRUE,FALSE,FALSE,FALSE)]
ind <- grep("rs",G)
G[ind] <- unlist(strsplit(G[ind],"[.]"))[c(TRUE,FALSE,FALSE,FALSE)]

gPP <- vector("list",K)

for(k in 1:K) {
tmp <- rownames(PPmarg[[k]])
rn <- character(length(tmp))
Grn <- character(length(tmp))
for(l in 1:length(tmp)) {
 tmp1 <- unlist(strsplit(tmp[l],"%"))
 sp <- unlist(strsplit(tmp1,"[.]"))[c(TRUE,FALSE,FALSE,FALSE)] 
 rn[l] <- paste(sp,collapse=".")
 
 gg <- NULL
 #print(sp)
 for(ll in 1:length(sp)) gg<- c(gg,G[grep(sp[ll],colnames(MPP))]) 
 #print(gg)
 Grn[l] <- paste(gg[order(gg)],collapse=".")
 }
 
rownames(PPmarg[[k]]) <- rn

 
 gPPmarg <- PPmarg[[k]]
 rownames(gPPmarg) <- Grn
 
 mods <- unique(rownames(gPPmarg))
 nm <- length(mods)
 gPP[[k]] <- matrix(0,ncol=length(shared),nrow=nm,dimnames=list(mods,shared))
 for(mm in 1:nm)  gPP[[k]][mm,] <- apply(matrix(gPPmarg[Grn==mods[mm],],ncol=length(shared),byrow=FALSE),2,sum)
  fname <- file.path(fdir,paste(alltraits[k],"-PP-by-group-shared-",traits,".txt",sep=""))
write.table(gPP[[k]],fname,row.names=TRUE,col.names=TRUE,quote=FALSE)


}

pname <- file.path(fdir,paste("PP-shared-",traits,".pdf",sep=""))
pdf(pname)
par(mfrow=c(1,K))
 for(k in 1:K) {
  MPP.plot.fn(t(as.matrix(PPmarg[[k]])),shared=shared,pthr=0.01)
  if(alltraits[k] != "GRAVES") title(alltraits[k])
  if(alltraits[k] == "GRAVES") title("AITD")
  }  
 dev.off()

pname <- file.path(fdir,paste("PP-by-group-shared-",traits,".pdf",sep=""))
pdf(pname)
par(mfrow=c(1,K))
 for(k in 1:K) {
 MPP.plot.fn(t(as.matrix(gPP[[k]])),shared=shared,pthr=0.01)
  if(alltraits[k] != "GRAVES") title(alltraits[k])
  if(alltraits[k] == "GRAVES") title("AITD")
 }      
 dev.off()
 