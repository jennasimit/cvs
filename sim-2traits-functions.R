library(data.table)
library(parallel)
library(annotSnpStats)
library(snpStats)
library(GUESSFM)

##
inv.logit.fn <-function(x) return(exp(x)/(1+exp(x)))
##

phen.gen.fn <-function(beta1=c(-2.3,.2,.2),beta2=c(-2.3,.2,.2),snpG,N0=100,N1=100,N2=100,causals1.ind,causals2.ind) {
# beta[1] = log(prev)
 N <-dim(snpG)[2] # number of indivs
 causal1 <- snpG[causals1.ind,]
 causal2 <- snpG[causals2.ind,]
 n0=0
 n1=0
 n2=0
 G0 <-NULL
 G1 <- NULL
 G2 <- NULL
 j=1
 while( (n0<=N0 | n1 <=N1 | n2 <=N2) & (j <= N) ) {
  p1 <-inv.logit.fn(beta1[1]+sum(beta1[-1]*causal1[,j]))
  p2 <-inv.logit.fn(beta2[1]+sum(beta2[-1]*causal2[,j]))
  u<-runif(2)
   if(u[1]>p1 & u[1]>p2) {
     n0<-n0 + 1
     G0 <-rbind(G0,snpG[,j])
   } else if(u[1]<=p1 & u[1]>p2) {
    n1<-n1 + 1
     G1 <-rbind(G1,snpG[,j])
   } else if(u[1]<=p2 & u[1]>p1) {
     n2<-n2 + 1
     G2 <-rbind(G2,snpG[,j])
    } else {
     if(u[2]<0.5 & n1<N1 ) { n1<-n1 + 1; G1 <-rbind(G1,snpG[,j]) 
     } else { n2<-n2 + 1; G2 <-rbind(G2,snpG[,j]) }
    }
 
j <- j+1
	}

G<-rbind(G0[1:N0,],G1[1:N1,],G2[1:N2,])
rownames(G)<-c(paste("control.",1:N0,sep=""),paste("case1.",1:N1,sep=""),paste("case2.",1:N2,sep=""))
colnames(G) <- rownames(snpG)
		
 return(list(G=G,y=c(rep(0,N0),rep(1,N1),rep(2,N2)) ))
  	  	     }
  	  	     
####

format.mod.fn <- function(k,out) {
 ind <- which(out[k,]==1)
 if(length(ind)>0) {
  mod <- paste(names(out[k,][ind]),sep="",collapse="%")
  } else {mod <- "0"}
  return(mod)
	}

#####

#' @title Calculate joint PP at tag SNP models and then expanded models
#' @param t1snp.data a SnpMatrix object for trait 1
#' @param t1pheno a phenotype vector for trait 1
#' @param t2snp.data a SnpMatrix object for trait 2
#' @param t2pheno a phenotype vector for trait 2
#' @param tags tag SNPs from common controls 
#' @param mppthr threshold for "best" SNPs MPP
#' @param mT1 minimum model size for trait 1
#' @param MT1 maximum model size for trait 1
#' @param mT2 minimum model size for trait 2
#' @param MT2 maximum model size for trait 2
#' @return a list with components pp (data.frame of joint PP, prior, logprior, logBF, and for each trait, model, logBF, model size) and logBF0, the offset term needed to calculate joint BFs
#' @export
sharedmPP.fn <- function(kappa,t1snp.data,t1pheno,t2snp.data,t2pheno,
			tags,mppthr,mT1,MT1,mT2,MT2,
			trait1="T1",trait2="T2",mydir) {

#' find "best" tag snps
 t1snps <- gfm.sel.snps.fn(snpG=t1snp.data,y=t1pheno,tags=tags,mppthr=mppthr,mydir)
 t2snps <- gfm.sel.snps.fn(snpG=t2snp.data,y=t2pheno,tags=tags,mppthr=mppthr,mydir)
 msnps <- union(t1snps[,"var"],t2snps[,"var"])
 s <- length(msnps)

#' generate all models from tag snps
 T1mod <- T1mods.fn(mT1,MT1,msnps)
 T2mod <- T1mods.fn(mT2,MT2,msnps)

 G1.mat <- as(t1snp.data,"numeric")
 data1 <- data.frame(Y=t1pheno,G1.mat)

 G2.mat <- as(t2snp.data,"numeric")
 data2 <- data.frame(Y=t2pheno,G2.mat)

 bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=T1mod$mod,family="binomial")[[1]] 
 bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=T2mod$mod,family="binomial")[[1]] 
 
 D1 <- data.frame(mod=T1mod$mod,size=T1mod$size,logbf=bf1$lBF)
 D2 <- data.frame(mod=T2mod$mod,size=T2mod$size,logbf=bf2$lBF)

 e1 <- expand.tags.bf(D1,tags)
 e2 <- expand.tags.bf(D2,tags)

 bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=e1$str,family="binomial")[[1]] 
 bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=e2$str,family="binomial")[[1]] 
 
 
 STR <- list(e1$str,e2$str)
 ABF <- list(bf1$lBF,bf2$lBF)
 pr1 <- prior.bin.fn(e1$size)
 pr2 <- prior.bin.fn(e2$size)
 pr <- list(pr1,pr2)
 p0 <- prior.bin.fn(0)

 mpp <- marginalpp(STR, ABF, pr, kappa, p0)

 pp1 <- t(mpp$shared.pp[[1]])
 rownames(pp1) <- mpp$STR[[1]]
 colnames(pp1) <- paste("pp",kappa,sep=".")
  
 pp2 <- t(mpp$shared.pp[[2]])
 rownames(pp2) <- mpp$STR[[2]]
 colnames(pp2) <- paste("pp",kappa,sep=".")
  
 return(list(pp1=pp1,pp2=pp2))
}
  	  	     
###

PP2collapse.fn <- function(PP) {

b <- data.frame(matrix(0,nrow=8,ncol=10))
b[,1] <- c("rs12722563","rs12722522","rs12722508","rs7909519","rs61839660","rs12722496","rs12722495","rs41295049")
b[,2] <- "rs2104286"
b[,3] <- c("rs11597367","rs35285258","rs11594656",rep("rs",5))
b[,4] <- c("rs3118475","rs62626317","rs41295055","rs62626325","rs41295105","rs56382813","rs41295079","rs")
b[,5] <- "rs6602437"
b[,6] <- c("rs41295121","rs41295159",rep("rs",6))
b[,7] <- "rs706779" # AITD
b[,8] <- "rs3118470" # AA
b[,9] <- "rs4147359" # UC
b[,10] <- "rs10795791" # RA



m <- length(PP)
check <- NULL
for(i in 1:m) {
  bsnp <- unlist(strsplit(as.character(names(PP)[i]),"%"))
  checki <- numeric(10)
  for(j in 1:10) checki[j] <- 1*(length(intersect(b[,j],bsnp))>0)
   check <- rbind(check,checki)		
  }
gnames <- c("A","B","C","D","E","F","AITD","AA","UC","RA")
colnames(check) <-  gnames
rownames(check)<-NULL

mod.names.fn <- function(k,check) {paste(colnames(check)[which(check[k,]>0)],collapse="-")}
k.mat <- matrix(1:m,ncol=1)
mod.names <- apply(k.mat,1,mod.names.fn,check)

 c2 <- combn(1:10,2,simplify=TRUE)
 c2names <- apply(combn(gnames,2,simplify=TRUE),2,name.fn <- function(cc) return(paste(cc,collapse=".")))
 c3 <- combn(1:10,3,simplify=TRUE)
 c3names <- apply(combn(gnames,3,simplify=TRUE),2,name.fn <- function(cc) return(paste(cc,collapse=".")))
 m2 <- matrix(0,nrow=length(c2names),ncol=10); for(i in 1:length(c2names)) { m2[i,c2[,i]] <- 1 }
 m3 <- matrix(0,nrow=length(c3names),ncol=10); for(i in 1:length(c3names)) { m3[i,c3[,i]] <- 1 }
 m1 <- diag(10)
 M <- rbind(rep(0,10),m1,m2,m3)
 c1names <- gnames
 rownames(M) <- c("m0",c1names,c2names,c3names)
 colnames(M) <-c("A","B","C","D","E","F","AITD","AA","UC","RA")


out <- data.frame(PP=PP,check)
names(out)[2:11] <- gnames

outPP <- as.data.frame(matrix(0,nrow = 1, ncol = dim(M)[1],dimnames = list(NULL, rownames(M))))


indr <- NULL
ind <- 1:m
r <- 0
indr <- NULL
while(length(ind)>0 & r<dim(outPP)[2])  {
  while(length(indr)==0) {
  r <- r+1; #print(c("new r",r)) 
  indr <- which(apply(out[,-1], 1, function(x) all(x == M[r,])) )
  rM <- rownames(M)[r] 
			}
 outPP[,rM] <- sum(PP[indr])
 ind <- setdiff(ind,indr)
 indr <- NULL
}

# write.table(outPP,fname1,quote=FALSE,row.names=FALSE,col.names=TRUE)


return(outPP)
}

##

prior.bin.fn <- function(nT1,s=100,mT1=3) {
    dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
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
    neighb <- neighb[order(neighb$logbf, decreasing = TRUE), ]
    neighb[, `:=`(rank, 1:nrow(neighb))]
    best <- as.data.frame(neighb)
    return(best)
}
