                                        #install.packages("googlesheets")
#install.packages("ggjoy")
library(googlesheets)
library(ggjoy)

(my_sheets <- gs_ls("Tables.xls"))
g <- gs_key("13etUMuUDcy4uLuT6hhCAkol2zw2G7sP31boDeI10SYs")

extracter <- function(data, rows,cols,method,sim) {
    d <- as.data.frame(data[rows,cols])
    colnames(d) <- c("model",paste0("N.",seq(1000,by=1000,length.out=ncol(d)-1)))
    d$method <- method
    d$sim <- sim
    d
}
library(data.table)
collapse <- function(x,wh) {
    xd <- as.data.table(x)
    for(w in wh)
        xd[model==w,model:="other"]
    xd <- xd[,value:=sum(value),by=c("model","method","sim","variable")]
    as.data.table(xd)
}
 
data <- gs_read(g,"B-sims")
B <- list(extracter(data,2:10,1:6,"Stochastic","B"),
          extracter(data,14:22,1:6,"Stepwise","B"))

data <- gs_read(g,"AD-sims")
data
which(letters %in% c("h","k","r","u"))
AD <- list(extracter(data,2:10,1:8,"Stochastic","A<D"),
           extracter(data,2:10,11:18,"Stochastic","A=D"),
           extracter(data,2:10,21:28,"Stochastic","A>D"),
           extracter(data,15:23,1:8,"Stepwise","A<D"),
           extracter(data,15:23,11:18,"Stepwise","A=D"),
           extracter(data,15:23,21:28,"Stepwise","A>D"))

library(magrittr)
library(reshape)
BAD <- c(B,AD) %>% lapply(., melt) %>% do.call("rbind",.)
BAD$sim %<>% factor(.,levels=c("B","A<D","A=D","A>D"))



data <- gs_read(g,"AC-sims")
AC <- list(extracter(data,2:8,1:8,"Stochastic","A<C"),
           extracter(data,2:8,11:18,"Stochastic","A=C"),
           extracter(data,2:8,21:28,"Stochastic","A>C"),
           extracter(data,13:19,1:8,"Stepwise","A<C"),
           extracter(data,13:19,11:18,"Stepwise","A=C"),
           extracter(data,13:19,21:28,"Stepwise","A>C"))

data <- gs_read(g,"sAITD-sims")
s <- list(extracter(data,2:8,1:6,"Stochastic","sAITD"),
          extracter(data,13:19,1:6,"Stepwise","sAITD"))

sAC <- c(s,AC) %>% lapply(., melt) %>% do.call("rbind",.)
sAC$sim %<>% factor(.,levels=c("sAITD","A<C","A=C","A>C"))

################################################################################

   
    
plotter <- function(B,what=c("AD","AC"),orient=c("h","v")) {
    what <- match.arg(what)
    x <- subset(B,!is.na(value))
    if(what=="AD") {
        del <- c("A.B.D","B.D","A.B")
        lev <- c("null","A","B","D","A.B","A.D","B.D","other")
    } else {
        del <-c("D","sRA","A.C.D","B.C","C.D","A.C.sAITD") 
        lev <- c("null","A","C","D","sAITD","sRA","A.C","B.C","C.D","A.C.D","A.C.sAITD","other")
    }
    x <- collapse(x,del)
    lev <- setdiff(lev,del)
    orient <- match.arg(orient)
    x$true <- x$model==gsub("<|>|=",".",x$sim)
    multi <- grepl("[=<>]",as.character(x$sim))
    ss1 <- substr(as.character(x$sim),1,1)
    ss3 <- substr(as.character(x$sim),3,3)
    x$partial <- multi==TRUE & (x$model==ss1 | x$model==ss3)
    x$model <- factor(as.character(x$model),levels=lev)
    x$n <- as.numeric(sub("N.","",x$variable))
    x$x <- as.numeric(x$model)
    xd <- as.data.table(subset(x,true==TRUE|partial==TRUE))
    xd <- xd[,.(minv=levels(variable)[min(as.numeric(variable))],
                maxv=levels(variable)[max(as.numeric(variable))],
                x=as.numeric(x)),
             by=c("method","sim","partial")]
    xl <- xu <- x
    xl$x <- xl$x - 0.3
    xu$x <- xu$x + 0.3
    xl0 <- xl
    xu0 <- xu
    xl0$value <- xu0$value <- 0
    xb <- rbind(xl0,x,xu0)
    ## x$group <- paste(x$variable,x$method)
    if(orient=="v") {
        ggplot(xb,aes(x=x,y=variable,height=value,group=variable,fill=model)) +
        geom_vline(aes(xintercept=x),data=xd[partial==FALSE,],col="darkblue") + #,linetype="dashed",size=1.5)+ 
        geom_vline(aes(xintercept=x),data=xd[partial==TRUE,],col="darkblue", linetype="dashed") + #,size=1.5)+ 
        geom_ridgeline_gradient(scale=1,stat="identity") +
        scale_x_continuous("model",breaks=seq_along(levels(x$model)),
                           labels=levels(x$model)) +
        scale_y_discrete("sample size",labels=sub("N.","",levels(x$variable))) +
       theme_joy() +
theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1))+ guides(fill = guide_legend(nrow = 1))
     } else {
        ggplot(x) +
        geom_rect(aes(xmin=minv,xmax=maxv,ymin=x,ymax=x+1),fill="yellow",data=xd[partial==FALSE,],col="yellow")+ 
        geom_rect(aes(xmin=minv,xmax=maxv,ymin=x,ymax=x+1),fill="lightyellow",data=xd[partial==TRUE,],col="lightyellow")+ 
        geom_ridgeline_gradient(aes(y=x,x=variable,height=value,group=x,fill=model),scale=1,stat="identity") +
        scale_y_continuous("model",breaks=seq_along(levels(x$model)),
                           labels=levels(x$model)) +
        scale_x_discrete("sample size",labels=sub("N.","",levels(x$variable))) +
 theme_joy() +
theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1))+ guides(fill = guide_legend(nrow = 1))
   }
}
plotter(BAD,what="AD",orient="v") + facet_grid(method ~ sim,scales="free")  
ggsave("joy-BAD.jpg",height=6,width=10,dpi=300)

plotter(sAC,what="AC",orient="v") + facet_grid(method ~ sim,scales="free")  
ggsave("joy-sAC.jpg",height=6,width=10,dpi=300)



library("googledrive")
#drive_find(pattern="IL2RA")
drive_upload("joy-BAD.jpg",path="IL2RA/Figures/")
drive_upload("joy-sAC.jpg",path="IL2RA/Figures/")

plotter(BAD,what="AD",orient="h") + facet_grid(method ~ sim,scales="free_x")  

plotter(C,what="AC",orient="h") + facet_grid(method ~ sim,scales="free")  
