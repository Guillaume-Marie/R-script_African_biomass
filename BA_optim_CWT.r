library("R2OpenBUGS")
library("coda")
library("lattice")
library("reshape2")
library("ggplot2")
library("quantreg")
library("RSQLite")
library("twosamples")
library("outliers")
library("mixdist")

# Create the function.
setprior <- function(cert,ftree_teta,brefp,fcp,LCTc){    
    a_MIX <- (ftree_teta-fcp$ftree_ESA_lw)/(fcp$ftree_ESA_up-fcp$ftree_ESA_lw)*(cert-2)+1
    b_MIX <- cert-a_MIX
    a_bare <- (fcp$fbare_LSCE-fcp$fbare_ESA_lw)/(fcp$fbare_ESA_up-fcp$fbare_ESA_lw)*(cert-2)+1
    b_bare <- cert-a_bare
    nlct <- length(LCTc)
    rs1 <- NULL; rs2 <- NULL; rs3 <- NULL; rs4 <-NULL
    for(i in 1:nlct){
        r1 <- rbeta(10000,a_MIX[i],b_MIX[i])*
       (fcp$ftree_ESA_up[i]-fcp$ftree_ESA_lw[i])+fcp$ftree_ESA_lw[i]
        r2 <- rbeta(10000,a_bare[i],b_bare[i])*
       (fcp$fbare_ESA_up[i]-fcp$fbare_ESA_lw[i])+fcp$fbare_ESA_lw[i]
        r3 <- rnorm(10000,brefp$Brefh_mu[i],brefp$Brefh_sigma[i])
        r4 <- rnorm(10000,brefp$Breft_mu[i],brefp$Breft_sigma[i])
        rs1 <- cbind(rs1,r1)
        rs2 <- cbind(rs2,r2)
        rs3 <- cbind(rs3,r3)
        rs4 <-cbind(rs4,r4)
    }   
    pr <- data.frame(fcp,a=a_MIX,b=b_MIX,PARAMETER="ftree_lc")
    prb <- data.frame(fcp,a=a_bare,b=b_bare,PARAMETER="fbare_lc")    
    rs1 <- data.frame(rs1)
    names(rs1) <- LCTc
    mrs1 <- melt(rs1)
    names(mrs1) <- c("LCT","val")
    mrs1$PARAMETER <- "ftree_lc"    
    rs2 <- data.frame(rs2)
    names(rs2) <- LCTc
    mrs2 <- melt(rs2)
    names(mrs2) <- c("LCT","val")
    mrs2$PARAMETER <- "fbare_lc"
    rs3 <- data.frame(rs3)
    names(rs3) <- LCTc
    mrs3 <- melt(rs3)
    names(mrs3) <- c("LCT","val")
    mrs3$PARAMETER <- "BMrefH"
    rs4 <- data.frame(rs4)
    names(rs4) <- LCTc
    mrs4 <- melt(rs4)
    names(mrs4) <- c("LCT","val")
    mrs4$PARAMETER <- "BMrefT"
    priors1 <- rbind(mrs1,mrs2,mrs3,mrs4)
    priors1$PARAMETER <- as.factor(priors1$PARAMETER)
    return(list(priors1,pr,prb))  
}

prior.frac <- read.table("prior_Fraction.csv",dec=",",sep=" ",header=TRUE)
model.struct <- "model_OPENBUGS.txt" 
biomass.map <- read.table("ybio.txt", header=TRUE) 

LCT <- c(10,11,30,40,50,60,61,62,100,110,120,130,150,153,160)
LCTc <- as.character(LCT)
LCTt <- c("h","h","n","n","t","t","t","t","n","n","h","h","h","h","t")
params <- c("ftree_lc","fherb_lc","fbare_lc","BMrefT","BMrefH","sigma","cert_w","cert_b")
nbiter <- 30000; ns <- 2000
nlct <- length(LCT)

j <- 0
final <- NULL; final2 <- NULL 
bioreft_up <- numeric(length(LCT)); bioreft_lw <- numeric(length(LCT))
biorefh_up <- numeric(length(LCT)); biorefh_lw <- numeric(length(LCT))
CIMAPul <- numeric(length(LCT)); CIMAPup <- numeric(length(LCT))
for(i in LCT){
   j <- j+1 
   ybio <- subset(biomass.map,LC_majority==i)
   qnt <- quantile(biomass.map$BM_majority, probs=c(.25, .75), na.rm = T)
   caps <- quantile(biomass.map$BM_majority, probs=c(.05, .95), na.rm = T)
   H <- 1.5 * IQR(biomass.map$BM_majority, na.rm = T)
   dd <- biomass.map[biomass.map$BM_majority > (qnt[1] - H),]
   ybio2 <- dd[dd$BM_majority < (qnt[2] + H),]
   ybio3 <- rm.outlier(ybio2$BM_majority)
   CIMAPup[j] <- quantile(ybio2$MAP_mean,0.75)
   CIMAPul[j] <- quantile(ybio2$MAP_mean,0.25)
   ulh <- getmode(ybio2$BM_majority)	
   uhh <- quantile(ybio2$BM_majority,probs=0.5) 
   uht <- quantile(ybio2$BM_majority,probs=0.95)
   ult <- quantile(ybio2$BM_majority,probs=0.85)
   x <- sample(1:dim(ybio2)[1],ns,replace=TRUE)
   ylc <- ybio2[x,] 
   bioreft_up[j] <- ifelse(uht<15,15,uht)
   bioreft_lw[j] <- ifelse(ult<9,9,ult)
   biorefh_up[j] <- ifelse(uhh>9,9,uhh)
   biorefh_up[j] <- ifelse(biorefh_up[j]<3,3, biorefh_up[j])
   biorefh_lw[j] <- ifelse(ulh>9,2,ulh)     
   if(j==1) final <- ylc["BM_majority"]
   else final <- cbind(final,ylc["BM_majority"])
   if(j==1) final2 <- ylc["MAP_mean"]
   else final2 <- cbind(final2,ylc["MAP_mean"])
}
final <- as.matrix(final)

Bref <- data.frame(LC_majority = LCT,
                   Brefh_mu = biorefh_up,
                   Breft_mu = bioreft_up,
                  )
Bref$Brefh_sigma <- Bref$Brefh_mu*0.15/4
Bref$Breft_sigma <- Bref$Breft_mu*0.15/4

for(i in 1:15) {
    print(paste("start for LCT",LCT[i]))
    data <- list(bb=as.vector(final[,i]), M = dim(final)[1],
              ftree_teta <- prior.frac$ftree_LSCE[i],         
              fbare_teta <- prior.frac$fbare_LSCE[i],         
              bfw <- Bref$Breft_mu[i],
              bfh <- Bref$Brefh_mu[i],
              bfw_s <- Bref$Breft_sigma[i],
              bfh_s <- Bref$Brefh_sigma[i],
              ESAup <- prior.frac$ftree_ESA_up[i],
              ESAlw <- prior.frac$ftree_ESA_lw[i],
              BAREup <- prior.frac$fbare_ESA_up[i],
              BARElw <- prior.frac$fbare_ESA_lw[i])

    assign(paste0("out.CWT",i),bugs(data,inits=NULL,params,model.struct,codaPkg=TRUE, 
                                    n.iter=nbiter,DIC=FALSE, n.chains=3,
                                    working.directory= paste0("/",LCT[i],"/"),
                                    restart=FALSE))
    print(paste("LCT",LCT[i],"....Done"))
}

coda.out15 <- read.bugs(out.CWT15)
dd15 <- gelman.diag(coda.out15)

listprior <- setprior(12,prior.frac$ftree_LSCE,Bref,prior.frac,LCTc)
pdf("Result_LCT_CWT_3.2_certRand_final_testtrue.pdf", width=36, height=20)
j <- 0; ddf <- NULL; ddfs <- NULL
for(i in LCT){
    j <- j+1
    out.CWT <- character(3)
    out.CWT[1] <- paste0("/certRand",i,"/CODAchain1.txt")
    out.CWT[2] <- paste0("/certRand",i,"/CODAchain2.txt")
    out.CWT[3] <- paste0("/certRand",i,"/CODAchain3.txt")
    coda.out <- read.bugs(out.CWT15)
    dd <- as.data.frame(summary(coda.out, q=c(0.025,0.975))[[2]])
    dds <- as.data.frame(summary(coda.out, q=c(0.025,0.975))[[1]])
    dd$LCT <- i; dds$LCT <- i
    ddf <- rbind(ddf,dd); ddfs <- rbind(ddf,dd)
    coda.melt0_CWT <- melt(as.data.frame(coda.out[[1]]),variable.name ="PARAMETER")
    coda.melt1_CWT <- melt(as.data.frame(coda.out[[2]]),variable.name ="PARAMETER")
    coda.melt2_CWT <- melt(as.data.frame(coda.out[[3]]),variable.name ="PARAMETER")
    mode0 <- getmode(as.vector(subset(coda.melt0_CWT,PARAMETER=="ftree_lc",value)[,1]))
    mode1 <- getmode(as.vector(subset(coda.melt1_CWT,PARAMETER=="ftree_lc",value)[,1]))
    mode2 <- getmode(as.vector(subset(coda.melt2_CWT,PARAMETER=="ftree_lc",value)[,1]))
    post.tetha <- mean(c(mode0, mode1, mode2))
    post.tetha0 <- data.frame(PARAMETER="ftree_lc",mode=post.tetha)                   
    print(ggplot(data=coda.melt0_CWT,aes(x=value))+geom_density()+
          geom_density(data=subset(coda.melt1_CWT, LCT==i),aes(x=value))+
          geom_density(data=subset(coda.melt2_CWT, LCT==i),aes(x=value))+
          geom_density(data=subset(listprior[[1]], LCT==i & PARAMETER=="BMrefT"),aes(x=val),lty=2)+
          geom_density(data=subset(listprior[[1]], LCT==i & PARAMETER=="BMrefH"),aes(x=val),lty=2)+          
          geom_vline(data=post.tetha0,aes(xintercept=mode),lty=2,size=1,col='blue')+        
          geom_vline(data=subset(listprior[[2]], LCT==i),aes(xintercept=ftree_ESA_lw),lty=2,size=1)+
          geom_vline(data=subset(listprior[[2]], LCT==i),aes(xintercept=ftree_ESA_up),lty=2,size=1)+
          geom_vline(data=subset(listprior[[3]], LCT==i),aes(xintercept=fbare_ESA_up),lty=2,size=1)+
          geom_vline(data=subset(listprior[[3]], LCT==i),aes(xintercept=fbare_ESA_lw),lty=2,size=1)+
          geom_vline(data=subset(listprior[[2]], LCT==i),aes(xintercept=ftree_LSCE),lty=2,size=1,col='red')+
          geom_vline(data=subset(listprior[[3]], LCT==i),aes(xintercept=fbare_LSCE),lty=2,size=1,col='red')+ 
          facet_wrap(PARAMETER~.,scales=c("free"))+ ggtitle(paste("LCT : ",i))
          )
}
dev.off()

write.table(ddf,"CI_posterior.csv", sep=";")
write.table(ddfs,"CI_posterior_normal.csv", sep=";")

