library("R2OpenBUGS")
library("coda")
library("lattice")
library("reshape2")0!
library("ggplot2")
library("quantreg")
library("RSQLite")
library("twosamples")
library("outliers")
library("mixdist")  
# Create the function.
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

theme_set(
    theme_classic(base_size = 30)
)

conn <- dbConnect(RSQLite::SQLite(), "/home/volarex/Images/temp_bm_02.gpkg")

BM_majority <- 0

dd1=dbGetQuery(conn,"SELECT LC_majority, MAP_mean, BM_majority, LC_variety FROM temp_bm_02")

dd2=dbGetQuery(conn,"SELECT LC_majority, MAP_mean, BM_majority FROM temp_bm_02 WHERE LC_variety==1 AND LC_majority IN (10,11,30,40,50,60,61,62,100,110,120,130,150,153,160,200)", params = BM_majority)

BMcount=tapply(dd2[,"BM_majority"],dd2[,"LC_majority"],length)
print(BMcount)


write.table(dd2,"/home/volarex/Images/summary_LC_1km.csv", sep=";")

ytot=read.table(paste0("/home/volarex/Images/summary_LC_1km.csv"),header=TRUE,sep=";")
ytot=na.omit(ytot)

# logical_and(A>=60,A<=62)*max(B,C)+logical_or(A==160,A==50)*C+logical_and(A>=10,A<=40)*B+logical_and(A>=100,A<160)*B+logical_and(A>160,A<=300)*B


ybio=subset(ytot, (BM_majority>0 & BM_majority<100 & (LC_majority==10 |LC_majority==11|LC_majority==30|LC_majority==40)))  
df <- mixgroup(ybio$BM_majority, breaks = c(seq(0,20,1),seq(21,100,5)))
(fitpro <- mix(mixdat=as.mixdata(df), mixpar=mixparam(mu=c(1,15,70), sigma=c(2,5,5)), dist="lnorm", emsteps = 3))  
plot(fitpro, main="Fit a Probability Distribution",xlim=c(0,100))
grid()  
legend("topright", lty=1, lwd=c(1, 1, 2), c("Original Distribution to be Fit", "Individual Fitted Distributions", "Fitted Distributions Combined"), col=c("blue", "red", rgb(0.2, 0.7, 0.2)), bg="white")


ybio=subset(ytot, (BM_majority>=0 & BM_majority<100 & (LC_majority==120)))  
df <- mixgroup(ybio$BM_majority, breaks= c(seq(0,20,1),seq(25,100,5)))
(fitpro <- mix(mixdat=as.mixdata(df), mixpar=mixparam(mu=c(1,10,80), sigma=c(5,5,5)), dist="lnorm", emsteps = 3))  
plot(fitpro, main="Fit a Probability Distribution",xlim=c(0,100))
grid()  
legend("topright", lty=1, lwd=c(1, 1, 2), c("Original Distribution to be Fit", "Individual Fitted Distributions", "Fitted Distributions Combined"), col=c("blue", "red", rgb(0.2, 0.7, 0.2)), bg="white")


ybio=subset(ytot, (BM_majority>0 & BM_majority<100 & LC_majority==60))  
df <- mixgroup(ybio$BM_majority)
(fitpro <- mix(mixdat=as.mixdata(df), mixpar=mixparam(mu=c(10,80), sigma=c(5,5)), dist="norm", emsteps = 3))  
plot(fitpro, main="Fit a Probability Distribution",xlim=c(0,100))
grid()  
legend("topright", lty=1, lwd=c(1, 1, 2), c("Original Distribution to be Fit", "Individual Fitted Distributions", "Fitted Distributions Combined"), col=c("blue", "red", rgb(0.2, 0.7, 0.2)), bg="white")



LCT=c(10,11,30,40,50,60,61,62,100,110,120,130,150,153,160)
LCTt=c("h","h","n","n","t","t","t","t","n","n","h","h","h","h","t")
LCTc=c("10","11","30","40","50","60","61","62","100","110","120","130","150","153","160")

nlct=length(LCT)
    ns=2000
    j=0
    final=NULL;final2=NULL 
    bioreft_up=numeric(length(LCT)); bioreft_lw=numeric(length(LCT))
    biorefh_up=numeric(length(LCT)); biorefh_lw=numeric(length(LCT))
    CIMAPul=numeric(length(LCT));CIMAPup=numeric(length(LCT))
    for(i in LCT){
      j=j+1 
      ybio=subset(ytot,LC_majority==i)
      qnt <- quantile(ybio$BM_majority, probs=c(.25, .75), na.rm = T)
      caps <- quantile(ybio$BM_majority, probs=c(.05, .95), na.rm = T)
      H <- 1.5 * IQR(ybio$BM_majority, na.rm = T)
      dd=ybio[ybio$BM_majority > (qnt[1] - H),]
      ybio2=dd[dd$BM_majority < (qnt[2] + H),]
      ybio3=rm.outlier(ybio2$BM_majority)
      CIMAPup[j] <- quantile(ybio2$MAP_mean,0.75)
      CIMAPul[j] <- quantile(ybio2$MAP_mean,0.25)
      ulh <- getmode(ybio2$BM_majority)	
      uhh <- quantile(ybio2$BM_majority,probs=0.5) 
      uht <- quantile(ybio2$BM_majority,probs=0.95)
      ult <- quantile(ybio2$BM_majority,probs=0.85)
      #ul <- getmode(ybio$BM_majority)
      x=sample(1:dim(ybio2)[1],ns,replace=TRUE)
      ylc=ybio2[x,] 
      bioreft_up[j]=ifelse(uht<15,15,uht)
      bioreft_lw[j]=ifelse(ult<9,9,ult)
      biorefh_up[j]=ifelse(uhh>9,9,uhh)
      biorefh_up[j]=ifelse(biorefh_up[j]<3,3, biorefh_up[j])
      biorefh_lw[j]=ifelse(ulh>9,2,ulh)     
      if(j==1) final=ylc["BM_majority"]
      else final=cbind(final,ylc["BM_majority"])
      if(j==1) final2=ylc["MAP_mean"]
      else final2=cbind(final2,ylc["MAP_mean"])
    }
    final=as.matrix(final)

summary(final)

Bref <- data.frame(LC_majority=LCT,Brefh_lw=biorefh_lw,Brefh_up=biorefh_up,Breft_lw=bioreft_lw,Breft_up=bioreft_up,CIMAPup,CIMAPul)

#Bref <-read.table(paste0("/home/volarex/Images/Bref.csv"),header=TRUE,sep=";")
Bref$Brefh_mu <- Bref$Brefh_up
Bref$Breft_mu <- Bref$Breft_up
Bref$Brefh_sigma <-Bref$Brefh_up*0.15/4
#*(CIMAPup-CIMAPul)/max((CIMAPup-CIMAPul))
Bref$Breft_sigma <-Bref$Breft_up*0.15/4
#*0.25*(CIMAPup-CIMAPul)/max((CIMAPup-CIMAPul))

write.table(Bref,"/home/volarex/Images/Bref.csv", sep=";",row.names = FALSE)

png(paste0("/home/volarex/Images/figure_2.3.4.png"), width=1950, height=1080)
ggplot(ytot, aes(BM_majority, after_stat(density))) + geom_freqpoly(bins=60,size=1)  + xlim(0,410)+
geom_vline(data=Bref,aes(xintercept=Brefh_up),lty=2,colour="#F8766D",size=1)+
geom_vline(data=Bref,aes(xintercept=Breft_up),lty=2,colour="#00BFC4",size=1)+
facet_wrap(~LC_majority,scales="free_y")+ xlab("Biomass concentration from CSBIO map (t/ha)")
dev.off()


png(paste0("/home/volarex/Images/figure_2.3.4_bis.png"), width=1950, height=1080, pointsize=40)
ggplot(ytot, aes(MAP_mean, after_stat(density))) + geom_freqpoly(bins=60,size=1)  + xlim(0,4200)+
    geom_vline(data=Bref,aes(xintercept=CIMAPup),lty=2,size=1, lty=2)+
    geom_vline(data=Bref,aes(xintercept=CIMAPul),lty=2,size=1, lty=2)+
facet_wrap(~LC_majority,scales="free_y")
dev.off()

prior <- read.table("/home/volarex/Images/prior_Fraction.csv",dec=",",sep=" ",header=TRUE)

setprior <- function(cert,ftree_teta,brefp,fcp,LCTc){
    
    a_MIX=(ftree_teta-fcp$ftree_ESA_lw)/(fcp$ftree_ESA_up-fcp$ftree_ESA_lw)*(cert-2)+1
    b_MIX=cert-a_MIX
    a_bare=(fcp$fbare_LSCE-fcp$fbare_ESA_lw)/(fcp$fbare_ESA_up-fcp$fbare_ESA_lw)*(cert-2)+1
    b_bare=cert-a_bare

    nlct <- length(LCTc)
    rs1=NULL; rs2=NULL; rs3=NULL; rs4=NULL
    for(i in 1:nlct){
        r1 <- rbeta(10000,a_MIX[i],b_MIX[i])*(fcp$ftree_ESA_up[i]-fcp$ftree_ESA_lw[i])+fcp$ftree_ESA_lw[i]
        r2 <- rbeta(10000,a_bare[i],b_bare[i])*(fcp$fbare_ESA_up[i]-fcp$fbare_ESA_lw[i])+fcp$fbare_ESA_lw[i]
        r3 <- rnorm(10000,brefp$Brefh_mu[i],brefp$Brefh_sigma[i])
        r4 <- rnorm(10000,brefp$Breft_mu[i],brefp$Breft_sigma[i])
        rs1=cbind(rs1,r1)
        rs2=cbind(rs2,r2)
        rs3=cbind(rs3,r3)
        rs4=cbind(rs4,r4)
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

png(paste0("/home/volarex/Images/figure_2.3.6.1.png"), width=1950, height=1080, pointsize=40)
ggplot(mrs3, aes(x=val)) + geom_density(size=1.5)+
    geom_vline(data=pr,aes(xintercept=mu),lty=1,colour="#00BFC4",size=1.5)+   
    geom_vline(data=pr,aes(xintercept=upper),lty=2,colour="#F8766D",size=1.5)+
    geom_vline(data=pr,aes(xintercept=lower),lty=2,colour="#F8766D",size=1.5)+
    facet_wrap(LCT~.,scales=c("free"))
dev.off()

png(paste0("/home/volarex/Images/figure_2.3.6.3.png"), width=1950, height=1080, pointsize=40)
ggplot(brs, aes(x=val)) + geom_density(size=1.5)+
    geom_vline(data=prb,aes(xintercept=mu),lty=1,colour="#00BFC4",size=1.5)+   
    geom_vline(data=prb,aes(xintercept=upper),lty=2,colour="#F8766D",size=1.5)+
    geom_vline(data=prb,aes(xintercept=lower),lty=2,colour="#F8766D",size=1.5)+
    facet_wrap(LCT~.,scales=c("free"))
dev.off()

for(i in c(6,12,50,100)){
    listprior <- setprior(i,prior$ftree_LSCE,Bref,prior,LCTc)
    model.file = "/home/volarex/Images/model_ftree_V2.txt" 
    params=c("ftree_lc","fherb_lc","fbare_lc","BMrefT","BMrefH","sigma")
    nbiter=50
    LCTlong=rep(sort(rep(LCT,nbiter/2)),length(params))
    data=list(bb=final, M = dim(final)[1], N = dim(final)[2],
              ftree_a = listprior[[2]]$a,
              ftree_b = listprior[[2]]$b,         
              fbare_a = listprior[[3]]$a,
              fbare_b = listprior[[3]]$b,         
              bfw = Bref$Breft_mu,
              bfh = Bref$Brefh_mu,
              bfw_s = Bref$Breft_sigma,
              bfh_s = Bref$Brefh_sigma,
              ESAup=prior$ftree_ESA_up,
              ESAlw=prior$ftree_ESA_lw,
              BAREup=prior$fbare_ESA_up,
              BARElw=prior$fbare_ESA_lw)
    
    assign(paste0("out.CWT",i),bugs(data,inits=NULL,params,model.file,codaPkg=TRUE, n.iter=nbiter,
                 DIC=FALSE, n.chains=3, working.directory=paste0("/home/volarex/Images/cert",i,"/"),restart=FALSE))    
}

coda.out6=read.bugs(out.CWT6)
dd6=gelman.diag(coda.out6)
print(dd6)

coda.out12=read.bugs(out.CWT12)
dd12=gelman.diag(coda.out12)
print(dd12)

coda.out50=read.bugs(out.CWT5)
dd50=gelman.diag(coda.out50)
print(dd50)

coda.out100=read.bugs(out.CWT100)
dd100=gelman.diag(coda.out100)
print(dd100)

densityplot(coda.out)
xyplot(out.coda)
dd=summary(coda.out, q=c(0.025,0.5,0.975))

coda.melt0_CWT <- melt(as.data.frame(coda.out6[[1]]))
fact <- gsub(pattern = "\\[|\\]",replacement = " ",x = coda.melt0_CWT$variable)
coda.melt0_CWT <- data.frame(colsplit(fact, " ", names = c("PARAMETER", "LCT")),fc=coda.melt0_CWT$value)
coda.melt0_CWT$LCT <- LCTlong

coda.melt1_CWT <- melt(as.data.frame(coda.out12[[1]]))
fact <- gsub(pattern = "\\[|\\]",replacement = " ",x = coda.melt1_CWT$variable)
coda.melt1_CWT <- data.frame(colsplit(fact, " ", names = c("PARAMETER", "LCT")),fc=coda.melt1_CWT$value)
coda.melt1_CWT$LCT <- LCTlong
    
coda.melt2_CWT <- melt(as.data.frame(coda.out100[[1]]))
fact <- gsub(pattern = "\\[|\\]",replacement = " ",x = coda.melt2_CWT$variable)
coda.melt2_CWT <- data.frame(colsplit(fact, " ", names = c("PARAMETER", "LCT")),fc=coda.melt2_CWT$value)
coda.melt2_CWT$LCT <- LCTlong

coda.melt3_CWT <- melt(as.data.frame(coda.out50[[1]]))
fact <- gsub(pattern = "\\[|\\]",replacement = " ",x = coda.melt3_CWT$variable)
coda.melt3_CWT <- data.frame(colsplit(fact, " ", names = c("PARAMETER", "LCT")),fc=coda.melt3_CWT$value)
coda.melt3_CWT$LCT <- LCTlong

pdf("/home/volarex/Images/Result_LCT_CWT_3.2_certR12-100.pdf", width=36, height=20)
for(i in LCT){
print(ggplot(data=subset(coda.melt1_CWT,LCT==i),aes(x=fc))+geom_density()+
      #geom_density(data=subset(listprior[[1]], LCT==i),aes(x=val),lty=2)+
      geom_density(data=subset(coda.melt2_CWT, LCT==i),aes(x=fc))+
      geom_density(data=subset(coda.melt3_CWT, LCT==i),aes(x=fc))+
      geom_density(data=subset(coda.melt0_CWT, LCT==i),aes(x=fc))+      
      geom_vline(data=subset(listprior[[2]], LCT==i),aes(xintercept=ftree_LSCE),lty=2,size=1)+     
      geom_vline(data=subset(listprior[[2]], LCT==i),aes(xintercept=ftree_ESA_lw),lty=2,size=1)+
      geom_vline(data=subset(listprior[[2]], LCT==i),aes(xintercept=ftree_ESA_up),lty=2,size=1)+
      geom_vline(data=subset(listprior[[3]], LCT==i),aes(xintercept=fbare_ESA_up),lty=2,size=1)+
      geom_vline(data=subset(listprior[[3]], LCT==i),aes(xintercept=fbare_ESA_lw),lty=2,size=1)+      
    facet_wrap(PARAMETER~.,scales=c("free"))+ ggtitle(paste("LCT : ",i)))
}
dev.off()


ftree_lc1=subset(coda.melt1_CWT,PARAMETER=="ftree_lc")
fbare_lc1=subset(coda.melt1_CWT,PARAMETER=="fbare_lc")

round(as.vector(tapply(ftree_lc1$fc,ftree_lc1$LCT,getmode)),3)
round(as.vector(tapply(fbare_lc1$fc,fbare_lc1$LCT,getmode)),3)

round(as.vector(tapply(ftree_lc1$fc,ftree_lc1$LCT,quantile,0.975)),3)
round(as.vector(tapply(ftree_lc1$fc,ftree_lc1$LCT,quantile,0.025)),3)

round(as.vector(tapply(fbare_lc1$fc,fbare_lc1$LCT,quantile,0.975)),3)
round(as.vector(tapply(fbare_lc1$fc,fbare_lc1$LCT,quantile,0.025)),3)


model.file = "/home/volarex/Images/model_ftree_V5.txt" 
params=c("ftree_lc","fherb_lc","fbare_lc","BMrefT","BMrefH","sigma","cert_w","cert_b")
nbiter=30000

for(i in 1:15) {
    print(paste("start for LCT",LCT[i]))
    data=list(bb=as.vector(final[,i]), M = dim(final)[1],
              ftree_teta = prior$ftree_LSCE[i],         
              fbare_teta = prior$fbare_LSCE[i],         
              bfw = Bref$Breft_mu[i],
              bfh = Bref$Brefh_mu[i],
              bfw_s = Bref$Breft_sigma[i],
              bfh_s = Bref$Brefh_sigma[i],
              ESAup=prior$ftree_ESA_up[i],
              ESAlw=prior$ftree_ESA_lw[i],
              BAREup=prior$fbare_ESA_up[i],
              BARElw=prior$fbare_ESA_lw[i])

    assign(paste0("out.CWT",i),bugs(data,inits=NULL,params,model.file,codaPkg=TRUE, n.iter=nbiter,
                                    DIC=FALSE, n.chains=3,
                                    working.directory= paste0("/home/volarex/Images/certRand",LCT[i],"/"),
                                    restart=FALSE))
    print(paste("LCT",LCT[i],"....Done"))
}


coda.out50=read.bugs(out.CWT5)
dd50=gelman.diag(coda.out50)
print(dd50)

coda.out50=read.bugs(out.CWT15)
dd50=gelman.diag(coda.out50)
print(dd50)

listprior <- setprior(12,prior$ftree_LSCE,Bref,prior,LCTc)
pdf("/home/volarex/Images/Result_LCT_CWT_3.2_certRand_final_testtrue.pdf", width=36, height=20)
j=0; ddf=NULL;ddfs=NULL
for(i in LCT){
    j=j+1
    #out.CWT=get(paste0("out.CWT",j))
    out.CWT <- character(3)
    out.CWT[1] <- paste0("/home/volarex/Images/certRand",i,"/CODAchain1.txt")
    out.CWT[2] <- paste0("/home/volarex/Images/certRand",i,"/CODAchain2.txt")
    out.CWT[3] <- paste0("/home/volarex/Images/certRand",i,"/CODAchain3.txt")
    coda.out=read.bugs(out.CWT)
    print(i)
    dd <- as.data.frame(summary(coda.out, q=c(0.025,0.975))[[2]])
    dds <- as.data.frame(summary(coda.out, q=c(0.025,0.975))[[1]])
    dd$LCT <- i
    dds$LCT <- i
    print(dd)
    ddf <- rbind(ddf,dd)
    ddfs <- rbind(ddf,dd)
    coda.melt0_CWT <- melt(as.data.frame(coda.out[[1]]),variable.name ="PARAMETER")
    coda.melt1_CWT <- melt(as.data.frame(coda.out[[2]]),variable.name ="PARAMETER")
    coda.melt2_CWT <- melt(as.data.frame(coda.out[[3]]),variable.name ="PARAMETER")
    mode0 <- getmode(as.vector(subset(coda.melt0_CWT,PARAMETER=="ftree_lc",value)[,1]))
    mode1 <- getmode(as.vector(subset(coda.melt1_CWT,PARAMETER=="ftree_lc",value)[,1]))
    mode2 <- getmode(as.vector(subset(coda.melt2_CWT,PARAMETER=="ftree_lc",value)[,1]))
    post.tetha <- mean(c(mode0, mode1, mode2))
    post.tetha0 <-data.frame(PARAMETER="ftree_lc",mode=post.tetha)                   
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
write.table(ddf,"/home/volarex/Images/CI_posterior.csv", sep=";")
write.table(ddfs,"/home/volarex/Images/CI_posterior_normal.csv", sep=";")
dev.off()
