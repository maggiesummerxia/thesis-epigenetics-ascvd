######################################################
#Install packages needed
######################################################

if (!requireNamespace("minfi", quietly = TRUE))
  BiocManager::install("minfi")
if (!requireNamespace("ENmix", quietly = TRUE))
  BiocManager::install("ENmix")
if (!requireNamespace("ENmix", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylationEPICmanifest")
if (!requireNamespace("ENmix", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")


######################################################
#Load libraries
######################################################

rm(list=ls())
library(dplyr)
library(minfi)
library(ENmix)

######################################################
#Edited ENmix's QCinfo function so that images are saved as PDF in the HPC (no jpeg/cairo support)
######################################################

QCinfo <- function(rgSet,detPthre=0.000001,detPtype="negative",nbthre=3,samplethre=0.05,CpGthre=0.05,
                   bisulthre=NULL,outlier=TRUE,distplot=TRUE)
{
  
  ##number of bead
  if(!is(rgSet, "rgDataSet") & !is(rgSet, "RGChannelSetExtended"))
    stop("[QCinfo] The input should be an object of 'rgDataSet' or 'RGChannelSetExtended'")
  
  if(is(rgSet, "rgDataSet")){
    cginfo=getCGinfo(rgSet)
    typeI <- cginfo[cginfo$Infinium_Design_Type %in% c("I","snpI"),]
    typeIred=typeI[typeI$Color_Channel=="Red",]
    typeIgrn=typeI[typeI$Color_Channel=="Grn",]
    typeII<-cginfo[cginfo$Infinium_Design_Type %in% c("II","snpII"),]
    locusNames=c(typeIred$Name,typeIgrn$Name,typeII$Name)
    ##detection P value
    detP<-calcdetP(rgSet,detPtype=detPtype)
    ctrls<-metadata(rgSet)$ictrl
  }else if(is(rgSet, "RGChannelSetExtended")){
    typeI <- getProbeInfo(rgSet, type = "I")
    typeII <- getProbeInfo(rgSet, type = "II")
    locusNames <- getManifestInfo(rgSet, "locusNames")
    ##detection P value
    detP<-detectionP(rgSet)
    ctrls<-getProbeInfo(rgSet,type="Control")
  }
  
  bc_I <- assays(rgSet)$NBeads[typeI$AddressA,]
  flag<-bc_I>assays(rgSet)$NBeads[typeI$AddressB,]
  bc_I[flag] <- assays(rgSet)$NBeads[typeI$AddressB,][flag];
  nbead <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                  dimnames = list(locusNames, colnames(rgSet)))
  nbead[typeI$Name,]<-bc_I
  nbead[typeII$Name,]<-assays(rgSet)$NBeads[typeII$AddressA,]
  rm(list=c("bc_I","flag"))
  
  if(!identical(rownames(detP),rownames(nbead))){
    idx=intersect(rownames(detP),rownames(nbead))
    detP=detP[idx,]
    nbead=nbead[idx,]
  }
  if(!identical(colnames(detP),colnames(nbead))){
    idx=intersect(colnames(detP),colnames(nbead))
    detP=detP[,idx]
    nbead=nbead[,idx]
  }
  
  ## bisulfite conversion internal controls
  ctrls=ctrls[ctrls$Address %in% rownames(rgSet),]
  ctrl_r <- assays(rgSet)$Red[ctrls$Address,]
  ctrl_g <- assays(rgSet)$Green[ctrls$Address,]
  cc=ctrls[(ctrls$Type %in% c("BISULFITE CONVERSION I")) & (ctrls$ExtendedType 
                                                            %in% c("BS Conversion I C1","BS Conversion I-C2","BS Conversion I-C3")),]
  I_green=colMeans(ctrl_g[cc$Address,])
  cc=ctrls[(ctrls$Type %in% c("BISULFITE CONVERSION I")) & (ctrls$ExtendedType
                                                            %in% c("BS Conversion I-C4","BS Conversion I-C5","BS Conversion I-C6")),]
  I_red=colMeans(ctrl_r[cc$Address,])
  cc=ctrls[ctrls$Type %in% c("BISULFITE CONVERSION II"),]
  II_red=colMeans(ctrl_r[cc$Address,])
  bisul=(I_green+I_red+II_red)/3
  
  #threshold of bisulfite conversion control intensity
  if(is.null(bisulthre)){bisulthre=mean(bisul,na.rm=TRUE)-3*sd(bisul,na.rm=TRUE)}
  
  ##low quality samples
  qcmat <- nbead<nbthre | detP>detPthre
  badValuePerSample <- apply(qcmat,2,sum)/nrow(qcmat)
  flag <- badValuePerSample > samplethre | bisul < bisulthre
  cat(sum(flag)," samples with percentage of low quality CpG value greater than ",
      samplethre, " or bisulfite intensity less than ", bisulthre, "\n")
  badsample=colnames(qcmat)[flag]
  
  ##low quality CpGs
  qcmat <- qcmat[,!flag]
  NbadValuePerCpG <- apply(qcmat,1,sum)
  badValuePerCpG <- NbadValuePerCpG/ncol(qcmat)
  flag2 <- badValuePerCpG>CpGthre & NbadValuePerCpG>1
  cat(sum(flag2)," CpGs with percentage of low quality value greater than ",
      CpGthre,"\n")
  badCpG <- rownames(qcmat)[flag2]
  qcmat=qcmat[!flag2,]
  
  #plotting quality scores
  cat("Ploting qc_sample.jpg ...")
  pdf(file="qc_sample.pdf")
  color=rep("black",length(badValuePerSample))
  color[flag]="red"
  plot(badValuePerSample,bisul,xlab="Percent of low quality data per sample",
       ylab="Average bisulfite conversion intensity",cex=1.5,col=color,
       main=paste(length(badsample)," samples were classified as low quality samples"))
  abline(h=bisulthre,lty=2,col="red")
  abline(v=samplethre,lty=2,col="red")
  dev.off()
  cat("Done\n")
  
  cat("Ploting qc_CpG.jpg ...")
  pdf(file="qc_CpG.pdf")
  par(mfrow=c(2,1))
  hist(badValuePerCpG,breaks=1000,xlab="Percent of low quality data per CpG",
       main=paste(length(badCpG)," CpGs were classified as low quality CpGs"))
  abline(v=CpGthre,lty=2,col="red")
  hist(badValuePerCpG,breaks=1000,xlim=c(0,0.1),
       xlab="Percent of low quality data per CpG",main="Zoom in view")
  abline(v=CpGthre,lty=2,col="red")
  dev.off()
  cat("Done\n")
  
  #Identifying outlier samples
  if(outlier)
  {
    cat("Identifying ourlier samples based on beta or total intensity values...\n")
    
    if(is(rgSet, "rgDataSet")){mdat=getmeth(rgSet)
    }else if(is(rgSet, "RGChannelSetExtended")){mdat=preprocessRaw(rgSet)}
    
    mdat=mdat[rownames(qcmat),]
    mdat=mdat[,colnames(qcmat)]
    #outliers based on total intensity values
    mu <- assays(mdat)$Meth+assays(mdat)$Unmeth
    mumean=apply(mu,2,mean,na.rm=TRUE)
    q2575 <- quantile(mumean, probs=c(0.25,0.75), na.rm=TRUE)
    qr <- q2575["75%"]-q2575["25%"]
    low=q2575["25%"] - 3*qr
    flag1=  mumean < low
    cat("After excluding low quality samples and CpGs\n")
    cat(sum(flag1)," samples are outliers based on averaged total intensity value","\n")
    
    #outliers in beta value distribution
    beta=getB(mdat, type="Illumina")
    qq=apply(beta,2,function(x) quantile(x, probs=c(0.25,0.5,0.75), na.rm=TRUE))
    q2575 <- apply(qq,1,function(x) quantile(x, probs=c(0.25,0.75), na.rm=TRUE))
    qr <- q2575["75%",]-q2575["25%",]
    up=q2575["75%",] + 3*qr
    low=q2575["25%",] - 3*qr
    flag=qq > up | qq < low
    flag=apply(flag,2,sum)>0
    cat(sum(flag)," samples are outliers in beta value distribution","\n")
    flag=flag | flag1
    badsample=c(badsample,colnames(beta)[flag])
    outlier_sample=colnames(beta)[flag]
    cat(sum(flag)," outlier samples were added into badsample list\n")
    if(ncol(qcmat)<15){
      cat("WARNING: Sample size may be too small to correctly identify outlier samples!\n")
      cat("RECOMMAND: set outlier=FALSE or double check total intensity and beta value 
     distribution plots to confirm\n")}
  }
  
  if(distplot)
  {
    
    if(is(rgSet, "rgDataSet")){mdat=getmeth(rgSet)
    }else if(is(rgSet, "RGChannelSetExtended")){mdat=preprocessRaw(rgSet)}
    
    beta=getB(mdat, type="Illumina")
    
    cat("Ploting freqpolygon_beta_beforeQC.jpg ...")
    pdf(file="freqpolygon_beta_beforeQC.jpg")
    color=rep("black",ncol(beta))
    color[colnames(beta) %in% badsample]="red"
    multifreqpoly(beta,cex.lab=1.4,cex.axis=1.5, col=color, legend="",
                  cex.main=1.5,main="Beta value distribution",
                  xlab="Methylation beta value")
    dev.off()
    cat("Done\n")
    
    beta=beta[!(rownames(beta) %in% badCpG),]
    beta=beta[,!(colnames(beta) %in% badsample)]
    cat("Ploting freqpolygon_beta_afterQC.jpg ...")
    pdf(file="freqpolygon_beta_afterQC.jpg")
    multifreqpoly(beta,cex.lab=1.4,cex.axis=1.5, col="black",legend="",
                  cex.main=1.5,main="Beta value distribution",
                  xlab="Methylation beta value")
    dev.off()
    cat("Done\n")
  }
  
  if(outlier)
  {list(detP=detP,nbead=nbead,bisul=bisul,badsample=badsample,badCpG=badCpG,
        outlier_sample=outlier_sample)
  }else{list(detP=detP,nbead=nbead,bisul=bisul,badsample=badsample,badCpG=badCpG)}
}
getB <- function(mdat,type="Illumina",offset=100)
{
  if(!is(mdat, "methDataSet") & !is(mdat, "MethylSet"))
  {stop("The input is not an object of methDataSet or MethylSet\n")}
  if(type=="Illumina"){offset=100}
  beta<-assays(mdat)$Meth/(assays(mdat)$Meth+assays(mdat)$Unmeth+offset)
  beta
}

######################################################
#Run QC! 
######################################################

#Processing of iDATs - as performed by Juan
#iDATs wd
setwd("/mnt/lustre/users/k1815935/FoodChallenge_data/iDATs")
#Import iDATs
basenames <- scan("./idat_list.csv","character")
#basenames <- basenames[c(1:75)] #(test the script with only a few samples)#
rgSet <- read.metharray(basenames=basenames, extended=TRUE, force=TRUE)
#Set wd for processed files
setwd("../iDATs_processing/")

#QC - as performed by Juan - does not work with the whole dataset on the local machine
#checks bisulfite intensity, detection p-value, nbeads, rates of missing data, intensity outliers, beta distribution outliers
qc <- QCinfo(rgSet, distplot=FALSE) # if needed: > BiocManager::install("IlluminaHumanMethylationEPICmanifest")
#checks median methylated and unmethylated signals, samples with medians below 10.5 wiil be excluded
mraw <- preprocessRaw(rgSet) #generates methylated and unmethylated signals
qc1 <- getQC(mraw) #checks median methylated and unmethylated signals, samples with medians below 10.5 will be excluded
badSampleCutoff <- 10.5
meds <- (qc1$mMed + qc1$uMed)/2
whichBad <- which((meds < badSampleCutoff))
badMedianIntensity <- rownames(qc1)[whichBad]
#List of samples and CpGs to be excluded
badsamples <- unique(c(qc$badsample, badMedianIntensity))
qc$badsample
badMedianIntensity
badsamples
badcpgs <- qc$badCpG
badcpgs

#Plot samples with median intensity below sample cut-off
pdf("badsamplecutoff.pdf")
plotQC(qc1, badSampleCutoff = 10.5)
dev.off()

#Sex predcition - as performed by Juan
gmraw <- mapToGenome(mraw)
predictedSex <- getSex(gmraw, cutoff = -2)$predictedSex
predictedSex

#Retrieve signals from SNP probes - as performed by Juan
snps <- getSnpBeta(rgSet)
rownames(snps)

#Normalisation, set to NA bad quality signals, exclusion of bad samples and CpGs - as performed by Juan
mset <- preprocessENmix(rgSet,bgParaEst="oob",dyeCorr=TRUE,nCores=20,QCinfo=qc)
mdat <- norm.quantile(mset,method="quantile1")
beta.rcp <- rcp(mdat,qcscore=qc)
detP <- qc$detP[!rownames(qc$detP) %in% qc$badCpG,!colnames(qc$detP) %in% qc$badsample]
nbead <- qc$nbead[!rownames(qc$nbead) %in% qc$badCpG,!colnames(qc$nbead) %in% qc$badsample]
beta.rcp[detP>0.000001 | nbead<3]=NA
beta.rcp <- beta.rcp[!rownames(beta.rcp) %in% badcpgs,!colnames(beta.rcp) %in% badsamples]
nrow(beta.rcp)
ncol(beta.rcp)
write.csv(beta.rcp,"DIMENSION_EPIC_beta_enmix_whole_blood.csv",quote=F)
sessionInfo()

#Done!
#RC 17/11/2020