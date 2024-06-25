BiocManager::install("CREAM")
#Primarily, it is used for calling clusters of cis-regulatory elements (COREs)
library(CREAM)
CREAM(in_path, WScutoff = 1.5, MinLength = 1000, peakNumMin = 2)
CREAM(system.file("extdata", "A549_Chr21.bed", package = "CREAM"),
      MinLength = 1000, peakNumMin = 2)
dat<-"/Users/user/Desktop/LacZ_MLL_peaks.bed"
dat1<-CREAM(dat,MinLength = 1000, peakNumMin = 2)
#using WScutoff=0.1 to output more LOCKs
dat1<-CREAM('LacZ_MEN1_MLL_peaks.broadPeak',WScutoff=0.01,MinLength = 1000, peakNumMin = 2)
dat1<-as.data.frame(dat1)
colnames(dat1)<-c('chr','start','end')
dat1$name<-paste0('LOCKs',1:nrow(dat1))
write.table(dat1,'Lacz_MLL_cream.bed',sep='\t',quote=F,col.names = F,row.names = F)

#ElementRecog is a function to identify COREs
ElementRecog(InputData, windowSize_Vec, peakNumMax, peakNumMin)

InputData <- read.table(system.file("extdata", "A549_Chr21.bed",package = "CREAM"), sep="\t")
colnames(InputData) <- c("chr", "start", "end")
MinLength <- 1000
if(nrow(InputData) < MinLength){
  stop(paste( "Number of functional regions is less than ", MinLength,
              ".", sep = "", collapse = ""))
}
peakNumMin <- 2
WScutoff <- 1.5
WindowVecFinal <- WindowVec(InputData, peakNumMin, WScutoff)
OutputList <- ElementRecog(InputData, WindowVecFinal,
                           (1+length(WindowVecFinal)), peakNumMin)
