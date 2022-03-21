## Klasifikace samplu (bedgraph z Agilent MethylSeq pomoci mnp klasifikatoru - jen single file
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
# VSTUPEM JE BEDGRAPH FILE
if (length(args)!=2) {
  stop("Three arguments must be supplied: called methylations .bedgraph.gz and intersected methylations .bedgraph .n", call.=FALSE)
}


##TEST
# setwd("/home/rj/ownCloud/WORK/mnp_ceitec_MetylaceKlasifikace/V2")
# args <- c("VK2056.bedGraph.gz","INTERSECTED.bed","VK2056_klasifikace.xlsx")


################################################################################

library(writexl)
library(mnp.v11b6)

load("10k_probes.Rdata")


bedgraph_file <- args[1]
intersected_file <- args[2]

################################################################################
# PART 1 : INTERSECTING WITH EXACTLY POSTIONED CpGs (those called by Bismark and overlapping with IlmnID MAPINFO positon )

# NEW DATA FRAME
cgID10k <-probesMETAGI

# BEDGRAPH FILE LOADING
bedgraphfile=gzfile(bedgraph_file,'rt')
bedgraph <- read.delim(bedgraphfile, header = FALSE, fill = TRUE,skip = 1)
headerC <- c("chromosome","start","stop", "methylationPercentage")
colnames(bedgraph) <- headerC
close(bedgraphfile)

# GENOME AVERAGE BETA VALUE
genome_average <- mean(bedgraph$methylationPercentage,na.rm=TRUE)

#BETA VALUES MINING FOR 10K ILLUMINA PROBES
for (i in 1:nrow(cgID10k))
{
  seqname<- paste("chr",cgID10k$CHR[i], sep="") #chromozome number
  tempC <- bedgraph[bedgraph$chromosome == seqname,] # probes v chromozomu
  x <- which(cgID10k$MAPINFO[i]==tempC$start) #z start a stop se rovnaji
  
  if (!length(x)==0){cgID10k$BetaValueCOV[i]<-tempC$methylationPercentage[x]}
  else{cgID10k$BetaValueCOV[i]<-NA}
}

print("Missing probes: ")
sum(is.na(cgID10k$BetaValueCOV))

################################################################################
# PART 2 : INTERSECTING WITH EXTENDED CpGs calls (+-3000bp to both ways from a IllmID)

# EXTENDED BEDGRAPH FILE LOADING
bedgraph <- read.delim(intersected_file, header = FALSE, fill = FALSE)
headerC <- c("chromosome1","start1","stop1", "IlmnID","chromosome2","start2","stop2","methylationPercentage")
colnames(bedgraph) <- headerC
# list of probes with NA beta values
missing_probes <- cgID10k[is.na(cgID10k$BetaValueCOV),"IlmnID"]
# cycle trhough missing_probes and get average value for this probe
for (i in 1:length(missing_probes))
{
  probe <- missing_probes[[i]]
  row <- which(cgID10k$IlmnID==probe)
  bedgraph[which(bedgraph$IlmnID==probe),c("methylationPercentage")]
  
  if( is.numeric( bedgraph[which(bedgraph$IlmnID==probe),c("methylationPercentage")] ) | length(bedgraph[which(bedgraph$IlmnID==probe),c("methylationPercentage")])>1 )
  {
    average <- mean(as.numeric ( bedgraph[which(bedgraph$IlmnID==probe),c("methylationPercentage")] ),  na.rm = TRUE) 
    cgID10k$BetaValueCOV[row] <- average
  }
}
print("Missing probes after extended CpG regions : ")
sum(is.na(cgID10k$BetaValueCOV))
################################################################################
# PART 3 : ASSIGN a constant beta value to missing probes
cgID10k[is.na(cgID10k$BetaValueCOV),"BetaValueCOV"] <- genome_average
print("Missing probes after imputing genome average methylation value : ")
sum(is.na(cgID10k$BetaValueCOV))

################################################################################
# BETA MATRIX do KLASIFIKATORU
betamatrix<-cgID10k
rownames(betamatrix) <- betamatrix$IlmnID
betamatrix$IlmnID <- NULL
betamatrix$CHR <- NULL
betamatrix$MAPINFO <- NULL
betamatrix$BEDstart <- NULL
betamatrix$BEDend <- NULL
betamatrix$Agilent <- NULL

# betamatrix[is.na(betamatrix)] <- 0
betamatrix <- as.matrix(betamatrix)
betamatrix<-betamatrix/100

#MNP CLASSIFIER
RGsetMETAGI <- MNPpredict_betas(betamatrix, calibrate=TRUE, type="prob", MCF=FALSE)

################################################################################

# TRANSFORM TO EXCEL SHEET
klasifikace<-as.data.frame(RGsetMETAGI)
klasifikace<-t(klasifikace)
klasifikace <- cbind(rownames(klasifikace), data.frame(klasifikace, row.names=NULL))
colnames(klasifikace)<-c("MethylationClass","ClassifierScore")
klasifikace<-klasifikace[order(klasifikace$ClassifierScore, decreasing=TRUE),]

sampleName<-strsplit(bedgraph_file, "[.]")[[1]][1]
write_xlsx(as.data.frame(klasifikace), sprintf("/home/data/%s_mnp_klasifikace.xlsx", sampleName))

message(sprintf("%s DONE", sampleName))

