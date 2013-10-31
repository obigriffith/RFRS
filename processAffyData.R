#Install necessary packages with Bioconductor
#Note: CustomCDF packages may not install on mac osx. 
#Recommend running this on a large memory linux box

#Load the appropriate libraries
library(affy)
library(gcrma)
library(hgu133ahsentrezgcdf) #cdfname="HGU133A_HS_ENTREZG"
library(hgu133ahsentrezgprobe)
library(hgu133ahsentrezg.db)
library(GEOquery)

#Set working directory
setwd("~/Rworkflow")

#Download CEL files for desired GSM ids
dir.create("CEL")
GSMlist=read.table(file="GSM_list.txt")
apply(GSMlist, 1, getGEOSuppFiles, makeDirectory = FALSE, baseDir = "CEL")
cels = list.files("CEL/", pattern = "CEL", ignore.case=TRUE)
sapply(paste("CEL", cels, sep="/"), gunzip)
cels = list.files("CEL/", pattern = "CEL", ignore.case=TRUE)

#Set CDF to use:
cdf="HGU133A_HS_ENTREZG"

#Read in the data
raw.data.ALL=ReadAffy(verbose=TRUE, celfile.path="CEL", cdfname=cdf)

##GCRMA normalization
data.gcrma.norm.ALL=gcrma(raw.data.ALL)

##Get the important stuff out of the data - the expression estimates for each array
gcrma.ALL=exprs(data.gcrma.norm.ALL)

#Remove control probes
gcrma.ALL=gcrma.ALL[1:12065,] #Remove Affy control probes, custom CDF

#Format to 5 decimal places
gcrma.ALL=format(gcrma.ALL, digits=5)

#Map probes to gene symbols
#To see all mappings for Entrez gene db associated with customCDF
ls("package:hgu133ahsentrezg.db") #customCDF
probes.ALL=row.names(gcrma.ALL)
symbol.ALL = unlist(mget(probes.ALL, hgu133ahsentrezgSYMBOL))
ID.ALL = unlist(mget(probes.ALL, hgu133ahsentrezgENTREZID))
gcrma.ALL=cbind(probes.ALL,ID.ALL,symbol.ALL,gcrma.ALL)

#Write GCRMA-normalized, mapped data to file
write.table(gcrma.ALL, file = "ALL_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
