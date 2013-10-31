library(randomForest)

#Set working directory and filenames for Input/output
setwd("C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/intellectual_property/appendices/RFRS/")

RF_model_file="RF_model_17gene_optimized"
datafile="patient_data.txt"
RelapseProbabilityPlotfile="RelapseProbabilityPlot.Rdata"
RelapseProbabilityFitfile="RelapseProbabilityFit.Rdata"
reportfile="patient_results.pdf"

#Load model file, choose (1) OR (2) and comment out the other (contains "rf_model" object)
RF_model_file="RF_model_17gene_optimized.Rdata" #1
#RF_model_file="RF_model_8gene_optimized.Rdata" #2
load(file=RF_model_file)

#Read in data (expecting a tab-delimited file with Gene Symbols as colnames and patient_id as rowname)
#	TOP2A	MAPKAPK2	SUPT4H1	FMOD	APOC1	FEN1	AURKA	TXNIP	EBP	CKS2	DTX4	SYNE2	DICER1	RACGAP1	AP1AR	NUP107	CCNB2
#GSM36893	7.0874	3.9958	7.6561	6.7689	10.268	8.8817	6.6811	8.3538	7.033	8.0512	6.0171	3.2419	6.272	10.0237	6.3404	8.9953	7.3143
patient_data=read.table(datafile, header = TRUE, row.names=1, na.strings = "NA", sep="\t")

#Run test data through forest
RF_predictions_response=predict(rf_model, patient_data, type="response")
RF_predictions_prob=predict(rf_model, patient_data, type="prob")
RFRS=RF_predictions_prob[,"Relapse"]

#Determine RFRS group according to previously determined thresholds
RF_risk_group=RF_predictions_prob[,"Relapse"]
RF_risk_group[RF_predictions_prob[,"Relapse"]<0.333]="low"
RF_risk_group[RF_predictions_prob[,"Relapse"]>=0.333 & RF_predictions_prob[,"Relapse"]<0.606]="int"
RF_risk_group[RF_predictions_prob[,"Relapse"]>=0.606]="high"

#Load Relapse probability plot, and loess fit to allow current patient to be plotted
load(file=RelapseProbabilityPlotfile)
load(file=RelapseProbabilityFitfile)
RelapseProb=predict(fit, RFRS)

pdf(file=reportfile)
replayPlot(RelapseProbabilityPlot)
points(x=RFRS, y=RelapseProb, pch=18, col="red",cex=2)
legend_text=c(paste("Patient: ",rownames(patient_data)),paste("RFRS =",round(RFRS,digits=4)),paste("risk group =",RF_risk_group),paste("Relapse prob. = ",round(RelapseProb,digits=1),"%",sep=""))
legend(x=0.6,y=11,legend=legend_text, bty="n",pch=c(18,NA,NA,NA),col=c("red",NA,NA,NA),pt.cex=2)
dev.off()


