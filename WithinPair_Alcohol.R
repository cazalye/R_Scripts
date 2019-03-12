#################################################################  
# 	 Within-Pair Analysis: Alcohol Phenotype		#
#################################################################


## To run on server ##
# cd [...path/file]
# grun -n [job name] -q [queue].q [path to R version] CMD BATCH [script name]




########################
##### START SCRIPT ##### 
########################
 
####################    Load Libraries etc    ###################     
library(minfi)
setwd("...") 
#################################################################

##############  Alcohol Pheno: Complete MZ Pairs  ###############  
pheno <- read.csv(file="pheno.csv", check.names=F)
dim(pheno) 
MZ <- pheno[which(pheno$Zygosity==1),]
dim(MZ) 
length(which(table(MZ$Family_ID)==2)) # how many are complete pairs?
MZ <-  MZ[which(MZ$Family_ID %in% names(which(table(MZ$Family_ID)==2))), ] # subset by complete pairs
dim(MZ)
#################################################################

###############    Methylation data for MZ Pairs   ##############  
load("methylation.RData")
beta <- methylation[, which(colnames(methylation) %in% MZ$Barcode)]  	# subset methylation by MZ pairs
beta <- t(beta[order(col.names(beta)),])				# order by barcode and transpose		
pheno <- MZ[which(MZ$Barcode %in% rownames(beta)),] 	# Subset MZ pheno by good QC samples 
row.names(pheno) <- pheno$Barcode
pheno <- pheno[order(row.names(pheno)),]
identical(row.names(pheno), row.names(beta))  # should be TRUE
################################################################

##################      Calculate Deltas     ###################  
pheno_beta <- merge(pheno, beta, by="row.names")  
pheno_beta <- pheno_beta[order(pheno_beta$Family_ID),]

# split into 2 dataframes 
A <- pheno_beta[seq(1,nrow(pheno),2),] 
B <- pheno_beta[-seq(1,nrow(pheno),2),]

# Calculate delta for variables of interest: CpGs + Alcohol + BMI + Smoking + CD8T + CD4T + NK + Bcell + Mono 
delta <- data.frame(A$Family_ID, A[10:12]-B[10:12], A[14:18]-B[14:18], A[30:ncol(pheno_beta)]-B[30:ncol(pheno_beta)]) 
colnames(delta)[1] <- "Family_ID"
rownames(delta) <- delta$Family_ID
################################################################

################   FUNCTION: Within Pair Model   ###############   
# while the EWAS model uses lmer from the lme4 package, here we do not adjust for family or zygosity so no need for a mixed model with random effects, instead just use lm from base R.

gDay_wp = function(i) {
	pheno$cpg = beta[,i]          
	cpg = colnames(beta)[i]
	result = data.frame()
	fit <- lm(cpg ~ Alcohol + BMI + Smoking + CD8T + CD4T + NK + Bcell + Mono, data=pheno)
	Estimate = summary(fit)$coefficients[2,1]
	SE= summary(fit)$coefficients[2,2]
	T= summary(fit)$coefficients[2,3]
	p = 2 * (1 - pnorm(abs(summary(fit)$coefficients[2,3]))) 
	Mean= mean(pheno$cpg, na.rm=TRUE) 		
	Median= median(pheno$cpg, na.rm=TRUE) 
 	result=c(cpg, Estimate, SE, T, p, Mean, Median)
}
################################################################

#################     Run Model  &  QQ-plot    #################
pheno <- delta[,1:11]
beta  <- delta[,12:ncol(delta)]

options(mc.cores=10) 

result = do.call(rbind, mclapply(1:NCOL(beta), gDay_wp))
result = as.data.frame(result)
colnames(result) = c("CpG", "Estimate", "SE", "t_stat", "p_val", "Mean", "Median")
result$Adjusted_p <- p.adjust(result$p_val, method="fdr") 
result <- result[order(result$Adjusted_p),]
row.names(result) <- result$CpG
write.csv(result, file=".csv", row.names=F)

png("qqPlot.png")
qq(result$p_val, main = "...")  # Draws a QQ-plot from p-values, should be a straight line
dev.off()
estlambda(result$p_val)   # ought to be around 1
################################################################



##### END SCRIPT ##### 
rm(list=ls())	   ###
###################### 