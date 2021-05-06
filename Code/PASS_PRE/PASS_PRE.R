library(here)
library(GSVA)
library(glmnet)
library(pROC)
library(ROCR) #probably disable later on
library(plotmo)
library(ggpubr)
source('code/FE_Model_Functions.R') #Feature Engineering & Modeling
source('code/Model_Development/Model_Devep_Fun.R')
####################################################################################
#### Load PASS-PRE Signature and get ssGSEA for each dataset
####################################################################################
#PASS-PRE Signature
load('output/Pathway_Singatures.Rdata')
PASS.PRE.Sigs <- prepare_sig(Pathway.Sigs$PASS_PRE)

#load Riaz Data
load('data/Riaz_data.RData')
Riaz.Pre.Samples <- split(Riaz_data$Samples,Riaz_data$Samples$PreOn)$Pre
Riaz.Pre.TPM <- Riaz_data$TPM[,match(Riaz.Pre.Samples$Sample, colnames(Riaz_data$TPM))]
Riaz.Pre.ssgsea <- gsva(expr = Riaz.Pre.TPM, gset.idx.list = PASS.PRE.Sigs, method='ssgsea')

#load Gide Data
load('data/Gide_data.RData')
Gide.Pre.Samples <- split(Gide_data$Samples,Gide_data$Samples$PREEDT)$PRE
Gide.Pre.TPM <- Gide_data$TPM[,match(Gide.Pre.Samples$Sample, colnames(Gide_data$TPM))]
Gide.Pre.ssgsea <- gsva(expr = Gide.Pre.TPM, gset.idx.list = PASS.PRE.Sigs, method='ssgsea')

#load Lee Data
load('data/Lee_data.RData')
Lee.Pre.Samples <- Lee_data$Pre_Samples
Lee.Pre.TPM <- Lee_data$TPM[,match(Lee.Pre.Samples$Pt_ID, colnames(Lee_data$TPM))]
Lee.Pre.ssgsea <- gsva(expr = as.matrix(Lee.Pre.TPM), gset.idx.list = PASS.PRE.Sigs, method='ssgsea')

#load MGH Data
load('data/MGH_PRE_data.RData')
MGH.Pre.Samples <- MGH_PRE_data$Samples
Batch14_PRE_ssgsea <-  gsva(expr = as.matrix(MGH_PRE_data$Batch14), gset.idx.list = PASS.PRE.Sigs, method='ssgsea')
Batch17_PRE_ssgsea <- gsva(expr = as.matrix(MGH_PRE_data$Batch17), gset.idx.list = PASS.PRE.Sigs, method='ssgsea')
SN0119610_PRE_ssgsea <- gsva(expr = as.matrix(MGH_PRE_data$SN0119610), gset.idx.list = PASS.PRE.Sigs, method='ssgsea')
SN0123099_PRE_ssgsea <- gsva(expr = as.matrix(MGH_PRE_data$SN0123099), gset.idx.list = PASS.PRE.Sigs, method='ssgsea')
SN0131794_PRE_ssgsea <- gsva(expr = as.matrix(MGH_PRE_data$SN0131794), gset.idx.list = PASS.PRE.Sigs, method='ssgsea')

MGH.Pre.ssgsea <- data.frame(cbind(Batch14_PRE_ssgsea, Batch17_PRE_ssgsea, SN0119610_PRE_ssgsea, SN0123099_PRE_ssgsea, SN0131794_PRE_ssgsea))
MGH.Pre.ssgsea <- MGH.Pre.ssgsea[,match(MGH.Pre.Samples$Sample, colnames(MGH.Pre.ssgsea))]
all(colnames(MGH.Pre.ssgsea) == MGH.Pre.Samples$Sample)

####################################################################################
#### Training, Validation & Testing Model
####################################################################################
#Train Model
set.seed(1028)
PASS.PRE.fit <- cv.glmnet(x=t(Riaz.Pre.ssgsea), y=as.factor(Riaz.Pre.Samples$Resp_NoResp), 
                          family="binomial",type.measure="auc",nfolds=3, alpha=0.5, offset = prior_correction(Riaz.Pre.Samples$Resp_NoResp, tau = 3/5))

# Validate on Riaz On ssGSEA
Riaz.PRE.prob <- predict(PASS.PRE.fit, newx=t(Riaz.Pre.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(Riaz.Pre.Samples$Resp_NoResp, tau = 3/5))
performance(prediction(predictions = Riaz.PRE.prob, labels = Riaz.Pre.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.7293907

# Test on Gide On ssGSEA
Gide.PRE.prob <- predict(PASS.PRE.fit, newx=t(Gide.Pre.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(Gide.Pre.Samples$Resp_NoResp, tau = 3/5)) #rep(0,72)
performance(prediction(predictions = Gide.PRE.prob, labels = Gide.Pre.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.6814815

# Test on Lee On ssGSEA
Lee.PRE.prob <- predict(PASS.PRE.fit, newx=t(Lee.Pre.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(Lee.Pre.Samples$Resp_NoResp,tau=3/5)) #rep(0,44)
roc(Lee.Pre.Samples$Resp_NoResp, Lee.PRE.prob)$auc #Area under the curve: 0.5496

# Test on MGH(Full Genes) On ssGSEA
MGH.PRE.prob <- predict(PASS.PRE.fit, newx = t(MGH.Pre.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(MGH.Pre.Samples$Resp_NoResp, tau = 3/5)) #rep(0, 19)
performance(prediction(predictions = MGH.PRE.prob, labels = MGH.Pre.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.6923077

#Plot Model
tiff(filename = 'output/PASS_PRE/PASS_PRE_Model_Plot.tiff',width=5, height=5,unit='in',compression = 'lzw',res=150)
plot(PASS.PRE.fit)
title('Elastic-Net Regularization CV AUC',line=2.5)
dev.off()

tiff(file="output/PASS_PRE/PASS_PRE_L3reg_plot_path.tiff",width=8, height=6.5,unit='in',compression = 'lzw',res=150)
plot_glmnet(glmnet(x=t(Riaz.Pre.ssgsea), y=as.factor(Riaz.Pre.Samples$Resp_NoResp), family="binomial", alpha=0.5),xvar='lambda',label=15)
abline(v=log(PASS.PRE.fit$lambda.1se), col="blue", lty=2)
text(x=-3, y=-20, label = 'log(lambda.1se)\n-1.850035\nNonzero Coeff 6',cex=0.8) #-1.850035
title('Elastic-Net Regularization solution paths',line=3)
dev.off()


coeft.mat <- as.matrix(coef(PASS.PRE.fit, s='lambda.1se'))
coeft.mat <- coeft.mat[c("(Intercept)",
                         "COMPLEMENT_CASCADE",
                         "REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
                         "BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS",
                         "PLASMA_LIPOPROTEIN_ASSEMBLY",
                         "INTERLEUKIN_2_FAMILY_SIGNALING",
                         "RA_BIOSYNTHESIS_PATHWAY"),]
coeft.mat <- round(coeft.mat,5)
names(coeft.mat)[1] <- 'Intercept'
#names(coeft.mat)[c(1,3)] <- c('Intercept', "REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT\n_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_\nFACTOR_BINDING_PROTEINS_IGFBPS")
coeft.df <- data.frame(cbind(names(coeft.mat),coeft.mat))
colnames(coeft.df) <- c("Pathway","Weight")
coeft.df

#coeft.df
tiff(file="output/PASS_PRE/PASS_PRE_Coefficient.tiff",width=15, height=6.5,unit='in',compression = 'lzw',res=150)
tab <- ggtexttable(coeft.df, rows = NULL, theme = ttheme("mOrange"))
tab %>%
  tab_add_title(text = "PASS-PRE model", face = "bold" ,size = 11, padding = unit(0.3, "line"))
dev.off()

####################################################################################
#### Output Result
####################################################################################
#Save SSGSEA Value
write.csv(Riaz.Pre.ssgsea, file='output/PASS_PRE/ssGSEA_Score/Riaz_PASS_PRE_ssgsea.csv')
write.csv(Gide.Pre.ssgsea, file='output/PASS_PRE/ssGSEA_Score/Gide_PASS_PRE_ssgsea.csv')
write.csv(Lee.Pre.ssgsea, file='output/PASS_PRE/ssGSEA_Score/Lee_PASS_PRE_ssgsea.csv')
write.csv(MGH.Pre.ssgsea, file='output/PASS_PRE/ssGSEA_Score/MGH_PASS_PRE_ssgsea.csv')

#Save Sig Score
PASS_PRE_ss <- list('Riaz' = Riaz.PRE.prob,
                   'Gide' = Gide.PRE.prob,
                   'Lee' = Lee.PRE.prob,
                   'MGH' = MGH.PRE.prob
                   )
save(PASS_PRE_ss, file='output/PASS_PRE/PASS_PRE_SigScore.Rdata')

#Save Model
save(PASS.PRE.fit, file='output/PASS_PRE/PASS_PRE_Model_Fit.Rdata')


