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
#### Load TimeANLS-PRE Signature and get ssGSEA for each dataset
####################################################################################
#TimeANLS-PRE Signature
load('output/Pathway_Singatures.Rdata')
TimeANLS.Sigs <- prepare_sig(Pathway.Sigs$TimeANLS_Sig)

#load Riaz Data (ipi-PROG-NAIVE Sample info & TPM)
load('data/Riaz_data.RData')
Riaz.Pre.Samples <- split(Riaz_data$Samples,Riaz_data$Samples$PreOn)$Pre
Riaz.Pre.TPM <- Riaz_data$TPM[,match(Riaz.Pre.Samples$Sample, colnames(Riaz_data$TPM))]
Riaz.TA.Pre.ssgsea <- gsva(expr = Riaz.Pre.TPM, gset.idx.list = TimeANLS.Sigs, method='ssgsea')

#load Gide Data
load('data/Gide_data.RData')
Gide.Pre.Samples <- split(Gide_data$Samples,Gide_data$Samples$PREEDT)$PRE
Gide.Pre.TPM <- Gide_data$TPM[,match(Gide.Pre.Samples$Sample, colnames(Gide_data$TPM))]
Gide.TA.Pre.ssgsea <- gsva(expr = Gide.Pre.TPM, gset.idx.list = TimeANLS.Sigs, method='ssgsea')

#load Lee Data
load('data/Lee_data.RData')
Lee.Pre.Samples <- Lee_data$Pre_Samples
Lee.Pre.TPM <- Lee_data$TPM[,match(Lee.Pre.Samples$Pt_ID, colnames(Lee_data$TPM))]
Lee.TA.Pre.ssgsea <- gsva(expr = as.matrix(Lee.Pre.TPM), gset.idx.list = TimeANLS.Sigs, method='ssgsea')

#load MGH Data
load('data/MGH_PRE_data.RData')
MGH.Pre.Samples <- MGH_PRE_data$Samples
Batch14_PRE_ssgsea <-  gsva(expr = as.matrix(MGH_PRE_data$Batch14), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
Batch17_PRE_ssgsea <- gsva(expr = as.matrix(MGH_PRE_data$Batch17), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
SN0119610_PRE_ssgsea <- gsva(expr = as.matrix(MGH_PRE_data$SN0119610), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
SN0123099_PRE_ssgsea <- gsva(expr = as.matrix(MGH_PRE_data$SN0123099), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
SN0131794_PRE_ssgsea <- gsva(expr = as.matrix(MGH_PRE_data$SN0131794), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
MGH.TA.Pre.ssgsea <- data.frame(cbind(Batch14_PRE_ssgsea, Batch17_PRE_ssgsea, SN0119610_PRE_ssgsea, SN0123099_PRE_ssgsea, SN0131794_PRE_ssgsea))
MGH.TA.Pre.ssgsea <- MGH.TA.Pre.ssgsea[,match(MGH.Pre.Samples$Sample, colnames(MGH.TA.Pre.ssgsea))]
all(colnames(MGH.TA.Pre.ssgsea) == MGH.Pre.Samples$Sample)

####################################################################################
#### Training, Validation & Testing Model
####################################################################################
#Train Model
set.seed(2678)
TimeANLS.PRE.fit <- cv.glmnet(x=t(Riaz.TA.Pre.ssgsea), y=as.factor(Riaz.Pre.Samples$Resp_NoResp), 
                              family="binomial",type.measure="auc",nfolds=3, alpha=0.5, offset = prior_correction(Riaz.Pre.Samples$Resp_NoResp, tau = 2/3))

# Validate on Riaz On ssGSEA
Riaz.TA.PRE.prob <- predict(TimeANLS.PRE.fit, newx=t(Riaz.TA.Pre.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(Riaz.Pre.Samples$Resp_NoResp, tau = 2/3))
performance(prediction(predictions = Riaz.TA.PRE.prob, labels = Riaz.Pre.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.8172043

###### Test on Gide Pre ssGSEA
Gide.TA.PRE.prob <- predict(TimeANLS.PRE.fit, newx=t(Gide.TA.Pre.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(Gide.Pre.Samples$Resp_NoResp, tau = 2/3))
performance(prediction(predictions = Gide.TA.PRE.prob, labels = Gide.Pre.Samples$Resp_NoResp),"auc")@y.values[[1]] #0.5967078

###### Test on Lee Pre ssGSEA
Lee.TA.PRE.prob <- predict(TimeANLS.PRE.fit, newx=t(Lee.TA.Pre.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(Lee.Pre.Samples$Resp_NoResp,tau=2/3))
performance(prediction(predictions = Lee.TA.PRE.prob, labels = Lee.Pre.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.4938017

###### Test on MGH(Full Genes) Pre ssGSEA
MGH.TA.PRE.prob <- predict(TimeANLS.PRE.fit, newx = t(MGH.TA.Pre.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(MGH.Pre.Samples$Resp_NoResp, tau = 2/3))
performance(prediction(predictions = MGH.TA.PRE.prob, labels = MGH.Pre.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.7564103

#Plot Model
tiff(filename = 'output/TimeANLS_PRE/TimeANLS_PRE_Model_Plot.tiff',width=5, height=5,unit='in',compression = 'lzw',res=150)
plot(TimeANLS.PRE.fit)
title('Elastic-Net Regularization CV AUC',line=2.5)
dev.off()

tiff(file="output/TimeANLS_PRE/TimeANLS_PRE_L3reg_plot_path.tiff",width=8, height=6.5,unit='in',compression = 'lzw',res=150)
plot_glmnet(glmnet(x=t(Riaz.TA.Pre.ssgsea), y=as.factor(Riaz.Pre.Samples$Resp_NoResp), family="binomial", alpha=0.5),xvar='lambda',label=15)
abline(v=log(TimeANLS.PRE.fit$lambda.1se), col="blue", lty=2)
text(x=-2.8, y=-20, label = 'log(lambda.1se)\n-3.787388\nNonzero Coeff 11',cex=0.8) #-3.787388
title('Elastic-Net Regularization solution paths',line=3)
dev.off()


coeft.mat <- as.matrix(coef(TimeANLS.PRE.fit, s='lambda.1se'))
coeft.mat <- coeft.mat[c("(Intercept)",
                         "INITIAL_TRIGGERING_OF_COMPLEMENT",
                         "SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_7ALPHA_HYDROXYCHOLESTEROL",
                         "RIPK1_MEDIATED_REGULATED_NECROSIS",
                         "CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
                         "PEROXISOMAL_LIPID_METABOLISM",
                         "OTHER_INTERLEUKIN_SIGNALING",
                         "PEROXISOMAL_PROTEIN_IMPORT",
                         "NICOTINATE_METABOLISM",
                         "ARACHIDONIC_ACID_METABOLISM",
                         "ZBP1_DAI_MEDIATED_INDUCTION_OF_TYPE_I_IFNS",
                         "ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS"),]
coeft.mat <- round(coeft.mat,5)
names(coeft.mat)[1] <- 'Intercept'

coeft.df <- data.frame(cbind(names(coeft.mat),coeft.mat))
colnames(coeft.df) <- c("Pathway","Weight")
coeft.df

#coeft.df
tiff(file="output/TimeANLS_PRE/TimeANLS_PRE_Coefficient.tiff",width=12, height=6.5,unit='in',compression = 'lzw',res=150)
tab <- ggtexttable(coeft.df, rows = NULL, theme = ttheme("mOrange"))
tab %>%
  tab_add_title(text = "TimeANLS-PRE model", face = "bold" ,size = 11, padding = unit(0.3, "line"))
dev.off()

####################################################################################
#### Output Result
####################################################################################
write.csv(Riaz.TA.Pre.ssgsea, file='output/TimeANLS_PRE/ssGSEA_Score/Riaz_TA_PRE_ssgsea.csv')
write.csv(Gide.TA.Pre.ssgsea, file='output/TimeANLS_PRE/ssGSEA_Score/Gide_TA_PRE_ssgsea.csv')
write.csv(Lee.TA.Pre.ssgsea, file='output/TimeANLS_PRE/ssGSEA_Score/Lee_TA_PRE_ssgsea.csv')
write.csv(MGH.TA.Pre.ssgsea, file='output/TimeANLS_PRE/ssGSEA_Score/MGH_TA_PRE_ssgsea.csv')

TA_PRE_ss <- list('Riaz' = Riaz.TA.PRE.prob,
                  'Gide' = Gide.TA.PRE.prob,
                  'Lee' = Lee.TA.PRE.prob,
                  'MGH' = MGH.TA.PRE.prob
                )
save(TA_PRE_ss, file='output/TimeANLS_PRE/TA_PRE_SigScore.Rdata')

#Save Model
save(TimeANLS.PRE.fit, file='output/TimeANLS_PRE/TA_PRE_Model_Fit.Rdata')

