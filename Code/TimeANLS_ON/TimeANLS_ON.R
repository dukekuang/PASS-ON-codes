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
#### Load TimeANLS-ON Signature and get ssGSEA for each dataset
####################################################################################
#TimeANLS-ON Signature
load('output/Pathway_Singatures.Rdata')
TimeANLS.Sigs <- prepare_sig(Pathway.Sigs$TimeANLS_Sig)

load('data/Riaz_data.RData')
Riaz.On.Samples <- split(Riaz_data$Samples,Riaz_data$Samples$PreOn)$On
Riaz.On.TPM <- Riaz_data$TPM[,match(Riaz.On.Samples$Sample, colnames(Riaz_data$TPM))]
Riaz.TA.On.ssgsea <- gsva(expr = Riaz.On.TPM, gset.idx.list = TimeANLS.Sigs, method='ssgsea')

#load Gide Data
load('data/Gide_data.RData')
Gide.On.Samples <- split(Gide_data$Samples,Gide_data$Samples$PREEDT)$EDT
Gide.On.TPM <- Gide_data$TPM[,match(Gide.On.Samples$Sample, colnames(Gide_data$TPM))]
Gide.TA.On.ssgsea <- gsva(expr = Gide.On.TPM, gset.idx.list = TimeANLS.Sigs, method='ssgsea')

#load Lee Data
load('data/Lee_data.RData')
Lee.On.Samples <- Lee_data$On_Samples
Lee.On.TPM <- Lee_data$TPM[,match(Lee.On.Samples$Pt_ID, colnames(Lee_data$TPM))]
Lee.TA.On.ssgsea <- gsva(expr = as.matrix(Lee.On.TPM), gset.idx.list = TimeANLS.Sigs, method='ssgsea')

#load MGH Genecode Data
load('data/MGH_ON_data.RData')
MGH.On.Samples <- MGH_ON_data$Samples
Batch14_ON_ssgsea <-  gsva(expr = as.matrix(MGH_ON_data$Batch14), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
Batch17_ON_ssgsea <- gsva(expr = as.matrix(MGH_ON_data$Batch17), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
SN0119610_ON_ssgsea <- gsva(expr = as.matrix(MGH_ON_data$SN0119610), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
SN0123099_ON_ssgsea <- gsva(expr = as.matrix(MGH_ON_data$SN0123099), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
SN0131794_ON_ssgsea <- gsva(expr = as.matrix(MGH_ON_data$SN0131794), gset.idx.list = TimeANLS.Sigs, method='ssgsea')
MGH.TA.On.ssgsea <- data.frame(cbind(Batch14_ON_ssgsea, Batch17_ON_ssgsea, SN0119610_ON_ssgsea, SN0123099_ON_ssgsea, SN0131794_ON_ssgsea))
MGH.TA.On.ssgsea <- MGH.TA.On.ssgsea[,match(MGH.On.Samples$Sample, colnames(MGH.TA.On.ssgsea))]

####################################################################################
#### Training, Validation & Testing Model
####################################################################################
#Train Model
set.seed(2678)
TimeANLS.On.fit <- cv.glmnet(x=t(Riaz.TA.On.ssgsea), y=as.factor(Riaz.On.Samples$Resp_NoResp), 
                             family="binomial",type.measure="auc",nfolds=3, alpha=0.5,
                             offset = prior_correction(Riaz.On.Samples$Resp_NoResp, tau=2/3)) #rep(0,54))

# Validate on Riaz On ssGSEA
Riaz.TA.On.prob <- predict(TimeANLS.On.fit, newx=t(Riaz.TA.On.ssgsea), type = "response", s="lambda.1se",newoffset = prior_correction(Riaz.On.Samples$Resp_NoResp,2/3))
roc(Riaz.On.Samples$Resp_NoResp, Riaz.TA.On.prob)$auc #Area under the curve: 0.8052
summary(Riaz.TA.On.prob)

# Test on Gide On ssGSEA
Gide.TA.On.prob <- predict(TimeANLS.On.fit, newx=t(Gide.TA.On.ssgsea), type = "response", s="lambda.1se",newoffset = prior_correction(Gide.On.Samples$Resp_NoResp, tau=2/3))
roc(Gide.On.Samples$Resp_NoResp, Gide.TA.On.prob)$auc #Area under the curve: 0.7532
summary(Gide.TA.On.prob)

# Test on Lee On ssGSEA
Lee.TA.On.prob <- predict(TimeANLS.On.fit, newx=t(Lee.TA.On.ssgsea), type = "response", s="lambda.1se",newoffset = prior_correction(Lee.On.Samples$Resp_NoResp_ResPROG, 2/3))
roc(Lee.On.Samples$Resp_NoResp_ResPROG, Lee.TA.On.prob)$auc #Area under the curve: 0.9138
summary(Lee.TA.On.prob)

# Test on MGH(Full Genes) On ssGSEA
MGH.TA.On.prob <- predict(TimeANLS.On.fit, newx = t(MGH.TA.On.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(MGH.On.Samples$Resp_NoResp, tau = 2/3))
roc(MGH.On.Samples$Resp_NoResp, MGH.TA.On.prob)$auc #Area under the curve: 0.8231
summary(MGH.TA.On.prob)

#Plot Model
tiff(filename = 'output/TimeANLS_ON/TimeANLS_ON_Model_Plot.tiff',width=5, height=5,unit='in',compression = 'lzw',res=150)
plot(TimeANLS.On.fit)
title('Elastic-Net Regularization CV AUC',line=2.5)
dev.off()

tiff(file="output/TimeANLS_ON/TimeANLS_ON_L3reg_plot_path.tiff",width=8, height=6.5,unit='in',compression = 'lzw',res=150)
plot_glmnet(glmnet(x=t(Riaz.TA.On.ssgsea), y=as.factor(Riaz.On.Samples$Resp_NoResp), family="binomial", alpha=0.5),xvar='lambda',label=15)
abline(v=log(TimeANLS.On.fit$lambda.1se), col="blue", lty=2)
text(x=-2.5, y=-100, label = 'log(lambda.1se)\n-1.020004\nNonzero Coeff 3',cex=0.8) #-1.020004
title('Elastic-Net Regularization solution paths',line=3)
dev.off()

coeft.mat <- as.matrix(coef(TimeANLS.On.fit, s='lambda.1se'))
coeft.mat <- coeft.mat[c("(Intercept)",
                         "OTHER_INTERLEUKIN_SIGNALING",
                         "ARACHIDONIC_ACID_METABOLISM",
                         "ZBP1_DAI_MEDIATED_INDUCTION_OF_TYPE_I_IFNS"),]
coeft.mat <- round(coeft.mat,5)
names(coeft.mat)[1] <- 'Intercept'
#names(coeft.mat)[c(1,3,5,12)] <- c('Intercept', 
#                                "SYNTHESIS_OF_BILE_ACIDS_AND\n_BILE_SALTS_VIA_7ALPHA_HYDROXYCHOLESTEROL",
#                                "CASPASE_ACTIVATION_VIA_EXTRINSIC\n_APOPTOTIC_SIGNALLING_PATHWAY",
#                                "ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR\n_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS")
coeft.df <- data.frame(cbind(names(coeft.mat),coeft.mat))
colnames(coeft.df) <- c("Pathway","Weight")
coeft.df

#coeft.df
tiff(file="output/TimeANLS_ON/TimeANLS_ON_Coefficient.tiff",width=12, height=6.5,unit='in',compression = 'lzw',res=150)
tab <- ggtexttable(coeft.df, rows = NULL, theme = ttheme("mOrange"))
tab %>%
  tab_add_title(text = "TimeANLS-ON model", face = "bold" ,size = 11, padding = unit(0.3, "line"))
dev.off()

####################################################################################
#### Output Result
####################################################################################
write.csv(Riaz.TA.On.ssgsea, file='output/TimeANLS_ON/ssGSEA_Score/Riaz_TimeANLS_ON_ssgsea.csv')
write.csv(Gide.TA.On.ssgsea, file='output/TimeANLS_ON/ssGSEA_Score/Gide_TimeANLS_ON_ssgsea.csv')
write.csv(Lee.TA.On.ssgsea, file='output/TimeANLS_ON/ssGSEA_Score/Lee_TimeANLS_ON_ssgsea.csv')
write.csv(MGH.TA.On.ssgsea, file='output/TimeANLS_ON/ssGSEA_Score/MGH_TimeANLS_ON_ssgsea.csv')

TA_ON_ss <- list('Riaz' = Riaz.TA.On.prob,
                   'Gide' = Gide.TA.On.prob,
                   'Lee' = Lee.TA.On.prob,
                   'MGH' = MGH.TA.On.prob
                 )
save(TA_ON_ss, file='output/TimeANLS_ON/TimeANLS_ON_SigScore.Rdata')

#Save PASS ON Model
save(TimeANLS.On.fit, file='output/TimeANLS_ON/TimeANLS_ON_Model_Fit.Rdata')








