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
#### Load PASS-ON Signature and get ssGSEA for each dataset
####################################################################################
#PASS-ON Signature
load('output/Pathway_Singatures.Rdata')
PASS.ON.Sigs <- prepare_sig(Pathway.Sigs$PASS_ON)

#PASS-ON ssGSEA
load('data/Riaz_data.RData')
Riaz.On.Samples <- split(Riaz_data$Samples,Riaz_data$Samples$PreOn)$On
Riaz.On.TPM <- Riaz_data$TPM[,match(Riaz.On.Samples$Sample, colnames(Riaz_data$TPM))]
Riaz.On.ssgsea <- gsva(expr = Riaz.On.TPM, gset.idx.list = PASS.ON.Sigs, method='ssgsea')

#load Gide Data
load('data/Gide_data.RData')
Gide.On.Samples <- split(Gide_data$Samples,Gide_data$Samples$PREEDT)$EDT
Gide.On.TPM <- Gide_data$TPM[,match(Gide.On.Samples$Sample, colnames(Gide_data$TPM))]
Gide.On.ssgsea <- gsva(expr = Gide.On.TPM, gset.idx.list = PASS.ON.Sigs, method='ssgsea')

#load Lee Data
load('data/Lee_data.RData')
Lee.On.Samples <- Lee_data$On_Samples
Lee.On.TPM <- Lee_data$TPM[,match(Lee.On.Samples$Pt_ID, colnames(Lee_data$TPM))]
Lee.On.ssgsea <- gsva(expr = as.matrix(Lee.On.TPM), gset.idx.list = PASS.ON.Sigs, method='ssgsea')

#load MGH Origin Data
load('data/MGH_ON_data.RData')
MGH.On.Samples <- MGH_ON_data$Samples
Batch14_ON_ssgsea <-  gsva(expr = as.matrix(MGH_ON_data$Batch14), gset.idx.list = PASS.ON.Sigs, method='ssgsea')
Batch17_ON_ssgsea <- gsva(expr = as.matrix(MGH_ON_data$Batch17), gset.idx.list = PASS.ON.Sigs, method='ssgsea')
SN0119610_ON_ssgsea <- gsva(expr = as.matrix(MGH_ON_data$SN0119610), gset.idx.list = PASS.ON.Sigs, method='ssgsea')
SN0123099_ON_ssgsea <- gsva(expr = as.matrix(MGH_ON_data$SN0123099), gset.idx.list = PASS.ON.Sigs, method='ssgsea')
SN0131794_ON_ssgsea <- gsva(expr = as.matrix(MGH_ON_data$SN0131794), gset.idx.list = PASS.ON.Sigs, method='ssgsea')
MGH.On.ssgsea <- data.frame(cbind(Batch14_ON_ssgsea, Batch17_ON_ssgsea, SN0119610_ON_ssgsea, SN0123099_ON_ssgsea, SN0131794_ON_ssgsea))
MGH.On.ssgsea <- MGH.On.ssgsea[,match(MGH.On.Samples$Sample, colnames(MGH.On.ssgsea))]
all(colnames(MGH.On.ssgsea) == MGH.On.Samples$Sample)

####################################################################################
#### Training, Validation & Testing Model
####################################################################################
#Train Model
set.seed(1028)
PASS.On.fit <- cv.glmnet(x=t(Riaz.On.ssgsea), y=as.factor(Riaz.On.Samples$Resp_NoResp), 
                         family="binomial",type.measure="auc",nfolds=3, alpha=0.6,
                         offset = prior_correction(Riaz.On.Samples$Resp_NoResp, tau=2/3)) #rep(0,54)

# Validate on Riaz On ssGSEA
Riaz.On.prob <- predict(PASS.On.fit, newx=t(Riaz.On.ssgsea), type = "response", s="lambda.1se",newoffset = prior_correction(Riaz.On.Samples$Resp_NoResp,2/3)) #prior_correction(Riaz.On.Samples$Resp_NoResp)
performance(prediction(predictions = Riaz.On.prob, labels = Riaz.On.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.8297258

# Test on Gide On ssGSEA
Gide.On.prob <- predict(PASS.On.fit, newx=t(Gide.On.ssgsea), type = "response", s="lambda.1se",newoffset = prior_correction(Gide.On.Samples$Resp_NoResp, tau=2/3)) #rep(0, 18)
performance(prediction(predictions = Gide.On.prob, labels = Gide.On.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.8831169

# Test on Lee On ssGSEA
Lee.On.prob <- predict(PASS.On.fit, newx=t(Lee.On.ssgsea), type = "response", s="lambda.1se",newoffset = prior_correction(Lee.On.Samples$Resp_NoResp_ResPROG, 2/3)) #rep(0,20)
performance(prediction(predictions = Lee.On.prob, labels = Lee.On.Samples$Resp_NoResp_ResPROG),"auc")@y.values[[1]] #[1] 0.8505747

# Test on MGH(Full Genes) On ssGSEA
MGH.On.prob <- predict(PASS.On.fit, newx = t(MGH.On.ssgsea), type = "response", s="lambda.1se", newoffset = prior_correction(MGH.On.Samples$Resp_NoResp, tau = 2/3)) #rep(0,32)
performance(prediction(predictions = MGH.On.prob, labels = MGH.On.Samples$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.8923077

#Plot Model
tiff(filename = 'output/PASS_ON/PASS_ON_Model_Plot.tiff',width=5, height=5,unit='in',compression = 'lzw',res=150)
plot(PASS.On.fit)
title('Elastic-Net Regularization CV AUC',line=2.5)
dev.off()

tiff(file="output/PASS_ON/PASS_ON_L3reg_plot_path.tiff",width=8, height=6.5,unit='in',compression = 'lzw',res=150)
plot_glmnet(glmnet(x=t(Riaz.On.ssgsea), y=as.factor(Riaz.On.Samples$Resp_NoResp), family="binomial", alpha=0.6),xvar='lambda',label=15)
abline(v=log(PASS.On.fit$lambda.1se), col="blue", lty=2)
text(x=-4, y=-100, label = 'log(lambda.1se)\n-2.216504\nNonzero Coeff 4',cex=0.8) #-2.216504
title('Elastic-Net Regularization solution paths',line=3)
dev.off()

coeft.mat <- as.matrix(coef(PASS.On.fit, s='lambda.1se'))
coeft.mat <- coeft.mat[c("(Intercept)","PEROXISOMAL_LIPID_METABOLISM","GENERATION_OF_SECOND_MESSENGER_MOLECULES","FATTY_ACID_METABOLISM","PD_1_SIGNALING"),]
coeft.mat <- round(coeft.mat,5)
names(coeft.mat)[1] <- 'Intercept'
coeft.df <- data.frame(cbind(names(coeft.mat),coeft.mat))
colnames(coeft.df) <- c("Pathway","Weight")
coeft.df

#coeft.df
tiff(file="output/PASS_ON/PASS_ON_Coefficient.tiff",width=6.5, height=4.5,unit='in',compression = 'lzw',res=150)
tab <- ggtexttable(coeft.df, rows = NULL, theme = ttheme("mOrange"))
tab %>%
  tab_add_title(text = "PASS-ON model", face = "bold" ,size = 11, padding = unit(0.3, "line"))
dev.off()

####################################################################################
#### Output Result
####################################################################################
write.csv(Riaz.On.ssgsea, file='output/PASS_ON/ssGSEA_Score/Riaz_PASS_ON_ssgsea.csv')
write.csv(Gide.On.ssgsea, file='output/PASS_ON/ssGSEA_Score/Gide_PASS_ON_ssgsea.csv')
write.csv(Lee.On.ssgsea, file='output/PASS_ON/ssGSEA_Score/Lee_PASS_ON_ssgsea.csv')
write.csv(MGH.On.ssgsea, file='output/PASS_ON/ssGSEA_Score/MGH_PASS_ON_ssgsea.csv')


PASS_ON_ss <- list('Riaz' = Riaz.On.prob,
                   'Gide' = Gide.On.prob,
                   'Lee' = Lee.On.prob,
                   'MGH' = MGH.On.prob
                  )

save(PASS_ON_ss, file='output/PASS_ON/PASS_ON_SigScore.Rdata')

#Save PASS ON Model
save(PASS.On.fit, file='output/PASS_ON/PASS_ON_Model_Fit.Rdata')

