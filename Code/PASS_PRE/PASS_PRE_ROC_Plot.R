library(here)
library(pROC)
library(ROCR) #probably disable later on
library(plotmo)
library(RColorBrewer)
colors4 <- c(brewer.pal(4, "Accent")[1:3],'skyblue3')
colors4
########################################
##### Load Sample Information
########################################
#Riaz
load('data/Riaz_data.RData')
#Gide
load('data/Gide_data.RData')
#Lee
load('data/Lee_data.RData')
#MGH
load('data/MGH_PRE_data.RData')

########################################
##### Load Sample Signature Score
########################################
load('output/PASS_PRE/PASS_PRE_SigScore.Rdata')

Riaz.res <- Riaz_data$Samples[match(rownames(PASS_PRE_ss$Riaz), Riaz_data$Samples$Sample),]
Gide.res <- Gide_data$Samples[match(rownames(PASS_PRE_ss$Gide), Gide_data$Samples$Sample),]
Lee.res <- Lee_data$Pre_Samples[match(rownames(PASS_PRE_ss$Lee), Lee_data$Pre_Samples$Pt_ID),]
MGH.res <- MGH_PRE_data$Samples[match(rownames(PASS_PRE_ss$MGH), MGH_PRE_data$Samples$Sample),]

########################################
##### ROC plot
########################################
#Riaz AUC & RIC
Riaz.pred <- prediction(PASS_PRE_ss$Riaz, Riaz.res$Resp_NoResp)
roc(Riaz.res$Resp_NoResp, PASS_PRE_ss$Riaz)$auc #0.7294
Riaz.ROC <- performance(Riaz.pred, measure = "tpr", x.measure = "fpr")

#Gide AUC & ROC
Gide.pred <- prediction(PASS_PRE_ss$Gide, Gide.res$Resp_NoResp)
roc(Gide.res$Resp_NoResp, PASS_PRE_ss$Gide)$auc #0.6815
Gide.ROC <- performance(Gide.pred, measure = "tpr", x.measure = "fpr")

#Lee AUC & ROC
Lee.pred <- prediction(PASS_PRE_ss$Lee, Lee.res$Resp_NoResp)
roc(Lee.res$Resp_NoResp, PASS_PRE_ss$Lee)$auc #0.5496
performance(prediction(predictions = PASS_PRE_ss$Lee, labels = Lee.res$Resp_NoResp),"auc")@y.values[[1]] #[1] 0.4504132
Lee.ROC <- performance(Lee.pred, measure = "tpr", x.measure = "fpr")

#MGH AUC & ROC
MGH.pred <- prediction(PASS_PRE_ss$MGH, MGH.res$Resp_NoResp)
roc(MGH.res$Resp_NoResp, PASS_PRE_ss$MGH)$auc #0.6923
MGH.ROC <- performance(MGH.pred, measure = "tpr", x.measure = "fpr")

####################################################################
######## ROC Plot
####################################################################
tiff(file="output/PASS_PRE/Riaz_PASS_PRE_ROC.tiff",width=5, height=5,unit='in',compression = 'lzw',res=150) #width=10
plot(Riaz.ROC, avg="vertical", lwd=3, col=colors4[1], spread.estimate="stderror",plotCI.lwd=3,main="ROC AUC",ylab='True Positive Rate')
abline(a=0, b= 1, lty=2, col='black')
legend(x="bottomright",'Riaz AUC 0.73',col=colors4[1],lwd=3,pt.cex=1,cex=0.65)
dev.off()


tiff(file="output/PASS_PRE/GiLeMg_PASS_PRE_ROC.tiff",width=5, height=5,unit='in',compression = 'lzw',res=150) #width=10
plot(Gide.ROC, avg="vertical", lwd=3, col=colors4[2], spread.estimate="stderror",plotCI.lwd=3,main="ROC AUC",ylab='True Positive Rate')
plot(Lee.ROC, avg="vertical", lwd=3, col=colors4[3], spread.estimate="stderror",plotCI.lwd=3,add=TRUE)
plot(MGH.ROC, avg="vertical", lwd=3, col=colors4[4], spread.estimate="stderror",plotCI.lwd=3,add=TRUE)
abline(a=0, b= 1, lty=2, col='black')
legend(x="bottomright",c('Gide AUC 0.69','Lee AUC 0.45','MGH AUC 0.69'),col=colors4[c(2,3,4)],lwd=3,pt.cex=1,cex=0.65)
dev.off()

####################################################################
#All test data combine together
all_test_sample_prob <- c(PASS_PRE_ss$Gide, PASS_PRE_ss$Lee, PASS_PRE_ss$MGH)
all_test_label <- c(Gide.res$Resp_NoResp, Lee.res$Resp_NoResp, MGH.res$Resp_NoResp)
table(all_test_label)
#No_Response    Response 
#         62          73 

roc(all_test_label, all_test_sample_prob)$auc #Area under the curve: 0.6469
var(roc(all_test_label, all_test_sample_prob)) #[1] 0.002417506
all.test.sample.prediction <- prediction(predictions = all_test_sample_prob, labels = all_test_label)
all.test.sample.ROC <- performance(all.test.sample.prediction, measure = "tpr", x.measure = "fpr")
all.test.sample.prec.cutff <- performance(all.test.sample.prediction, measure = "prec", x.measure = "cutoff") 
all.test.sample.recall.cutff <- performance(all.test.sample.prediction, measure = "rec", x.measure = "cutoff")

tiff(file="output/PASS_PRE/All_Test_Sample_PASS_PRE_ROC.tiff",width=5, height=5,unit='in',compression = 'lzw',res=150) #width=10
plot(all.test.sample.ROC, avg="vertical", lwd=3, col=colors4[1], spread.estimate="stderror",plotCI.lwd=3,main="ROC AUC",ylab='True Positive Rate')
abline(a=0, b= 1, lty=2, col='black')
legend(x="bottomright",'All Pre-Trt. Test Samples AUC 0.65',col=colors4[1],lwd=3,pt.cex=1,cex=0.65)
dev.off()

####################################################################
#Accuracy
cal_metrics <- function(label, pred){
  roc.p=pROC::roc(label, pred)
  if (roc.p$auc>0.5){
    cutoff=roc.p$thresholds[which.max(roc.p$sensitivities+roc.p$specificities)]
    sensitivity=roc.p$sensitivities[which.max(roc.p$sensitivities+roc.p$specificities)]
    specificity=roc.p$specificities[which.max(roc.p$sensitivities+roc.p$specificities)]
    df=data.frame(type='positive classification',
                  auc=round(roc.p$auc,3),cutoff=cutoff,
                  sensitivity=sensitivity,specificity=specificity)
    return(df)
  }
  else{
    cutoff=roc.p$thresholds[which.min(roc.p$sensitivities+roc.p$specificities)]
    sensitivity=roc.p$sensitivities[which.min(roc.p$sensitivities+roc.p$specificities)]
    specificity=roc.p$specificities[which.min(roc.p$sensitivities+roc.p$specificities)]
    df=data.frame(type='negative classification',
                  auc=1-round(roc.p$auc,3),cutoff=cutoff,
                  sensitivity=1-sensitivity,specificity=1-specificity)
    return(df)
  }
}

Riaz.cutoff.df <- cal_metrics(label = Riaz.res$Resp_NoResp, pred = PASS_PRE_ss$Riaz)
Riaz.cutoff.df
#                     type  auc    cutoff   sensitivity specificity
#1 positive classification 0.729 0.3834639   0.5555556   0.8709677

all.test.sample.result.df <- data.frame(Prob = c(PASS_PRE_ss$Gide, PASS_PRE_ss$Lee, PASS_PRE_ss$MGH), 
                                        Label =  c(Gide.res$Resp_NoResp, Lee.res$Resp_NoResp, MGH.res$Resp_NoResp)
                                        )

colnames(all.test.sample.result.df) <- c('Prob', 'Label')
all.test.sample.result.df$Best_Cutoff <- ifelse(all.test.sample.result.df$Prob >= Riaz.cutoff.df$cutoff, "Response", "No_Response")
table(all.test.sample.result.df$Label,all.test.sample.result.df$Best_Cutoff)
#            No_Response Response
#No_Response          11       51
#Response              4       69
(11+69)/(11+51+4+69) #[1] 0.5925926


#Output Sigscore
all(rownames(PASS_PRE_ss$Riaz) == Riaz.res$Sample)
Riaz.res$Sig_Score <- PASS_PRE_ss$Riaz
write.csv(Riaz.res,file='output/PASS_PRE/SigScore/Riaz_PASS_PRE_SigScore.csv',row.names=F)

all(rownames(PASS_PRE_ss$Gide) == Gide.res$Sample)
Gide.res$Sig_Score <- PASS_PRE_ss$Gide
write.csv(Gide.res,file='output/PASS_PRE/SigScore/Gide_PASS_PRE_SigScore.csv',row.names=F)

all(rownames(PASS_PRE_ss$Lee) == Lee.res$Pt_ID)
Lee.res$Sig_Score <- PASS_PRE_ss$Lee
write.csv(Lee.res,file='output/PASS_PRE/SigScore/Lee_PASS_PRE_SigScore.csv',row.names=F)

all(rownames(PASS_PRE_ss$MGH) == MGH.res$Sample)
MGH.res$Sig_Score <- PASS_PRE_ss$MGH
write.csv(MGH.res,file='output/PASS_PRE/SigScore/MGH_PASS_PRE_SigScore.csv',row.names=F)
