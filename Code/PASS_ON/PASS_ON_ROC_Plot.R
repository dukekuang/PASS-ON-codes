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
load('data/MGH_ON_data.RData')

########################################
##### Load Sample Signature Score
########################################
load('output/PASS_ON/PASS_ON_SigScore.Rdata')

Riaz.res <- Riaz_data$Samples[match(rownames(PASS_ON_ss$Riaz), Riaz_data$Samples$Sample),]
Gide.res <- Gide_data$Samples[match(rownames(PASS_ON_ss$Gide), Gide_data$Samples$Sample),]
Lee.res <- Lee_data$On_Samples[match(rownames(PASS_ON_ss$Lee), Lee_data$On_Samples$Pt_ID),]
MGH.res <- MGH_ON_data$Samples[match(rownames(PASS_ON_ss$MGH), MGH_ON_data$Samples$Sample),]

########################################
##### ROC plot
########################################
#Riaz AUC & RIC
Riaz.pred <- prediction(PASS_ON_ss$Riaz, Riaz.res$Resp_NoResp)
roc(Riaz.res$Resp_NoResp, PASS_ON_ss$Riaz)$auc #Area under the curve: 0.8297
var(roc(Riaz.res$Resp_NoResp, PASS_ON_ss$Riaz)) #[1] 0.004073928
Riaz.ROC <- performance(Riaz.pred, measure = "tpr", x.measure = "fpr")

#Gide AUC & ROC
Gide.pred <- prediction(PASS_ON_ss$Gide, Gide.res$Resp_NoResp)
roc(Gide.res$Resp_NoResp, PASS_ON_ss$Gide)$auc #Area under the curve: 0.8831
var(roc(Gide.res$Resp_NoResp, PASS_ON_ss$Gide)) #[1] 0.007724743
Gide.ROC <- performance(Gide.pred, measure = "tpr", x.measure = "fpr")

#Lee AUC & ROC
Lee.pred <- prediction(PASS_ON_ss$Lee, Lee.res$Resp_NoResp_ResPROG)
roc(Lee.res$Resp_NoResp_ResPROG, PASS_ON_ss$Lee)$auc #Area under the curve: 0.8506
var(roc(Lee.res$Resp_NoResp_ResPROG, PASS_ON_ss$Lee)) #[1] 0.008420625
Lee.ROC <- performance(Lee.pred, measure = "tpr", x.measure = "fpr")

#MGH AUC & ROC
MGH.pred <- prediction(PASS_ON_ss$MGH, MGH.res$Resp_NoResp)
roc(MGH.res$Resp_NoResp, PASS_ON_ss$MGH)$auc #Area under the curve: 0.8923
var(roc(MGH.res$Resp_NoResp, PASS_ON_ss$MGH)) #[1] 0.005608284
MGH.ROC <- performance(MGH.pred, measure = "tpr", x.measure = "fpr")

####################################################################
########3 ROC Plot
####################################################################
tiff(file="output/PASS_ON/Riaz_PASS_ON_ROC.tiff",width=5, height=5,unit='in',compression = 'lzw',res=150) #width=10
plot(Riaz.ROC, avg="vertical", lwd=3, col=colors4[1], spread.estimate="stderror",plotCI.lwd=3,main="ROC AUC",ylab='True Positive Rate')
abline(a=0, b= 1, lty=2, col='black')
legend(x="bottomright",'Riaz AUC 0.83',col=colors4[1],lwd=3,pt.cex=1,cex=0.65)
dev.off()


tiff(file="output/PASS_ON/GiLeMg_PASS_ON_ROC.tiff",width=5, height=5,unit='in',compression = 'lzw',res=150) #width=10
plot(Gide.ROC, avg="vertical", lwd=3, col=colors4[2], spread.estimate="stderror",plotCI.lwd=3,main="ROC AUC",ylab='True Positive Rate')
plot(Lee.ROC, avg="vertical", lwd=3, col=colors4[3], spread.estimate="stderror",plotCI.lwd=3,add=TRUE)
plot(MGH.ROC, avg="vertical", lwd=3, col=colors4[4], spread.estimate="stderror",plotCI.lwd=3,add=TRUE)
abline(a=0, b= 1, lty=2, col='black')
legend(x="bottomright",c('Gide AUC 0.88','Lee AUC 0.85','MGH AUC 0.89'),col=colors4[c(2,3,4)],lwd=3,pt.cex=1,cex=0.65)
dev.off()

####################################################################
#All data combine together
all_sample_prob <- c(PASS_ON_ss$Riaz, PASS_ON_ss$Gide, PASS_ON_ss$Lee, PASS_ON_ss$MGH)
all_label <- c(Riaz.res$Resp_NoResp, Gide.res$Resp_NoResp, Lee.res$Resp_NoResp_ResPROG, MGH.res$Resp_NoResp)
table(all_label)
#No_Response    Response 
#95          43 

roc(all_label, all_sample_prob)$auc #Area under the curve: 0.8135

####################################################################
#All test data combine together
all.test.sample.prob <- c(PASS_ON_ss$Gide, PASS_ON_ss$Lee, PASS_ON_ss$MGH)
all.test.label <- c(Gide.res$Resp_NoResp, Lee.res$Resp_NoResp_ResPROG, MGH.res$Resp_NoResp)
table(all.test.label)
#No_Response    Response 
#         62          22 

roc(all.test.label, all.test.sample.prob)$auc #Area under the curve: 0.8798
var(roc(all.test.label, all.test.sample.prob)) #[1] 0.001612903
all.test.sample.prediction <- prediction(predictions = all.test.sample.prob, labels = all.test.label)
all.test.sample.ROC <- performance(all.test.sample.prediction, measure = "tpr", x.measure = "fpr")
all.test.sample.prec.cutff <- performance(all.test.sample.prediction, measure = "prec", x.measure = "cutoff") 
all.test.sample.recall.cutff <- performance(all.test.sample.prediction, measure = "rec", x.measure = "cutoff")

tiff(file="output/PASS_ON/All_test_Sample_PASS_ON_ROC.tiff",width=5, height=5,unit='in',compression = 'lzw',res=150) #width=10
plot(all.test.sample.ROC, avg="vertical", lwd=3, col=colors4[1], spread.estimate="stderror",plotCI.lwd=3,main="ROC AUC",ylab='True Positive Rate')
abline(a=0, b= 1, lty=2, col='black')
legend(x="bottomright",'All ON-Trt. Test Samples AUC 0.88',col=colors4[1],lwd=3,pt.cex=1,cex=0.65)
dev.off()

####################################################################
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

Riaz.cutoff.df <- cal_metrics(label = Riaz.res$Resp_NoResp, pred = PASS_ON_ss$Riaz)
Riaz.cutoff.df
#                   type  auc    cutoff sensitivity specificity
#positive classification 0.83 0.4247767   0.7619048   0.8787879

#All Test Sample ACC
all.test.sample.result.df <- data.frame(Prob = all.test.sample.prob, 
                                        Label =  all.test.label
                                        )
colnames(all.test.sample.result.df) <- c('Prob', 'Label')
all.test.sample.result.df$Best_Cutoff <- ifelse(all.test.sample.result.df$Prob >= Riaz.cutoff.df$cutoff, "Response", "No_Response")
table(all.test.sample.result.df$Label,all.test.sample.result.df$Best_Cutoff)
#            No_Response Response
#No_Response          55        7
#Response              8       14
print('Accuracy')
(55+14)/(55+7+8+14) #[1] 0.8214286

#Output Sigscore
all(rownames(PASS_ON_ss$Riaz) == Riaz.res$Sample)
Riaz.res$Sig_Score <- PASS_ON_ss$Riaz
write.csv(Riaz.res,file='output/PASS_ON/SigScore/Riaz_PASS_ON_SigScore.csv',row.names=F)

all(rownames(PASS_ON_ss$Gide) == Gide.res$Sample)
Gide.res$Sig_Score <- PASS_ON_ss$Gide
write.csv(Gide.res,file='output/PASS_ON/SigScore/Gide_PASS_ON_SigScore.csv',row.names=F)

all(rownames(PASS_ON_ss$Lee) == Lee.res$Pt_ID)
Lee.res$Sig_Score <- PASS_ON_ss$Lee
write.csv(Lee.res,file='output/PASS_ON/SigScore/Lee_PASS_ON_SigScore.csv',row.names=F)

all(rownames(PASS_ON_ss$MGH) == MGH.res$Sample)
MGH.res$Sig_Score <- PASS_ON_ss$MGH
write.csv(MGH.res,file='output/PASS_ON/SigScore/MGH_PASS__SigScore.csv',row.names=F)


