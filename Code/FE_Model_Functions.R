library(data.table)

####################################################################
#Prepare signaure to ssGSEA
prepare_sig <- function(sig){
  leadingEdge.list <- sig$leadingEdge
  names(leadingEdge.list) <- sig$pathway
  return(leadingEdge.list)
}

####################################################################
#Cost-Sensitive learning
prior_correction <- function(data_label , tau=0.5){
  tau_ratio = (1- tau)/(tau)
  label_porp <- data.frame(table(data_label))
  label_ratio <- (label_porp$Freq[2]/sum(label_porp$Freq))/(1-label_porp$Freq[2]/sum(label_porp$Freq))
  correction.value <- log(tau_ratio * label_ratio)
  return(rep(correction.value, length(data_label)))
}

####################################################################
#Model Signature
#sig_score <- function(prob){
  #prob.scale <-  scale(prob,center=T,scale=T)
  #prob.scale2 <- ifelse(prob.scale < -2, -2, ifelse(prob.scale >2,2,prob.scale))
  #prob.scale2 <- ifelse(prob.scale < -2, -2, ifelse(prob.scale >2,2,prob.scale))
  #return(prob.scale2)
#}


