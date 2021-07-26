##############################################################################################################
model.fit.dev <- function(data,labels, model.offset, alpha, nfolds=3, seeds=1028){
  set.seed(seeds)
  #Train
  model.fit <- cv.glmnet(x=t(data$riaz), y=as.factor(labels$riaz), 
                         family="binomial",type.measure="auc",nfolds=nfolds, alpha=alpha,
                         offset = model.offset$train)
  #Validate Riaz
  Riaz.On.prob <- predict(model.fit, newx=t(data$riaz), type = "response", s="lambda.1se", newoffset = model.offset$riaz_test)
  #Test Gide
  Gide.On.prob <- predict(model.fit, newx=t(data$gide), type = "response", s="lambda.1se", newoffset = model.offset$gide_test)
  #Test Lee
  Lee.On.prob <- predict(model.fit, newx=t(data$lee), type = "response", s="lambda.1se", newoffset = model.offset$lee_test)
  #Test MGH
  MGH.On.prob <- predict(model.fit, newx = t(data$MGH), type = "response", s="lambda.1se", newoffset = model.offset$MGH_test)
  #Test MGH origin
  #MGH.orgin.On.prob <- predict(model.fit, newx = t(data$MGH_origin), type = "response", s="lambda.1se", newoffset = model.offset$MGH_origin_test)
  
  #Output
  #prob.list <- list('Riaz' = Riaz.On.prob, 'Gide' = Gide.On.prob, 'Lee' = Lee.On.prob, 'MGH' = MGH.On.prob, 'MGH_origin' = MGH.orgin.On.prob)
  prob.list <- list('Riaz' = Riaz.On.prob, 'Gide' = Gide.On.prob, 'Lee' = Lee.On.prob, 'MGH' = MGH.On.prob)
  
  #model.prob <- list('model' = model.fit, 'prob.list' = prob.list)
  #return(model.prob)
  return(prob.list)
}

##############################################################################################################
Model_Performace <- function(prob, labels, MGH_origin=F){
  Riaz <- performance(prediction(predictions = prob$Riaz, labels = labels$riaz),"auc")@y.values[[1]]
  Gide <- performance(prediction(predictions = prob$Gide, labels = labels$gide),"auc")@y.values[[1]]
  Lee <- performance(prediction(predictions = prob$Lee, labels = labels$lee),"auc")@y.values[[1]]
  if(MGH_origin == F){
    MGH <- performance(prediction(predictions = prob$MGH, labels = labels$MGH),"auc")@y.values[[1]]
  }else{
    MGH <- performance(prediction(predictions = prob$MGH_origin, labels = labels$MGH_Origin),"auc")@y.values[[1]]
  }
  auc.df <- data.frame('AUC' = c(Riaz,Gide,Lee,MGH), row.names=c('Riaz','Gide','Lee','MGH'))
  return(auc.df)
}
##############################################################################################################


