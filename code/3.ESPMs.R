#----------------------------------------------------Predict model--------------------------------------------#
#This model requires the use of python in the R language
#Calling python
library(reticulate)
use_condaenv("/home/anaconda3/envs/zzw_envs/",required=T)
#neural networks APl
library(keras)
library(pROC)
library(caret)
library(ROCR)

#callback function for neural networks
earlystopping <- keras::callback_early_stopping(monitor = 'val_loss',
                                                mode = 'min' ,
                                                patience =10,
                                                verbose = 0)

checkpoint    <- keras::callback_model_checkpoint("./best_weights.hdf5",
                                                  monitor = 'val_loss',
                                                  mode='min',
                                                  save_best_only=T,
                                                  verbose = 0)
callback_list <- list(earlystopping,checkpoint)
#neural network parameters
train_ANN_model <- function(trainData, trainLabels, epochs = 100, dropout = 0.4, activation_func = "selu", validation_split = 0.1){
  ind <- !is.na(trainLabels)
  x_train <- trainData[ind, ]
  y_train <- factor(trainLabels)
  levels <- levels(y_train)
  class <- length(levels)
  y_train <- keras::to_categorical(as.numeric(y_train[ind]) - 1, class)
  k_clear_session()
  init_methods <- "glorot_uniform"
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = 1024, activation = activation_func, input_shape = ncol(x_train), kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_gaussian_dropout(rate = dropout) %>%
    layer_dense(units = 512, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_gaussian_dropout(rate = dropout) %>%
    layer_dense(units = 256, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_gaussian_dropout(rate = dropout) %>%
    layer_dense(units = 128, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_gaussian_dropout(rate = dropout) %>%
    layer_dense(units = 64, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_gaussian_dropout(rate = dropout) %>%
    layer_dense(units = 16, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_gaussian_dropout(rate = dropout) %>%
    layer_dense(units = 8, activation = activation_func, kernel_initializer = init_methods) %>%
    layer_batch_normalization() %>%
    layer_gaussian_dropout(rate = dropout) %>%
    layer_dense(units = class, activation = 'softmax')
  model %>% compile(
    loss = "binary_crossentropy",
    optimizer_adam(lr=0.001, beta_1=0.9, beta_2=0.999, decay = 1e-06),
    metrics = c('AUC')
  )
  
  history <- model %>% keras::fit(
    x_train, y_train,
    epochs = epochs, batch_size = 256,
    view_metrics = F,
    callbacks=callback_list,
    validation_split = validation_split,
    verbose = 0
  )
  model %>% compile(
    loss = "binary_crossentropy",
    optimizer = optimizer_sgd(lr = 1e-04, momentum = 0.9, decay = 1e-07),
    metrics = c('AUC')
  )
  
  history <- model %>% keras::fit(
    x_train, y_train,
    epochs = epochs, batch_size = 256,
    view_metrics = F,
    validation_split = validation_split,
    callbacks=callback_list,
    verbose = 0
  )
  list(classifier = model, levels = levels)
}
get_predict_prob <- function(ANNmodel, newData){
  
  pre <- predict(ANNmodel$classifier, newData)
  colnames(pre) <- ANNmodel$levels
  pre
}
'---------------------------------EVO-------------------------------------'

#Training evolutionary models
cancer=c("HNSC","LIHC","LUAD","LUSC")
for (m in 1:length(cancer)){ 
evo<-paste("./Input_features/",cancer[m],"/",cancer[m],"_exp_best_feature_mean.csv",sep="")
nonevo<-paste("./Input_features/",cancer[m],"/",cancer[m],"_exp_non*_best_feature_mean.csv",sep="")
tcga<- read.table(evo,sep = ',',row.names  = 1,header = T,check.names = F)
set.seed(123)
folds <-createMultiFolds(y=tcga$class,k=5,times = 20)
gt<-as.data.frame(array(,dim=c(100,6)))
colnames(gt)<-c("AUC","accuracy","precision","recall",'specificity',"F1")
for (i in 1:100) {
  tcga_train<-tcga[folds[[i]],]
  tcga_test<-tcga[-folds[[i]],]
  #labels<-tcga_train[,"classby3"]
  labels<-tcga_train[,"class"]
  tcga_train<-tcga_train[,-ncol(tcga_train)]
  fs<-as.matrix(tcga_train)
  #new_labels<-tcga_test[,"classby3"]
  new_labels<-tcga_test[,"class"]
  tcga_test<-tcga_test[,-ncol(tcga_test)]
  new_fs<-as.matrix(tcga_test)
  tensorflow::set_random_seed(3690)
  ANN_model <- train_ANN_model(fs, labels)
  pred<-get_predict_prob(ANN_model,new_fs)
  #Model performance evaluation metrics
  gt[i,1]<-as.numeric(auc(new_labels,pred[,2]))
  pred.class <- as.integer(pred[,2] > 0.5)
  cft <-confusionMatrix(as.factor(pred.class), as.factor(new_labels), positive = "1")
  tp <- cft$table[2, 2]
  tn <- cft$table[1, 1]
  fp <- cft$table[2, 1]
  fn <- cft$table[1, 2]
  ACC <-(tp+tn)/(tp+tn+fp+fn)
  precision<- tp/(tp+fp)
  recall<- tp/(tp+fn)
  specificity <- tn/(tn + fp)
  F1=2*(precision*recall)/(precision+recall)
  gt[i,2] <- ACC
  gt[i,3] <- precision
  gt[i,4] <- recall
  gt[i,5] <- specificity
  gt[i,6] <- F1
  
}
#write.table(gt,paste(cancer[m],"/",cancer[m],"_5_results.txt",sep=""),row.names = F,quote=F,sep="\t")
################-------------------------------NON-EVO---------------------------------##################
#Training non-evolutionary models
nonevos <- Sys.glob(nonevo)
result<-as.data.frame(matrix("",5,2))
gt<-as.data.frame(array(,dim=c(500,6)))
colnames(gt)<-c("AUC","accuracy","precision","recall",'specificity',"F1")
set.seed(123)
for (j in seq_along(nonevos)){
  tcga<-read.table(nonevos[j], sep = ',', header = TRUE,check.names = F,row.names = 1)
  result[j,1]<-j
  print(j)
    for (i in 1:100) {
      tcga_train<-tcga[folds[[i]],]
      tcga_test<-tcga[-folds[[i]],]
      #labels<-tcga_train[,"classby3"]
      labels<-tcga_train[,"class"]
      tcga_train<-tcga_train[,-ncol(tcga_train)]
      fs<-as.matrix(tcga_train)
      #new_labels<-tcga_test[,"classby3"]
      new_labels<-tcga_test[,"class"]
      tcga_test<-tcga_test[,-ncol(tcga_test)]
      new_fs<-as.matrix(tcga_test)
      tensorflow::set_random_seed(3690)
      ANN_model <- train_ANN_model(fs, labels)
      pred<-get_predict_prob(ANN_model,new_fs)
      gt[100*(j-1)+i,1]<-as.numeric(auc(new_labels,pred[,2]))
      pred.class <- as.integer(pred[,2] > 0.5)
      cft <-confusionMatrix(as.factor(pred.class), as.factor(new_labels), positive = "1")
      tp <- cft$table[2, 2]
      tn <- cft$table[1, 1]
      fp <- cft$table[2, 1]
      fn <- cft$table[1, 2]
      ACC <-(tp+tn)/(tp+tn+fp+fn)
      precision<- tp/(tp+fp)
      recall<- tp/(tp+fn)
      specificity <- tn/(tn + fp)
      F1=2*(precision*recall)/(precision+recall)
      gt[100*(j-1)+i,2] <- ACC
      gt[100*(j-1)+i,3] <- precision
      gt[100*(j-1)+i,4] <- recall
      gt[100*(j-1)+i,5] <- specificity
      gt[100*(j-1)+i,6] <- F1
    }
    result[j,2]<-mean(gt$AUC[(1+100*(j-1)):(100*j)])}
#write.table(gt,paste(cancer[m],"/",cancer[m],"_non_5_results.txt",sep=""),sep="\t",row.names = F,quote=F)
}
