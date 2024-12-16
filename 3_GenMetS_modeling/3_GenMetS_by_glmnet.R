###########################
#GenMetS by glmnet
##########################
rm(list=ls())
library(dplyr)
library(data.table)
library(caret)
library(glmnet)
library(glmnetUtils)

N=31409
myfile="s934_g434_data_for_training_testing_SNP31409.rds"

myraw<-readRDS(file = myfile)
dim(myraw)
cat("columns in myraw\n")
print(unique(substr(colnames(myraw), 1, 2)))

idx<-which(myraw$cohort %in% c("gusto", "spresto"))
table(myraw$MetS_score, myraw$cohort) 

SNPidx<- which(grepl("^rs", colnames(myraw)))
metsidx<-which(colnames(myraw)=="MetS_score")
Featureidx <- c(SNPidx) 


##########################################################
#data:
#train=(spresto + gusto)*80%
#test =(spresto + gusto)*20%
###########################################################
set.seed(1234)

spresto_data <- myraw %>% filter(cohort=="spresto")
gusto_data   <- myraw %>% filter(cohort == "gusto")
tmptrain     <- rbind(spresto_data, gusto_data)
x_train      <- tmptrain %>% dplyr::select(all_of(Featureidx)) %>% data.frame() %>%as.matrix() ; 
dim(x_train)
y_train      <- tmptrain$MetS_score; length(y_train)
print(dim(x_train)) 

inTraining <- createDataPartition(y_train, p = .80, list = FALSE)
x_train80 <- x_train[inTraining, ]
y_train80 <- y_train[inTraining]

x_test20 <- x_train[-inTraining, ]
y_test20 <- y_train[-inTraining]

write.table(tmptrain[-inTraining, c(1,2)],  "result/test_data_subjectID.txt", row.names=F, sep="\t", quote=F)

n_train <- dim(x_train80)[1]
K=10
flds <- createFolds(y_train80 , k = K, list = TRUE, returnTrain = FALSE)
foldids = rep(1,length(y_train80 ))
for(I in 2:K){ foldids[unlist(flds[I])] = I }


###############################
##find the best alpha
###############################
outfile="result/GenMetS_cva_training_model.RData"
if(!file.exists(outfile)){
set.seed(1234)
cva_result<- cva.glmnet(x_train80 , y_train80 , 
                  foldid = foldids, 
                  alpha = seq(0.1, 1, len = 21)^3,
                  family = "gaussian", 
                  sparse=TRUE,
                  keep = TRUE, 
                  intercept = FALSE) # standarize=TRUE/FALSE does not influence the result.
save(cva_result, file=outfile)
}
load(outfile)

###############report cva_result
  n<- length(cva_result$alpha)
  myresult <- matrix(rep(NA, 7*n), ncol=7)
  
  for(idx in c(1:n)){
    cvfit  <- cva_result$modlist[[idx]]
    myalpha<- cva_result$alpha[idx]
    mylambda<- cvfit$lambda.min 
    tmp    <- which(cvfit$lambda==cvfit$lambda.min) #the best fit
    mynzero<- cvfit$nzero[tmp]
    mycvm  <- cvfit$cvm[tmp]
    mycvsd <- cvfit$cvsd[tmp]
    myr2   <- cvfit$glmnet$dev.ratio[tmp]
     
    myresult[idx, ] <- c(idx, myalpha, mylambda,mynzero,mycvm, mycvsd, myr2)
  }
  colnames(myresult) <- c("idx", "alpha", "lambda", "nzero", "cvm", "cvsd", "fit_r2")

  
  #######################
myresult<- myresult %>% data.frame()
print(myresult)
myresultfile =  "result/GenMetS_training_result.txt"
write.table(myresult, myresultfile, row.names=F, sep="\t", quote=F)


modelidx<- which( myresult$cvm==min(myresult$cvm[1:n]) )
cat("idx=", modelidx, "\n")
print(myresult[modelidx, ])
cvfit<- cva_result$modlist[[modelidx]]

y_test20_pred <- predict(cvfit, s = cvfit$lambda.min, newx = x_test20)
print(R2(y_test20_pred, y_test20))

tmptrain     <- rbind(spresto_data, gusto_data)
tmptest <- tmptrain[-inTraining, ]
tmptest$GenMetS<- y_test20_pred
write.table(tmptest, "result/test20_result.txt", row.names=F, sep="\t", quote=F)


mycoef <- data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
mycoef$id<- rownames(mycoef)
selected_featureidx<- which(abs(mycoef$s1) >0.00000)
cat("length of coef: ", length(selected_featureidx), "\n")
df <- mycoef[selected_featureidx, ]
featurefile=paste0("result/GenMetS_feature", length(selected_featureidx), ".txt")
cat("featurefile=",featurefile, "\n")
write.table(df, featurefile, sep="\t", row.names=F, quote=F)

