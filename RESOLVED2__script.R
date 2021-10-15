setwd("...")

#install.packages('ggplot2')
library(ggplot2)
#install.packages(c("ggalt","ggfortify","ggpubr","ggthemes","ztable"))
library(ggalt)
#install.packages("stringi")
library(stringi)
library(ggfortify)
library(dplyr)
library(tidyr)
library(scales)
#install.packages("ggpubr")
library(ggpubr)
library('scales')     # for scale_y_continuous(label = percent)
#install.packages("ggthemes")
library('ggthemes')   # for scale_fill_few('medium')
library(survival)
library(cowplot)
#install.packages("survminer")
library(survminer)
library(prodlim)
#install.packages("pec")
library(pec)
#install.packages("ModelMetrics")
#install.packages("RcppRoll")
#install.packages("ddalpha")
#install.packages("kernlab")
#install.packages("ipred")
#install.packages("tidyselect")
library(caret)
library(MASS)
#install.packages('gtools')
#install.packages('caTools')
#install.packages('pROC')
library(pROC)
#install.packages('glmnet')
library(glmnet)

### load and prepare data - merge drugbank pharmacological data and pubmed. 

FDA_clean <- read.csv2("F:/actif E&R IGR DITEP/RESOLVED²/SCRIPTS_WD/data_final/FDA2018_by_drugbank2606.latest_clean.txt", 
                       sep = "\t", header = TRUE, dec = ",")
DB_cat <- read.csv2("F:/actif E&R IGR DITEP/RESOLVED²/SCRIPTS_WD/data_final/DRUGBANK/drugbank_mining_category.latest.txt", 
                    sep = "\t", header = TRUE)
DB_target <- read.csv2("F:/actif E&R IGR DITEP/RESOLVED²/SCRIPTS_WD/data_final/DRUGBANK/drugbank_mining_targets.latest.txt", 
                       sep = "\t", header = TRUE)

FDA_clean_ok <- FDA_clean %>% filter(DRUG_TO_KEEP_FOR_ANALYSIS==1)
FDA_clean_ok <- left_join(FDA_clean_ok, DB_cat, by=c("COMMON_DRUGBANK_ALIAS" = "X"))
FDA_clean_ok <- left_join(FDA_clean_ok, DB_target, by=c("COMMON_DRUGBANK_ALIAS" = "X"))

# select columns of interest
colnames(FDA_clean_ok[,1:30])
FDA_clean_ok_drugbank <- FDA_clean_ok[,c(1,17,18,27,31,33,37,39,42,43,47,49,58:ncol(FDA_clean_ok))]
colnames(FDA_clean_ok_drugbank[,1:30])
# keep only rows without NA
FDA_clean_ok_drugbank_noNA <- FDA_clean_ok_drugbank %>% filter(Clinical_activity_detected_as.following_recurrent_responder_OR_clinical_activity_OR_prolonged_stability != "NA")
FDA_clean_ok_drugbank_noNA <- FDA_clean_ok_drugbank_noNA %>% filter(DLT_identified_or_MTD_reached != "NA")
# check
sum(is.na(FDA_clean_ok_drugbank_noNA))
FDA_clean_ok_drugbank <-FDA_clean_ok_drugbank_noNA
nrow(FDA_clean_ok_drugbank)
dim(FDA_clean_ok_drugbank)

# remove unrelevant categories
drops <- c("category.Antineoplastic.Agents", 
           "category.Myelosuppressive.Agents", 
           "category.Antineoplastic.and.Immunomodulating.Agents", 
           "category.Pharmacologic.Actions",             
           "category.Therapeutic.Uses", 
           "category.Chemical.Actions.and.Uses",
           "category.Immunosuppressive.Agents")
FDA_clean_ok_drugbank <- FDA_clean_ok_drugbank[ , !(names(FDA_clean_ok_drugbank) %in% drops)]
names(FDA_clean_ok_drugbank[,1:15])
dim(FDA_clean_ok_drugbank)

# COMPUTE FOLLOW-UP
library(prodlim)
quantile(prodlim(Hist(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN,FDA_APPROVAL)~1
                 ,data=FDA_clean_ok_drugbank,reverse=TRUE))
range(FDA_clean_ok_drugbank$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN)

summary(FDA_clean_ok_drugbank$FDA_APPROVAL)

sum(FDA_clean_ok_drugbank$FDA_APPROVAL==1)/nrow(FDA_clean_ok_drugbank)

View(FDA_clean_ok_drugbank[,1:5])

summary(survfit(Surv(FDA_clean_ok_drugbank$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
             FDA_clean_ok_drugbank$FDA_APPROVAL) ~ 1), times=c(6))


# Split the data in training and test sets 70:30
## Randomly sample cases to create independent training and test data, which have been used for initial model training and test
partition = createDataPartition(FDA_clean_ok_drugbank[,'FDA_APPROVAL'], 
                                times = 1, p = 0.7, list = FALSE)
training = FDA_clean_ok_drugbank[partition,] # Create the training sample
dim(training)
# save original training set that will be used for model fitting
write.table(training, "training_set.txt", sep = "\t")

test = FDA_clean_ok_drugbank[-partition,] # Create the test sample
dim(test)
# save original test set that will be used for model evaluation
write.table(test, "test_set.txt", sep = "\t")

# for further verifications: load directly the training and test sets used for model fit
test <- read.csv2("RESOLVED2 - 62 - SUPP DATA 2 test_set.txt" , sep="\t", header=TRUE, dec = ".")
training <- read.csv2("RESOLVED2 - 61 - SUPP DATA 1 training_set.txt", sep="\t", header=TRUE, dec = ".")
test <- test[,2:1416]
training <- training[,2:1416]

dim(test)
dim(training)

View(training[1:10,])
sum(is.na(test))
sum(is.na(training))
FDA_clean_ok$OLDEST_DATE_OF_PUBLICATION_manually_cured

# check date distribution in test and training

FDA_clean_ok_date_of_pub <- FDA_clean_ok[,c("COMMON_DRUGBANK_ALIAS","OLDEST_DATE_OF_PUBLICATION_manually_cured")]

test_date <- left_join(test, FDA_clean_ok_date_of_pub, by="COMMON_DRUGBANK_ALIAS")
training_date <- left_join(training, FDA_clean_ok_date_of_pub, by="COMMON_DRUGBANK_ALIAS")

test_date$date <- as.Date(test_date$OLDEST_DATE_OF_PUBLICATION_manually_cured, "%d/%m/%Y")
test_date$YEAR <- as.Date(cut(test_date$date,
                        breaks = "year"))

training_date$date <- as.Date(training_date$OLDEST_DATE_OF_PUBLICATION_manually_cured, "%d/%m/%Y")
training_date$YEAR <- as.Date(cut(training_date$date,
                              breaks = "year"))

library("ggplot2")
library("scales")

p_test <- ggplot(test_date, aes(YEAR, frequency(YEAR))) + 
  geom_histogram(stat="identity") +
  theme_bw() + xlab("years") + ylab("number of publication per year") + ggtitle("test set") 

p_training <- ggplot(training_date, aes(YEAR, frequency(YEAR))) + 
  geom_histogram(stat="identity") +
  theme_bw() + xlab("years") + ylab("number of publication per year") + ggtitle("training set") 

splots  <- list()
splots[[1]] <- p_test
splots[[2]] <- p_training
plots <- ggarrange(p_test, p_training, ncol = 2, nrow = 1)
ggsave('RESOLVED2 - training and test date distribution.tif',
       plots,width=10, height=5, 
       unit='in', dpi=600, device = "jpeg")

# 1/ compute the RESOLVED model
# SELECT K BEST FEATURES using cross validation 
# put variables into a matrix
variables_matrix = as.matrix(training[,5:ncol(training)])
# compute the surv function
y <- Surv(training$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, training$FDA_APPROVAL)
# set the folds
fold=rep(1:100,length.out=nrow(training))

# Fit Cox regression using lasso with 100-fold cross-validation, in order to better features
cox.lasso=cv.glmnet(variables_matrix
                    ,y
                    ,family="cox",
                    foldid=fold, 
                    alpha = 1)


# Plot cross-validation curve
tiff('RESOLVED2 - 22 - SUPP  Figure 3 - cross-validation curve.tif', 
     width=12, 
     height=10, unit='in', res=600, compression = "jpeg")
plot(cox.lasso)
dev.off()

# plot estimates as a function of lambda:
# = L1 normalisation
tiff('RESOLVED2 - 23 - SUPP  Figure 4 - Distribution of features coefficients according to lambda values.tif', 
    width=12, 
    height=8, unit='in', res=600, compression = "jpeg")
plot(cox.lasso$glmnet.fit,xvar="lambda")
dev.off()

# Value of lambda that gives the minimum of the cross-validation curve
cox.lasso$lambda.min
# minimum value of the cross-validation curve
min(cox.lasso$cvm)

# fit model with lambda min identified by CV.
model_final = glmnet(variables_matrix
                    ,y
                    ,family="cox"
                    ,lambda = cox.lasso$lambda.min
                    ,alpha = 1)

# Get estimated positive coefficients:
coefficients=coef(model_final, s=model_final$lambda.min)
coefficients_pos=  coefficients[coefficients[,1]>0,]

model_final_coef <- as.data.frame(coefficients_pos)
model_final_coef$category <- rownames(model_final_coef)
model_final_coef$category <- gsub("category.", "", model_final_coef$category)
model_final_coef$category <- gsub("\\.", " ", model_final_coef$category)

ggplot(model_final_coef, aes(x = reorder(category, -coefficients_pos), 
                           y = coefficients_pos, 
                           na.rm=TRUE)) +
  geom_histogram(stat="identity") +
  ggtitle("Distribution of Beta coefficient from Lasso penalized Cox model", subtitle = NULL) +
  labs(x = "Drug features", y = "Coefficient") +
  theme(axis.text.x = element_text(model_final_coef$category), 
        panel.background = element_blank() )+
  coord_flip()+ 

ggsave('RESOLVED2 - 12 - Figure 1 - model.tif', 
       width=12, height=6, unit='in', dpi=600,
         device = "tiff", compression = "jpeg")

# compute a MANUAL model
# get positive coefficients and put it in a dataframe
coefficients <- coef(model_final, s=model_final$lambda.min)
coefficients_nonzero =  coefficients[coefficients[,1]>0,]
fit_coef <- as.data.frame(coefficients_nonzero)

# save the model for further use
write.table(coefficients_nonzero, "RESOLVED2 - 31 - SUPP Table 1 - Lasso-penalized Beta coefficients.txt", sep = "\t")

# for further validation steps: load the model
coefficients_nonzero <- read.csv2("RESOLVED2 - 31 - SUPP Table 1 - Lasso-penalized Beta coefficients.txt", sep="\t", header=TRUE, dec = ".")
dim(coefficients_nonzero)
fit_coef <- as.data.frame(coefficients_nonzero)

# get features names
variables <- rownames(fit_coef)

# prepare test set with only selected features and prepare a matrix
test_select <- test %>% dplyr::select(variables)
test_mat <- as.matrix(test_select)
coefficients_nonzero <- as.matrix(coefficients_nonzero)

# compute predicted probabilities
test$probs <- exp(test_mat%*%coefficients_nonzero)

# observe probs distribution
plot(sort(test$probs, decreasing = TRUE))
hist(test$probs,breaks = 100, xlim = c(0,50))

# compute basic cindex
install.packages("survcomp")
library(survcomp)
ConcIndex<-concordance.index(x=test$probs, 
                             surv.time=Surv(test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                            test$FDA_APPROVAL), 
                             surv.event=test$FDA_APPROVAL, 
                             method="noether")
# print
output <- data.frame(C.index=ConcIndex$c.index,
                     Lower=ConcIndex$lower,
                     Higher=ConcIndex$upper,
                     pvalue=ConcIndex$p.value)
output$a <- "["; output$b <- "-"; output$c <- "]"
output <- unite(output, 'C-index [95% Confidence interval]', c(C.index, a , Lower, b, Higher , c), sep = "", remove = T)
print(output)
write.table(output, "cindex.txt", sep = "\t")

# compute better c-index IPCW
library(pec)
install.packages("htmlTable")
library(htmlTable)
install.packages("sandwich")
library(sandwich)
library(survival)

# compute cindex
variables <- read.csv2("RESOLVED2 - 31 - SUPP Table 1 - Lasso-penalized Beta coefficients.txt", sep="\t", header=TRUE)
features <- rownames(variables)
test_features = as.data.frame(test[,c(features,"DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN", "FDA_APPROVAL")])

# inverse of the probability of censoring weigthed estimate of the concordance probability to adjust for right censoring
cindex <- pec::cindex(as.matrix(1/test$probs),
            formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                           FDA_APPROVAL)~.,
            data=test,
            cens.model="marginal")
cindex

# now check whether c-index is stable regardless of subset of the test set depending on year of publication

variables <- read.csv2("RESOLVED2 - 31 - SUPP Table 1 - Lasso-penalized Beta coefficients.txt", sep="\t", header=TRUE)
features <- rownames(variables)

# select columns of interest in test set
test_features = as.data.frame(test[,c(features,"DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN", "FDA_APPROVAL","probs")])
test_features$COMMON_DRUGBANK_ALIAS <- test$COMMON_DRUGBANK_ALIAS

# get dates
FDA_clean_ok_date_of_pub <- FDA_clean_ok[,c("COMMON_DRUGBANK_ALIAS","OLDEST_DATE_OF_PUBLICATION_manually_cured")]
test_features_date <- left_join(test_features, FDA_clean_ok_date_of_pub, by="COMMON_DRUGBANK_ALIAS")

test_features_date$date <- as.Date(test_features_date$OLDEST_DATE_OF_PUBLICATION_manually_cured, "%d/%m/%Y")
test_features_date$YEAR <- as.Date(cut(test_features_date$date,
                              breaks = "year"))

# now cut test_features_date by dates - 2 cuts
test_2011to2017 <- test_features_date %>% filter(YEAR>"2010-01-01")
test_2003to2011 <- test_features_date %>% filter(YEAR>"2002-01-01"&YEAR<"2011-01-01")
test_1970to2003 <- test_features_date %>% filter(YEAR<"2003-01-01")

# reload split test sets with RESOLVED features
test_2011to2017_features = as.data.frame(test_2011to2017[,c(features,"DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN", "FDA_APPROVAL")])
test_2003to2011_features = as.data.frame(test_2003to2011[,c(features,"DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN", "FDA_APPROVAL")])
test_1970to2003_features = as.data.frame(test_1970to2003[,c(features,"DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN", "FDA_APPROVAL")])

str(test_2010to2017_features)

# compute IPCW for each test set
cindex_2011to2017 <- pec::cindex(as.matrix(1/test_2011to2017$probs),
                      formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                     FDA_APPROVAL)~.,
                      data=test_2011to2017_features,
                      cens.model="marginal")

cindex_2003to2011 <- pec::cindex(as.matrix(1/test_2003to2011$probs),
                                 formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                                FDA_APPROVAL)~.,
                                 data=test_2003to2011_features,
                                 cens.model="marginal")

cindex_1970to2003 <- pec::cindex(as.matrix(1/test_1970to2003$probs),
                                 formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                                FDA_APPROVAL)~.,
                                 data=test_1970to2003_features,
                                 cens.model="marginal")
cindex_2011to2017
cindex_2003to2011
cindex_1970to2003


##### 2/ as a comparison, we fitted another model based on conventionnal data considered in the go/no go decision: 
# clinical activity, complete response and toxicity
# select features & put variables into a matrix
# this model will be named "EffTOX" model
variables_matrix_alter = as.matrix(training[,c("Clinical_activity_detected_as.following_recurrent_responder_OR_clinical_activity_OR_prolonged_stability",
                                         "Complete_response_reported",
                                         "DLT_identified_or_MTD_reached")])
# compute the surv function
y <- Surv(training$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, training$FDA_APPROVAL)
fold=rep(1:100,length.out=nrow(training))

# compute model fitting
cox.lasso=cv.glmnet(variables_matrix_alter
                    ,y
                    ,family="cox",
                    foldid=fold, 
                    alpha = 1)

plot(cox.lasso$glmnet.fit,xvar="lambda")
plot(cox.lasso)

cox.lasso$lambda.min

model_final_alter = glmnet(variables_matrix_alter
                     ,y
                     ,family="cox"
                    ,lambda = cox.lasso$lambda.min
                    ,alpha = 1
                     )

coefficients=coef.glmnet(model_final_alter)
coefficients_pos=coefficients[coefficients[,1]!=0,]

write.table(as.data.frame(as.matrix(coefficients_pos)), "efficacity x tox model coefficients.txt", sep = "\t")

model_final_alter_coef <- as.data.frame(as.matrix(coefficients_pos))

model_final_alter_coef$category <- rownames(model_final_alter_coef)
model_final_alter_coef$category <- gsub("category.", "", model_final_alter_coef$category)
model_final_alter_coef$category <- gsub("\\.", " ", model_final_alter_coef$category)
model_final_alter_coef$category <- gsub("Clinical_activity_detected_as following_recurrent_responder_OR_clinical_activity_OR_prolonged_stability"
                                        , "Clinical activity", model_final_alter_coef$category)
model_final_alter_coef$category <- gsub("_", " ", model_final_alter_coef$category)
model_final_alter_coef$Coefficients <- model_final_alter_coef$s0

ggplot(model_final_alter_coef, aes(x = reorder(category, -Coefficients), 
                             y = Coefficients, 
                             na.rm=TRUE)) +
  geom_histogram(stat="identity") +
  ggtitle("Distribution of Beta coefficient from Cox model with conventional features", subtitle = NULL) +
  labs(x = "Features", y = "Coefficient") +
  theme(axis.text.x = element_text(model_final_alter_coef$category), 
        panel.background = element_blank() )+
  coord_flip()+ 
  
  ggsave('RESOLVED2 - 24 - SUPP  Figure 7 - conventional features coefficients efftox.tif', 
         width=12, height=6, unit='in', dpi=600,
         device = "tiff",
         compression = "jpeg")

# compute a MANUAL model
# for further validation steps: load the model
coefficients_nonzero <- read.csv2("efficacity x tox model coefficients.txt", sep="\t", header=TRUE, dec = ".")
fit_coef <- as.data.frame(coefficients_nonzero)
# get features names
variables <- rownames(fit_coef)
# prepare test set with only selected features and prepare a matrix
test_select <- test %>% dplyr::select(variables)
test_mat <- as.matrix(test_select)
coefficients_nonzero <- as.matrix(coefficients_nonzero)
# compute predicted probabilities
test$probs_efftox <- exp(test_mat%*%coefficients_nonzero)
# observe probs distribution
plot(sort(test$probs_efftox, decreasing = TRUE))
hist(test$probs_efftox,breaks = 100)
# compute basic cindex
library(survcomp)
ConcIndex_effTox<-concordance.index(x=test$probs_efftox, 
                             surv.time=Surv(test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                            test$FDA_APPROVAL), 
                             surv.event=test$FDA_APPROVAL, 
                             method="noether")
ConcIndex_effTox$c.index
# print
output_efftox <- data.frame(C.index=ConcIndex_effTox$c.index,
                     Lower=ConcIndex_effTox$lower,
                     Higher=ConcIndex_effTox$upper,
                     pvalue=ConcIndex_effTox$p.value)
output$a <- "["; output$b <- "-"; output$c <- "]"
output <- unite(output, 'C-index [95% Confidence interval]', c(C.index, a , Lower, b, Higher , c), sep = "", remove = T)
print(output)
write.table(output, "model de reference eff x tox/cindex_from_effXtox_model.txt", sep = "\t")

# compute better c-index IPCW
library(pec)
#install.packages("htmlTable")
library(htmlTable)
#install.packages("sandwich")
library(sandwich)

# compute cindex
variables <- read.csv2("efficacity x tox model coefficients.txt", sep="\t", header=TRUE)
features <- rownames(variables)
test_features = as.data.frame(test[,c(features,"DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN", "FDA_APPROVAL")])
# inverse of the probability of censoring weigthed estimate of the concordance probability to adjust for right censoring
cindex_IPCW_effTox <- pec::cindex(as.matrix(1/test$probs_efftox),
                      formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                     FDA_APPROVAL)~.,
                      data=test,
                      cens.model="marginal"
           #           splitMethod = "BootCv"
                      )
cindex_IPCW_effTox

# compute AUC for both model: COMPUTE TIME DEPENDANT AUROC 
library(survivalROC)
surv_roc_resolved_3 = survivalROC(Stime= test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN,  
                     status= test$FDA_APPROVAL,      
                     marker = test$probs,     
                     predict.time = 3, 
                     cut.values = test$probs,
                     method="KM")
surv_roc_resolved_3$AUC

surv_roc_resolved_6 = survivalROC(Stime= test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN,  
                                  status= test$FDA_APPROVAL,      
                                  marker = test$probs,     
                                  predict.time = 6, 
                                  cut.values = test$probs,
                                  method="KM")
surv_roc_resolved_6$AUC

surv_roc_efftox_3 = survivalROC(Stime= test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN,  
                                  status= test$FDA_APPROVAL,      
                                  marker = test$probs_efftox,     
                                  predict.time = 3, 
                                  cut.values = test$probs_efftox,
                                  method="KM")
surv_roc_efftox_3$AUC

surv_roc_efftox_6 = survivalROC(Stime= test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN,  
                                  status= test$FDA_APPROVAL,      
                                  marker = test$probs_efftox,     
                                  predict.time = 6, 
                                  cut.values = test$probs_efftox,
                                  method="KM")
surv_roc_efftox_6$AUC


tiff('RESOLVED2 - 24 - SUPP  Figure 6 - ROC curves.tif',width=8, height=8, 
     unit='in', res = 600, compression = "jpeg")
par(mfrow=c(2,2))

plot1 <- plot(surv_roc_resolved_3$FP, surv_roc_resolved_3$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
              xlab=paste( "FP", "", "
AUC = ",round(surv_roc_resolved_3$AUC,3)), 
              ylab="TP",main="A. RESOLVED2 - 3 years")+abline(0,1)

plot3 <- plot(surv_roc_efftox_3$FP, surv_roc_efftox_3$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
              xlab=paste( "FP", "", "
AUC = ",round(surv_roc_efftox_3$AUC,3)), 
              ylab="TP",main="B. Efficacy & toxicity model - 3 years")+ abline(0,1)

plot2 <- plot(surv_roc_resolved_6$FP, surv_roc_resolved_6$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
              xlab=paste( "FP", "", "
AUC = ",round(surv_roc_resolved_6$AUC,3)), 
              ylab="TP",main="C. RESOLVED2 - 6 years")+ abline(0,1)

plot4 <- plot(surv_roc_efftox_6$FP, surv_roc_efftox_6$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
              xlab=paste( "FP", "", "
AUC = ",round(surv_roc_efftox_6$AUC,3)), 
              ylab="TP",main="D. Efficacy & toxicity model - 6 years")+abline(0,1)
dev.off()


# 3/ third control: compare the cox regression approach to a model base on RF
test <- read.csv2("test_set.txt" , sep="\t", header=TRUE, dec = ".")
training <- read.csv2("training_set.txt", sep="\t", header=TRUE, dec = ".")
test <- test[,2:1416]
training <- training[,2:1416]

dim(FDA_clean_ok_drugbank)

colnames(FDA_clean_ok_drugbank[,1:10])
FDA_clean_ok_drugbank_ok <- FDA_clean_ok_drugbank[,3:1415]
colnames(FDA_clean_ok_drugbank_ok[,1:10])
#install.packages("randomForestSRC")
library(randomForestSRC)

# data format
# DATA: a data frame or matrice with covariables and outcome values in survival format (OS,STATUS)

# run the rf ; https://cran.r-project.org/web/packages/survxai/vignettes/How_to_compare_models_with_survxai.html
# the ~. means that you use all remainig data of your DATA matrice to compute the model (except OS and STATUS for sure)
# check further hyperparmaters that you can tune in the information: https://cran.r-project.org/web/packages/randomForestSRC/randomForestSRC.pdf
# details: https://kogalur.github.io/randomForestSRC/theory.html

train <- sample(1:nrow(FDA_clean_ok_drugbank_ok), 
                round(nrow(FDA_clean_ok_drugbank_ok) * 0.70))

RF_RESOLVED2 <- rfsrc(Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN,
                          FDA_APPROVAL) ~ ., data = FDA_clean_ok_drugbank_ok[train, ], 
              ntree = 500, importance = T)

print(RF_RESOLVED2)

# check relative importance of the varibales and errors
vimp(RF_RESOLVED2)

# Strong variables have minimal depth less that or equal to the threshold.
v.max <- max.subtree(RF_RESOLVED2)
print(v.max$threshold)

# This corresponds to the top set of variables above the threshold.
print(v.max$topvars)

# check model
tiff('RESOLVED2 - 25 - SUPP RESULTS FIG 1 - RF PLOT.tif', 
     width=12, 
     height=12, unit='in', res=600 #compression = "jpeg"
     )
plot.survival(RF_RESOLVED2, plots.one.page = TRUE,
              show.plots = TRUE, 
              span = "cv", cens.model = "rfsrc")
dev.off()

# identify features selected in the model
tiff('RESOLVED2 - 25 - SUPP RESULTS FIG 2 - RF FEATURES.tif', 
     width=20, 
     height=10, unit='in', res=600)
plot.rfsrc(RF_RESOLVED2, plots.one.page = FALSE)
dev.off()

RF_RESOLVED2_predictions <- predict(RF_RESOLVED2, FDA_clean_ok_drugbank_ok[-train , ])

print(RF_RESOLVED2_predictions)

cindex_IPCW_RF_RESOLVED2 <- pec::cindex(RF_RESOLVED2,
                                  formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                                 FDA_APPROVAL)~.,
                                  data=FDA_clean_ok_drugbank_ok[-train , ],
                                  cens.model="marginal")
cindex_IPCW_RF_RESOLVED2


########## the next part of the script is aiming to calibrate and validate the classification model
########## from the predicted probabilities.
# cut-off/threshold calibration is performed on the training set using log-rank test
# cut-off/threshold validation is performed on the test set to estimate the predictive performance.

# compute a function to annotate each drug as predicted approved or unapproved
score_model = function(df, threshold){
  df$score = ifelse(df$probs < threshold, 'predicted non-approved', 'predicted approved')
  df
  }
# example
test = score_model(test, 4.553632)
# check output
test[1:20, c('FDA_APPROVAL','probs','score')]

########### threshold calibration using logrank: LR test's pvalue for each threshold

# first compute predicted probabilities on training set
coefficients_nonzero <- read.csv2("RESOLVED2 - 31 - SUPP Table 1 - Lasso-penalized Beta coefficients.txt", sep="\t", header=TRUE, dec = ".")
fit_coef <- as.data.frame(coefficients_nonzero)
# get features names
variables <- rownames(fit_coef)
# prepare test set with only selected features and prepare a matrix
training_select <- training %>% dplyr::select(variables)
training_mat <- as.matrix(training_select)
coefficients_nonzero <- as.matrix(coefficients_nonzero)
# compute predicted probabilities
training$probs <- exp(training_mat%*%coefficients_nonzero)

# observe probs distribution
plot(sort(training$probs, decreasing = TRUE))
max(training$probs)
hist(training$probs,breaks = 100, xlim = c(0,50))

# create a sequence of threshold to test with exponential distribution to match with probs distribution
lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

log_sequence <- lseq(from = min(training$probs)+0.1, to = 20, length.out = 50)
hist(log_sequence)  

data <- training

output_threshold <- data.frame(matrix(ncol = 3, nrow = length(log_sequence)))
colnames(output_threshold) <- c("threshold","p","chisq")
output_threshold$threshold <- log_sequence

for (i in 1:length(log_sequence)) {
  
  data = score_model(data, output_threshold$threshold[i])
  
  FDA_FS_threshold_i <- survdiff(Surv(data$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
                                    data$FDA_APPROVAL) ~ data$score)
  
  output_threshold[i,2] <- pchisq(FDA_FS_threshold_i$chisq, length(FDA_FS_threshold_i$n)-1, lower.tail = FALSE)
  output_threshold[i,3] <- FDA_FS_threshold_i$chisq
}

plot(sort(log(output_threshold[1:1000,2]), decreasing = FALSE), xlab = "", ylab = "log(pvalue)")
plot(sort(output_threshold[1:1000,3], decreasing = FALSE), xlab = "", ylab = "chisq")
plot(output_threshold[1:1000,3] , log(output_threshold[1:1000,2]), 
     ylab = "log(pvalue)", xlab = "chisq")

# identify minimal pvalue
sorted_by_p <- output_threshold[order(output_threshold$p),]
optimal_threshold = sorted_by_p$threshold[1]
View(sorted_by_p)


tiff('RESOLVED2 - 24 - SUPP  Figure 5 - threshold calibration using log-rank test pvalue.tif', 
     width=12, 
     height=10, unit='in', res=600, compression = "jpeg")
plot(output_threshold[1:1000,1] , log(output_threshold[1:1000,2]), 
     xlab = "threshold", ylab = "log(pvalue)")+ title("Minimal pvalue for threshold = 4.553632")
dev.off()

plot(output_threshold[1:1000,1] , output_threshold[1:1000,3], 
     xlab = "threshold", ylab = "chisq")


# visualize kaplan meier curve on training set
training = score_model(training, 4.553632)
data <- training

FDA_FS_by_prediction <- survfit(Surv(data$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
                                     data$FDA_APPROVAL) ~ data$score)

summary(survfit(Surv(data$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
                     data$FDA_APPROVAL) ~ data$score),times=c(6))
library(survminer)
library(ggplot2)
FDA_FS_plot_training <- ggsurvplot(
  FDA_FS_by_prediction,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  risk.table.height = 0.2,
  legend.labs = c("predicted approved", "predicted non-approved"),
  legend.title = "",
  legend = "none",
  palette = c("dark grey", "black"),
  #linetype = c(1,2,3),
  data = data,  # data used to fit survival curves. 
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  break.time.by = 2,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  ,title = "A"
  ,xlab = "Years"
  ,ylab = "FDA approval free survival (%)"
  , xlim = c(0,36)
  , surv.median.line =c("hv")
  , tables.col = "strata"
  # , pval.coord = c(240, 0.25)
) 

splots  <- list()
splots[[1]] <- FDA_FS_plot_training
plots <- arrange_ggsurvplots(splots, print = FALSE, title = NA, ncol = 1, nrow = 1)
ggsave('RESOLVED2 - 15 - Figure 4a. FDA_FS~prediction on training set.tif',plots,width=14, height=8, 
       unit='in', dpi=600, device = "jpeg")

# then used the calibrated threshold with the model to classify drugs in the test set and visualize 
# compute probabilities on the test set
coefficients_nonzero <- read.csv2("RESOLVED2 - 31 - SUPP Table 1 - Lasso-penalized Beta coefficients.txt", sep="\t", header=TRUE, dec = ".")
fit_coef <- as.data.frame(coefficients_nonzero)
# get features names
variables <- rownames(fit_coef)
# prepare test set with only selected features and prepare a matrix
test_select <- test %>% dplyr::select(variables)
test_mat <- as.matrix(test_select)
coefficients_nonzero <- as.matrix(coefficients_nonzero)
# compute predicted probabilities
test$probs <- exp(test_mat%*%coefficients_nonzero)
# observe probs distribution
plot(sort(test$probs, decreasing = TRUE))
hist(test$probs,breaks = 100, xlim = c(0,50))
# classify according to the optimal threshold calibrated in the training set
optimal_threshold = 4.553632
test = score_model(test, optimal_threshold)
# visualize
data <- test

sum(test$score=="predicted non-approved")/nrow(test)

FDA_FS_by_prediction <- survfit(Surv(data$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
                                     data$FDA_APPROVAL) ~ data$score)

summary(survfit(Surv(data$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
                     data$FDA_APPROVAL) ~ data$score),times=c(6))

FDA_FS_plot_test <- ggsurvplot(
  FDA_FS_by_prediction,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  risk.table.height = 0.2,
  legend.labs = c("predicted approved", "predicted non-approved"),
  legend.title = "",
  legend = "none",
  palette = c("dark grey", "black"),
  #linetype = c(1,2,3),
  data = data,  # data used to fit survival curves. 
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  break.time.by = 2,     # break X axis in time intervals by 500.
  ggtheme = theme_classic(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = TRUE # show bars instead of names in text annotations
  ,title = "B"
  ,xlab = "Years"
  ,ylab = "FDA approval free survival (%)"
  #  , xlim = c(2,20)
  , surv.median.line =c("hv")
  , tables.col = "strata"
  # , pval.coord = c(240, 0.25)
)


splots  <- list()
splots[[1]] <- FDA_FS_plot_test
plots <- arrange_ggsurvplots(splots, print = FALSE, title = NA, ncol = 1, nrow = 1)
ggsave('RESOLVED2 - 15 - Figure 4b. FDA_FS~prediction on test set.tif',plots,width=14, height=8, 
       unit='in', dpi=600, device = "jpeg")


### save plot in Figure 4A/B
splots  <- list()
splots[[1]] <- FDA_FS_plot_training
splots[[2]] <- FDA_FS_plot_test
plots <- arrange_ggsurvplots(splots, print = FALSE, title = NA, ncol = 2, nrow = 1)
ggsave('RESOLVED2 - 14 - Figure 3 - FDA_FS_V2.tif',
       plots,width=20, height=7, 
       unit='in', dpi=600, device = "jpeg")

## perform coxph function in the test set to compute HR of relative likelihood of approval for predicted drugs
test$score <- factor(test$score, levels = c("predicted non-approved","predicted approved"))

summary(coxph(Surv(test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
     test$FDA_APPROVAL) ~ test$score))

# check PH assumption
mod <- coxph(Surv(test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
                  test$FDA_APPROVAL) ~ test$score)
plot(cox.zph(mod))





###################### test on other drugs ######################
setwd("F:/actif E&R IGR DITEP/RESOLVED²/SCRIPTS_WD/data_final/")
test <- read.csv2("RESOLVED2 - 33 - SUPP Table 4 - cases - test CM.txt", sep="\t", header=TRUE, dec = ".")

#View(test)
coefficients_nonzero <- read.csv2("RESOLVED2 - 31 - SUPP Table 1 - Lasso-penalized Beta coefficients.txt", sep="\t", header=TRUE, dec = ".")
fit_coef <- as.data.frame(coefficients_nonzero)
# get features names
variables <- rownames(fit_coef)
# prepare test set with only selected features and prepare a matrix
library(dplyr)
test_select <- test %>% dplyr::select(variables)
test_mat <- as.matrix(test_select)
coefficients_nonzero <- as.matrix(coefficients_nonzero)
# compute predicted probabilities
test$probs <- exp(test_mat%*%coefficients_nonzero)
# observe probs distribution
plot(sort(test$probs, decreasing = TRUE))
hist(test$probs,breaks = 100, xlim = c(0,50))
# classify according to the optimal threshold calibrated in the training set
optimal_threshold = 4.553632

score_model = function(df, threshold){
  df$score = ifelse(df$probs < threshold, 'predicted non-approved', 'predicted approved')
  df
}

test = score_model(test, optimal_threshold)
# output
test[,c("COMMON_DRUGBANK_ALIAS","score", "probs")]
write.table(test, "predictions_on_test.txt", sep = "\t")


# compare mean score of FP and TP
test_with_prediction <- read.csv2("predictions_on_test.txt", sep="\t", header=TRUE, dec = ".")
View(test_with_prediction)

test_with_prediction <- test_with_prediction %>% filter(FPTP != "")

t.test(test_with_prediction$probs ~ test_with_prediction$FPTP)
