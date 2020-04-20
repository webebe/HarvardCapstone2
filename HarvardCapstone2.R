if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(dslabs)) install.packages("dslabs", repos = "http://cran.us.r-project.org")
if(!require(HistData)) install.packages("HistData", repos = "http://cran.us.r-project.org")
if(!require(ggpubr)) install.packages("ggpubr", repos = "http://cran.us.r-project.org")
if(!require(broom)) install.packages("broom", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org")
if(!require(purrr)) install.packages("purrr", repos = "http://cran.us.r-project.org")
if(!require(pdftools)) install.packages("pdftools", repos = "http://cran.us.r-project.org")
if(!require(matrixStats)) install.packages("matrixStats", repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require(randomForest)) install.packages("randomForest", repos = "http://cran.us.r-project.org")
if(!require(Rborist)) install.packages("Rborist", repos = "http://cran.us.r-project.org")
if(!requirey(rpart)) install.packages("rpart", repos = "http://cran.us.r-project.org")
if(!require(rpart.plot)) install.packages("rpart.plot", repos = "http://cran.us.r-project.org")
if(!require(gam)) install.packages("gam", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org")
if(!require(RColorBrewer)) install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(party)) install.packages("party", repos = "http://cran.us.r-project.org")
library(party)
library(tidyverse)
library(dslabs)
library(HistData)
library(ggpubr)
library(broom)
library(caret)
library(lubridate)
library(purrr)
library(pdftools)
library(matrixStats)
library(dplyr)
library(randomForest)
library(Rborist)
library(rpart)
library(rpart.plot)
library(gam)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(gam)
library(data.table)

#Read Dataset
url_file<-"https://raw.githubusercontent.com/webebe/HarvardCapstone2/master/chronic_kidney_disease.txt"
chronic_kidney_disease<-read.csv(url(url_file))

#Many observations include NA values in several variables, therefor the dataset was cleaned.
#####Data cleaning
dat<-chronic_kidney_disease
#clean dataset: factorize variable where possible; adjust obviously wrong entries;
dat<-dat[-which(dat$age=="notckd"),]
dat$age<-as.numeric(dat$age)

dat$sg<-factor(dat$sg)
dat$su<-factor(dat$su)
dat$pcv<-as.numeric(dat$pcv)
dat$wbcc<-as.numeric(dat$wbcc)
dat$rbcc<-as.numeric(dat$rbcc)

ind<-dat$dm=="\tno"
dat$dm[ind]="no"
ind2<-dat$dm==" yes"
dat$dm[ind2]="yes"
ind<-dat$dm=="\tyes"
dat$dm[ind]="yes"
dat$dm<-factor(dat$dm)

ind<-dat$cad=="\tno"
dat$cad[ind]="no"
ind2<-dat$cad=="\tyes"
dat$cad[ind2]="yes"
dat$cad<-factor(dat$cad)
ind<-dat$appet=="no"  
dat$appet[ind]="poor"
dat$appet<-factor(dat$appet)
ind<-dat$pe=="good"
dat$pe[ind]="no"
dat$pe<-factor(dat$pe)

ind<-dat$class=="ckd\t"
dat$class[ind]="ckd"
ind<-dat$class=="no"
dat$class[ind]="notckd"
dat$class<-factor(dat$class)

##dealing with missing data
x<-seq(1,nrow(dat))
missingValues_mean<-sapply(x, function(x) mean(is.na(dat[x,])))
boxplot(missingValues_mean, main = "Missing values mean")
summary(missingValues_mean)
#decision: Ignore observations with more than 16% missing values ->The 0.75% quantile
missing_index<-which(missingValues_mean>quantile(missingValues_mean,0.75))
dat_clean<-dat[-missing_index,]

#Using Mean Substitution to handle missing values in continous variables, median in categorical vars
#Use the mean values of the uncleaned dataset in order to minimize variations in the analysis dataset
#define a Mode function for replacing missing values in factor variables with the median value
Mode <- function (x, na.rm) {
  xtab <- table(x)
  xmode <- names(which(xtab == max(xtab)))
  if (length(xmode) > 1) xmode <- ">1 mode"
  return(xmode)
}

#age
dat_clean[which(is.na(dat_clean$age)),"age"]<-round(mean(dat$age,na.rm = TRUE))
#bp
dat_clean[which(is.na(dat_clean$bp)),"bp"]<-round(mean(dat$bp,na.rm = TRUE))
#sg
dat_clean[is.na(dat_clean[,"sg"]),"sg"] <- Mode(dat[,"sg"], na.rm = TRUE)
#su
dat_clean[is.na(dat_clean[,"su"]),"su"] <- Mode(dat[,"su"], na.rm = TRUE)


#pc column was removed due to missing values and no way of knowing how to substitute these values
pc_index<-which(colnames(dat_clean)=="pc")
dat_clean<-dat_clean[,-pc_index]

#pcc  -->also remove
pcc_index<-which(colnames(dat_clean)=="pcc")
dat_clean<-dat_clean[,-pcc_index]

#ba --<also remove
ba_index<-which(colnames(dat_clean)=="ba")
dat_clean<-dat_clean[,-ba_index]

#bgr
dat_clean[which(is.na(dat_clean$bgr)),"bgr"]<-round(mean(dat$bgr,na.rm = TRUE))
#bu
dat_clean[which(is.na(dat_clean$bu)),"bu"]<-round(mean(dat$bu,na.rm = TRUE))
#sc
dat_clean[which(is.na(dat_clean$sc)),"sc"]<-round(mean(dat$sc,na.rm = TRUE))
#sod
dat_clean[which(is.na(dat_clean$sod)),"sod"]<-round(mean(dat$sod,na.rm = TRUE))
#pot
dat_clean[which(is.na(dat_clean$pot)),"pot"]<-round(mean(dat$pot,na.rm = TRUE))
#hemo
dat_clean[which(is.na(dat_clean$hemo)),"hemo"]<-round(mean(dat$hemo,na.rm = TRUE))
#pcv
dat_clean[which(is.na(dat_clean$pcv)),"pcv"]<-round(mean(dat$pcv,na.rm = TRUE))
#wbcc
dat_clean[which(is.na(dat_clean$wbcc)),"wbcc"]<-round(mean(dat$wbcc,na.rm = TRUE))
#rbcc
dat_clean[which(is.na(dat_clean$rbcc)),"rbcc"]<-round(mean(dat$rbcc,na.rm = TRUE))

#rbc ->calcualte from rbcc
range(which(dat$rbc=="normal"))
range(which(dat$rbc=="abnormal"))
mean(which(dat$rbc=="normal"))
mean(which(dat$rbc=="abnormal"))
#fill NA values in rbc with 'normal' if closer to the mean of other normal rbcc values and with 'abnormal' 
dat_clean$rbc<-ifelse(is.na(dat_clean$rbc),
                ifelse(abs(dat_clean$rbcc-mean(which(dat$rbc=="normal")))>abs(dat_clean$rbcc-mean(which(dat$rbc=="abnormal"))),
                       "2","3"),
                dat_clean$rbc)

dat_clean$rbc
dat_clean$rbc<-ifelse(dat_clean$rbc=="2","abnormal","normal")%>%factor()

#After manually investigating the dataset, some observations where removed due to lack of input in more than 2 variables
which(is.na(dat_clean$htn))
which(is.na(dat_clean$dm))
which(is.na(dat_clean$cad))
which(is.na(dat_clean$appet))
which(is.na(dat_clean$pe))
which(is.na(dat_clean$ane))
which(is.na(dat_clean$ba))
table(dat_clean$ba)

#factor varibles again to get rid of unused factors
dat_clean$htn<-factor(dat_clean$htn)
dat_clean<-dat_clean[-which(dat_clean$dm==""),]
dat_clean$dm<-factor(dat_clean$dm)
dat_clean$ane<-factor(dat_clean$ane)

dat_clean<-dat_clean[-which(is.na(dat_clean$htn)),]
dat_clean<-dat_clean[-which(is.na(dat_clean$pe)),]

#After cleaning the dataset we end up with a set of 300 Observations of 21 variables.
abbreviation<-c("age","bp","sg","su","rbc","bgr","bu","sc","sod","pot","hemo","pcv","wbcc","rbcc","htn","dm",
                "cad","appet","pe","ane","class")

vars<-c("age","blood pressure","specific gravity","sugar","red blood cells",
        "blood glucose random","blood urea","serum creatinine","sodium","potassium","hemoglobin","packed cell volume",
        "white blood cell count","red blood cell count","hypertension","diabetes mellitus","coronary artery disease",
        "appetite","pedal edema","anemia","class")

tabeldat <- data_frame(abreviation=abbreviation, description = vars )
tabeldat %>% knitr::kable()

####Methods

#The class variable contains the outcome we want to predict using the predictors.

#summarize and visualize data
head(dat_clean)
#age and CKD
dat_clean%>%ggplot(aes(x=age))+geom_histogram(aes(color=class,fill=class))
#bp and CKD
dat%>%ggplot(aes(y=bp))+geom_point(aes(color=class,fill=class),stat="count")
#SG and CKD
dat_clean%>%ggplot(aes(x=sg))+geom_histogram(aes(color=class,fill=class),stat="count")
#al and CKD
dat_clean%>%ggplot(aes(x=al))+geom_histogram(aes(color=class,fill=class),stat="count")
#RBC and CKD
dat_clean%>%ggplot(aes(x=rbc))+geom_histogram(aes(color=class,fill=class),stat="count")
#Bgr and CKD
dat_clean%>%ggplot(aes(y=bgr))+geom_point(aes(color=class,fill=class),stat="count")
#bu and CKD
dat_clean%>%ggplot(aes(y=bu))+geom_point(aes(color=class,fill=class),stat="count")
#sc and CKD
dat_clean%>%ggplot(aes(y=sc))+geom_point(aes(color=class,fill=class),stat="count")
#sod and CKD
dat_clean%>%ggplot(aes(y=sod))+geom_point(aes(color=class,fill=class),stat="count")
#pot  and CKD
dat_clean%>%ggplot(aes(y=pot))+geom_point(aes(color=class,fill=class),stat="count")
#hemo and CKD
dat_clean%>%ggplot(aes(y=hemo))+geom_point(aes(color=class,fill=class),stat="count")
#pcv  and CKD
dat_clean%>%ggplot(aes(y=pcv))+geom_point(aes(color=class,fill=class),stat="count")
#wbcc and CKD
dat_clean%>%ggplot(aes(y=wbcc))+geom_point(aes(color=class,fill=class),stat="count")
#rbcc and CKD
dat_clean%>%ggplot(aes(y=rbcc))+geom_point(aes(color=class,fill=class),stat="count")
#htn and CKD
dat_clean%>%ggplot(aes(x=htn))+geom_histogram(aes(color=class,fill=class),stat="count")
#dm  and CKD
dat_clean%>%ggplot(aes(x=dm))+geom_histogram(aes(color=class,fill=class),stat="count")
#cad  and CKD
dat_clean%>%ggplot(aes(x=cad))+geom_histogram(aes(color=class,fill=class),stat="count")
#appet  and CKD
dat_clean%>%ggplot(aes(x=appet))+geom_histogram(aes(color=class,fill=class),stat="count")
#pe  and CKD
dat_clean%>%ggplot(aes(x=pe))+geom_histogram(aes(color=class,fill=class),stat="count")
#ane  and CKD
dat_clean%>%ggplot(aes(x=ane))+geom_histogram(aes(color=class,fill=class),stat="count")

#From the visualisation above we see that there are several variables indicating predictive power to determine the class
#The dataset will be split into a test and a training dataset. Algorithms will be trained on the train dataset and then
#applied to the test set. Performance of diffent model will be compared by accuracy, specificity and sensitivity.

#Create Test and training dataset
test_index<-createDataPartition(dat_clean$class,times=1,p=0.2,list = FALSE)
test_set<-slice(dat_clean[test_index,])
train_set<-slice(dat_clean[-test_index,])
train_x<-train_set[,-21]
test_x<-test_set[,-21]
train_y<-train_set[,21]
test_y<-test_set[,21]
#Check if Class distribution if even in the two subsets
mean(test_set$class=="ckd")
mean(train_set$class=="ckd")

###Logistic Regression
fit_logreg<-glm(as.numeric(train_set$class=="notckd")~.,data = train_set, family="binomial")
y_hat_logreg<-predict(fit_logreg,newdata = test_set,type="response")
y_hat_logreg<-ifelse(y_hat_logreg>0.5,"notckd","ckd")%>%factor()
acc_logreg<-confusionMatrix(y_hat_logreg,test_set$class)$overall["Accuracy"]

###Knn
fit_knn <- knn3(train_y~.,data=train_set,k=5)
y_hat_knn<-predict(fit_knn, newdata = test_set, type="class")
acc_Knn<-confusionMatrix(y_hat_knn,reference = test_y)$overall["Accuracy"]

###LDA
fit_lda<-train(class~.,method= "lda", data= train_set)
y_hat_lda<-predict(fit_lda,newdata = test_set)
acc_lda<-confusionMatrix(y_hat_lda,test_set$class)$overall["Accuracy"]

###loess
fit_loess <- train(train_x,train_y,method="gamLoess")
y_hat_loess<-predict(fit_loess,test_x)
acc_loess<-confusionMatrix(y_hat_loess,test_y)$overall["Accuracy"]

###Random Forrest
set.seed(9,sample.kind = "Rounding")
fit_rf<-train(train_x,train_y,method="rf",tuneGrid = data.frame(mtry = seq(3,9,2)), importance=TRUE)
plot(fit_rf)
y_hat_rf<-predict(fit_rf,test_x)
acc_rf<-confusionMatrix(y_hat_rf,test_y)$overall["Accuracy"]
#Variable importance helps with interpretation of archived results
varImp(fit_rf)

##Visualize Forest plot
x <- ctree(class ~ ., data=train_set)
plot(x, type="simple")

###Random forest with cross validation to choose tuning parameter
fit_rf2<-train(class~.,method="Rborist",tuneGrid=data.frame(predFixed=2,minNode=c(3,21)),data=train_set)
y_hat_rf2<-predict(fit_rf2,test_set)
acc_rf2<-confusionMatrix(y_hat_rf2,test_y)$overall["Accuracy"]

###Results
classification_results <- data_frame(Method = "Logistic Regression", Overall_Accuracy = acc_logreg)
classification_results <- bind_rows(classification_results,
                          data_frame(Method="KNN",
                                     Overall_Accuracy = acc_Knn ))
classification_results <- bind_rows(classification_results,
                                    data_frame(Method="LDA",
                                               Overall_Accuracy = acc_lda ))
classification_results <- bind_rows(classification_results,
                                    data_frame(Method="Loess",
                                               Overall_Accuracy = acc_loess ))
classification_results <- bind_rows(classification_results,
                                    data_frame(Method="Random Forest",
                                               Overall_Accuracy = acc_rf ))
classification_results <- bind_rows(classification_results,
                                    data_frame(Method="Random Forest Cross Validated",
                                               Overall_Accuracy = acc_rf2 ))
classification_results %>% knitr::kable()
