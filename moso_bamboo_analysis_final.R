library(dplyr)       # for data manipulation (dplyr) 
library(broom)       # for making model summary tidy
library(visreg)      # for potting logodds and probability 
library(margins)     # to calculate Average Marginal Effects
library(rcompanion)  # to calculate pseudo R2
library(ROCR)        # to compute and plot Reciever Opering Curve
library(yardstick)   # to estimate the accuracy of model
# load the diabetes dataset
#getwd()
#setwd("E:\Dipesh_HP_office_i3\nfolgio\R programming")
#dat<-read.csv("Book.csv")

getwd()
setwd("E:/Dipesh_HP_office_i3/nfolgio/R programming/Data")
dat <- read.csv("Book_1.csv")
head(dat)



# converting node as factor
dat$node <-as.factor(dat$node)

# See the data strcuture 
glimpse(dat)

## Plot
plot(sort(dat$germination))



# Train and split
#whole data set split into 80% train and 20% test data set
# Total number of rows in the Diabetes data frame
n <- NROW(dat$germination)
n
# Number of rows for the training set (80% of the dataset)
n_train <- round(0.80 * n)

# Create a vector of indices which is an 80% random sample
set.seed(123)
train_indices <- sample(1:n, n_train)
# Subset the germination data frame to training indices only
train <- dat[train_indices, ]
# Exclude the training indices to create the test set
test <- dat[-train_indices, ]
paste("train sample size: ", dim(train)[1])
paste("test sample size: ", dim(test)[1])

#Fitting a binary logistic regression
model_logi <- glm(germination ~ node + treatment, data = dat, family = "binomial")
#Model summary
summary(model_logi)

# Pseudo R_squared values and Likelyhood ratio test
nagelkerke(model_logi)

# Probabilities of germination wrt node
visreg(model_logi, "node", scale="response", rug=2, xlab="node level",
       ylab="P(germination)")

# Probabilities of germination wrt treatment
visreg(model_logi, "treatment", scale="response", rug=2, xlab="treatment level",
       ylab="P(germination)")

#obtaining odds log coefficient
(exp(coef(model_logi))) #obtaining odds ratio

tidy(model_logi, exponentiate = TRUE, conf.level = 0.95) #odds ratio

# Calculate average marginal effect
effects_logit_dia = margins(model_logi)
# Summary of marginal effect
summary(effects_logit_dia)

plot(effects_logit_dia) #plotting marginal effects
# predict the test dataset

# model evaluation
#Confusion matrix and classification accuracy
pred <- predict(model_logi,test, type = "response") 
predicted <- round(pred) # round of the value; >0.5 will convert to 
# 1 else 0
# Creating a contigency table
tab <- table(Predicted = predicted, Reference = test$germination)
tab

# Creating a dataframe of observed and predicted data
test$germination <- as.factor(test$germination)
act_pred <- data.frame(observed = test$germination, predicted = factor(predicted))
# Calculating Accuracy
accuracy_est <- accuracy(act_pred, observed, predicted)
#rlang::last_error()
#print(accuracy_est)
accuracy_est
# Note: the accuracy( ) function gives rounded estimate i.e. 0.764

##cutts off value vs. accuracy
pred.rocr <- prediction(pred, test$germination)
eval <- performance(pred.rocr,"acc")
plot(eval)

## finding threshold 
#res <- predict(model_logi,train, type = "response")
#ROCRpred <- prediction(res, train$germination)
#ROCRperf <- performance(ROCRpred, "tpr", "fpr")
#plot(ROCRperf, colorize = TRUE, print.cutoffs.at = seq(0.1, by = 0.1))



# Identifying the best cutoff value that maximizes accuracy
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max] #y.values are accuracy 
#measures
cut <- slot(eval, "x.values")[[1]][max] #x.values are cutoff 
#measures
print(c(Accuracy = acc, Cutoff = cut))

# Pseudo R_squared values and Likelyhood ratio test
nagelkerke(model_logi)

library(yardstick)
# Creating a actual/observed vs predicted dataframe
act_pred <- data.frame(observed = test$germination, predicted =  
                         factor(predicted))


# Calculating precision, recall and F1_score

prec <- precision(act_pred, observed, predicted)
rec <- recall(act_pred, observed, predicted)
F1_score <- f_meas(act_pred, observed, predicted) #called f_measure
print(prec)
print(rec)
print(F1_score)

## Area under reciever operating curve
perf_rocr <- performance(pred.rocr, measure = "auc",
                         x.measure = "cutoff")
perf_rocr@y.values[[1]] <- round(perf_rocr@y.values[[1]], digits = 
                                   4)
perf.tpr.fpr.rocr <- performance(pred.rocr, "tpr", "fpr")
plot(perf.tpr.fpr.rocr, colorize=T, 
     main = paste("AUC:", (perf_rocr@y.values)))
abline(a = 0, b = 1)

