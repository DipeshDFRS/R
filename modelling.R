#Day 1 8-100
#Day 2 101-350
#Day 3 351-470
#Day 4 480-700
#Day 5 700-
#######################################################################
# Scientific Data Analysis and Report Writing Training
# Forest Researcha and Training Centre
# 2076/5/20-26
# Resource Person: Prakash Aryal, Goldengate International College
######################################################################
#Day I
#Installation
#Install R Cran (RGUI)
#Overview of R
#Simple mathematics e.g 2+2 enter
# Install R Studio
# Perform Mathematics
8^(1/3)
625^(1/4)
1000^(1/7)
round (10.35365456476,3)

######################################################################
# Getting help
?anova
??anova
demo(graphics)
getwd()
setwd("C:/Users/hp/Desktop/R_Training")
#####################################################################
#create data
#check what are vector and scaler data?
# create  vector data
a<-1:20
a
# check class and structure
class(a)
str(a)

b<-21:40
class(b)
str(b)

plot(a,b)


#################################################
#c-concatenate (link)
a1<-c(1,2,4,6,8,9,20)
b1<-c(2,4,6,8,9,6,14)
c1<-c("a","b","c","d","e","f","g")
class(c1)
str(c1)

plot(a1,c1)

c11<-as.numeric(c1)
class(c11)

plot(a1,c11)
##we cannot plot such vectors where one is numeric and another is character, though
# the class is changed from character to numeric
str(c11)  # no value inside

mean(a1)
max(a1)
min(a1)
range(a1) # min max
var(a1)   #variance
sum(a1)
median(a1)
sort(a1)
rank(a1)
sd(a1)


#####################################################################
# Forestry data example
# create two data: one DBH and another Height
DBH<- c(10,14,13,18,17,16,19,20,22,19)
H<-c(5,5.5,6,6.7,7,6.1,6.4,7.2,8.3,6.9)
#Now lets see the plot
plot(DBH,H)


#for generating a model, install package lm.beta
mod1<-lm(DBH~H)
summary(mod1)
plot(DBH,H)
abline(mod1, lwd=2, xlab="DBH", ylab="Height", main="model1")



#obtaining and installing packages
# packages
# install
#citing packages
citation("lattice")
citation("ggplot2")

#see how many/ which objects are active in R environment
ls()
objects()

#remove objects
rm(b1)
#
###################################################################
#DAY 2
#atomic and list vector
#dataframe
#matrix
#vector functions
#desciptive statistics
#normality checking of data
#set:creation, union and intersection
#import data from other sources in R
###################################################################
c<-"black"
c
str(c)
typeof(c)
length(c)

m<-"5000"
m
str(m)
typeof(m)
length(m)

## vector ## atomic and list ###
vec.at<-c(1,2,3,6,7,8,9,13)
class(vec.at)
str(vec.at)
typeof(vec.at)
dim(vec.at)

vec.li<-c("Ram", "Shyam", c(1:3))
class(vec.li)
str(vec.li)
typeof(vec.li)


attributes(vec.li)
is.list(vec.li)
as.list(vec.li)

# make list
vec.li<-list("Ram", "Shyam", c(1:3))
class(vec.li)
str(vec.li)
typeof(vec.li)

#now became list
is.list(vec.li)


## dataframe
# creat a dataframe with id, x and y

dat <- data.frame(id = c("plot1", "plot2", "plot3", "plot4", "plot5"),
                  x = 11:15,
                  y = 21:25)
dat
dim(dat)
head(dat)
class(dat)
str(dat)
typeof(dat)

## matrix
# creat a matrix 
mat1<-matrix(1:16, nrow = 4, byrow = T)
head(mat1)
attributes(mat1)

mat2<-matrix(1:100, nrow = 10, byrow = F)
head(mat2)
attributes(mat2)


############################################
##vector functions
x<-rnorm(20,5,1.5)
x
mean(x)

# set seed
# fixes the algorithm, if a same number of seed is set from all devices,
# result will be the same for all
# set.seed(seed, kind = NULL, normal.kind = NULL, sample.kind = NULL)
set.seed(1234)
x<-rnorm(20,5,1.5)
x
mean(x)


## Descriptive statistics
# create x and y
x<-rnorm(20,5,1.5)
x
mean(x)

y<-rnorm(20,5,0.5)
y
mean(y)

#descriptive statistics
#max, min, mean, range, median, mode, etc.
max(x)
max(y)
summary(x)
summary(y)
sum(x)
sum(y)
sd(x)
sd(y)

library(prettyR)
cor(x,y)
cor(y,x)

################################################################3
#check normality of data

#normality test 1
#histogram
#quick function 
#use histogram
#####################################################
hist(x)
hist(y)

#par(mfcol = c(1,1))

z<-rnorm(5000, 250, 5.0)
z
hist(z)

z1<-cbind(x,y)
plot(z1)

# normality test 2
# density plot test
#####################################################

plot(density(x))

hist(log(x))
hist(log(y))
hist(log(z))
hist(log(z1))

hist(sqrt(x))
plot(density(sqrt(x)))

## normality test 3
## shapiro.test(x)    
# to test Shapiro-Wilk normality test
######################################################
shapiro.test(x)

##################################################################
## normality test 4
## bartlett.test    
bartlett.test(x)
#did not response
#get help
?bartlett.test
#bartlett test is used for individual number and class category
# two class data
# e.g category (colour) and number  (count of groups)
plot(count ~ spray, data = InsectSprays)
bartlett.test(InsectSprays$count, InsectSprays$spray)
bartlett.test(count ~ spray, data = InsectSprays)

#view data from R help example
InsectSprays$count
InsectSprays$spray

### end of normality test
######################################################################
x<-rnorm(20,5,1.5)
y<-rnorm(20,5,0.5)

length(x)
length(y)
z<-cbind(x,y)
# z<-rbind(x,y)
class(z)

z<-data.frame(x,y)
class(z)


colMeans(z)
colSums(z)
rowMeans(z)
rowSums(z)
plot(z)
plot(z, col="Red")
plot(z, col="Red", line=1)


# set
#create, union and intersection

setA<-c("A", "B", "C")
setB<-c("A", "M", "N")
union(setA, setB)
intersect(setA, setB)

#############################################################
##############################################################
#Import data in R
#until now, we created the data ourself in R
# very often we mioght need to import data from excel, text, csv, and or web
#clipboard<-temporary storage in a computer

#First copy the data from excel sheet, only data

dat<-read.table("clipboard", header=T)
names(dat)
head(dat)
class(dat)
str(dat)
summary(dat)  #gives summary of all variables of the dataframe
summary(dat$DBHmax)  #gives summary of only variable of interest
#check normality of DBHmax

hist(dat$DBHmax)
shapiro.test(dat$DBHmax)

####################################################################
# next way to import data
# importing the csv file (data) from file location
dat1<-read.csv("C://Users//hp//Desktop//R_Training//Data_share//Data.csv")
names(dat1)
head(dat1)
summary(dat1)

mod1<-lm(dat1$DBHmax~dat1$Heightmax)

summary(mod1)
plot(dat1$DBHmax, dat1$Heightmax, xlab="Tree DBH max", ylab="Tree Height max", 
     main="Scatter plot DBH and Height")

abline(mod1, lwd=3, col="red")
########################################################################
########################################################################
#DAY III
#Data exploration
#Hypothesis test
# import data.csv file
# merge
data<-read.csv("C://Users//hp//Desktop//R_Training//Data_share//Data.csv")
names(data)
head(data)
tail(data)
str(data)
summary(data)
plot(data)

plot(data, panel=panel.smooth)

##############################################################
# one sample t-test
# assumption:data come from normal distribution
# assumption:data from independent random variables

hist(data$TRmean.dbh)
shapiro.test(data$TRmean.dbh)
t.test(data$TRmean.dbh)

hist(data$DBHmax)
shapiro.test(data$DBHmax)
t.test(data$DBHmax)

hist(log(data$TRmean.dbh))      # log transformation to check normality
shapiro.test(log(data$TRmean.dbh))   # p value<0.05, then data is not normal
t.test(data$TRmean.dbh)

hist(log(data$Trmean.ht))
shapiro.test(log(data$Trmean.ht))
t.test(data$Trmean.ht)      # two-tailed t-test

t.test(data$Trmean.ht, alt="g", conf.level = 0.99)      # one-tailed t-test

t.test(data$Trmean.ht, alt="l", conf.level = 0.99)      # one-tailed t-test

boxplot (data$Trmean.ht)
boxplot (log(data$Trmean.ht))

boxplot(data$TRmean.dbh)
boxplot (log(data$TRmean.dbh))
##################################################
# NON-PARAMETRIC T-TEST
# wilcoxon signed-rank test
# smaller sample size (but >6)
# for non-normal data

wilcox.test(data$TRmean.dbh)  # mu (pop mean if available only)
wilcox.test(data$TRmean.dbh, mu=30.39146)  # mu (pop mean if available only)

#####################################################
#two-sample t-test
#####################################################

data.icimod<-read.table("clipboard", header=T)
names(data.icimod)
head(data.icimod)
str(data.icimod)

boxplot(data.icimod$AGTB.t.per.ha~data.icimod$Strata)
# OR
plot(data.icimod$AGTB.t.per.ha~data.icimod$Strata)
#OR
plot(data.icimod$Strata, data.icimod$AGTB.t.per.ha)

#bartlett test checks the homogeinity of variances
bartlett.test(data.icimod$AGTB.t.per.ha~data.icimod$Strata)

# Fligner-Killeen test of homogeneity of variances
fligner.test(data.icimod$AGTB.t.per.ha~data.icimod$Strata)


shapiro.test(data.icimod$AGTB.t.per.ha[data.icimod$Year=="2012"])
shapiro.test(data.icimod$AGTB.t.per.ha[data.icimod$Year=="2014"])

par(mfcol = c(1,2))

hist(data.icimod$AGTB.t.per.ha[data.icimod$Year=="2012"])
hist(data.icimod$AGTB.t.per.ha[data.icimod$Year=="2014"])

m12<-mean(data.icimod$AGTB.t.per.ha[data.icimod$Year=="2012"])
m12
m14<-mean(data.icimod$AGTB.t.per.ha[data.icimod$Year=="2014"])
m14



#performing t-test 
# Welch Two Sample t-test
# for independent data
t.test((data.icimod$AGTB.t.per.ha[data.icimod$Year=="2012"]),
      (data.icimod$AGTB.t.per.ha[data.icimod$Year=="2014"]))

t.test((data.icimod$AGTB.t.per.ha[data.icimod$Year=="2012"]),
       (data.icimod$AGTB.t.per.ha[data.icimod$Year=="2014"]),
       var.equal = TRUE)

# above one is for independent data
# thus for dependent data like this, we need to do paired t-test
#check variance of data first,
var(data.icimod$AGTB.t.per.ha[data.icimod$Year=="2012"])
var(data.icimod$AGTB.t.per.ha[data.icimod$Year=="2014"])

#paired t-test
# formula: t.test(x,y,paired=T)

t.test((data.icimod$AGTB.t.per.ha[data.icimod$Year=="2012"]),
       (data.icimod$AGTB.t.per.ha[data.icimod$Year=="2014"]), 
       paired = T)


t.test(data.icimod$AGTB.t.per.ha~data.icimod$Strata, 
       paired=T) # for independent data

#merge two dataframe
# merge two data frames by ID
merge.data <- merge(data,dat1,by="Code")
#####################################################################

### Day Iv

#####################################################################
# Binomial test to compare two proportion
## Formula: prop.test (c(N.male, N.female),c(N.disease1, N.disease2))
# e.g
prop.test(c(100,110), c(250, 500))

# chisq test for association and goodness of fit
# chisq test is not for cause and effect relationship

# import data "bird_behaviour from day_III_data
bird<-read.table("clipboard", header=T)
str (bird)
class(bird)

#create a table of "bird"

table (bird)
ftable(bird)

table(bird$Activity, bird$Time) # makes a table with Time in rows

#plots
plot(table)
plot(bird$Activity, bird$Time)

plot(table(bird$Activity, bird$Time))  #general plot

cdplot (bird$Activity, bird$Time)  #continuous density plot
cdplot (bird$Time, bird$Activity)  #continuous density plot


# plot(bird$Time, bird$Activity)
#par(mfcol=c(1,1))

#chi square test
chisq.test (bird$Activity, bird$Time)
chisq.test (table(bird$Activity, bird$Time))

# chisq test for goodness of fit can be used in this way
kisql<-c(24,25,27) # data of diseased bat in three forest  
chisq.test(kisql)
# simulated chisq test


#chisq.test, simulate()
# Fisher's exact test can be used if the table is 2 X 2
# Fisher's Exact Test for Count Data
fisher.test(table(bird$Activity, bird$Time))


#########################################################
str(bird) # "bird" is a dataframe
t1<-table(bird)
t1  # t1 is a table from "bird"
str(t1)
t2<-as.data.frame(t1) # converted to dataframe again
t2
str(t2)
plot(t2)

#######################################################
#create a new matrix
D<-matrix(c(23, 7, 18, 12), nrow = 2, byrow = T)
D
# row names and col names
rownames (D)<- c("Chitwan", "Bardiya")
colnames (D)<- c("Road.accident", "Poaching")
D

##################################################
#Correlation and linear regression
#correlation
#formula
# for parametric test
# data having normal distribution
# cor(x,y)
# cor.test(x,y)
cor(data$TRmean.dbh, data$Trmean.ht)
cor.test(data$TRmean.dbh, data$Trmean.ht)
cor.test(data$TRmean.dbh, data$Trmean.ht, simulate.p.value=T,B=2000)
plot(data$TRmean.dbh, data$Trmean.ht)

# Rank correlation coefficient
# for non-parametric data
# data having non normal distribution
# cor.test(x,y, method="spearman")
cor.test(data$TRmean.dbh, data$Trmean.ht, method="spearman")


# Kendall's rank correlation tau
# when there is a tie (e.g. 2 or more items with the highest score)
# Kendall correlation coeff.- to deal with data samples with tied ranks
# (when spearman suffers!)
# use kendall test
# Formula= cor.test(x, y, method="kendall")
cor.test(data$TRmean.dbh, data$Trmean.ht, method="kendall")

#kendall test gives lower value corr. coeff. but is a robust method
# for smaller sample size like 10-15, use kendall test
# for tie case (spearman gets suffered), use kendall test

str(data)
    

# # when data is missing, use (na.action="na.omit") or (complete.obs)


# Regression, linear
# prediction of a dependent variable (response) with the help of independent var.
# e.g y=per year C stock (dependent); x=DBH
# non of the model is correct or incorrect,
# models are useful, thus we use them
# observation of data

plot(data$TRmean.dbh, data$Trmean.ht)
# OR
plot(data$Trmean.ht~data$TRmean.dbh)
m0<-lm(data$Trmean.ht~data$TRmean.dbh)
m0 #view model
summary(m0)
abline(m0)
abline(m0, xlab="Mean DBH", ylab="Mean Height",
       main="DBH-Height model of Urban Trees in Kathmandu",
       col="red", lwd=2.5,)

#model is used but is it useful?
# Diagnostics:

par(mfrow=c(2,2))
plot(m0)

# Above model diagnostic shows
# This (linear) model does not seem good for this modeling
# Lets go for polynomial regression (still linear)
# residuals
# m0$residuals

res<-resid(m0)
res
absres<-abs(res)
levene<-lm(abres~data$TRmean.dbh)
anova(levene)    # check p value (p<0.05?)

shapiro.test(res)
# OR
shapiro.test(resid(m0))

hist(m0$residuals)

###############################################
# next model
m1<-lm(data$Trmean.ht~data$TRmean.dbh+I(data$TRmean.dbh^2))
summary(m1)
plot(m1)
plot(data$Trmean.ht~data$TRmean.dbh+I(data$TRmean.dbh^2))

# higher level of polynomial model
###############################################
m2<-lm(data$Trmean.ht~data$TRmean.dbh+I(data$TRmean.dbh^2)+I(data$TRmean.dbh^3))
summary(m2)
plot(m2)
plot((data$Trmean.ht~data$TRmean.dbh+I(data$TRmean.dbh^2)+I(data$TRmean.dbh^3)))

###################################################
# higher level of polynomial model
###############################################
m3<-lm(data$Trmean.ht~data$TRmean.dbh+I(data$TRmean.dbh^2)+
         I(data$TRmean.dbh^3)+I(data$TRmean.dbh^4))
summary(m3)
plot(m3)
plot(data$Trmean.ht~data$TRmean.dbh+I(data$TRmean.dbh^2)+
       I(data$TRmean.dbh^3)+I(data$TRmean.dbh^4))

###################################################
# anova can be used to compare the models
anova(m0,m1,m2,m3)

# decision: none of the models can better perform
# use of DBH only can cover hardly 54-55% of the prediction
# Thus extra parameters are to be used in addition to DBH to obtain better model
# Multiple regression model

names(data)
m4<-lm(data$Trmean.ht~data$TRmean.dbh+data$Canopy.cover+
         data$Ground.cover+data$TRichness++data$Trabundance)
summary(m4)
par(mfrow=c(2,2))
plot(m4)
# only dbh and Ground.cover (-ve) are significant

# try next model with only significant variables

m5<-lm(data$Trmean.ht~data$TRmean.dbh+data$Ground.cover)
summary(m5)
par(mfrow=c(2,2))
plot(m5)

#########################################
# machine given model using step
m5<-step(m4)
summary(m5)

# copied machine built model
m6<-lm(data$Trmean.ht ~ data$TRmean.dbh + data$Ground.cover)
summary(m6)
par(mfrow=c(2,2))
plot(m6)
coefficients (m6)

################################################
#complex model
m7<-lm(data$Trmean.ht~data$TRmean.dbh*data$Canopy.cover*
         data$Ground.cover*data$TRichness*data$Trabundance)
summary(m7)
par(mfrow=c(2,2))
plot(m7)
#step
m8<-step(m7)
summary(m8)


############################################################################
# Day V
###############################################################
# Anova (F-test) (ratio-test) compares means of three or more samples
# Anova: single response (dependent) variable e.g. height of invasive species
#         dependent variable is Numeric; Independent variable is 
#         categorical e.g. site class
#         unlike Regression where independent var is also numeric
# Two-way Anova: If there are two categorical independent variables
# ANACOVA: If there are two independent variables; 1 categorical & 1 numeric
######################################################################
# Manova
# Manova: multiple response (dependent) variables e.g. height & density
#of invasive species
#############################################################################
# aov (y~x)
?anova
?manova
###########################################
data.an<-read.table("clipboard", header = T)
head(data.an)
str(data.an)
plot(data.an)
data.an

# first check the data
tapply(data.an$AGTB.t.per.ha, data.an$Grazing.level, var)  
# gives variance of AGTB at all grazing levels

bartlett.test(data.an$AGTB.t.per.ha, data.an$Grazing.level, var)
fligner.test(data.an$AGTB.t.per.ha, data.an$Grazing.level, var)

boxplot(data.an$AGTB.t.per.ha ~ data.an$Grazing.level)





# now go for one-way anova
m0<-aov(data.an$AGTB.t.per.ha ~ data.an$Grazing.level)
summary (m0)
# Diagnostics
par(mfrow=c(2,2))
plot(m0)
par(mfrow=c(2,3))
plot(m0,which=(1:6))

summary.lm(m0) # gives detailed summary of anova "m0"

# pre-hoc test : before field data
# post-hoc test: after data 
#post-hoc: TukeyHSD(m0)
TukeyHSD(m0)
plot (TukeyHSD(m0))
# seems due to difference of Low-High density plots
pairwise.t.test(data.an$AGTB.t.per.ha, data.an$Grazing.level) # alternative of Tukey
# but is not as good as TukeyHSD

# Conclusion: data not good for ANOVA (anova suffers) not normality in data

m1<-aov(data.an$AGTB.t.per.ha ~ data.an$District)
summary (m1)

m2<-aov(data.an$AGTB.t.per.ha ~ data.an$Strata)
summary (m2)

m3<-aov(data.an$AGTB.t.per.ha~data.an$Grazing.level+data.an$District)
summary (m3)
summary.lm (m3)

plot(m3)
par(mfrow=c(2,3))
plot(m3,which=(1:6))
TukeyHSD(m3)

m5<-aov(data.an$AGTB.t.per.ha~data.an$Grazing.level*data.an$District)
summary (m5)
summary.lm (m5)

m4<-step(m3)

m6<-step(m5)

# last line of code is formula for the best model (minimal adequate model)
M.final<-aov(data.an$AGTB.t.per.ha ~ data.an$Grazing.level)
plot(M.final)


# Option: Kruskal-Wallis (non-parametric) test
kruskal.test((data.an$AGTB.t.per.ha ~ data.an$Grazing.level))

# result: not significant (Grazing has no role to predict Biomass)

#################################################################
# ANCOVA
data.anc<-read.table("clipboard", header = T)
# grazing=predictor (categorical)
# root=predictor (numeric/continuous)
# fruit=
#ANCOVA=Reg+anova
# get your data ready, explore
par(mfcol=c(1,1))
plot(data.anc$Root,data.anc$Fruit, pch=56, col=c("blue", "red")
     [as.numeric(data.anc$Grazing)])
levels(data.anc$Grazing)

# abline (lm(Fruit~Root))

abline(lm(data.anc$Fruit[data.anc$Grazing=="Grazed"]~data.anc$Root
          [data.anc$Grazing=="Grazed"]), col="blue")

abline(lm(data.anc$Fruit[data.anc$Grazing=="Ungrazed"]~data.anc$Root
          [data.anc$Grazing=="Ungrazed"]), col="red")

tapply(data.anc$Fruit, data.anc$Grazing, mean)

#legend(locator(2), legend=c("Grazed", "Ungrazed"), pch=c(20,20),
 #      col=c("blue", "red"), title = "Grazing Effects on Roots and Fruits")


##########################
# Day VI
###########
#########
data.anc<-read.table("clipboard", header = T)

#OR 
file.choose()
data.ancova<- read.table("C:\\Users\\hp\\Desktop\\R_Training\\anco.txt", 
                         header = T)
str(data.ancova)
head(data.ancova)
summary(data.ancova)
tapply(data.ancova$Fruit, data.ancova$Grazing, mean)

x<-data.ancova
###############################################################
plot(x)
plot(x$Root,x$Fruit, pch=56, col=c("blue", "red")
     [as.numeric(x$Grazing)])
levels(x$Grazing)

# abline (lm(Fruit~Root))

abline(lm(x$Fruit[x$Grazing=="Grazed"]~x$Root
          [x$Grazing=="Grazed"]), col="blue")

abline(lm(x$Fruit[x$Grazing=="Ungrazed"]~x$Root
          [x$Grazing=="Ungrazed"]), col="red")

tapply(x$Fruit, x$Grazing, mean)

#legend(locator(2), legend=c("Grazed", "Ungrazed"), pch=c(20,20),
 #     col=c("blue", "red"), title = "Grazing Effects on Roots and Fruits")

# check normality
# script not final!!!!
hist(x$Fruit, x$Grazing)

shapiro.test(x$Fruit, x$Grazing)

bartlett.test(x$Fruit, x$Root)
fligner.test(x$Fruit, x$Grazing)

boxplot(data.an$AGTB.t.per.ha ~ data.an$Grazing.level)

M1<-lm(x$Fruit~x$Root+I(x$Grazing))
# script not final
#upto here!!!!!
#############################################################
ancova<-lm(x$Fruit~x$Grazing*x$Root)
summary(ancova)
anova(ancova)

#####
ancova2<-update(ancova, ~.-x$Grazing:x$Root)
anova(ancova, ancova2)

ancova3<-update(ancova2, ~.-x$Grazing)
anova(ancova2, ancova3)  # 3is a bad model
summary(ancova2)
anova(ancova2)
step(ancova)   # you can opt but be careful
