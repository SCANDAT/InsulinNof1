library(rio)
library(ggplot2)
library(lme4)
library(nlme)
library(lmerTest)

#This works with R studio, use alternative approach or do manually if working with other IDE
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Import SAS data
dataset <- import('../Data/studydataset.sas7bdat')
head(dataset)

#First investigate autocorrelation
typical_block_data <- subset(dataset, Block == 1)
acf(typical_block_data$value)
pacf(typical_block_data$value)

acf(dataset$value)
pacf(dataset$value)

#Run AR(1) model
model <- lme(value ~ insulintype, 
             random = ~ 1 | person, 
             correlation = corAR1(form = ~ local_time | person), 
             data = dataset)
summary(model)

#construct confidence intervals
intervals(model, which="fixed")

#Store for later use 
save(model,file="LMEAR1 model.rds")


#Test also AR(2) model (no big difference)
model2 <- lme(value ~ insulintype, 
             random = ~ 1 | person, 
             correlation = corARMA(p=2, form= ~ local_time | person), 
             data = dataset)
summary(model2)

#construct confidence intervals
intervals(model2, which="fixed")

#Store for later use 
save(model2,file="LMEAR2 model.rds")
