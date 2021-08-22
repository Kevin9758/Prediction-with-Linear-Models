


ctrain <- read.csv ("protein-train.csv")
pfull <- lm(accuracy ~ .,data = ctrain)
par ( mfrow = c (2 ,2) )
plot(pfull$fitted.values, pfull$residuals, xlab = "Fitted Values", 
     ylab = "Residuals")

plot(1: nrow (ctrain), pfull$residuals , xlab = " Index ", ylab =
       " Residuals ")

qqnorm(pfull$residuals)
qqline(pfull$residuals, col="blue", lwd = 2)
hist(pfull$residuals)
fullrsq <- summary(pfull)$r.squared




# Methods






# Cook's distance
load("premov")

cd <- cooks.distance(premov) # Cook's distance

plot(cd, ylab = "Cook's distance")
abline(h = 0.5, col = "red", lty = 2)
h <- unname(which(cd > 0.5))
# Cook's distance greater than 0.5 is generally considered large
# Thus, there are no outliers 
# Other rules of thumb are 4/n, and 4/n-k-1

# otrain <- ptrain[-h,]
# new <- lm(accuracy ~ .,data = ptrain)
# summary(new)





load("ptrain")
set.seed(20688307)
N <- nrow(ptrain) 
trainInd <- sample(1: N ,round(N*0.8) ,replace = F) 
# 80/20 train/validation split
trainSet <- ptrain[trainInd ,]
validSet <- ptrain[-trainInd ,]

load("m_aicf")
load("m_bicf")
load("m_aicb")
load("m_bicb")
load("m_aich")
load("m_bich")
load("m_bestA")
load("m_bestB")
load("m_bestC")




load("full")

predfull <- predict(full, newdata = validSet) 
fullv <- sqrt(mean((validSet$accuracy - predfull)^2)) #RMSE for validation
fullt <- sqrt(mean(full$residuals^2)) # RMSE for train



# AIC-forward selection

predaf <- predict(m_aicf, newdata = validSet)
aicfv <- sqrt(mean((validSet$accuracy - predaf)^2)) 
aicft <- sqrt(mean(m_aicf$residuals ^2)) 

faic <- AIC(m_aicf)

# AIC-backward selection

predab <- predict(m_aicb, newdata = validSet) 
aicbv <- sqrt(mean((validSet$accuracy - predab)^2)) 
aicbt <- sqrt(mean(m_aicb$residuals^2))

baic <- AIC(m_aicb)

# AIC-forward/backward selection

predah <- predict(m_aich, newdata = validSet)
aichv <- sqrt(mean((validSet$accuracy - predah )^2))
aicht <- sqrt(mean(m_aich$residuals ^2))

haic <- AIC(m_aich)

# BIC-forward selection

predbf <- predict (m_bicf , newdata = validSet)
bicfv <- sqrt(mean((validSet$accuracy - predbf )^2)) 
bicft <- sqrt(mean(m_bicf$residuals^2)) 

fbic <- BIC(m_bicf)

# BIC-backward selection

predbb <- predict(m_bicb, newdata = validSet) 
bicbv <- sqrt(mean((validSet$accuracy - predbb)^2)) 
bicbt <- sqrt(mean(m_bicb$residuals^2)) 

bbic <- BIC(m_bicb)


# BIC-forward/backward selection

predbh <- predict(m_bich, newdata = validSet) 
bichv <- sqrt(mean((validSet$accuracy - predbh)^2))
bicht <- sqrt(mean(m_bich$residuals ^2))

hbic <- BIC(m_bich)



# ICM RMSE - we get varying RMSE's since different orderings for ICM
# Sometimes we get a higher RMSE, sometimes lower. So we can cross validate to
# check further.

# ICM - no penalty

bestA <- predict(bestA_model, newdata = validSet) 
bestAv <- sqrt(mean((validSet$accuracy - bestA)^2))
bestAt <- sqrt(mean(bestA_model$residuals ^2))

bestaic <- AIC(bestA_model)


# ICM -BIC penalty

bestB <- predict(bestB_model, newdata = validSet) 
bestBv <- sqrt(mean((validSet$accuracy - bestB)^2))
bestBt <- sqrt(mean(bestB_model$residuals ^2))

bestbic <- BIC(bestB_model)

# ICM -2 * BIC penalty

bestC <- predict(bestC_model, newdata = validSet) 
bestCv <- sqrt(mean((validSet$accuracy - bestC)^2))
bestCt <- sqrt(mean(bestC_model$residuals ^2))







load("RMSE1")
load("RMSE2")
load("RMSE3")
load("RMSE4")

RMSPE1 <- mean(RMSE1)
RMSPE2 <- mean(RMSE2)
RMSPE3 <- mean(RMSE3)
RMSPE4 <- mean(RMSE4)



load("mfinal")
finalrsq <- summary(mfinal)$r.squared



# We check the model for linear regression assumptions again


par ( mfrow = c (2 ,2) )

 plot(mfinal$fitted.values, mfinal$residuals, xlab = "Fitted Values", 
      ylab = "Residuals")
 plot (1: nrow (ptrain), mfinal$residuals , xlab = " Index ", ylab =
         " Residuals ")
 
 qqnorm(mfinal$residuals)
 qqline(mfinal$residuals, col="blue", lwd = 2)
 hist(mfinal$residuals)

MSPE <- RMSPE1^2
 








# Detect and remove multicollinearity / commented since I saved after running

library("faraway")

ptrain <- read.csv ("protein-train.csv")

i <- 0
pred <- c()

# full model
pfull <- lm(accuracy ~ .,data = ptrain)

# vector of VIF's of predictors
# npred <- vif(pfull)
# vifpred <- unname(npred)

# copy of full model
ctrain <- ptrain[]

load("ptrain")

# # while the largest VIF >= 10
# while(max(vifpred) >= 10)
# {
#   
#   # find the index of the max VIF 
#   j <- which.max(vifpred) + 1 
#   ptrain <- ptrain[,-j] 
#   # vifpred only includes the predictors (685) while ptrain has the
#   # response as well (686), so we add 1 to j
#   
#   
#   premov <- lm(accuracy ~ .,data = ptrain) # fit new model without predictor 
#   # with largest VIF
#   
#   npred <- vif(premov)
#   vifpred <- unname(npred) # find new vector of VIF's and loop until
#   # every VIF is less than 10
#   
# 
#   i <- i + 1 # count of how many predictors removed
#   
# }
load("premov")




# after removing multicollinearity, we assess model linear regression 
# assumptions

#summary(premov)
plot(premov$fitted.values, premov$residuals, xlab = "Fitted Values", 
     ylab = "Residuals")
plot (1: nrow (ptrain), premov$residuals , xlab = " Index ", ylab =
        " Residuals ")

qqnorm(premov$residuals)
qqline(premov$residuals, col="blue", lwd = 2)
hist(premov$residuals)

library(MASS)

# quantities for individual observations
sres <- studres(premov) # studentized residuals
lev <- hatvalues(premov) # leverage
cd <- cooks.distance(premov) # Cook's distance


# plot of studentized residuals vs fitted values
plot(premov$fitted.values, sres, xlab = "Fitted values", 
     ylab = "Studentized residuals")
# plots look the same
abline(h = c(3,-3), col = "red", lty=2) # lie within 3 standard deviations
# of 0
which(abs(sres) > 3) # find the observations that don't lie within 3
# standard deviations
#  23  226  527  846 1196 1324 1344 1686 1733
qqnorm(sres)
qqline(sres, col = "blue", lwd = 2)


#leverage

plot(lev, ylab = "Leverage")
abline(h = 2 * mean(lev), col = "red", lty = 2) # leverages twice than average
# leverage are considered high
which(lev > 2 * mean(lev)) # 1635


# Cook's distance

plot(cd, ylab = "Cook's distance")
abline(h = 0.5, col = "red", lty = 2)
h <- unname(which(cd > 0.5))
# Cook's distance greater than 0.5 is generally considered large
# Thus, there are no outliers 
# Other rules of thumb are 4/n, and 4/n-k-1

otrain <- ptrain[-h,]
new <- lm(accuracy ~ .,data = ptrain)
# summary(new)





# Check model regression assumptions after outliers removed

# plot(new$fitted.values, new$residuals, xlab = "Fitted Values", 
#      ylab = "Residuals")
# qqnorm(new$residuals)
# qqline(new$residuals, col="blue", lwd = 2)
# hist(new$residuals)



# First we test out different model selection methods with a single 
# train/validation split.

library(MASS)

N <- nrow(ptrain)
set.seed(20688307)

trainInd <- sample(1: N ,round(N*0.8) ,replace = F) 
# 80/20 train/validation split
trainSet <- ptrain[trainInd ,]
validSet <- ptrain[-trainInd ,]

full <- lm(accuracy ~ ., data = trainSet)
empty <- lm(accuracy ~ 1, data = trainSet)

# stepwise forward with aic

load("m_aicf")

#m_aicf <- stepAIC(object = empty, scope = list(upper = full, lower = empty), direction = #"forward")



# stepwise backward with aic

load("m_aicb")

#m_aicb <- stepAIC(object = full, scope = list(upper = full, lower = empty), direction = #"backward")



# stepwise forward/backward with aic

load("m_aich")

#m_aich <- stepAIC(object = empty, scope = list(upper = full, lower = empty), direction
#= "both")



# stepwise forward with bic

load("m_bicf")

#m_bicf <- stepAIC(object = empty, scope = list(upper = full, lower = empty), direction = #"forward", k =log(nrow(trainSet)))



# stepwise backward with bic

load("m_bicb")

#m_bicb <- stepAIC(object = full, scope = list(upper = full, lower = empty), direction = "backward", #k=log(nrow(trainSet)))



#stepwise forward-backward with bic

load("m_bich")

#m_bich <- stepAIC(object = empty, scope = list(upper = full, lower = empty), direction = "both", k #=log(nrow(trainSet)))



# We calculate the RMSE on train and validation for each model


# Full model (after removing multicollinearity)

predfull <- predict(full, newdata = validSet) 
fullv <- sqrt(mean((validSet$accuracy - predfull)^2)) #RMSE for validation
fullt <- sqrt(mean(full$residuals^2)) # RMSE for train



# AIC-forward selection

predaf <- predict(m_aicf, newdata = validSet)
aicfv <- sqrt(mean((validSet$accuracy - predaf)^2)) 
aicft <- sqrt(mean(m_aicf$residuals ^2)) 

faic <- AIC(m_aicf)

# AIC-backward selection

predab <- predict(m_aicb, newdata = validSet) 
aicbv <- sqrt(mean((validSet$accuracy - predab)^2)) 
aicbt <- sqrt(mean(m_aicb$residuals^2))

baic <- AIC(m_aicb)

# AIC-forward/backward selection

predah <- predict(m_aich, newdata = validSet)
aichv <- sqrt(mean((validSet$accuracy - predah )^2))
aicht <- sqrt(mean(m_aich$residuals ^2))

haic <- AIC(m_aich)

# BIC-forward selection

predbf <- predict (m_bicf , newdata = validSet)
bicfv <- sqrt(mean((validSet$accuracy - predbf )^2)) 
bicft <- sqrt(mean(m_bicf$residuals^2)) 

fbic <- BIC(m_bicf)

# BIC-backward selection

predbb <- predict(m_bicb, newdata = validSet) 
bicbv <- sqrt(mean((validSet$accuracy - predbb)^2)) 
bicbt <- sqrt(mean(m_bicb$residuals^2)) 

bbic <- BIC(m_bicb)


# BIC-forward/backward selection

predbh <- predict(m_bich, newdata = validSet) 
bichv <- sqrt(mean((validSet$accuracy - predbh)^2))
bicht <- sqrt(mean(m_bich$residuals ^2))

hbic <- BIC(m_bich)



#ICM with no penalty (AIC)
library(MASS)

N <- nrow(ptrain)

load("m_bestA")

#set.seed(20688307)

# trainInd <- sample(1: N ,round(N*0.8) ,replace = F) 
# # 80/20 train/validation split
# trainSet <- ptrain[trainInd ,]
# validSet <- ptrain[-trainInd ,]
# 
# full <- lm(accuracy ~ ., data = trainSet)
# empty <- lm(accuracy ~ 1, data = trainSet)
# 
# varlist = c()
# varnames = names(trainSet)
# n = nrow(trainSet)
# varoder <- sample(1:ncol(trainSet)) #random order of variables
# minCrit = Inf
# noChange = F
# 
# while(!noChange) {
#   noChange = T
#   for (i in varoder) {
#     if (i == 1)
#       next
#     
#     if (i %in% varlist & length(varlist) > 1) {
#       index = c(1, varlist[varlist != i])
#     trainVars = trainSet[, index]
#     
#       fit = lm(accuracy ~ ., data = trainVars)
#       
#       if (AIC(fit) < minCrit) {
#           minCrit = AIC (fit)
#           varlist = varlist[varlist != i]
#           #print(paste0(" Criterion : ", round(minCrit, 1) , ", variables : ",
#             #           paste0 (varnames[varlist], collapse = " ")))
#           bestA_model = fit
#           noChange = F
#       }
#       
#     }  else if  (!i %in% varlist) {
#       index = c(1, varlist, i)
#       trainVars = trainSet[, index]
#     
#       fit = lm(accuracy ~ ., data = trainVars)
#     
#       if (AIC(fit) < minCrit) {
#           minCrit = AIC(fit)
#           varlist = c(varlist, i)
#          # print(paste0(" Criterion : ", round(minCrit, 1), ", variables : ",
#                    #  paste0(varnames[varlist], collapse = " ")))
#   
#           bestA_model = fit
#           noChange = F
#        }
#     } 
#   }
# }


# ICM with BIC penalty
library(MASS)

N <- nrow(ptrain)
#set.seed(20688307) different AIC each time without seed

load("m_bestB")

# trainInd <- sample(1: N ,round(N*0.8) ,replace = F) 
# # 80/20 train/validation split
# trainSet <- ptrain[trainInd ,]
# validSet <- ptrain[-trainInd ,]
# 
# full <- lm(accuracy ~ ., data = trainSet)
# empty <- lm(accuracy ~ 1, data = trainSet)
# 
# pen <- log (nrow (trainSet))
# varlist = c()
# varnames = names(trainSet)
# n = nrow(trainSet)
# varoder <- sample(1:ncol(trainSet)) #random order of variables
# minCrit = Inf
# noChange = F
# 
# while(!noChange) {
#   noChange = T
#   for (i in varoder) {
#     if (i == 1)
#       next
#     
#     if (i %in% varlist & length(varlist) > 1) {
#       index = c(1, varlist[varlist != i])
#     trainVars = trainSet[, index]
#     
#       fit = lm(accuracy ~ ., data = trainVars)
#       
#       if (AIC(fit, k=pen) < minCrit) {
#           minCrit = AIC (fit, k = pen)
#           varlist = varlist[varlist!= i]
#         #  print(paste0(" Criterion : ", round(minCrit, 1) , ", variables : ",
#         #               paste0 (varnames[varlist], collapse = " ")))
#           bestB_model = fit
#           noChange = F
#       }
#       
#     }  else if  (!i %in% varlist) {
#       index = c(1, varlist, i)
#       trainVars = trainSet[, index]
#     
#       fit = lm(accuracy ~ ., data = trainVars)
#     
#       if (AIC(fit, k = pen) < minCrit) {
#           minCrit = AIC(fit, k = pen)
#           varlist = c(varlist, i)
#         #  print(paste0(" Criterion : ", round(minCrit, 1), ", variables : ",
#                 #     paste0(varnames[varlist], collapse = " ")))
#   
#           bestB_model = fit
#           noChange = F
#        }
#     } 
#   }
# }




# ICM with 2*BIC penalty
library(MASS)

N <- nrow(ptrain)
#set.seed(20688307) different AIC each time without seed

load("m_bestC")

# trainInd <- sample(1: N ,round(N*0.8) ,replace = F)
# # 80/20 train/validation split
# trainSet <- ptrain[trainInd ,]
# validSet <- ptrain[-trainInd ,]
# 
# full <- lm(accuracy ~ ., data = trainSet)
# empty <- lm(accuracy ~ 1, data = trainSet)
# 
# pen <- 2 * log (nrow (trainSet))
# varlist = c()
# varnames = names(trainSet)
# n = nrow(trainSet)
# varoder <- sample(1:ncol(trainSet)) #random order of variables
# minCrit = Inf
# noChange = F
# 
# while(!noChange) {
#   noChange = T
#   for (i in varoder) {
#     if (i == 1)
#       next
# 
#     if (i %in% varlist & length(varlist) > 1) {
#       index = c(1, varlist[varlist != i])
#     trainVars = trainSet[, index]
# 
#       fit = lm(accuracy ~ ., data = trainVars)
# 
#       if (AIC(fit, k=pen) < minCrit) {
#           minCrit = AIC (fit, k = pen)
#           varlist = varlist[varlist!= i]
#         #  print(paste0(" Criterion : ", round(minCrit, 1) , ", variables : ",
#         #               paste0 (varnames[varlist], collapse = " ")))
#           bestC_model = fit
#           noChange = F
#       }
# 
#     }  else if  (!i %in% varlist) {
#       index = c(1, varlist, i)
#       trainVars = trainSet[, index]
# 
#       fit = lm(accuracy ~ ., data = trainVars)
# 
#       if (AIC(fit, k = pen) < minCrit) {
#           minCrit = AIC(fit, k = pen)
#           varlist = c(varlist, i)
#         #  print(paste0(" Criterion : ", round(minCrit, 1), ", variables : ",
#                 #     paste0(varnames[varlist], collapse = " ")))
# 
#           bestC_model = fit
#           noChange = F
#        }
#     }
#   }
# }




# ICM RMSE - we get varying RMSE's since different orderings for ICM
# Sometimes we get a higher RMSE, sometimes lower. So we can cross validate to
# check further.

# ICM - no penalty

bestA <- predict(bestA_model, newdata = validSet) 
bestAv <- sqrt(mean((validSet$accuracy - bestA)^2))
bestAt <- sqrt(mean(bestA_model$residuals ^2))

bestaic <- AIC(bestA_model)


# ICM -BIC penalty

bestB <- predict(bestB_model, newdata = validSet) 
bestBv <- sqrt(mean((validSet$accuracy - bestB)^2))
bestBt <- sqrt(mean(bestB_model$residuals ^2))

bestbic <- BIC(bestB_model)

# ICM -2 * BIC penalty

bestC <- predict(bestC_model, newdata = validSet) 
bestCv <- sqrt(mean((validSet$accuracy - bestC)^2))
bestCt <- sqrt(mean(bestC_model$residuals ^2))




load("RMSE1")
load("RMSE2")
load("RMSE3")
load("RMSE4")

RMSPE1 <- mean(RMSE1)
RMSPE2 <- mean(RMSE2)
RMSPE3 <- mean(RMSE3)
RMSPE4 <- mean(RMSE4)


# Cross Validation to determine the best selection method
# library (MASS)
# N <- nrow(ptrain)
# set.seed(20688307)
# 
# K<-5
# 
# validSetSplits <- sample((1:N)%%K + 1)
# RMSE1 <- c()
# RMSE2 <- c()
# RMSE3 <- c()
# RMSE4 <- c()
# 
# for (p in 1:K) {
#   validSet <- ptrain[validSetSplits == p ,]
#   trainSet <- ptrain[validSetSplits != p ,]
# 
#   full <- lm (accuracy ~ . , data = trainSet)
#   empty <- lm (accuracy ~ 1 , data = trainSet)
# 
# 
# m1 <- stepAIC(object = empty, scope = list(upper = full, lower = empty), direction = "forward", trace = #0)
# 
# pred1 <- predict( m1, newdata = validSet)
# RMSE1[p] <- sqrt(mean((validSet$accuracy - pred1)^2))
# 
# 
# 
# varlist = c()
# varnames = names(trainSet)
# n = nrow(trainSet)
# varoder <- sample(1:ncol(trainSet)) #random order of variables
# minCrit = Inf
# noChange = F
# 
# while(!noChange) {
#   noChange = T
#   for (i in varoder) {
#     if (i == 1)
#       next
# 
#     if (i %in% varlist & length(varlist) > 1) {
#       index = c(1, varlist[varlist != i])
#     trainVars = trainSet[, index]
# 
#       fit = lm(accuracy ~ ., data = trainVars)
# 
#       if (AIC(fit) < minCrit) {
#           minCrit = AIC (fit)
#           varlist = varlist[varlist != i]
#           #print(paste0(" Criterion : ", round(minCrit, 1) , ", variables : ",
#                        #paste0 (varnames[varlist], collapse = " ")))
#           m2 = fit
#           noChange = F
#       }
# 
#     }  else if  (!i %in% varlist) {
#       index = c(1, varlist, i)
#       trainVars = trainSet[, index]
# 
#       fit = lm(accuracy ~ ., data = trainVars)
# 
#       if (AIC(fit) < minCrit) {
#           minCrit = AIC(fit)
#           varlist = c(varlist, i)
#          # print(paste0(" Criterion : ", round(minCrit, 1), ", variables : ",
#                   #   paste0(varnames[varlist], collapse = " ")))
# 
#            m2 = fit
#            noChange = F
#         }
#      }
#    }
# }
# 
#  pred2 <- predict(m2 ,newdata = validSet)
#  RMSE2[p] <- sqrt(mean((validSet$accuracy - pred2)^2))
# 
# 
# 
# 
#  m3 <- stepAIC(object = empty, scope = list(upper = full, lower = empty), direction = "forward", trace = 0, k =log(nrow(trainSet)))
#  pred3 <- predict(m3 , newdata = validSet )
#  RMSE3[p] <- sqrt(mean(( validSet$accuracy - pred3)^2))
# 
# 
# 
# 
# pen <- log (nrow (trainSet))
# varlist = c()
# varnames = names(trainSet)
# n = nrow(trainSet)
# varoder <- sample(1:ncol(trainSet)) #random order of variables
# minCrit = Inf
# noChange = F
# 
# while(!noChange) {
#   noChange = T
#   for (z in varoder) {
#     if (z == 1)
#       next
# 
#     if (z %in% varlist & length(varlist) > 1) {
#       index = c(1, varlist[varlist != z])
#     trainVars = trainSet[, index]
# 
#       fit = lm(accuracy ~ ., data = trainVars)
# 
#       if (AIC(fit, k=pen) < minCrit) {
#           minCrit = AIC (fit, k = pen)
#           varlist = varlist[varlist!= z]
#           #print(paste0(" Criterion : ", round(minCrit, 1) , ", variables : ",
#                        #paste0 (varnames[varlist], collapse = " ")))
#           m4 = fit
#           noChange = F
#       }
# 
#     }  else if  (!z %in% varlist) {
#       index = c(1, varlist, z)
#       trainVars = trainSet[, index]
# 
#       fit = lm(accuracy ~ ., data = trainVars)
# 
#       if (AIC(fit, k = pen) < minCrit) {
#           minCrit = AIC(fit, k = pen)
#           varlist = c(varlist, z)
#          # print(paste0(" Criterion : ", round(minCrit, 1), ", variables : ",
#                      #paste0(varnames[varlist], collapse = " ")))
# 
#           m4 = fit
#           noChange = F
#        }
#     }
#   }
# }
# 
# pred4 <- predict(m4, newdata = validSet)
# RMSE4[p] <- sqrt(mean((validSet$accuracy - pred4)^2))
# 
#                   
# }              
#                   
#                   

# We apply selected procedure to the full dataset

load("mfinal")

# full <- lm(accuracy ~ ., data = ptrain)
# empty <- lm(accuracy ~ 1, data = ptrain)
# 
# mfinal <- stepAIC(object = empty, scope = list(upper = full, lower = empty), direction = "forward", trace = 0)

mfinal$coefficients


# We check the model for linear regression assumptions again


#summary(premov)
par ( mfrow = c (2 ,2) )

plot(mfinal$fitted.values, mfinal$residuals, xlab = "Fitted Values", 
     ylab = "Residuals")
plot (1: nrow (ptrain), mfinal$residuals , xlab = " Index ", ylab =
        " Residuals ")

qqnorm(mfinal$residuals)
qqline(mfinal$residuals, col="blue", lwd = 2)
hist(mfinal$residuals)



# use final model to predict response for protein-test.csv

ptest <- read.csv ("protein-test.csv")


predfinal <- predict(mfinal, newdata = ptest)

writeLines(as.character(predfinal), "mypreds1.txt.txt")













