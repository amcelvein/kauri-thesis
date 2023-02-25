install.packages("randomForest")
library(randomForest)

#import and clean data
setwd('C://Users/amce702/OneDrive - The University of Auckland/Documents/R_Annie')

#for all points
RFdata <- read.csv('modeling_points_206.csv',stringsAsFactors=T)
RFdata <- RFdata[, -1] #remove cluster column

#for all points without soil
RFdata <- read.csv('no_soil_modeling_points_206.csv',stringsAsFactors=T)
RFdata <- RFdata[, -1] #remove cluster column
RFdata <- RFdata[, -6] #remove forest floor depth

#for RF_GWR
RFdata <- subset(RFdata, select=c(C..,kauri.circumference, forest.floor.depth, elevation, slope.deg, dist.track))

#make sure all variables are numeric
RFdata$C.. <- as.numeric(RFdata$C..)
RFdata$kauri.circumference<- as.numeric(RFdata$kauri.circumference)
RFdata$dist.rivers <- as.numeric(RFdata$dist.rivers)
RFdata$slope.deg <- as.numeric(RFdata$slope.deg)
RFdata$dist.road <- as.numeric(RFdata$dist.road)
RFdata$topo.flow.dir <- as.numeric(RFdata$topo.flow.dir)
RFdata$PK.presence <- as.numeric(RFdata$PK.presence)
RFdata$North.presence <- as.numeric(RFdata$North.presence)
RFdata$HW.presence <- as.numeric(RFdata$HW.presence)
RFdata$PA.presence <- as.numeric(RFdata$PA.presence)

#####################################################
#####training on 2/3 of 206 and predicting###########
#####################################################

##split 75%
library(caTools)
set.seed(2)
sample = sample.split(RFdata$C.., SplitRatio = .75)
train = subset(RFdata, sample == TRUE)
test  = subset(RFdata, sample == FALSE)

##find optimal number of variables selected at each split (mtry)
##MAKE SURE YOU SET TRAIN[-X] TO WHATEVER C IS
mtry <- tuneRF(train[-5],train$C..
               , ntreeTry=5000
               , stepFactor=1.5
               ,improve=0.001
               , trace=TRUE
               , plot=FALSE)

best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(mtry)
print(best.m)
tr = data.frame(mtry)
tr$RMSE = sqrt(tr[,2])

##create mtry rf
set.seed(2)
mtry_rf <- randomForest(C..~.
                        , data=train
                        , mtry=best.m
                        , ntree=5000
                        , importance=T)

## use tuneRF to quickly check the range of mtrys
##then use caret cross-validation to check
library(mlbench)
library(caret)
set.seed(2)
control <- trainControl(method="cv", number=5)
#control <- trainControl(method="boot") #if using bootstrapping method
tunegrid <- expand.grid(.mtry=c(7:22)) ##replace this with range from mtry

#create caret rf model
set.seed(2)
caret_rf <- train(C..~.
                   , data=train
                   , method="rf"
                   , metric="RMSE"
                   , tuneGrid=tunegrid
                   , ntree = 100
                   , trControl=control)

#comparing mtry and caret
ggplot(tr,aes(x=mtry,y=RMSE))+
  stat_summary(fun=mean,geom="line")

ggplot(caret_rf,aes(x=mtry,y=RMSE))+
  stat_summary(fun=mean,geom="line")

##test rf models, use best moving forward
print(caret_rf) #slightly better for all 206
print(mtry_rf)

######using caret moving forward
set.seed(2)
rf_class <- predict(caret_rf, newdata = test)
predictions <- cbind(data.frame(train_preds=rf_class),tested=(test$C..))

#scatterplot with caret results
plot(predictions$train_preds
     , predictions$tested
     , xlab='predicted'
     , ylab='actuals'
     , main='%C predicted vs actuals'
     , sub= 'dataset: 206 points, GWR varaibles')
abline(a = 0,                                     
       b = 1,
       col = "red",
       lwd = 2)
#text(28,17,"1:1 relationship", col='red')

########IMPORTANCE PLOTS####################
#variable importance using caret
library(caret)
set.seed(2)
varImp(caret_rf, top=10)

#export as 500x300 image
plot(varImp(caret_rf), top=10, main='Variable Importance Plot: RF_all_but_soil', ylab='Top 10 Variables')

#partial dependence plots
library(pdp)
library(randomForest)

#export as 9.5x7 pdf
imp <- importance(caret_rf$finalModel)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
op <- par(mfrow=c(3, 3))
for (i in seq_along(impvar)) {
  partialPlot(caret_rf$finalModel, train, impvar[i], xlab=impvar[i],
              main=paste("Partial Dependence on", impvar[i]))
}
par(op)

############RESULTS#########
#all results
rf_results <- test
rf_results$predictedC <- predictions$train_preds
rf_results$globalRes <- rf_results$C.. - rf_results$predictedC
rf_results$stGlobalRes <- (((rf_results$globalRes) - mean(rf_results$globalRes)) / sd(rf_results$globalRes))
summary(rf_results)
write.csv(rf_results, "C://Users/amce702/OneDrive - The University of Auckland/Documents/R_Annie/RF_all_but_soil.csv")

#residuals
shapiro.test(rf_results$stGlobalRes)
#RSS
sum(rf_results$globalRes^2)
#metrics
library(Metrics)
rmse(predictions$tested, predictions$train_preds)
mse(predictions$tested, predictions$train_preds)
mae(predictions$tested, predictions$train_preds)
mape(predictions$tested, predictions$train_preds)

#morans
library(ape)
rfpoint.dist <- as.matrix(dist(cbind(rf_results$Lat, rf_results$Long)))
rfpoint.dist.inv <- 1/rfpoint.dist
diag(rfpoint.dist.inv) <- 0
rfmorans<-Moran.I(rf_results$stGlobalRes, rfpoint.dist.inv)

######Adding in 553 dataset#######################
setwd('C://Users/amce702/OneDrive - The University of Auckland/Documents/GWR')

#import 552 data from CSV and project
pointdata552 <- read.csv('552_expansion.csv',stringsAsFactors=F)
pointdata552 <- pointdata552[, -1] ##remove cluster column

library(sp)
pointdata552 <- SpatialPointsDataFrame(pointdata552[,25:26], pointdata552)
NZ2000 <- CRS("+init=EPSG:2193")
pointdata552 <- spTransform(pointdata552, NZ2000)
proj4string(pointdata552) <- CRS("+proj=longlat +ellps=WGS84")

#pointdata552 <- subset(pointdata552, select=c(C..,kauri.circumference, forest.floor.depth, elevation, slope.deg, dist.track, Lat, Long))
pointdata552 <- subset(pointdata552, select=-c(forest.floor.depth))

#cleaning data type
pointdata552$C..<- as.numeric(pointdata552$C..)
pointdata552$kauri.circumference<- as.numeric(pointdata552$kauri.circumference)
pointdata552$dist.rivers <- as.numeric(pointdata552$dist.rivers)
pointdata552$slope.deg <- as.numeric(pointdata552$slope.deg)
pointdata552$dist.road <- as.numeric(pointdata552$dist.road)
pointdata552$topo.flow.dir <- as.numeric(pointdata552$topo.flow.dir)
pointdata552$PK.presence <- as.numeric(pointdata552$PK.presence)
pointdata552$North.presence <- as.numeric(pointdata552$North.presence)
pointdata552$HW.presence <- as.numeric(pointdata552$HW.presence)
pointdata552$PA.presence <- as.numeric(pointdata552$PA.presence)

#new predictions & export
rf_class552 <- predict(caret_rf, newdata = pointdata552)
predictions552 <- cbind(data.frame(train_preds=rf_class552),tested=(pointdata552$C..))
rf_results552 <- pointdata552
rf_results552$predictedC <- predictions552$train_preds
summary(rf_results552)
write.csv(rf_results552, "C://Users/amce702/OneDrive - The University of Auckland/Documents/R_Annie/RF_all_but_soil552.csv")
