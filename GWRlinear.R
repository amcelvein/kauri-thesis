##preprocessing 

library(GWmodel)
library(sp)
library(rgdal)
library(tmap)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(SciViews) #for shapiro.test
library(MASS) #for stepwise
library(ape) #for morans
library(AICcmodavg)
library(raster) #for export
library(ggplot2)

setwd('C://Users/amce702/OneDrive - The University of Auckland/Documents/GWR')

##import data from Shapefile
pointdata <- readOGR('.', 'modeling_points_arc206', stringsAsFactors=F)
parkextent <- readOGR('.', 'parkextent_selection_merged')
pointdata <- pointdata[, -1] #remove annoying column Arc adds
pointdata <- pointdata[, -1] #remove cluster column
pointdata <- pointdata[, -4] #remove PA presence column
names(pointdata)

##confirming shapefile data is projected correctly
#proj4string(pointdata) <- CRS("+proj=longlat +ellps=WGS84")
#proj4string(parkextent) <- CRS("+proj=longlat +ellps=WGS84")
NZ2000 <- CRS("+init=EPSG:2193")
pointdata <- spTransform(pointdata, NZ2000)
parkextent <- spTransform(parkextent, NZ2000)

#########create map of C###########
tm_shape(parkextent)+ tm_borders("black", lwd = .5) +
  tm_shape(pointdata) + tm_symbols(col = "C__", scale = .5) +
  tm_legend(legend.outside = TRUE, legend.outside.size = .4)


#####################################
#######GWR  model####################
#####################################

summary(pointdata)
DeVar <- "C__"
InDeVars <- c( "PK_presenc", "HW_presenc", "N__", "C_N_ratio", "pH", "bulk_densi"
               , "forest_flo", "canopy_hei", "kauri_circ", "dist_ocean", "dist_river", "dist_track", "dist_road" 
               , "dist_clean", "dist_detec", "TCT_bright", "TCT_wetnes", "TCT_greenn", "NDVI", "elevation" 
               , "slope_deg", "aspect_eas", "aspect_nor", "topo_wetne", "topo_flow_", "topo_rugge", "Long"      
               , "Lat")
points_subset = subset(pointdata, select = -c(North_pres,N__, C_N_ratio, pH, bulk_densi,TCT_bright,TCT_wetnes,TCT_greenn,topo_rugge, topo_wetne, topo_flow_,dist_clean, Lat, Long), na.action=na.omit)

InDeVars_subset = c("PK_presenc", "HW_presenc"
                    , "forest_flo", "canopy_hei", "kauri_circ", "dist_ocean", "dist_river", "dist_track", "dist_road" 
                    , "dist_detec", "NDVI", "elevation" 
                    , "slope_deg", "aspect_eas", "aspect_nor")

##distance matrix
distMat <- gw.dist(dp.locat=coordinates(points_subset))

##set optimal bandwidth
optimalBW <- bw.gwr(C__ ~ .
                    , data=points_subset
                    , longlat = TRUE
                    , approach="AICc"
                    , kernel="bisquare"
                    , adaptive=TRUE
                    , dMat=distMat
                    )

#model selection
modelSel <- gwr.model.selection(DeVar
                                , InDeVars_subset
                                , data=points_subset
                                , kernel="bisquare"
                                , adaptive=TRUE
                                , bw=optimalBW
                                , dMat=distMat
                                , longlat= TRUE)

#Visualise model selection in a circular plot
sortedModels <- model.sort.gwr(modelSel, numVars <- length(InDeVars_subset), ruler.vector = modelSel[[2]][,2])
modelList <- sortedModels[[1]]

pdf('gwrmodelselection.pdf', height=9, width=9)
gwr.model.view(DeVar, InDeVars_subset, model.list=modelList)
dev.off()

##AICc
n <- length(InDeVars_subset)
AICcList <- sortedModels[[2]][,2]
indices <- rep(n,n)
for (i in 2:n) {
  indices[i]=indices[i-1]+((n-i)+1)
}
AICcBestModelValues <- AICcList[indices]
BestModels <- sortedModels[[1]][indices]
BestModels[n]

variablesAsAdded <- c("kauri_circ", "forest_flo" ,"elevation", "slope_deg"
                      , "dist_track","dist_road", "dist_ocean" 
                      , "HW_presenc", "PK_presenc","dist_river", "dist_detec"
                      , "aspect_nor", "canopy_hei", "aspect_eas", "NDVI"
                      )

plot(cbind(1:15,AICcBestModelValues), col = "black", pch = 20, lty = 5, main = "AICc optimisation", ylab = "AICc", type = "b", axes=FALSE)
par(las=2)
axis(1, at=1:15, labels=variablesAsAdded)
axis(2, at=NULL, labels=TRUE)

#AICc difference of 3 or less (in this research- first 5 vars)
AICcDifference <- AICcBestModelValues[1:(n-1)]-AICcBestModelValues[2:n]

#run GWR with vars
gwrmodel <- gwr.basic(C__~kauri_circ+forest_flo+elevation+slope_deg+dist_track
                      , data=points_subset
                      , bw=optimalBW
                      , kernel="bisquare"
                      , adaptive=TRUE
                      , dMat=distMat
                      ) 

#create results table
results_GWR <- as.data.frame(gwrmodel$SDF)
GWRresults <- cbind(points_subset, as.matrix(results_GWR))

#####################################
#######linear  model#################
#####################################

#create linear model using GWR variables
GWR_var_linear_model <- lm(C__~kauri_circ+forest_flo+elevation+slope_deg+dist_track, data=points_subset)
summary(GWR_var_linear_model)

##predicted values
linear_results <- points_subset
linear_results$predictedC <- 18.387276 + 0.039659 *points_subset$kauri_circ + 0.323770*points_subset$forest_flo + 0.029234*points_subset$elevation - 0.956677*points_subset$slope_deg - 0.008786*points_subset$dist_track 
linear_results$globalRes <- points_subset$C__ - linear_results$predictedC
linear_results$stGlobalRes <- (((linear_results$globalRes) - mean(linear_results$globalRes)) / sd(linear_results$globalRes))
summary(linear_results)

##############552 version
#import 552 datapoints for linear
pointdata552 <- readOGR('.', 'modeling_points_arc552', stringsAsFactors=F)
pointdata552 <- pointdata552[, -1] #remove annoying column Arc adds
pointdata552 <- pointdata552[, -1] #remove cluster column
pointdata552 <- pointdata552[, -4] #remove PA presence column

#project
NZ2000 <- CRS("+init=EPSG:2193")
pointdata552 <- spTransform(pointdata552, NZ2000)

#make sure everything is numeric
pointdata552$slope_deg <- as.numeric(pointdata552$slope_deg)

#create subset of vars with no multicollinearity 
points_subset552 = subset(pointdata552, select = -c(North_pres,TCT_bright,TCT_wetnes,TCT_greenn,topo_rugge, topo_wetne, topo_flow_,dist_clean, Lat, Long), na.action=na.omit)

#create model
GWR_var_linear_model552<-lm(C__~kauri_circ+forest_flo+elevation+slope_deg+dist_track, data=points_subset)
summary(GWR_var_linear_model552)

#results & export
linear_results552 <- points_subset552
linear_results552$predictedC <- 18.387276 + 0.039659 *points_subset552$kauri_circ + 0.323770*points_subset552$forest_flo + 0.029234*points_subset552$elevation - 0.956677*points_subset552$slope_deg - 0.008786*points_subset552$dist_track 
summary(linear_results552)
shapefile(linear_results552, filename='linear_552.csv', overwrite=TRUE)

#########################################################################################
###compare GWR to linear global model############################################
############################################################################################

##map R-squared stat
tm_shape(parkextent)+ tm_borders("black", lwd = .5) +
tm_shape(GWRresults) + tm_symbols(col="Local_R2", palette="Greens", n=7, scale=.5) +
tm_legend(legend.outside = TRUE, legend.outside.size = .4)

##global moran's for the residuals
#double checked in Arc
GWRresults <- GWRresults[, -41] #remove coordx.1
GWRresults <- GWRresults[, -41] #remove coordx.2
shapefile(GWRresults, filename='GWRresults.csv')
shapefile(linear_results, filename='linear_results.shp')

#GWR morans
#first generate a distance matrix, then take inverse of the matrix values and replace the diagonal entries with zero
# p-value > 0.05, so the data is not autocorrelated
GWRpoint.dist <- as.matrix(dist(cbind(GWRresults$coords.x1, GWRresults$coords.x2)))
GWRpoint.dist.inv <- 1/GWRpoint.dist
diag(GWRpoint.dist.inv) <- 0
GWRmorans<- Moran.I(GWRresults$Stud_residual, GWRpoint.dist.inv)

#Linear morans
# p-value < 0.05, so the data is autocorrelated
# observed > expected, positively autocorrelated
linearpoint.dist <- as.matrix(dist(cbind(linear_results$coords.x1, linear_results$coords.x2)))
linearpoint.dist.inv <- 1/linearpoint.dist
diag(linearpoint.dist.inv) <- 0
Linearmorans<-Moran.I(linear_results$stGlobalRes, linearpoint.dist.inv)

#are residuals normally distributed?
shapiro.test(linear_results$stGlobalRes) #p=.00168, not normal
shapiro.test(GWRresults$Stud_residual) #p=.003876, not normal

#QQ plot of residuals
ggplot() +
  geom_qq(aes(sample = linear_results$stGlobalRes)) +
  geom_abline(color = "red") +
  coord_fixed()

ggplot() +
  geom_qq(aes(sample = GWRresults$Stud_residual)) +
  geom_abline(color = "red") +
  coord_fixed()

#histogram of residuals
GWRresultsdataframe <- data.frame(GWRresults)
ggplot(GWRresultsdataframe, aes(x=Stud_residual)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 

linearresultsdataframe <- data.frame(linear_results)
ggplot(linearresultsdataframe, aes(x=stGlobalRes)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")

########map local and global residuals

##creating map cutoffs
minT <- floor(min(points$stGlobalRes))
maxT <- ceiling(max(points$stGlobalRes))
diffT <- (maxT-minT)/3
breaksT <- seq(minT,maxT,diffT)
breaksresiduals <- c(-2.58, -1.96, 0, 1.96, 2.58)

##mapping global residuals from linear
tm_shape(parkextent)+ tm_borders("black", lwd = .5) +
tm_shape(linear_results) + tm_symbols(col="stGlobalRes", style="fixed"
, breaks=breaksresiduals, palette="RdYlGn", scale=.5) +
tm_legend(legend.outside = TRUE, legend.outside.size = .4)

#local residuals from GWR
tm_shape(parkextent)+ tm_borders("black", lwd = .5) +
tm_shape(GWRresults) + tm_symbols(col="Stud_residual", style="fixed", breaks=breaksresiduals, palette="RdYlGn", scale=.5) +
tm_legend(legend.outside = TRUE, legend.outside.size = .4)
summary(GWRresults)

#####betas and t values
whereNonSig <- which(GWRresults$PctUnem1_TV>-1.96 & mapGWR$PctUnem1_TV<1.96)
mapGWR$PctUnem1_beta_sig <- mapGWR$PctUnem1_beta
mapGWR$PctUnem1_beta_sig[whereNonSig] <- NA

##AICC
models1<- list(GWR_var_linear_model)
mod.names1<-c('linear model')
aictab(cand.set = models1, modnames = mod.names1)
GWRadjustedT<-gwr.t.adjust(gwrmodel)

##more comparison metrics
library(Metrics)
rmse(GWRresults$y, GWRresults$yhat)
mse(GWRresults$y, GWRresults$yhat)
mae(GWRresults$y, GWRresults$yhat)
mape(GWRresults$y, GWRresults$yhat)
rse(GWRresults$y, GWRresults$yhat)

rmse(linear_results$C__, linear_results$predictedC)
mse(linear_results$C__, linear_results$predictedC)
mae(linear_results$C__, linear_results$predictedC)
mape(linear_results$C__, linear_results$predictedC)
rse(linear_results$C__, linear_results$predictedC)