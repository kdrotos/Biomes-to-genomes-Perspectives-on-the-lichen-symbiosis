# KDrotos PhD Chapter 2 #
# First pass at a maxent model to examine global lichen coverage #

# last updated 2023-03-13

# SETUP ----

setwd("C:/Users/Katherine/Documents/Guelph/PhD/THESIS/Chapter 2 - Global lichen coverage/MaxentModels/C. rangiferina maxent run/Cr rerun output")

# most of this is pulled from the MaxEnt tutorial, which has some R stuff at the end

# MaxEnt packages etc

install.packages("ROCR", dependencies=TRUE)
install.packages("vcd", dependencies=TRUE)
install.packages("raster")


library(ROCR)
library(vcd)
library(boot)
library(raster)
library(rgdal) # this will be retired in 2023, move to GDAL and PROJ

library(ggplot2)
library(exactextractr)

# read in data (produced by MaxEnt software)

Cr_presence <- read.csv("Cladonia_rangiferina_samplePredictions.csv")
Cr_avgs <- read.csv("Cladonia_rangiferina_sampleAverages.csv")
Cr_background <- read.csv("Cladonia_rangiferina_backgroundPredictions.csv")

Cr_map <- raster("Cladonia_rangiferina.asc")

# ROC ------
# this is in prep to do an ROC, which we don't need to do this time (it's part of MaxEnt software's default output)
Cr_pp <- Cr_presence$Cloglog.prediction # get the prediction column
testCrpp <- Cr_pp[Cr_presence$Test.or.train=="test"] # select only test points
trainCrpp <- Cr_pp[Cr_presence$Test.or.train=="train"] # select only train points

Crbb <- Cr_background$Cloglog

# ROCR requires a specific format for the prediction values:

combined <- c(testCrpp, Crbb) # combine into a single vector
label <- c(rep(1,length(testCrpp)),rep(0,length(Crbb))) # labels: 1=present, 0=random
pred <- prediction(combined, label) # labeled predictions
perf <- performance(pred, "tpr", "fpr") # True / false positives, for ROC curve
plot(perf, colorize=TRUE) # Show the ROC curve
performance(pred, "auc")@y.values[[1]] # Calculate the AUC


# Working from the map -----

hist(Cr_map)

# cropping the raster

cropbox1 <- drawExtent() # allows clicking to define top-left and bottom-right boundaries of cropped box

Canada_crop <- crop(Cr_map, cropbox1)
plot(Canada_crop)

hist(Canada_crop, maxpixels=3600000)

# manually cropping to just Canada and Alaska

cropbox2 <-c(-169.875, -50, 42, 83)
Canada_crop2 <- crop(Cr_map, cropbox2)
plot(Canada_crop2)

# histogram of frequency vs. values
Can_Cr_hist1 <- hist(Canada_crop2, maxpixels=3600000, col="lightblue", xlim=c(0,1), xlab="Occurrence probability")

# make the map nice

plot(Canada_crop2, legend = TRUE, legend.args = list(text="Occurrence prediction probability"))

# going to export it without the legend title, it's overlapping and not necessary at this time

plot(Canada_crop2, legend = TRUE)

# summary for Canada_crop2

#class      : RasterLayer 
#dimensions : 984, 2877, 2830968  (nrow, ncol, ncell)
#resolution : 0.04166667, 0.04166667  (x, y)
#extent     : -169.875, -50, 42, 83  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +no_defs 
#source     : memory
#names      : Cladonia_rangiferina 
#values     : 0.00562724, 0.741056  (min, max)



# estimating total cover -----

cellStats(Canada_crop2, 'sum')
# returned = 506,676.8

# check if ocean (values of 0) are being included

summary(getValues(Canada_crop2))

# they are not! so: 
# total pixels - NA pixels = polygon pixels
2830968 - 1390031 = 1440937

# then can divide the sum of values by ncell to yield value/pixel

506676.8/1440937
# returned = 0.3516301
