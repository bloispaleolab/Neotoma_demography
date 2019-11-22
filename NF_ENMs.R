########### CODE FOR ENM ANALYSES ####################################
#load in packages 
library(rJava)
library(ENMeval)
library(dismo)
library(spThin)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)
library(RColorBrewer)

#Load in combined locality data
neotoma_fuscipes <- read.csv("Data/Neo_fus_sorted.csv", header=T)


#load in locs and thin 
NF_north <- read.csv("Data/north.csv")

NF_south <- read.csv("Data/south.csv")

#spatial filter in spThin
NF_thin <- thin(loc.data = NF_north, 
                lat.col = "latitude", long.col = "longitude", 
                spec.col = "id", 
                thin.par = 5, reps = 100, 
                locs.thinned.list.return = TRUE, 
                write.files = TRUE, 
                max.files = 5,
                out.dir = "Data/thin/", out.base = "NF_north_f", 
                write.log.file = TRUE,
                log.file = "Data/thin/NF_thin_north_log.txt" )

#load in spatial thin dataset
north_F <- read.csv("Data/thin/NF_north_f_thin1.csv")


NF_thin <- thin(loc.data = NF_south, 
                lat.col = "latitude", long.col = "longitude", 
                spec.col = "id", 
                thin.par = 5, reps = 100, 
                locs.thinned.list.return = TRUE, 
                write.files = TRUE, 
                max.files = 5,
                out.dir = "Data/thin/", out.base = "NF_south_f", 
                write.log.file = TRUE,
                log.file = "Data/thin/NF_thin_south_log.txt" )

#load in spatial thin dataset
south_F <- read.csv("Data/thin/NF_south_f_thin1.csv")

#put north and south datasets together 
coords <- list (north_F[2:3], south_F[2:3])

#Load in env data 
env_w <- list.files("Data/lorenz/0BP/", pattern='tif', full.names=TRUE)
env_w
predictors <- stack(env_w)
plot(predictors[[1]])
##Generating study region 
## function to create a minimum convex polygon (written by Jamie Kass 2015) 
simpleMCP <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}

# Generate MCP  
MCP <- lapply(coords, simpleMCP)
plot(MCP[[1]], add=T)
########BELOW HAS THE BUFFER TO CHANGE########
# Generate a buffered MCP (here 3.0 degrees)
MCP_buff <- lapply(MCP, gBuffer, width=3)
plot(MCP_buff[[1]], add=T)
#Clipping layers
layers <- list()
for (i in 1:length(MCP_buff))  {
  layers[[i]] <- crop(predictors, MCP_buff[[i]])
}
env <- list()
for (i in 1:length(MCP_buff))    {
  env[[i]] <- mask(layers[[i]], MCP_buff[[i]])
}

#ENMs in ENMeval
NF_pops <- list()
for (i in 1:length(coords))    {
  NF_pops[[i]] <- ENMevaluate(coords[[i]], env[[i]], mod.name="maxent.jar",
                              tune.args = list(fc = c("L", "LQ", "H", "LQH"), rm = seq(0.5, 6, 0.5)), 
                              partitions = "jackknife", overlap=F, updateProgress = F, numCores = 4)
}

## Look at eval. stats
neotoma_north_results <-NF_pops[[1]]@results
neotoma_north_results
write.csv(neotoma_north_results, "../output/ENMs/neotoma_north_results.csv", quote=FALSE, row.names=FALSE)
neotoma_south_results <-NF_pops[[2]]@results
neotoma_south_results
write.csv(neotoma_south_results, "../output/ENMs/neotoma_south_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for N. fuscipes northern population 
## Using the top perfroming model (for now)  ##
#projected to SR; #NF_seq_north best settings were Hinge; RM = 1 
args=c("noaddsamplestobackground","noautofeature", "noproduct","nolinear","noquadratic","nothreshold","betamultiplier=1.0", "responsecurves=true") 
pred.args_N <- c("outputformat=Cloglog", "doclamp=TRUE")

north_model <- maxent(env[[1]], north_F[2:3], args = args)
north_predict <- predict(env[[1]], north_model, args=pred.args_N) 
plot(north_predict)        

#NF_seq_south best settings were Linear, Quadratic; RM = 1.5
args=c("noaddsamplestobackground","noautofeature", "noproduct","nohinge","nothreshold","betamultiplier=1.5", "responsecurves=true") 
pred.args_S <- c("outputformat=Cloglog", "doclamp=TRUE")

south_model <- maxent(env[[2]], south_F[2:3], args = args)
south_predict <- predict(env[[2]], south_model, args=pred.args_S) 
plot(south_predict)  

#load in LGM layers
#load in variables 
env_lgm <- list.files("Data/lorenz/21K/", pattern='tif', full.names=TRUE)
env_lgm
predictors_lgm <- stack(env_lgm)

#generate MCP for all locs 
MCP_LGM <- simpleMCP(neotoma_fuscipes[3:2])

########BELOW HAS THE BUFFER TO CHANGE########
# Generate a buffered MCP (here 3.0 degrees)
MCP_buff_lgm <- gBuffer(MCP_LGM, width=3)
plot(MCP_buff[[1]], add=T)

#Clipping layers to predict to a larger area 
layers <-  crop(predictors_lgm, MCP_buff_lgm)
predict_lgm <- mask(layers, MCP_buff_lgm)
plot(predict_lgm)

pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#LGM predictions
north_predict_lgm <- predict(predict_lgm, north_model, args=pred.args) 
south_predict_lgm <- predict(predict_lgm, south_model, args=pred.args) 
plot(north_predict_lgm)
plot(south_predict_lgm)

#add in range expansion data 
N <- data.frame(39.17, -121.39)
S <- data.frame(37.9, -122.18)

####Load in US shapefile
US<- readOGR("Data/states_21basic/states.shp")

####Load in US lakes
lakes <- readOGR(dsn = "Data/50m-rivers-lake-centerlines/", layer = "50m-rivers-lake-centerlines")
####Load in paleo ice/lakes
load("Data/IceSheets.RData")

#plots current
pdf("../Figures/enms.pdf")
par(mfrow=c(2,2))
#north
plot(north_predict, col=rev(terrain.colors(10)), ylim=c(25, 50), xlim=c(-100, -80), breaks= seq(0, 1, by = .1))
plot(lakes, col="blue", add=T);
plot(US, ylim=c(25, 50), xlim=c(-100, -80), add=T)
points(north_F[2:3], pch=20)
#south
plot(south_predict,col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1))
plot(lakes, col="blue", add=T)
plot(US, ylim=c(25, 50), xlim=c(-100, -80), add=T)
points(south_F[2:3], pch=20)
#LGM; north
plot(north_predict_lgm, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1))
plot(iceList[[43]], col="lightblue", add=T)
plot(US, ylim=c(25, 50), xlim=c(-100, -80), add=T)
points(N[2:1], pch = 25, col="red")
#south
plot(south_predict_lgm, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1))
plot(iceList[[43]], col="lightblue", add=T)
plot(US, ylim=c(25, 50), xlim=c(-100, -80), add=T)
points(S[2:1], pch = 25, col="red")
dev.off()