########### CODE FOR IBD ANALYSES ####################################
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)
library(spThin)

#Load in vcf file for %85 of data 
wd <- ".."
setwd(wd)
nf_85 <- read.vcfR("Data/NF_85.recode.vcf")
nf_50 <- read.vcfR("Data/NF_50.recode.vcf")

#load in pop information 
nf_pop <- read.csv("Data/NF_groups.csv")
### convert to genlight
NF_85.genlight <- vcfR2genlight(nf_85, n.cores=2)
NF_50.genlight <- vcfR2genlight(nf_50, n.cores=2)

#add pop names
pop(NF_85.genlight)<- nf_pop[,2]
pop(NF_50.genlight)<- nf_pop[,2]

#generate list of genlight objects 
NF_SNPs <- list(NF_85.genlight, NF_50.genlight)

#Calculation and visualization of Nei’s distances (using lapply)
### Calculate Nei's distances between individuals/pops
NF.ind <- lapply(NF_SNPs, stamppNeisD, pop = FALSE) # Nei's 1972 distance between indivs
NF.pop <- lapply(NF_SNPs, stamppNeisD, pop = TRUE) # Nei's 1972 distance between pops


NF_85.genlight@ploidy <- as.integer(ploidy(NF_85.genlight))
NF_50.genlight@ploidy <- as.integer(ploidy(NF_50.genlight))

### Isolation by distance
coords <- read.csv ("Data/Neo_fus_sorted.csv") # tab-separated file for all pops
xy.coords.only<- subset(coords, select=c("DEC_LAT","DEC_LONG"))
Dgeo <- dist(xy.coords.only)

## Calculate distance in km 
DistMat <- rdist.earth(x1 = coords[2:3], miles = FALSE)
DistMat <- as.dist(DistMat)

# create the dist objects used in analyses below
colnames(NF.ind[[1]]) <- rownames(NF.ind[[1]])
NF.85.ind.dist<-as.dist(NF.ind[[1]], diag=T)
attr(NF.85.ind.dist, "Labels")<-rownames(NF.ind[[1]]) # name the rows of a matrix

colnames(NF.ind[[2]]) <- rownames(NF.ind[[2]])
NF.50.ind.dist<-as.dist(NF.ind[[2]], diag=T)
attr(NF.50.ind.dist, "Labels")<-rownames(NF.ind[[2]]) # name the rows of a matrix

#generate list 
NF_ind <- list(NF.85.ind.dist, NF.50.ind.dist)
#test IBD 
IBD_all <- lapply(NF_ind, mantel.randtest, Dgeo, nrepet=10000)
IBD_all
plot(Dgeo,NF.85.ind.dist, pch=20,cex=.5)
abline(lm(NF.85.ind.dist~Dgeo))

dis <- list(NF.85.ind.dist, NF.50.ind.dist)
#plot
pdf("Figures/NF_IBD_plot.pdf")
par(mfrow = c(2, 1))
for (i in 1:length(dis)) 
  plot(DistMat, dis[[i]], pch=20,cex=.5, xlab="Geographic Distance (km)",
       ylab="Genetic Distance", abline(lm(dis[[i]]~DistMat)), ylim=c(0,.4),
       cex.main=.9)
dev.off()

#############################
##IBD for northern A pops 
###plot AFS per one pop
NF_85.genlight.sep <- seppop(NF_85.genlight, drop=TRUE)
NF_50.genlight.sep <- seppop(NF_50.genlight, drop=TRUE)
#check if it worked 
NF_85.genlight.sep$North_A		
### Calculate Nei's distances between individuals/pops for 85
NF.85_north_A <- stamppNeisD(NF_85.genlight.sep$North_A, pop = FALSE) # Nei's 1972 distance between indivs
colnames(NF.85_north_A) <- rownames(NF.85_north_A)
NF.85.dist_a<-as.dist(NF.85_north_A, diag=T)
attr(NF.85.dist_a, "Labels")<-rownames(NF.85_north_A) # name the rows of a matrix

North_A <- read.csv("Data/NF_north_a.csv")
North_A_geo <- merge(North_A, coords, by.x = "Ind", by.y = "ID")
xy.coords.only<- subset(North_A_geo, select=c("DEC_LAT","DEC_LONG"))
Dgeo_a <- dist(xy.coords.only)

## Calculate distance in km 
DistMat_a <- rdist.earth(x1 = North_A_geo[3:4], miles = FALSE)
DistMat_a <- as.dist(DistMat_a)

#test IBD
IBD_85.A <- mantel.randtest(DistMat_a, NF.85.dist_a)
IBD_85.A
plot(DistMat_a, NF.85.dist_a, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(NF.85.dist_a~DistMat_a))

#####IBD for northern B pops
NF_85.genlight.sep$North_B
### Calculate Nei's distances between individuals/pops
NF.85_north_B <- stamppNeisD(NF_85.genlight.sep$North_B, pop = FALSE) # Nei's 1972 distance between indivs
colnames(NF.85_north_B) <- rownames(NF.85_north_B)
NF.85.dist_b<-as.dist(NF.85_north_B, diag=T)
attr(NF.85.dist_b, "Labels")<-rownames(NF.85_north_B) # name the rows of a matrix

North_B <- read.csv("Data/NF_north_b.csv")
North_B_geo <- merge(North_B, coords, by.x = "Ind", by.y = "ID")
xy.coords.only<- subset(North_B_geo, select=c("DEC_LAT","DEC_LONG"))
Dgeo_b <- dist(xy.coords.only)

## Calculate distance in km 
DistMat_b <- rdist.earth(x1 = North_B_geo[3:4], miles = FALSE)
DistMat_b <- as.dist(DistMat_b)

#test IBD
IBD_85.B <- mantel.randtest(DistMat_b, NF.85.dist_b)
IBD_85.B
plot(DistMat_b, NF.85.dist_b, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(NF.85.dist_b~DistMat_b))

#####IBD for southern pops
NF_85.genlight.sep$South
### Calculate Nei's distances between individuals/pops
NF.85_South <- stamppNeisD(NF_85.genlight.sep$South, pop = FALSE) # Nei's 1972 distance between indivs
colnames(NF.85_South) <- rownames(NF.85_South)
NF.85.dist_c<-as.dist(NF.85_South, diag=T)
attr(NF.85.dist_c, "Labels")<-rownames(NF.85_South) # name the rows of a matrix

South <- read.csv("GIS/NF_south.csv")
South_geo <- merge(South, coords, by.x = "Ind", by.y = "ID")
xy.coords.only<- subset(South_geo, select=c("DEC_LAT","DEC_LONG"))
Dgeo_c <- dist(xy.coords.only)

## Calculate distance in km 
DistMat_c <- rdist.earth(x1 = South_geo[3:4], miles = FALSE)
DistMat_c <- as.dist(DistMat_c)

#test IBD
IBD_85.C <- mantel.randtest(DistMat_c, NF.85.dist_c)
IBD_85.C
plot(DistMat_c, NF.85.dist_c, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(NF.85.dist_c~DistMat_c))

#####IBD for north pops
NF_85.genlight.2 <- NF_85.genlight
#add pop names
pops <- read.csv("Data/NF_2groups.csv")
pop(NF_85.genlight.2)<- pops[,2]
NF_85.genlight.2.sep <- seppop(NF_85.genlight.2, drop=TRUE)   #separate genlights per population
NF_85.genlight.2.sep$North
### Calculate Nei's distances between individuals/pops
NF_85.D.ind_N <- stamppNeisD(NF_85.genlight.2.sep$North, pop = FALSE) # Nei's 1972 distance between indivs
colnames(NF_85.D.ind_N) <- rownames(NF_85.D.ind_N)
NF_85.D.ind.dist_N<-as.dist(NF_85.D.ind_N, diag=T)
attr(NF_85.D.ind.dist_N, "Labels")<-rownames(NF_85.D.ind_N) # name the rows of a matrix

NF_N <- read.csv("Data/NF_north.csv")
NF_N_geo <- merge(NF_N, coords, by.x = "Ind", by.y = "ID")
xy.coords.only<- subset(NF_N_geo, select=c("DEC_LAT","DEC_LONG"))
Dgeo_N <- dist(xy.coords.only)

## Calculate distance in km 
DistMat_N <- rdist.earth(x1 = NF_N_geo[3:4], miles = FALSE)
DistMat_N <- as.dist(DistMat_N)

#test IBD
IBD <- mantel.randtest(DistMat_N, NF_85.D.ind.dist_N)
IBD
plot(DistMat_N, NF_85.D.ind.dist_N, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(NF_85.D.ind.dist_N~DistMat_N))


##save IBD plots 
gen_dist <- list(NF_85.D.ind.dist_N, NF.85.dist_a, NF.85.dist_b, NF.85.dist_c)
geo_dist <- list(DistMat_N, DistMat_a, DistMat_b, DistMat_c)
#plot
pdf("Figures/NF.85_IBD_plot.pdf")
par(mfrow = c(2, 2))
for (i in 1:length(gen_dist)) 
  plot(geo_dist[[i]], gen_dist[[i]], pch=20,cex=.5, xlab="Geographic Distance (km)",
       ylab="Genetic Distance", abline(lm(gen_dist[[i]]~geo_dist[[i]])), ylim=c(0,.15),
       cex.main=.9)
dev.off()
