########### CODE FOR PCA/DAPC ANALYSES ####################################
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)

#Load in vcf file for %85 threshold of data 
wd <- "."
setwd(wd)
nf_85 <- read.vcfR("Data/NF_85.recode.vcf")
head(nf_85)               #check the vcf object
nf_85@fix[1:10,1:5]       #check


#quick check read depth distribution per individual
pdf("Figures/read_depth_per_ind_85.pdf")
dp <- extract.gt(nf_85, element='DP', as.numeric=TRUE)
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)",
        las=2, cex=0.4, cex.axis=0.5)
dev.off()


#load in pop information 
nf_pop <- read.csv("Data/NF_groups.csv")
### convert to genlight
NF_85.genlight <- vcfR2genlight(nf_85, n.cores=2)
#add pop names
pop(NF_85.genlight)<- nf_pop[,2]

# check the genlight
NF_85.genlight                        # check the basic info on the genlightobject
indNames(NF_85.genlight)              # check individual names
as.matrix(NF_85.genlight)[1:16,1:10]  # see tiny bit of the data
pop(NF_85.genlight)                   # population assignment
# look at the total data matrix (0,1,2; white = missing data)
glPlot (NF_85.genlight)  # takes some time
# N missing SNPs per sample
x_85 <- summary(t(as.matrix(NF_85.genlight)))
#write.table(x_85[7,], file = "PCA/missing.persample.txt", sep = "\t")  # NAs, if present, are in seventh row of summary

##PCA
toRemove_85 <- is.na(glMean(NF_85.genlight, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove_85) # position of entirely non-typed loci
b_85 <- NF_85.genlight[, !toRemove_85]
pca.85 <- glPca(b_85, nf=300, n.cores=4) # this should work
#pca.85 <- glPca(NF_85.genlight, nf=300, n.cores=4)     # retain first 300 axes (for later use in find.clusters); slow function

# proportion of explained variance by first three axes
pca.85$eig[1]/sum(pca.85$eig) # proportion of variation explained by 1st axis
pca.85$eig[2]/sum(pca.85$eig) # proportion of variation explained by 2nd axis
pca.85$eig[3]/sum(pca.85$eig) # proportion of variation explained by 3rd axis
# save fig
pdf ("Figures/PCA_NF_85.pdf", width=14, height=7)
col <- c("blue", "green", "orange")
g1 <- s.class(pca.85$scores, pop(NF_85.genlight), xax=1, yax=2,
              col=transp(col,.6),
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T,
              pgrid.draw =F, plot = FALSE)
g2 <- s.label (pca.85$scores, xax=1, yax=2, ppoints.col = "red", plabels =
                 list(box = list(draw = FALSE),
                      optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#DAPC
grp_85 <- find.clusters(NF_85.genlight, max.n.clust=10, glPca = pca.85, perc.pca =
                          100, n.iter=1e6, n.start=1000)
dapc_85 <- dapc(NF_85.genlight, grp_85$grp, glPca = pca.85)

groups <- data.frame(grp_85$grp)
col <- c("orange", "green", "blue")
scatter(dapc_85, col = col)

pdf ("Figures/NF_85_DPAC.pdf", width=14, height=7)
#par(mfrow = c(1, 2))
scatter(dapc_85, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=col, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("South", "North_B", "North_A"))
dev.off()


#Load in vcf file for %50 threshold of data 
nf_50 <- read.vcfR("Data/NF_50.recode.vcf")
head(nf_50)               #check the vcf object
nf_50@fix[1:10,1:5]       #check


#quick check read depth distribution per individual
pdf("Figures/read_depth_per_ind_50.pdf")
dp <- extract.gt(nf_50, element='DP', as.numeric=TRUE)
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)",
        las=2, cex=0.4, cex.axis=0.5)
dev.off()

### convert to genlight
NF_50.genlight <- vcfR2genlight(nf_50, n.cores=2)
#add pop names
pop(NF_50.genlight)<- nf_pop[,2]

# check the genlight
NF_50.genlight                        # check the basic info on the genlightobject
indNames(NF_50.genlight)              # check individual names
as.matrix(NF_50.genlight)[1:16,1:10]  # see tiny bit of the data
pop(NF_50.genlight)                   # population assignment
# look at the total data matrix (0,1,2; white = missing data)
glPlot (NF_50.genlight)  # takes some time
# N missing SNPs per sample
x_50 <- summary(t(as.matrix(NF_50.genlight)))
#write.table(x_85[7,], file = "PCA/missing.persample.txt", sep = "\t")  # NAs, if present, are in seventh row of summary

##PCA
toRemove_50 <- is.na(glMean(NF_50.genlight, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove_50) # position of entirely non-typed loci
b_50 <- NF_50.genlight[, !toRemove_50]
pca.50 <- glPca(b_50, nf=300, n.cores=4) # this should work
#pca.85 <- glPca(NF_85.genlight, nf=300, n.cores=4)     # retain first 300 axes (for later use in find.clusters); slow function

# proportion of explained variance by first three axes
pca.50$eig[1]/sum(pca.50$eig) # proportion of variation explained by 1st axis
pca.50$eig[2]/sum(pca.50$eig) # proportion of variation explained by 2nd axis
pca.50$eig[3]/sum(pca.50$eig) # proportion of variation explained by 3rd axis
# save fig
pdf ("Figures/PCA_NF_50.pdf", width=14, height=7)
col <- c("blue", "green", "orange")
g1 <- s.class(pca.50$scores, pop(NF_50.genlight), xax=1, yax=2,
              col=transp(col,.6),
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T,
              pgrid.draw =F, plot = FALSE)
g2 <- s.label (pca.50$scores, xax=1, yax=2, ppoints.col = "red", plabels =
                 list(box = list(draw = FALSE),
                      optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#DAPC
grp_50 <- find.clusters(NF_50.genlight, max.n.clust=10, glPca = pca.50, perc.pca =
                          100, n.iter=1e6, n.start=1000)
dapc_50 <- dapc(NF_50.genlight, grp_50$grp, glPca = pca.50)

groups <- data.frame(grp_50$grp)
col <- c("green", "orange", "blue")
scatter(dapc_50, col = col)

pdf ("Figures/NF_50_DPAC.pdf", width=14, height=7)
#par(mfrow = c(1, 2))
scatter(dapc_50, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=col, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("North_B", "South", "North_A"))
dev.off()

