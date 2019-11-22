########### CODE FOR Range Expansion ANALYSES ####################################
#installing Range expansion package 
source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
library(snpStats)

devtools::install_github("BenjaminPeter/rangeExpansion", ref="package")
library(rangeExpansion)
library(sp)
wd <- ".."
setwd(wd)


############ running it on the different populations 
#South
region <- list("REGION_1")
###raw data SNAPP
snp.file <- "Data/south/NF_85.snapp"
coords.file <- "Data/south/locs_south.csv" 

ploidy <- 2 #diploid individuals

#outgroup_columns <- NULL 
raw_data <- load.data.snapp(snp.file, coords.file, sep=',', ploidy=ploidy) 

#we calculate the population-level data from individual data and calculate all pairwise statistics:
pops <- make.pop(raw_data, ploidy)
psi_groups <- get.all.psi(pops)

#Find origin 
results_S <- run.regions(region=region, pop=pops, psi=psi_groups, xlen=10, ylen=20)

summary(results_S)
pdf("plots/NF_grps_refugia.pdf")
plot(results_S, add.map=T, add.samples=F, add.sample.het = F)
dev.off()


############ running it on the different populations 
#North
region <- list("REGION_1")
###raw data SNAPP
snp.file.N <- "Data/north/NF_north_all.snapp"
coords.file.N <- "Data/north/locs_north.csv" 

ploidy <- 2 #diploid individuals


#outgroup_columns <- NULL 
raw_data.N <- load.data.snapp(snp.file.N, coords.file.N, sep=',', ploidy=ploidy) 

#we calculate the population-level data from individual data and calculate all pairwise statistics:
pops.N <- make.pop(raw_data.N, ploidy)
psi.N <- get.all.psi(pops.N)

#Find origin 
results_N <- run.regions(region=region, pop=pops.N, psi=psi.N, xlen=10, ylen=20)

summary(results_N)

pdf("NF_range_expan.pdf")
par(mfrow=c(1,2))
plot(results_S, add.map=T, add.samples=T, add.sample.het = T)
plot(results_N, add.map=T, add.samples=T, add.sample.het = T)
dev.off()

