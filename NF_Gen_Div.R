########## CODE FOR GENETIC DIVERISTY ANALYSES ####################################
###USing thetamater to calculate genetic diversity 
library(devtools)
install_github("radamsRHA/ThetaMater")
library(ThetaMater) # Load package ThetaMater
library(MCMCpack) # Load dependancy phybase
library(ape)
library(phangorn)
wd <- ".."
setwd(wd)

NF_data_S <- Read.AllelesFile(alleles.file = "Data/NF_south.loci")
NF_data_N.B <- Read.AllelesFile(alleles.file = "Data/NF_northB.loci")
NF_data_N.A <- Read.AllelesFile(alleles.file = "Data/NF_northA.loci")

# Let's look at the data
NF_data_S$k.vec # number of segregating sites
NF_data_S$l.vec # locus lengths
NF_data_S$n.vec # number of samples
NF_data_S$c.vec # number of observations (i.e., sum(NF_data$c.vec) = number of loci)


NF.MCMC.S <- ThetaMater.M1(k.vec = NF_data_S[[1]], l.vec = NF_data_S[[3]], n.vec = NF_data_S[[2]], 
                           c.vec = NF_data_S[[4]], ngens = 1000000, burnin = 100000, theta.shape = 50, theta.scale = .0001, thin = 10)


NF.MCMC.Na <- ThetaMater.M1(k.vec = NF_data_N.A[[1]], l.vec = NF_data_N.A[[3]], n.vec = NF_data_N.A[[2]], 
                            c.vec = NF_data_N.A[[4]], ngens = 1000000, burnin = 100000, theta.shape = 50, theta.scale = .0001, thin = 10)


NF.MCMC.Nb <- ThetaMater.M1(k.vec = NF_data_N.B[[1]], l.vec = NF_data_N.B[[3]], n.vec = NF_data_N.B[[2]], 
                            c.vec = NF_data_N.B[[4]], ngens = 1000000, burnin = 100000, theta.shape = 50, theta.scale = .0001, thin = 10)

mean(NF.MCMC.S) 
mean(NF.MCMC.Na)
mean(NF.MCMC.Nb)

# Boxplot of theta for each population  
South <- NF.MCMC.S
North_A <- NF.MCMC.Na
North_B <- NF.MCMC.Nb

NF.MCMC <- data.frame(North_A, North_B, South)

pdf("Figures/NF_gen_div.pdf")
boxplot(NF.MCMC, main="Genetic Diversity, θ", 
        xlab="Populations", ylab="Genetic Diversity")

dev.off()
