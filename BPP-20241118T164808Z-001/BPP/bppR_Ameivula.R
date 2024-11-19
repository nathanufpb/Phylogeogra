#devtools::install_github("dosreislab/bppr")
library(bppr)
library(ape)

setwd("/Users/felipemedeiros/Dropbox/Disciplinas UFPB/Filogeografia/Filogeografia_Felipe_Camura/Dados_alunos/Dia4/BPP/BPP_7/demografia/results/")

#import Ameivula MCMC file
files <- list.files(pattern="mcmc")
mcmc_concat <- do.call(`rbind`,lapply(files, read.table, header=T))

summary(mcmc_concat)

time <- bppr::msc2time.r(mcmc_concat, u.m = 0.001,
                         u.sd = 9e-6, g.m = 1.5, g.sd = 0.5)

# posterior means:
apply(time, 2, mean) -> mean.IM

# posterior HPDs:
coda::HPDinterval(coda::as.mcmc(time)) -> HPD.IM

#save in a table
mean.IM <- as.data.frame(mean.IM)
HPD.IM <- as.data.frame(HPD.IM)

bpp_demography <- data.frame(mean.IM,HPD.IM)
write.csv(bpp_demography, "bpp_demography.csv")

Ameivula_tree <- read.tree("../../species_tree")
#sceloporus.tree2 <- read.nexus("../tree.nexus")
ape::plot.phylo(Ameivula_tree)

mcmc2densitree(Ameivula_tree, time, "t_", col = "darkgreen",thin=0.025, alpha=0.025,pfract = 0.15,y.offset = 0.1)
title(xlab="Divergence time (Ma)")


#import Ameivula MCMC file (MSC-M model)
setwd("/Users/felipemedeiros/Dropbox/Disciplinas UFPB/Filogeografia/Filogeografia_Felipe_Camura/Dados_alunos/Dia4/BPP/BPP_7/demografia_migracao2/results/")
  
files.mig <- list.files(pattern="mcmc")
mcmc_concat.mig <- do.call(`rbind`,lapply(files.mig, read.table, header=T))

summary(mcmc_concat.mig)

time.mig <- bppr::msc2time.r(mcmc_concat.mig, u.m = 0.001,
                         u.sd = 9e-6, g.m = 1.5, g.sd = 0.5)
# posterior means:
apply(time.mig, 2, mean) -> mean.mig

# posterior HPDs:
coda::HPDinterval(coda::as.mcmc(time.mig)) -> HPD.mig

#save in a table
mean.mig <- as.data.frame(mean.mig) 
HPD.mig <- as.data.frame(HPD.mig)

bpp_demography_mig <- data.frame(mean.mig,HPD.mig)
write.csv(bpp_demography_mig, "bpp_demography_mig.csv")

mcmc2densitree(Ameivula_tree, time.mig, "t_", col = "darkgreen",thin=0.09, alpha=0.15,pfract = 0.15,y.offset = 0.1)
title(xlab="Divergence time (Ma)")

