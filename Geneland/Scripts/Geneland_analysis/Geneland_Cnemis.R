# Carrega a biblioteca Geneland
library(Geneland)

# Define o diretório principal onde serão criadas pastas e salvos os resultados
maindir <- "/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/Scripts/Geneland_analysis/"
setwd(maindir)  # Define o diretório de trabalho

# Carrega os arquivos de entrada
# Carrega a matriz genotípica (nuclear)
geno.temp <- read.table("Proceratophrys_boiei_matrix_sorted_nuclear_no_individual_column_NEW.txt", na.string="NA")
geno <- as.matrix(geno.temp)  # Converte para matriz

# Carrega as coordenadas geográficas dos indivíduos
coord.temp <- read.table("/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/Scripts/Proceratophrys_boiei_haplotype_matrix_sorted_nuclear_with_individuals_NEW_geo.txt", na.string='NA')
coord <- as.matrix(coord.temp)

# Carrega os dados mitocondriais
mt.temp <- read.table('/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/Scripts/Geneland_analysis/filtered_mitochondrial_data_no_header.txt', na.string='NA')
mt <- as.matrix(mt.temp)

###se eu tiver dado morfológicos tem como usar
#morf.temp <- read.table('morfologia.txt', na.string='NA')
#morf <- as.matrix(morf.temp)


# Parâmetros para a análise MCMC
chain <- 1000000  # Número de iterações
freq <- 1000      # Frequência de amostragem
burnin <- 200     # Número de iterações descartadas no início
pop <- 5          # Número máximo de populações esperadas
all.freq.model <- 'Uncorrelated'  # Modelo de frequências alélicas (correlacionado ou não)

# Diretório para resultados de análises sem admixagem
noadmix.dir <- paste(maindir, 'noadmix', '/', sep="")

# Cria o diretório para resultados de análises sem admixagem, se não existir
if (file.exists(noadmix.dir)){
  setwd(file.path(noadmix.dir))
} else {
  dir.create(file.path(noadmix.dir))
  setwd(file.path(noadmix.dir))
} # Define o diretório de trabalho

# Executa a análise MCMC sem admixagem
MCMC(coordinates=coord,
     geno.dip.codom=geno,
     geno.hap=mt,  # Dados mitocondriais; comente esta linha se não quiser usá-los
     varnpop=TRUE,
     npopmax=pop,
     spatial=TRUE,
     freq.model=all.freq.model,
     nit=chain,
     thinning=freq,
     path.mcmc="./")

# Processa os resultados após a análise MCMC
PostProcessChain(coordinates=coord,
                 path.mcmc="./",
                 nxdom=200,
                 nydom=200,
                 burnin=burnin)

#### plot num of pops
Plotnpop(path.mcmc="./", burnin=200,)

# Plota o número de populações inferidas
Plotnpop(path.mcmc="./", burnin=burnin, printit=TRUE, file="Number_of_Clusters.pdf", format="pdf")

# Gera um mapa mostrando os clusters
PosteriorMode(coordinates=coord, path.mcmc="./", file="map.pdf")

# Calcula estatísticas F a partir dos genótipos
Fstat.output(genotypes=geno, path.mcmc="./")

# Cria um mapa de probabilidade de pertencimento aos clusters
PlotTessellation(coord, path.mcmc='./', printit=TRUE, path='./')


#######Loop for multiple runs#######################################################################

maindir <- "/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/Scripts/Geneland_analysis/" #### dir 


nrun <- 5        ### number of runs
burnin <- 200     ### burnin
chain <- 500000  ### number of steps in chain
freq <-      500 ### sampling frequency
pop <- 5          ### max number of pops 

for(irun in 1:nrun)
{
  ## define path to MCMC directory
  path.mcmc <- paste(maindir,'10_Run',irun,"/",sep="")
  
  ### create directory
  if (file.exists(path.mcmc)){
    setwd(file.path(path.mcmc))
  } else {
    dir.create(file.path(path.mcmc))
    setwd(file.path(path.mcmc))
  }
  
  
  system(paste(path.mcmc))
  MCMC(coordinates=coord,
       geno.dip.codom=geno,
       #geno.hap=mt,
       varnpop=TRUE,
       npopmax=pop,
       spatial=TRUE,
       freq.model="Uncorrelated",
       nit=chain,
       thinning=freq,
       path.mcmc=path.mcmc)
  
  
  PostProcessChain(coordinates=coord,
                   path.mcmc=path.mcmc,
                   nxdom=200,
                   nydom=200,
                   burnin=burnin)
}


#### log posterior density calculation ###
lpd <- rep(NA,nrun)

for (irun in 1:nrun) {
  path.mcmc <- paste(maindir, "10_Run", irun, "/", sep="") # Ajuste para "10_Run"
  path.lpd <- paste(path.mcmc, "log.posterior.density.txt", sep="")
  
  if (file.exists(path.lpd)) { # Verifica se o arquivo existe
    lpd[irun] <- mean(scan(path.lpd)[-(1:burnin)])
  } else {
    warning(paste("Arquivo não encontrado:", path.lpd))
    lpd[irun] <- NA
  }
}


order(lpd,decreasing=TRUE)
lpd



##### plot number of pops

for(irun in 1:nrun)
{
  path.mcmc <- paste(maindir,"10_Run",irun,"/",sep="")
  setwd(file.path(path.mcmc))
  Plotnpop(path.mcmc=path.mcmc, burnin=200)
}

for(irun in 1:nrun)
{
  path.mcmc <- paste(maindir,"10_Run",irun,"/",sep="")
  setwd(file.path(path.mcmc))
  PosteriorMode(coordinates=coord, path.mcmc="./", file="map.pdf")
}








###Admixture

geno <- as.matrix(as.numeric(geno))
geno[is.na(geno)] <- NA

geno <- as.matrix(apply(geno, 2, as.numeric))

HZ(coordinates=coord,   
   geno.dip.codom=geno,
   path.mcmc.noadm=maindir,"10_Run",
   nit=20000,
   thinning=10,
   path.mcmc.adm="/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/Scripts/Geneland_analysis/")


admix.dir <- "C:/Users/Marcelo/Dropbox/Pseudis tocantins/Analises/GENELAND/nuclear/admix"
noadmix.dir <- "C:/Users/Marcelo/Dropbox/Pseudis tocantins/Analises/GENELAND/nuclear/1/"

if (file.exists(admix.dir)){
  setwd(file.path(admix.dir))
} else {
  dir.create(file.path(admix.dir))
  setwd(file.path(admix.dir))
}



HZ(coordinates=coord,
   geno.dip.codom=geno,
   path.mcmc.noadm=noadmix.dir,
   a.init=1,
   b.init=1000,
   a.max=10,
   estimate.a=FALSE,
   estimate.b=TRUE,
   nit=chain,
   thinning=freq,
   path.mcmc.adm=admix.dir)

