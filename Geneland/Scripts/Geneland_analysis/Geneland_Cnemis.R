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

# Parâmetros para a análise MCMC
chain <- 1000000  # Número de iterações
freq <- 1000      # Frequência de amostragem
burnin <- 200     # Número de iterações descartadas no início
pop <- 4          # Número máximo de populações esperadas
all.freq.model <- 'Uncorrelated'  # Modelo de frequências alélicas (correlacionado ou não)

# Diretório para resultados de análises sem admixagem
noadmix.dir <- paste(maindir, 'noadmix', '/', sep="")

# Cria o diretório para resultados de análises sem admixagem, se não existir
if (!file.exists(noadmix.dir)) {
  dir.create(file.path(noadmix.dir))
}
setwd(file.path(noadmix.dir))  # Define o diretório de trabalho

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

# Plota o número de populações inferidas
Plotnpop(path.mcmc="./", burnin=burnin, printit=TRUE, file="Number_of_Clusters.pdf", format="pdf")

# Gera um mapa mostrando os clusters
PosteriorMode(coordinates=coord, path.mcmc="./", file="map.pdf")

# Calcula estatísticas F a partir dos genótipos
Fstat.output(genotypes=geno, path.mcmc="./")

# Cria um mapa de probabilidade de pertencimento aos clusters
PlotTessellation(coord, path.mcmc='./', printit=TRUE, path='./')

# Configuração para executar múltiplas análises (loop)
nrun <- 10        # Número de execuções
chain <- 5000000  # Número de iterações
freq <- 5000      # Frequência de amostragem

for (irun in 1:nrun) {
  # Define o caminho para os resultados de cada execução
  path.mcmc <- paste(maindir, '10_Run', irun, "/", sep="")
  
  # Cria o diretório para cada execução, se não existir
  if (!file.exists(path.mcmc)) {
    dir.create(file.path(path.mcmc))
  }
  setwd(file.path(path.mcmc))  # Define o diretório de trabalho
  
  # Executa a análise MCMC
  MCMC(coordinates=coord,
       geno.dip.codom=geno,
       varnpop=TRUE,
       npopmax=pop,
       spatial=TRUE,
       freq.model="Uncorrelated",
       nit=chain,
       thinning=freq,
       path.mcmc=path.mcmc)
  
  # Processa os resultados
  PostProcessChain(coordinates=coord, path.mcmc=path.mcmc, nxdom=200, nydom=200, burnin=burnin)
}

# Calcula a densidade posterior para cada execução
lpd <- rep(NA, nrun)
for (irun in 1:nrun) {
  path.mcmc <- paste(maindir, '/', irun, "/", sep="")
  path.lpd <- paste(path.mcmc, "log.posterior.density.txt", sep="")
  lpd[irun] <- mean(scan(path.lpd)[-(1:burnin)])
}
order(lpd, decreasing=TRUE)  # Ordena os resultados por densidade posterior
lpd

# Gera gráficos para as execuções
for (irun in 1:nrun) {
  path.mcmc <- paste(maindir, '/', irun, "/", sep="")
  setwd(file.path(path.mcmc))
  Plotnpop(path.mcmc=path.mcmc, burnin=burnin)
  PosteriorMode(coordinates=coord, path.mcmc="./", file="map.pdf")
}

# Executa análises considerando admixagem
admix.dir <- "C:/Users/Marcelo/Dropbox/Pseudis tocantins/Analises/GENELAND/nuclear/admix"
noadmix.dir <- "C:/Users/Marcelo/Dropbox/Pseudis tocantins/Analises/GENELAND/nuclear/1/"

# Cria diretório para admixagem, se não existir
if (!file.exists(admix.dir)) {
  dir.create(file.path(admix.dir))
}
setwd(file.path(admix.dir))  # Define o diretório de trabalho

# Executa análise com admixagem
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
