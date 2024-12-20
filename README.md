# Filogeografia

## Sumário

1. [Download e limpeza dos dados](#1-download-e-limpeza-dos-dados)
2. [Alinhamento e Faseamento](#2-alinhamento-e-faseamento)
3. [Seleção de Modelos](#3-seleção-de-modelos)
4. [Árvores de Genes](#4-árvores-de-genes)
5. [Inferência da Estrutura Populacional](#5-inferência-da-estrutura-populacional)
6. [Estimativas genéticas, Análise de variância molecular (AMOVA) e Redes de haplótipos](#6-estimativas-genéticas-análise-de-variância-molecular-amova-e-redes-de-haplótipos)
7. [Demografia Histórica e Difusão espaço temporal](#7-demografia-histórica-e-difusão-espaço-temporal)

## Objetivo:
Este relatório, elaborado para a disciplina de Filogeografia, tem como objetivo descrever as etapas envolvidas na aplicação prática de ferramentas de bioinformática em análises filogeográficas. Dados de sequenciamento genético multilocus obtidos através do trabalho ***(CARNAVAL; BATES, 2007)***, foram utilizados para realizar as seguintes análises: construção de árvores filogenéticas, delimitação de linhagens, avaliação da diversidade genética, análise da estrutura populacional e reconstrução da história demográfica da nossa espécie modelo.

## 1. Download e limpeza dos dados.

Um script personalizado em python("***Search_ncbi_and_clear_data.ipynb***") foi utilizado para baixar as sequências de DNA no formato *.fasta através de um intervalo de código de acesso do GenBank (“https://www.ncbi.nlm.nih.gov/), e limpar os cabeçalhos de modo que a identificação de cada voucher fosse preservado e as sequências estivessem prontas para realizar o alinhamento separadamente para cada sequência, o Script está disponível no Github (https://github.com/nathanufpb/Phylogeogra.git). Os dados incluem os genes mitocondriais ***NADH dehydrogenase 2*** (***ND2***) e o ***Citocromo b*** (***Cytb***), além do gene nuclear ***proto-oncogene MYC fator de transcrição*** ***bHLH (c-myc)***.

## 2. Alinhamento e Faseamento.

Realizamos o alinhamento de cada sequência separadamente utilizando o algoritmo ***Muscle*** (EDGAR, 2004) através do software ***MEGA 11*** (TAMURA; STECHER; KUMAR, 2021) Figura 1. Utilizamos o algoritmo PHASE (STEPHENS; DONNELLY, 2003; STEPHENS; SMITH; DONNELLY, 2001) para separar as fases haplotípicas do gene nuclear ***c-myc***  utilizando a 5ª versão do software ***DnaSP*** (LIBRADO; ROZAS, 2009), obtendo o conjunto de diferentes alelos herdados juntos em um mesmo cromossomo (Figura 2). Um corte nas phases com similaridades abaixo de 60% foi realizado na sequência de alinhamento do gene ***c-myc*** utilizando o script ***python remove_low_phases.ipynb*** escrito no formato Júpiter Notebook contido no diretório do github, restaram 80 sequências (fases diferentes) no total.
Em contrapartida, os genes mitocondriais ***ND2*** e ***Cytb*** não foram submetidos a esse processo, uma vez que sua herança exclusivamente materna torna o faseamento desnecessário.

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/Alignament_mega11.png)

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/Phase_nuc.png)

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/Phase_nuc.png)

## 3. Seleção de Modelos

Utilizamos o ***Jmodeltest*** (DARRIBA et al., 2012; POSADA, 2008) para testar quais modelos de substituição nucleotídica se adequam melhor aos nossos dados de sequências alinhadas (Figura 4). Os resultados são mostrados na Tabela 1.

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/table_models.png)

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/Jmodeltest_analysis.jpg)

## 4. Árvores de Genes

Utilizamos uma abordagem bayesiana para calcular árvores filogenéticas através do software  ***BEAST2*** (BOUCKAERT et al., 2014). Primeiramente as sequências foram carregadas através do ***BEAUTI2*** onde foram especificados cada modelo para cada sequência (Figura 5). Em seguida a análise seguiu utilizando o ***BEAST2*** (Figura 6).

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/BEAST_analysis.jpg)

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/BEAUTI_insert_on_BEAST.png)

O programa ***Tracer*** versão 1.7.2 (RAMBAUT et al., 2018) foi utilizado para verificar o status e o resultado da análise bayesiana incorporando o arquivo de *.log gerado pelo output do ***BEAST2*** (Figura 7).

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/Tracer_result.png)

Em seguida, o software ***Treeannotator*** (BOUCKAERT et al., 2014) juntamente com o ***BEAST2*** foi utilizado para gerar uma árvore de consenso das árvores geradas (Figura 8). Após o consenso, a árvore foi editada para melhor visualização utilizando o ***Figtree*** versão 1.4 (Figura 9).

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/Treeannotator.png)

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/images/Trees_results_beast.jpg)

## 5. Inferência da Estrutura Populacional

### GeneLand

Geneland é um pacote escrito em linguagem R que permite fazer uso de genótipos multilocus individuais georreferenciados para a inferência do número de populações e da localização espacial de descontinuidades genéticas entre essas populações. As principais suposições do método são: ***(i)*** o número de populações é desconhecido e todos os valores são considerados a priori igualmente prováveis, ***(ii)*** as populações são espalhadas por áreas dadas por uma união de alguns polígonos de localização desconhecida no domínio espacial, ***(iii)*** o equilíbrio de ***Hardy-Weinberg*** é assumido dentro de cada população e ***(iv)*** as frequências de alelos em cada população são desconhecidas e tratadas como variáveis aleatórias seguindo o chamado modelo de Dirichlet ou modelo de Falush [(GUILLOT; MORTIER; ESTOUP, 2005)](https://www.zotero.org/google-docs/?oUPvZg).

Utilizamos o ***Geneland*** para verificar se existe estrutura populacional na espécie Procerathrophys boiei coletada ao longo do nordeste brasileiro sob brejos de altitude por Carnaval 2007.

![1732239873345](image/README/1732239873345.png)

A análise do ***Geneland*** começou criando pastas para os genes nucleares phaseados e separados em 2 exons e 1 intron, e uma pasta para os dois outros genes mitocondriais.

os scripts em python dentro da pasta Geneland/Scripts (***fasta2haplotype_matrix_mitocondrial.py*** e ***fasta2haplotype_matrix_nuclear_gen.py***)
foram utilizados para transformar os dados .fas em uma matriz haplotípica.

o script "***Configure_txt_file_and_create_coordenate_file.ipynb***" foi utilizado para criar um novo arquivo da matriz haplotípica dos genes mitocondriais contendo apenas as linhas que possuiam o mesmo código do gene nuclear(***42 genes***), e na construção de um arquivo contendo as coordenadas de cada linha mapeada para os códigos do arquivo nuclear.

o novo arquivo dos genes mitocondriais com o mesmo número de amostras e o mesmo código do arquivo com genes nucleares:

mt.temp = ***filtered_mitochondrial_data_no_header.txt***

geno.temp = ***Proceratophrys_boiei_matrix_sorted_nuclear_no_individual_column_NEW.txt***

coord.temp = ***Proceratophrys_boiei_haplotype_matrix_sorted_nuclear_with_individuals_NEW_geo.txt***

maindir = /home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/Scripts/Geneland_analysis/   ***#defina um diretório***

Com esses arquivos foi possível executar o Geneland através do Rscript "***Geneland_Cnemis.R***"

![1732239609054](image/README/1732239609054.png)

![1732239630984](image/README/1732239630984.png)

![1732239646866](image/README/1732239646866.png)

O Programa estimou três populações evolutivamente separadas entre PE e CE e entre PB e PE mas as probabilidades permaneceram baixas

## 6. Estimativas genéticas, Análise de variância molecular (AMOVA) e Redes de haplótipos


   Foi estimada a diversidade genética do Procerathrophys boiei utilizando o script (/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/DnaSP/Script_Python_statistics***Estatistica_Sumaria_final.py***)

   ![1732646575278](image/README/1732646575278.png)

   O resultado está em: /home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/DnaSP/Script_Python_statistics/***Tabela_formatada_1.pdf***


#### Rede de Haplótipos:
   O script ***fas2phylip.ipynb*** foi utilizado para converter as sequências de alinhamento do formato .fas para phylip .phy. Um arquivo de Traits utilizado para inserção no Popart foi criado através do script ***Create_traits_file.ipynb*** onde foram testadas a hipótese nula de que as populações de águas azul, Tapacurá e Maranguape (***norte***) não compartilham haplótipos com as populações de Brejo Madre de Deus, Brejo dos Cavalos e Frei Caneca (***sul***). O software ***PopArt*** foi utilizado para criar o mapa de haplótipos utilizando o algoritmo ***Median Joining Network***. As sequências de alinhamento no formato .phy e o arquivo de Traits para criar a rede de haplótipos para os três pedaços do gene ***Cmyc*** (***2 exons*** e ***1 Intron***) e para os dois genes mitocondriais (***NAD2*** e ***CYTB***).

   ![plot_haplotipo_MITO](images/MITOS_haplotip_haplotype_network.jpg).

   ![plot_haplotipo_NUC](images/Cmyc_haplotip_haplotype_network.jpg).




## 7. Demografia Histórica e Difusão espaço temporal

#### Mismatch distriutions (Análise da Distribuição das Diferenças ):

Foi utilizado o ***DnaSP*** para realizar a análise de Distribuição das Diferenças.

![1732721210226](image/README/1732721210226.png)


Os genes mitocondriais do ***Procerathropys boiei*** mostram uma possível expansão populacional antiga com diversos eventos de redução e aumento nas diferênças par-a-par. Porém, os genes Nucleares mostram que a expansão populacional mas com uma pequena redução recente, isso pode ser um reflexo do conservadorismo do gene nuclear.


#### Skylineplot

utilizamos os softwares ***BEASTX*** e o ***TRACER*** para verificar se existe concordância das análises de expansão populacional através da análise de Skylinegrid do gene mitocondrial CYTB utilizando o modelo de evolução do relógio molecular como Strict.

![skyline_plot](images/SKYLINE_PLOT.jpg)


### Reconstrução espaço-temporal

Utilizamos o ***BEAUTI***, ***BEAST***, ***TREEANNOTATOR*** e o ***SPREAD*** (Versão 1.0.7rc, 2013) para realizar a reconstrução temporal do gene ***Cytb***. Para isso foi necessário um arquivo dos pontos com as coordenadas e um arquivo contendo as sequências alinhadas. o ***BEAUTI*** e o ***BEAST*** foi utilizado para aplicar o método bayesiano para estimar a história filogenética, o ***TREEANNOTATOR*** foi utilizado para estimar a árvore de máxima credibilidade criada pela análise bayesiana que em seguida foi carregada no ***SPREAD*** para gerar o arquivo ***.KML*** para visualização no ***Google Earth Pro*** (versão 7.3.6.9345).
 

#

A animação a seguir mostra a difusão espaço-temporal do gene ***Cytb***  na espécies Procerathrophys boiei ao longo do tempo, visualizada no Google Earth Pro.

![Animação_da_Filogeografia](images/tenor.gif)


