# Filogeografia
 repository for phlogeographic analysis


A análise do Geneland começou criando pastas para os genes nucleares phaseados e separados em 2 exons e 1 intron, e uma pasta para os dois outros genes mitocondriais.

os scripts em python dentro da pasta Geneland/Scripts (fasta2haplotype_matrix_mitocondrial.py e fasta2haplotype_matrix_nuclear_gen.py)
foram utilizados para transformar os dados .fas em uma matriz haplotípica.

o script "Configure_txt_file_and_create_coordenate_file.ipynb" foi utilizado para criar um novo arquivo da matriz haplotípica dos genes mitocondriais contendo apenas as linhas que possuiam o mesmo código do gene nuclear(42 genes), e na construção de um arquivo contendo as coordenadas de cada linha mapeada para os códigos do arquivo nuclear.

o arquivo novo dos genes mitocondriais com o mesmo número de amostras e o mesmo código do arquivo com genes nucleares:

mt.temp = filtered_mitochondrial_data_no_header.txt

geno.temp = Proceratophrys_boiei_matrix_sorted_nuclear_no_individual_column_NEW.txt

coord.temp = Proceratophrys_boiei_haplotype_matrix_sorted_nuclear_with_individuals_NEW_geo.txt

maindir = /home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/Scripts/Geneland_analysis/

com esses arquivos foi possível executar o Geneland através dp Rscript "Geneland_Cnemis.R"

![Geneland Analysis Plot](https://github.com/nathanufpb/Phylogeogra/blob/main/Geneland/Scripts/Geneland_analysis/noadmix/Rplot01.png)