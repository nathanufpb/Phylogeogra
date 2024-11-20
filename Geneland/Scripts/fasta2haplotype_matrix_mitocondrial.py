import os
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

# Função para processar cada arquivo FASTA e gerar a matriz de haplótipos para o gene
def process_fasta_file(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    haplotype_dict = {}  # Dicionário para contar haplótipos
    haplotype_matrix = []
    individual_map = {}  # Dicionário para mapear o ID do indivíduo aos alelos

    haplotype_counter = 1  # Contador para gerar números de haplótipos únicos

    # Agrupar sequências haploides por indivíduo
    for seq_record in sequences:
        individual_id = seq_record.id  # ID do indivíduo (não há sufixo 'a' ou 'b' aqui)
        allele = str(seq_record.seq)  # Converter a sequência em string
        
        # Verificar se o alelo já foi registrado no dicionário
        if allele not in haplotype_dict:
            haplotype_dict[allele] = haplotype_counter  # Atribui um número único ao haplótipo
            haplotype_counter += 1  # Incrementa o contador de haplótipos

        # Mapear o haplótipo para o indivíduo
        individual_map[individual_id] = haplotype_dict[allele]  # Armazenar o número do haplótipo

    # Preencher a matriz final com os dados
    for individual_id, haplotype in individual_map.items():
        haplotype_matrix.append([individual_id, haplotype])

    return haplotype_matrix

# Função para processar todos os arquivos FASTA na pasta e gerar a matriz concatenada
def process_all_fasta_files(folder_path):
    all_haplotypes = []  # Lista para armazenar todas as matrizes de haplótipos
    gene_names = []  # Lista para armazenar os nomes dos genes
    all_individuals = set()  # Conjunto para armazenar todos os indivíduos únicos

    # Iterar por todos os arquivos FASTA na pasta
    for fasta_file in os.listdir(folder_path):
        if fasta_file.endswith(".fasta") or fasta_file.endswith(".fas"):
            gene_name = os.path.splitext(fasta_file)[0]  # Obtém o nome do gene sem a extensão
            gene_names.append(gene_name)
            print(f"Processando o arquivo: {fasta_file}")
            
            # Caminho completo do arquivo
            fasta_file_path = os.path.join(folder_path, fasta_file)
            
            # Processar o arquivo e gerar a matriz de haplótipos
            haplotype_matrix = process_fasta_file(fasta_file_path)

            if haplotype_matrix:
                all_haplotypes.append(haplotype_matrix)
                all_individuals.update(row[0] for row in haplotype_matrix)

    # Criar o DataFrame final com todos os indivíduos
    final_data = []
    for individual in sorted(all_individuals):
        row = [individual]
        for haplotype_matrix in all_haplotypes:
            # Procurar os haplótipos do indivíduo para cada gene
            for record in haplotype_matrix:
                if record[0] == individual:
                    row.append(record[1])  # Adicionar o haplótipo do gene
                    break
            else:
                # Se o indivíduo não tiver dados para o gene, adicionar -999
                row.append(-999)  # Adicionar -999 para dados faltantes
        final_data.append(row)

    # Criar as colunas com base nos nomes dos genes
    columns = ["Individual"] + [f"{gene}" for gene in gene_names]
    df = pd.DataFrame(final_data, columns=columns)

    # Ordenar o DataFrame pela coluna 'Individual' em ordem alfabética
    df = df.sort_values(by="Individual", ascending=True)

    return df

# Caminho para a pasta contendo os arquivos FASTA
folder_path = "/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/mito_files"  # Substitua pelo caminho correto

# Processar todos os arquivos FASTA na pasta e gerar o DataFrame final
final_df = process_all_fasta_files(folder_path)

# Exibir o DataFrame final
print(final_df.head())

# 1. Gerar o primeiro output (com cabeçalho e coluna de indivíduos)
final_df.to_csv("Proceratophrys_boiei_haplotype_matrix_sorted_mitocondrial_with_individuals.txt", sep="\t", index=False)

# 2. Gerar o segundo output (sem coluna de indivíduos e sem cabeçalho)
final_df.drop(columns=["Individual"]).to_csv("Proceratophrys_boiei_haplotype_matrix_sorted_mitocondrial_no_individual_column.txt", sep="\t", index=False, header=False)
