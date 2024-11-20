## Script que irá transformar alinhamentos em formato fasta em matriz de haplótipos no formato de input do Geneland
## Fornecer ao script o número de caracteres do sufixo ao final nome de cada sequência denotando os alelos
## Por exemplo: >IND1_1 (2 caracteres: "_1", denotando o primeiro alelo do IND1); ou >IND1a (1 caracter: "a", denotando o primeiro alelo do IND1)
## O Script precisa identificar o nome do individuo em comum aos dois aldos (sem o sufixo determinando cada fase) para criar a matriz final de forma correta

import os
import sys
import subprocess

# Função para verificar e instalar pacotes
def check_and_install_packages(packages):
    for package in packages:
        try:
            __import__(package)
        except ImportError:
            print(f"Pacote {package} não encontrado. Instalando...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        else:
            print(f"Pacote {package} já está instalado.")

# Verificar e instalar pacotes necessários
required_packages = ["Bio", "pandas"]
check_and_install_packages(required_packages)

# Importar pacotes após a verificação
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

def detect_suffixes(sequences, suffix_length):
    """Detecta os sufixos únicos com base no comprimento fornecido."""
    suffixes = set(seq.id[-suffix_length:] for seq in sequences)
    if len(suffixes) != 2:
        raise ValueError(f"Esperados exatamente 2 sufixos únicos, mas encontrados: {suffixes}")
    return list(suffixes)

def process_fasta_file(fasta_file, suffix_length, allele_suffixes):
    """Processa um arquivo FASTA e gera a matriz de haplótipos."""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    haplotype_dict = {}
    haplotype_matrix = []
    individual_map = defaultdict(lambda: [None, None])
    haplotype_counter = 1

    for seq_record in sequences:
        individual_id = seq_record.id[:-suffix_length]
        allele = str(seq_record.seq)

        if allele not in haplotype_dict:
            haplotype_dict[allele] = haplotype_counter
            haplotype_counter += 1

        if seq_record.id.endswith(allele_suffixes[0]):
            individual_map[individual_id][0] = haplotype_dict[allele]
        elif seq_record.id.endswith(allele_suffixes[1]):
            individual_map[individual_id][1] = haplotype_dict[allele]

    for individual_id, (haplotype_a, haplotype_b) in individual_map.items():
        haplotype_matrix.append([
            individual_id,
            haplotype_a if haplotype_a is not None else -999,
            haplotype_b if haplotype_b is not None else -999
        ])

    return haplotype_matrix

def process_all_fasta_files(folder_path, suffix_length):
    """Processa todos os arquivos FASTA na pasta e gera o DataFrame final."""
    all_haplotypes = []
    all_individuals = set()
    allele_suffixes = None
    gene_names = []  # Lista para armazenar nomes de genes (baseados nos arquivos)

    for fasta_file in os.listdir(folder_path):
        if fasta_file.endswith((".fasta", ".fas")):
            print(f"Processando o arquivo: {fasta_file}")
            fasta_file_path = os.path.join(folder_path, fasta_file)
            sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

            if allele_suffixes is None:
                allele_suffixes = detect_suffixes(sequences, suffix_length)
                print(f"Sufixos detectados: {allele_suffixes}")

            haplotype_matrix = process_fasta_file(fasta_file_path, suffix_length, allele_suffixes)
            all_haplotypes.append(haplotype_matrix)

            # Extrair o nome do gene do nome do arquivo
            gene_name = os.path.splitext(fasta_file)[0]
            gene_names.append(gene_name)

            for row in haplotype_matrix:
                all_individuals.add(row[0])

    final_data = []
    for individual in sorted(all_individuals):
        row = [individual]
        for haplotype_matrix in all_haplotypes:
            for record in haplotype_matrix:
                if record[0] == individual:
                    row.extend([int(record[1]), int(record[2])])
                    break
            else:
                row.extend([-999, -999])
        final_data.append(row)

    columns = ["Individual"]
    for gene in gene_names:
        columns.extend([f"{gene}_1", f"{gene}_2"])
    
    df = pd.DataFrame(final_data, columns=columns)
    
    return df

# Caminho da pasta e tamanho do sufixo especificado pelo usuário
folder_path = "/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Geneland/nuc_files"  # Substitua pelo caminho correto
suffix_length = int(input("Digite o número de caracteres do sufixo: "))

final_df = process_all_fasta_files(folder_path, suffix_length)

print(final_df.head())

# Salvar os resultados
final_df.to_csv("Proceratophrys_boiei_haplotype_matrix_sorted_nuclear_with_individuals_NEW.txt", sep="\t", index=False)
final_df.drop(columns=["Individual"]).to_csv("Proceratophrys_boiei_matrix_sorted_nuclear_no_individual_column_NEW.txt", sep="\t", index=False, header=False)
