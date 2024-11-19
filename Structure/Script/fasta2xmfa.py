import subprocess
import sys

try:
    import Bio
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])

import os
from Bio import SeqIO

# Função para processar os arquivos FASTA da pasta e gerar o arquivo para xmfa2struct
def generate_xmfa2struct_input(folder_path, output_file):
    all_sequences = {}

    # Listar todos os arquivos FASTA na pasta
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith(".fas")]
    
    # Processar todos os arquivos FASTA na pasta
    for fasta_file in fasta_files:
        gene_name = os.path.basename(fasta_file).split('.')[0]  # Nome do gene baseado no nome do arquivo
        fasta_file_path = os.path.join(folder_path, fasta_file)
        sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

        # Criar um dicionário para armazenar as sequências do gene por indivíduo
        for seq in sequences:
            if seq.id not in all_sequences:
                all_sequences[seq.id] = {}

            all_sequences[seq.id][gene_name] = str(seq.seq)

    # Criar o arquivo de saída com o formato para xmfa2struct
    with open(output_file, 'w') as out_file:
        for gene in fasta_files:
            gene_name = os.path.basename(gene).split('.')[0]

            # Escrever o nome do gene no cabeçalho
            out_file.write(f"#{gene_name}\n")

            # Obter o tamanho total do gene a partir do primeiro arquivo FASTA (todos têm o mesmo tamanho)
            gene_length = len(list(SeqIO.parse(os.path.join(folder_path, gene), "fasta"))[0].seq)

            # Ordenar os indivíduos por nome
            sorted_individuals = sorted(all_sequences.keys())

            # Escrever as sequências de cada indivíduo
            for individual in sorted_individuals:
                # Garantir que todos os genes sejam representados, mesmo com dados ausentes
                gene_sequence = all_sequences[individual].get(gene_name, None)
                
                if gene_sequence is None:
                    # Se o gene não está presente para o indivíduo, preencher com "-" do tamanho total do gene
                    gene_sequence = '-' * gene_length

                # Escrever o cabeçalho de cada sequência com o sinal ">"
                out_file.write(f">{individual}\n{gene_sequence}\n")
            
            # Final do bloco de gene
            out_file.write("=\n\n")

# Caminho para a pasta contendo os arquivos FASTA
folder_path = "/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Structure/Fasta_files/"  # Substitua pelo caminho correto

# Caminho de saída do arquivo para xmfa2struct
output_file = "/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/Structure/Result_2xmfa/xmfa_input.fas"

# Gerar o arquivo de entrada para xmfa2struct
generate_xmfa2struct_input(folder_path, output_file)
