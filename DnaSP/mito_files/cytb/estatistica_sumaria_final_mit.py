## Organizar um diretório contendo subpastas com os alinhamentos por população (o nome do subdiretório deve ser o nome de suas populações)
## Os genes mitocondriais devem ter o prefixo "mt_" e nuclear o prefixo "nu_". Exemplo: mt_16S.fasta; nu_RAG1.fasta

import os
import sys
import subprocess
from Bio import SeqIO
import numpy as np
import pandas as pd
from importlib.metadata import version, PackageNotFoundError
import itertools  # Garantindo a importação de itertools

# Função para garantir que os pacotes necessários estão instalados
def install_packages(packages):
    missing = []
    for package in packages:
        try:
            version(package)
        except PackageNotFoundError:
            missing.append(package)
    
    if missing:
        print(f"Instalando pacotes ausentes: {', '.join(missing)}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", *missing])
    else:
        print("Todos os pacotes necessários já estão instalados.")

# Pacotes necessários
required_packages = ["biopython", "numpy", "pandas"]
install_packages(required_packages)

def __check(sequences: list) -> None:
    """
    Check for valid sequence proportions (1 < n, sample length).
    Ensures all sequences have the same length.
    """
    if len(sequences) < 2:
        raise ValueError("At least 2 sequences are required.")
    
    

def calculate_summary_statistics(fasta_file, population):
    # Determinar o tipo de gene baseado no prefixo
    loci_name = os.path.splitext(os.path.basename(fasta_file))[0]
    loci_name = loci_name.replace("mt_", "").replace("nu_", "")

    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if len(sequences) == 0:
        raise ValueError(f"O arquivo {fasta_file} está vazio ou mal formatado.")

    # Diferenciar ploidia por tipo de gene
    if fasta_file.startswith("mt_"):
        ploidy = 1  # Mitocondrial
    elif fasta_file.startswith("nu_"):
        ploidy = 2  # Nuclear
    else:
        raise ValueError(f"Arquivo {fasta_file} não tem o prefixo esperado ('mt_' ou 'nu_').")

    if ploidy == 2:
        grouped_sequences = group_alleles(sequences)
        if not grouped_sequences:
            raise ValueError(f"Nenhum indivíduo válido encontrado em {fasta_file}.")
        alignment = np.array([list(allele) for alleles in grouped_sequences.values() for allele in alleles])
        sample_size = alignment.shape[0]
    else:
        alignment = np.array([list(str(seq.seq)) for seq in sequences])
        sample_size = alignment.shape[0]

    nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}
    alignment_numeric = np.vectorize(nucleotide_map.get)(alignment)

    fragment_length = alignment_numeric.shape[1]

    # Cálculos compartilhados
    unique_haplotypes = np.unique(alignment, axis=0)
    haplotype_count = unique_haplotypes.shape[0]

    haplotype_counts = np.unique(alignment, axis=0, return_counts=True)[1]
    haplotype_diversity = calculate_haplotype_diversity(haplotype_counts, sample_size)

    total_differences, total_pairs = calculate_pairwise_differences(["".join(seq) for seq in alignment])
    pi = total_differences / (total_pairs * fragment_length)

    segregating_sites = count_segregating_sites(alignment_numeric)

    # Ajustar tamanho da amostra para ploidia
    adjusted_sample_size = sample_size // 2 if ploidy == 2 else sample_size

    # Cálculos específicos
    theta_w = calculate_theta(segregating_sites, adjusted_sample_size, fragment_length)
    theta_w_from_s = calculate_theta_from_s(segregating_sites, adjusted_sample_size)
    tajima_d_value = tajimas_d(sequences)

    return {
        "Loci": loci_name,
        "Population": population,
        "Length (bp)": fragment_length,
        "Sample Size (N)": sample_size,
        "Segregating Sites (S)": segregating_sites,
        "Number of Haplotypes (H)": haplotype_count,
        "Haplotypes Diversity (Hd)": f"{haplotype_diversity:.4f}",
        "Nucleotide Diversity (pi)": f"{pi:.5f}",
        "Tajima's D": f"{tajima_d_value:.4f}",
        "Theta": f"{theta_w:.4f}",
        "Theta-W (from S)": f"{theta_w_from_s:.4f}",
        "Average Nucleotide Differences": f"{total_differences / total_pairs:.4f}"
    }

def process_all_files(folder_path):
    summary_data = []
    population = os.path.basename(folder_path)
    
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith(".fasta")]
    if not fasta_files:
        print(f"Aviso: Nenhum arquivo FASTA encontrado na pasta '{folder_path}'")
        return pd.DataFrame()

    for fasta_file in fasta_files:
        file_path = os.path.join(folder_path, fasta_file)
        try:
            stats = calculate_summary_statistics(file_path, population)
            summary_data.append(stats)
        except ValueError as e:
            print(f"Aviso: {e}")

    return pd.DataFrame(summary_data)

def generate_summary_for_all_populations(base_folder):
    all_data = []
    for population_folder in os.listdir(base_folder):
        folder_path = os.path.join(base_folder, population_folder)
        if os.path.isdir(folder_path):
            population_df = process_all_files(folder_path)
            if not population_df.empty:
                all_data.append(population_df)
        else:
            print(f"Ignorando '{population_folder}' pois não é uma pasta.")
    
    if not all_data:
        raise ValueError("Nenhum dado foi encontrado para processar.")
    
    final_df = pd.concat(all_data, ignore_index=True)

    mito_genes = final_df[final_df['Loci'].str.contains('mt', case=False)]
    nuclear_genes = final_df[~final_df['Loci'].str.contains('mt', case=False)].sort_values(by='Loci')

    ordered_df = pd.concat([mito_genes, nuclear_genes])
    ordered_df = ordered_df.sort_values(by=['Loci', 'Population'])

    return ordered_df

# Caminho base contendo subpastas de populações
base_folder = "/home/nathan/Documents/Doutorado_diversidade_genética/Phylogeogra/DnaSP/mito_files/cytb"

try:
    summary_df = generate_summary_for_all_populations(base_folder)
    summary_df.to_csv("summary_statistics_output.csv", index=False)
    print("As estatísticas sumárias estão salvas em summary_statistics_output.csv")
except Exception as e:
    print(f"Erro: {e}")