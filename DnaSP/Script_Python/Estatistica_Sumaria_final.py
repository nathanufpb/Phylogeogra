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
    
    # Verifica se todas as sequências têm o mesmo comprimento
    seq_len = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_len:
            raise ValueError(f"All sequences must have the same length. Sequence {seq.id} is not of length {seq_len}.")

def count_segregating_sites(alignment):
    """Conta os sítios segregantes corretamente, ignorando gaps."""
    segregating_sites = 0
    for col in range(alignment.shape[1]):
        column_values = set(alignment[:, col])
        column_values.discard(4)  # Remove gaps (representados como 4)
        if len(column_values) > 1:  # Sítios variáveis
            segregating_sites += 1
    return segregating_sites

def __segregating(sequences: list) -> int:
    """
    Counts the number of segregating sites.
    For each position in the sequences, check if the character differs
    from the rest, indicating a segregating site.
    """
    seg_sites = 0
    # For each position in sequence
    for i in range(len(sequences[0])):
        s = sequences[0][i]  # Take the char of first sequence as reference
        for seq in sequences:  # For each other sequence
            if seq[i] != s:  # Segregating site if different
                seg_sites += 1  # Add 1 if not equal
                break  # And stop comparing this position
    return seg_sites

def group_alleles(sequences):
    """Agrupa sequências por indivíduo com base nos sufixos (_1, _2 ou -1, -2)."""
    grouped = {}
    for seq in sequences:
        individual_id = seq.id.rsplit('_', 1)[0].rsplit('-', 1)[0]
        if individual_id not in grouped:
            grouped[individual_id] = []
        grouped[individual_id].append(str(seq.seq))

    filtered_grouped = {ind: alleles for ind, alleles in grouped.items() if len(alleles) == 2}

    if len(filtered_grouped) < len(grouped):
        print(f"Aviso: Alguns indivíduos têm ploidia incorreta e foram ignorados.")
    return filtered_grouped

def calculate_haplotype_diversity(haplotype_counts, sample_size):
    """Calcula a diversidade haplotípica."""
    if sample_size <= 1:
        return 0  # Diversidade haplotípica não definida para amostras pequenas
    haplotype_frequencies = haplotype_counts / sample_size
    haplotype_diversity = (sample_size / (sample_size - 1)) * (1 - np.sum(haplotype_frequencies ** 2))
    return haplotype_diversity

def calculate_pairwise_differences(sequences):
    """Calcula as diferenças nucleotídicas entre todos os pares de sequências."""
    total_differences = 0
    total_pairs = 0
    
    for seq1, seq2 in itertools.combinations(sequences, 2):
        differences = sum(nt1 != nt2 for nt1, nt2 in zip(seq1, seq2) if nt1 != '-' and nt2 != '-')
        total_differences += differences
        total_pairs += 1
    
    return total_differences, total_pairs

def pi_estimator(sequences: list, safe=True) -> float:
    """
    Computes Θ_π, the Pi estimator.
    Θ_π = Number of pairwise differences / binomial(n, 2)

    Parameters
    ----------
        sequences: List[str]
            List of sequences.
        safe: bool, default: True
            Check if sequences have the same length.

    Returns
    -------
        Θ_π: float
            Pi estimator value.
    """
    if safe:
        __check(sequences)

    pairwise = itertools.combinations(sequences, 2)
    cs = [sum([not charA == charB for charA, charB in zip(seqA, seqB)]) for seqA, seqB in pairwise]
    n = len(sequences)
    binomial = ((n - 1) * n) / 2  # Binomial(n, 2)

    return sum(cs) / binomial

def __harmonic(n: int) -> float:
    """
    Computes the n-1th harmonic number.

    Parameters
    ----------
        n: int

    Returns
    -------
        h: float
            N-th harmonic number
    """
    return sum([1 / i for i in range(1, n)])

def calculate_theta_from_s(segregating_sites, sample_size):
    """Calcula Theta-W diretamente a partir dos sítios segregantes."""
    a1 = sum(1 / i for i in range(1, sample_size))
    return segregating_sites / a1


def calculate_theta(segregating_sites, sample_size, fragment_length): # Theta-W normalizado por fragment_length
    a1 = sum(1 / i for i in range(1, sample_size))
    return (segregating_sites / a1) / fragment_length


def tajimas_d(sequences: list) -> float:
    """
    Computes Tajima's D for a list of sequences.
    See: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203831/
    and https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-
    genomics-fall-2005/study-materials/tajimad1.pdf

    Parameters
    ----------
        sequences: List[str]
            List of sequences.
    Returns
    -------
        Θ_D: float
            Tajima's D value.
    """
    __check(sequences)

    seg_sites = __segregating(sequences)  # Number of segregating sites

    # Prevent division by 0
    if seg_sites == 0:
        return 0

    theta_pi = pi_estimator(sequences, safe=False)  # Pi Estimator

    num_seq = len(sequences)  # Number of sequences
    harmonic = __harmonic(num_seq)  # N-1th harmonic number

    harmonic = sum([1 / i for i in range(1, num_seq)])  # Ref 3, aka. a1
    a2 = sum([1 / (i**2) for i in range(1, num_seq)])  # Ref 4

    b1 = (num_seq + 1) / (3 * (num_seq - 1))  # Ref 8
    b2 = (2 * (num_seq**2 + num_seq + 3)) / (9 * num_seq * (num_seq - 1))  # Ref 9

    c1 = b1 - 1 / harmonic
    c2 = b2 - ((num_seq + 2) / (harmonic * num_seq)) + (a2 / (harmonic**2))

    e1 = c1 / harmonic
    e2 = c2 / ((harmonic**2) + a2)

    delta_Theta = theta_pi - (seg_sites / harmonic)  
    tD = delta_Theta / (((e1 * seg_sites) + (e2 * seg_sites * (seg_sites - 1))) ** 0.5)  # Ref 27
    return float(tD)

def calculate_summary_statistics(fasta_file, population):
    loci_name = os.path.splitext(os.path.basename(fasta_file))[0]
    loci_name = loci_name.replace("mt_", "").replace("nu_", "")

    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if len(sequences) == 0:
        raise ValueError(f"O arquivo {fasta_file} está vazio ou mal formatado.")

    ploidy = 1 if "mt" in fasta_file.lower() else 2

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

    unique_haplotypes = np.unique(alignment, axis=0)
    haplotype_count = unique_haplotypes.shape[0]

    haplotype_counts = np.unique(alignment, axis=0, return_counts=True)[1]
    haplotype_diversity = calculate_haplotype_diversity(haplotype_counts, sample_size)

    total_differences, total_pairs = calculate_pairwise_differences(["".join(seq) for seq in alignment])
    pi = total_differences / (total_pairs * fragment_length)

    segregating_sites = count_segregating_sites(alignment_numeric)
    
    # Ajuste do sample_size para ploidia diploide
    adjusted_sample_size = sample_size // 2 if ploidy == 2 else sample_size

    # Cálculo de Theta-W e D de Tajima
    theta_w = calculate_theta(segregating_sites, adjusted_sample_size, fragment_length)

    # Cálculo adicional: Theta-W a partir dos sítios segregantes
    theta_w_from_s = calculate_theta_from_s(segregating_sites, adjusted_sample_size)

    # Cálculo de Tajima's D
    tajima_d_value = tajimas_d(sequences)  # Chamando a função para calcular Tajima's D

    return {
        "Loci": loci_name,
        "Population": population,
        "Length (bp)": fragment_length,
        "Sample Size (N)": sample_size,
        "Segregating Sites (S)": segregating_sites,
        "Number of Haplotypes (H)": haplotype_count,
        "Haplotypes Diversity (Hd)": f"{haplotype_diversity:.4f}",
        "Nucleotide Diversity (pi)": f"{pi:.5f}",
        "Tajima's D": f"{tajima_d_value:.4f}",  # Exibindo o valor de Tajima's D
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
base_folder = "/Users/felipemedeiros/Dropbox/Disciplinas UFPB/Filogeografia/Filogeografia_Felipe_Camura/Dados_alunos/Dia3/DnaSP_estimativas genéticas/DNAsp"

try:
    summary_df = generate_summary_for_all_populations(base_folder)
    summary_df.to_csv("summary_statistics_output.csv", index=False)
    print("As estatísticas sumárias estão salvas em summary_statistics_output.csv")
except Exception as e:
    print(f"Erro: {e}")
