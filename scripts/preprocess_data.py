# Clean and transform data for analysis and model training
from pathlib import Path
import pandas as pd 

# Read the CSV file (metadata) 
metadata = pd.read_csv("~/data/processed_data/metadata/metadata.csv")
print(metadata)

# Filter data based on quality metrics
def filter_quality(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df['gcPercent'].between(30, 70)]
    df = df[df['numberOfContigs'] < 500]
    df = df[df['assemblyLevel'].isin(['Complete Genome', 'Chromosome'])]  # accepted levels
    return df

# Keep complete and non-redundant samples
def deduplicate_samples(df: pd.DataFrame) -> pd.DataFrame:
    df = df.drop_duplicates(subset=['key'])
    df = df.dropna(subset=['resistance_label'])  # enlever les échantillons sans label
    return df

# Convert FASTA format into analyzable formats for ML
def convert_fasta():
    pass

# sequence into numerical vectors nucleic acids or amino acids
def fasta_to_vectors():
    pass

# One-hot encoding:
def encode_sequences():
    pass

# k-mer encoding:
def kmer_encoding():
    pass

# embedding-based encoding:
def embedding_encoding():
    pass

# generate tab datasets (X.npy, y.csv) for model training
def generate_datasets():
    pass

# save everything into processed_data/features/
def save_processed_data():
    pass