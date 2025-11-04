# Clean and transform data for analysis and model training
import pandas as pd 

# Read the CSV file (metadata) 
metadata = pd.read_csv("~/data/processed_data/metadata/metadata.csv")
print(metadata)

# Compute relevant statistics from metadata
def compute_quality_stats(df: pd.DataFrame) -> pd.DataFrame:
    df['gc_percent'] = (df['gc_content'] / df['genome_length']) * 100
    return df

# Compute contig count
def compute_contig_stats(df: pd.DataFrame) -> pd.DataFrame:
    df['contig_count'] = df['num_contigs']
    return df

# Compute assembly level
def compute_assembly_level(df: pd.DataFrame) -> pd.DataFrame:
    level_map = {1: 'Complete Genome', 2: 'Chromosome', 3: 'Scaffold', 4: 'Contig'}
    df['assembly_level'] = df['assembly_level_code'].map(level_map)
    return df

# Filter data based on quality metrics
def filter_quality(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df['gc_percent'].between(30, 70)]  # garder GC% entre 30 et 70
    df = df[df['contig_count'] < 500]          # garder contig_count < 500
    df = df[df['assembly_level'].isin(['Complete Genome', 'Chromosome'])]  # niveaux acceptés
    return df

# Keep complete and non-redundant samples
def deduplicate_samples():

# Convert FASTA format into analyzable formats for ML
def convert_fasta():

# sequence into numerical vectors nucleic acids or amino acids
def fasta_to_vectors():


# One-hot encoding:
def encode_sequences():

# k-mer encoding:
def kmer_encoding():

# embedding-based encoding:
def embedding_encoding():

# generate tab datasets (X.npy, y.csv) for model training
def generate_datasets():

# save everything into processed_data/features/
def save_processed_data():
