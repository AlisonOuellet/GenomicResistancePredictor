# Clean and transform data for analysis and model training
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from typing import List, Dict, Iterable, Tuple, Optional

try:
    from sklearn.feature_extraction.text import TfidfVectorizer
    from sklearn.decomposition import TruncatedSVD
    _HAS_SK = True
except Exception:
    _HAS_SK = False

# Utilities
DNA_ALPHABET = "ACGTN" # nucleic acids + N (unknown)
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY" # 20 standard amino acids

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[1]   # scripts/ -> racine du repo
DATA = PROJECT_ROOT / "data"
PROC = DATA / "processed_data"
META_D = PROC / "metadata"
GENOME_D = PROC / "genomes"
CDS_D = PROC / "cds"
FEAT_D = PROC / "features"
META_CSV = META_D / "metadata.csv"

for p in [DATA, PROC, META_D, GENOME_D, CDS_D, FEAT_D]:
    p.mkdir(parents=True, exist_ok=True)

print(f"[INFO] DATA   = {DATA.resolve()}")
print(f"[INFO] FEAT_D = {FEAT_D.resolve()}")
print(f"[INFO] META_CSV = {META_CSV.resolve()} (exists={META_CSV.exists()})")


# Helpers
def _clean_dna(seq: str) -> str:
    s = seq.upper().replace("U", "T")
    return "".join(ch if ch in {"A", "C", "G", "T", "N"} else "N" for ch in s)

def _clean_aa(seq: str) -> str:
    s = seq.upper().replace("U", "C")
    return "".join(ch if ch in set(AA_ALPHABET) else "X" for ch in s)

def _iter_kmers(seq: str, k: int) -> Iterable[str]:
    for i in range(0, len(seq) - k + 1):
        yield seq[i:i+k]

def _kmer_counts(seq: str, k: int, is_dna: bool = True) -> Dict[str, int]:
    seq = _clean_dna(seq) if is_dna else _clean_aa(seq)
    skip_char = "N" if is_dna else "X"
    bag: Dict[str, int] = {}
    for kmer in _iter_kmers(seq, k):
        if skip_char in kmer:
            continue
        bag[kmer] = bag.get(kmer, 0) + 1
    return bag

def _vectorize_counts(counts_list: List[Dict[str, int]]) -> Tuple[np.ndarray, List[str]]:
    vocab: Dict[str, int] = {}
    for d in counts_list:
        for k in d:
            if k not in vocab:
                vocab[k] = len(vocab)
    V = len(vocab)
    X = np.zeros((len(counts_list), V), dtype=np.float32)
    for i, d in enumerate(counts_list):
        for k, c in d.items():
            X[i, vocab[k]] = c
    # Normalisation L1
    denom = X.sum(axis=1, keepdims=True) + 1e-8
    X = X / denom
    inv_vocab = [None] * V
    for k, j in vocab.items():
        inv_vocab[j] = k
    return X, inv_vocab

def _one_hot_sequence(seq: str, alphabet: str, max_len: int) -> np.ndarray:
    seq = seq.upper()
    idx = {ch: i for i, ch in enumerate(alphabet)}
    L = min(len(seq), max_len)
    X = np.zeros((max_len, len(alphabet)), dtype=np.float32)
    for i in range(L):
        j = idx.get(seq[i])
        if j is not None:
            X[i, j] = 1.0
    return X


# Read the CSV file (metadata) 
def read_csv(metadata_path: Path) -> pd.DataFrame:
    df = pd.read_csv(metadata_path)
    df.columns = [c.strip() for c in df.columns]
    return df

# Filter data based on quality metrics
def filter_quality(df: pd.DataFrame) -> pd.DataFrame:
    if "genome_size_bp" in df.columns:
        df = df[pd.to_numeric(df["genome_size_bp"], errors="coerce").between(4_000_000, 6_500_000, inclusive="both")]
    if "gcPercent" in df.columns:
        df = df[pd.to_numeric(df["gcPercent"], errors="coerce").between(45.0, 55.0, inclusive="both")]
    if "numberOfContigs" in df.columns:
        df = df[pd.to_numeric(df["numberOfContigs"], errors="coerce") <= 200]
    if "assemblyLevel" in df.columns:
        keep_levels = {"complete genome", "chromosome", "scaffold"}
        df = df[df["assemblyLevel"].astype(str).str.lower().isin(keep_levels)]
    def _exists_or_empty(p):
        p = str(p or "")
        return (p == "") or Path(p).exists()
    if "genomic_fasta" in df.columns and "cds_fasta" in df.columns:
        df = df[df["genomic_fasta"].apply(_exists_or_empty) & df["cds_fasta"].apply(_exists_or_empty)]
    return df

# Keep complete and non-redundant samples
def deduplicate_samples(df: pd.DataFrame) -> pd.DataFrame:
    # priorité: GCF puis moins de contigs puis niveau d’assemblage
    rank_map = {"complete genome": 3, "chromosome": 2, "scaffold": 1, "contig": 0}
    def _score(row):
        src = 1 if str(row.get("chosen_source", "")).upper() == "GCF" else 0
        contigs = -float(row.get("numberOfContigs", 1e9) or 1e9)  # moins c'est mieux
        level = rank_map.get(str(row.get("assemblyLevel", "")).lower(), -1)
        return (src, contigs, level)
    df = df.copy()
    df["__score"] = df.apply(_score, axis=1)

    if "biosample" in df.columns and df["biosample"].notna().any():
        df = df.sort_values("__score", ascending=False).drop_duplicates(subset=["biosample"], keep="first")
    elif "key" in df.columns:
        df = df.sort_values("__score", ascending=False).drop_duplicates(subset=["key"], keep="first")
    else:
        df = df.sort_values("__score", ascending=False).drop_duplicates(subset=["accession"], keep="first")

    return df.drop(columns=["__score"])

# Convert FASTA format into analyzable formats for ML
def convert_fasta(metadata_csv: Path = META_CSV,
                  concat_contigs: bool = True,
                  out_path: Path = FEAT_D / "sequences_dna.parquet") -> Path:
    df = pd.read_csv(metadata_csv, dtype=str)
    rows = []
    for _, r in df.iterrows():
        acc = r.get("accession")
        fp = r.get("genomic_fasta")
        if not fp or not Path(fp).exists():
            continue
        seqs = []
        for i, rec in enumerate(SeqIO.parse(fp, "fasta"), start=1):
            s = _clean_dna(str(rec.seq))
            if concat_contigs:
                seqs.append(s)
            else:
                rows.append({"accession": acc, "segment_id": f"contig_{i}", "length": len(s), "sequence": s})
        if concat_contigs and seqs:
            joined = "".join(seqs)
            rows.append({"accession": acc, "segment_id": "genome", "length": len(joined), "sequence": joined})
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(out_path, index=False)
    print(f"[OK] sequences écrites -> {out_path} ({len(rows)} lignes)")
    return out_path

# sequence into numerical vectors nucleic acids or amino acids
def fasta_to_vectors(parquet_path: Path,
                     kind: str = "dna",
                     onehot_len: int = 10_000,
                     k: int = 5) -> Dict[str, np.ndarray]:
    """
    Charge un parquet (accession, sequence) et calcule:
      - onehot: [N, L, |alphabet|]
      - kmer_counts: [N, |V|] (normalisé L1)
      - kmer_docs: liste de strings de k-mers (pour TF-IDF)
    """
    df = pd.read_parquet(parquet_path)
    if "sequence" not in df.columns or "accession" not in df.columns:
        raise ValueError("Parquet doit contenir 'accession' et 'sequence'.")

    is_dna = (kind.lower() == "dna")
    alphabet = DNA_ALPHABET if is_dna else AA_ALPHABET
    seqs = df["sequence"].astype(str).tolist()
    accs = df["accession"].astype(str).tolist()

    # One-hot (clipping/padding)
    if is_dna:
        seqs_clean = [_clean_dna(s) for s in seqs]
    else:
        seqs_clean = [_clean_aa(s) for s in seqs]

    onehots = np.stack([_one_hot_sequence(s, alphabet, onehot_len) for s in seqs_clean], axis=0)

    # k-mer counts
    counts_list = [_kmer_counts(s, k=k, is_dna=is_dna) for s in seqs_clean]
    X_kmer, vocab = _vectorize_counts(counts_list)

    # docs of kmers (space-separated)
    docs = []
    skip_char = "N" if is_dna else "X"
    for s in seqs_clean:
        toks = [km for km in _iter_kmers(s, k) if skip_char not in km]
        docs.append(" ".join(toks))

    return {
        "accessions": np.array(accs),
        "onehot": onehots.astype(np.float32),      # [N, L, |alphabet|]
        "kmer": X_kmer.astype(np.float32),         # [N, |V|]
        "kmer_vocab": np.array(vocab, dtype=object),
        "kmer_docs": np.array(docs, dtype=object)  # pour TF-IDF ultérieur
    }

# One-hot encoding:
def encode_sequences(parquet_path: Path,
                     kind: str = "dna",
                     onehot_len: int = 10_000,
                     out_name: Optional[str] = None) -> Path:
    vecs = fasta_to_vectors(parquet_path, kind=kind, onehot_len=onehot_len, k=4)
    if out_name is None:
        out_name = f"onehot_{kind}_L{onehot_len}.npz"
    out_path = FEAT_D / out_name
    np.savez_compressed(out_path, X=vecs["onehot"], accessions=vecs["accessions"])
    print(f"[OK] one-hot {kind} → {out_path} shape={vecs['onehot'].shape}")
    return out_path

# k-mer encoding:
def kmer_encoding(parquet_path: Path,
                  kind: str = "dna",
                  k: int = 5,
                  out_name: Optional[str] = None) -> Path:
    vecs = fasta_to_vectors(parquet_path, kind=kind, onehot_len=1, k=k)  # onehot_len=1 pour éviter gros calcul inutile
    if out_name is None:
        out_name = f"kmer_{kind}_k{k}.npz"
    out_path = FEAT_D / out_name
    np.savez_compressed(out_path,
                        X=vecs["kmer"],
                        vocab=vecs["kmer_vocab"],
                        accessions=vecs["accessions"])
    print(f"[OK] k-mer {kind} (k={k}) → {out_path} shape={vecs['kmer'].shape} |V|={len(vecs['kmer_vocab'])}")
    return out_path

# embedding-based encoding:
def embedding_encoding(parquet_path: Path,
                       kind: str = "dna",
                       k: int = 6,
                       svd_dim: int = 256,
                       tfidf_max_features: Optional[int] = None,
                       tfidf_min_df: int = 1,
                       tfidf_max_df: float = 1.0,
                       out_name: Optional[str] = None) -> Path:
    if not _HAS_SK:
        raise RuntimeError("scikit-learn n'est pas disponible pour TF-IDF+SVD.")
    vecs = fasta_to_vectors(parquet_path, kind=kind, onehot_len=1, k=k)
    docs = vecs["kmer_docs"].tolist()

    vectorizer = TfidfVectorizer(
        analyzer="word",
        token_pattern=r"[^ ]+",
        lowercase=False,
        min_df=tfidf_min_df,
        max_df=tfidf_max_df,
        max_features=tfidf_max_features
    )
    X = vectorizer.fit_transform(docs)
    if X.shape[1] <= 1:
        # rien à factoriser
        emb = X.toarray().astype(np.float32)
    else:
        svd = TruncatedSVD(n_components=min(svd_dim, X.shape[1]-1), random_state=42)
        emb = svd.fit_transform(X).astype(np.float32)

    if out_name is None:
        out_name = f"embed_{kind}_k{k}_svd{svd_dim}.npz"
    out_path = FEAT_D / out_name
    np.savez_compressed(out_path, X=emb, accessions=vecs["accessions"])
    # Sauvegarde des artefacts (vocabulaire + composantes) pour la repro
    vocab_json = FEAT_D / f"embed_{kind}_k{k}_tfidf_vocab.json"
    import json
    with open(vocab_json, "w", encoding="utf-8") as f:
        json.dump(vectorizer.vocabulary_, f)
    if X.shape[1] > 1:
        svd_npz = FEAT_D / f"embed_{kind}_k{k}_svd.npz"
        np.savez_compressed(svd_npz,
                            components=getattr(svd, "components_", None),
                            explained_variance=getattr(svd, "explained_variance_", None))
    print(f"[OK] embedding {kind} (k={k}, svd={svd_dim}) → {out_path} shape={emb.shape}")
    return out_path

# generate tab datasets (X.npy, y.csv) for model training
def generate_datasets(kind: str = "dna") -> Dict[str, Path]:
    """
    Construit rapidement un petit paquet de sorties tabulaires:
      - sequences_dna.parquet (ou AA si tu veux étendre)
      - onehot_{kind}.npz
      - kmer_{kind}.npz
      - embed_{kind}.npz
      - samples_index.csv (mapping accessions → meta)
    """
    # 1) Charger/filtrer/dédupliquer metadata.csv
    df = read_csv(META_CSV)
    df_f = filter_quality(df)
    df_f = deduplicate_samples(df_f)
    # export index pour traçabilité
    idx_csv = FEAT_D / "samples_index.csv"
    cols_keep = [c for c in ["accession", "organism", "strain", "biosample",
                             "assemblyLevel", "numberOfContigs", "gcPercent",
                             "genomic_fasta", "cds_fasta"] if c in df_f.columns]
    df_f[cols_keep].to_csv(idx_csv, index=False)

    # 2) Convertir FASTA → parquet séquences (DNA)
    parquet = FEAT_D / "sequences_dna.parquet"
    convert_fasta(metadata_csv=META_CSV, concat_contigs=True, out_path=parquet)

    # 3) One-hot (optionnel, gros)
    onehot = encode_sequences(parquet, kind=kind, onehot_len=10_000)

    # 4) k-mer
    kmer = kmer_encoding(parquet, kind=kind, k=5)

    # 5) Embedding (TF-IDF+SVD)
    embed = None
    try:
        embed = embedding_encoding(parquet, kind=kind, k=6, svd_dim=256,
                                   tfidf_max_features=200_000, tfidf_min_df=1, tfidf_max_df=1.0)
    except RuntimeError as e:
        print(f"[WARN] Embedding sauté: {e}")

    # 6) y.csv (placeholder): par défaut, pas d'étiquettes => on met juste accession
    y_csv = FEAT_D / "y.csv"
    pd.DataFrame({"accession": df_f["accession"].tolist()}).to_csv(y_csv, index=False)

    return {
        "index_csv": idx_csv,
        "parquet": parquet,
        "onehot": onehot,
        "kmer": kmer,
        "embed": embed if embed else Path(""),
        "y": y_csv
    }

# save everything into processed_data/features/
def save_processed_data(kind: str = "dna") -> None:
    """
    Fonction « orchestrateur » : lance un pipeline raisonnable par défaut.
    """
    print("[INFO] Génération des jeux de caractéristiques…")
    out = generate_datasets(kind=kind)
    print("[DONE] Artéfacts écrits dans:", FEAT_D)
    for k, v in out.items():
        print(f"  - {k}: {v}")

if __name__ == "__main__":
    save_processed_data("dna")
