import os, re, glob, csv, json, shutil
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(__file__))
RAW  = os.path.join(ROOT, "data", "raw_data")
PROC = os.path.join(ROOT, "data", "processed_data")
GENOMES = os.path.join(PROC, "genomes")
CDS     = os.path.join(PROC, "cds")
META    = os.path.join(PROC, "metadata")
os.makedirs(GENOMES, exist_ok=True)
os.makedirs(CDS, exist_ok=True)
os.makedirs(META, exist_ok=True)

def parse_accession(fname):
    """Retourne ('GCF'/'GCA', 'XXXXXXXXX.Y') si trouvé, sinon (None, None)."""
    m = re.search(r'(GC[AF])_(\d+\.\d+)', fname)
    if not m:
        return None, None
    return m.group(1), m.group(2)

def acc_root(acc_with_ver: str) -> str:
    return acc_with_ver.split('.')[0]

genomic_files = glob.glob(os.path.join(RAW, "**", "*_genomic.fna"), recursive=True)
cds_files     = glob.glob(os.path.join(RAW, "**", "cds_from_genomic_*.fna"), recursive=True)
jsonl_files   = glob.glob(os.path.join(RAW, "**", "assembly_data_report_*.jsonl"), recursive=True)

items = defaultdict(lambda: {"GCF": {"genomic": None, "cds": None},
                             "GCA": {"genomic": None, "cds": None}})

for fp in genomic_files:
    src, acc = parse_accession(os.path.basename(fp))
    if not src: 
        continue
    items[acc_root(acc)][src]["genomic"] = fp

for fp in cds_files:
    src, acc = parse_accession(os.path.basename(fp))
    if not src: 
        continue
    items[acc_root(acc)][src]["cds"] = fp

meta_by_key = {}
for jf in jsonl_files:
    with open(jf, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            acc = rec.get("accession")
            if not acc:
                continue
            src, numver = parse_accession(acc)
            if not src:
                continue
            meta_by_key[acc_root(numver)] = rec

rows = []
for key, bundle in items.items():
    chosen_src = "GCF" if (bundle["GCF"]["genomic"] or bundle["GCF"]["cds"]) else "GCA"
    chosen = bundle[chosen_src]
    if not (chosen["genomic"] or chosen["cds"]):
        continue

    acc_with_ver = None
    for fp in (chosen["genomic"], chosen["cds"]):
        if fp:
            _, acc = parse_accession(os.path.basename(fp))
            acc_with_ver = f"{chosen_src}_{acc}"
            break
    if not acc_with_ver:
        acc_with_ver = f"{chosen_src}_{key}.X"
    acc_nover = acc_with_ver.split("_")[1].split(".")[0]

    out_genomic = ""
    out_cds = ""
    if chosen["genomic"]:
        out_genomic = os.path.join(GENOMES, f"{chosen_src}_{acc_nover}_genomic.fna")
        shutil.copy2(chosen["genomic"], out_genomic)
    if chosen["cds"]:
        out_cds = os.path.join(CDS, f"{chosen_src}_{acc_nover}_cds.fna")
        shutil.copy2(chosen["cds"], out_cds)

    rec = meta_by_key.get(key, {})
    organism = rec.get("organism", {}).get("organismName", "")
    strain = rec.get("organism", {}).get("infraspecificNames", {}).get("strain", "")
    biosample = rec.get("biosample", {}).get("accession", "")
    level = rec.get("assemblyInfo", {}).get("assemblyLevel", "")
    contigs = rec.get("assemblyInfo", {}).get("numberOfContigs", "")
    size = rec.get("annotationInfo", {}).get("stats", {}).get("totalSequenceLength", "")
    gc = rec.get("annotationInfo", {}).get("stats", {}).get("gcPercent", "")

    rows.append({
        "accession": acc_with_ver,
        "key": key,
        "chosen_source": chosen_src,
        "genomic_fasta": out_genomic,
        "cds_fasta": out_cds,
        "organism": organism,
        "strain": strain,
        "biosample": biosample,
        "assemblyLevel": level,
        "numberOfContigs": contigs,
        "genome_size_bp": size,
        "gcPercent": gc,
    })

meta_csv = os.path.join(META, "metadata.csv")
if rows:
    with open(meta_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"[OK] {len(rows)} assemblies normalisés")
    print(f"→ {meta_csv}")
    print(f"→ genomes: {GENOMES}")
    print(f"→ cds:     {CDS}")
else:
    print("[WARN] Aucun assembly détecté. Vérifie data/raw_data/")
