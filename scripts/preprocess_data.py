# Nettoye et transforme es données pour l'analyse

# lit metadata.csv et filtre les génomes de mauvaise qualité (assembly_level, contig_count, GC%, etc)
# convertir les fichiers FASTA en formats analytiques
# séquences --> vecteurs d'acides nucléiques ou d'acides aminés
# encodage One-Hot, k-mer, embeddings protéiques, etc
# génère un dataset tabulaire ou nompy (X.npy, y.csv) prêt pour le ML
# sauvegarde tout dans processed_data/features/
