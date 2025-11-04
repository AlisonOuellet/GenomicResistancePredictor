# train the model to predict resistance
import pandas as pd

# Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]  # repo root
DATA = PROJECT_ROOT / "data"
PROC = DATA / "processed_data"
FEAT_D = PROC / "features"
MODEL_D = PROC / "models"
MODEL_D.mkdir(parents=True, exist_ok=True)

# load feature windows (X.npy) and labels (y.csv) from processed_data/features/
def load_data():
    """
    Charger les données d'entraînement à partir de processed_data/features/.

    Ce que cette fonction doit faire :
    1. Charger les features (par ex. 'kmer_dna_k5.npz', 'embed_dna_k6_svd256.npz' ou 'onehot_dna_L10000.npz')
       → Contient X (features) et accessions.
    2. Charger le fichier 'y.csv' qui contient les labels (résistance / non-résistance).
       → Doit contenir au minimum : 'accession' et une colonne cible ('label', 'resistance', etc.)
    3. Aligner les features X et les labels y selon les accessions (pour garder le même ordre).
    4. Retourner X (numpy array), y (numpy array) et éventuellement des métadonnées (colonnes, shape, classes).

    Exemple attendu :
    >>> X, y = load_data()
    >>> X.shape
    (120, 512)  # 120 échantillons, 512 features
    >>> np.unique(y)
    array([0, 1])
    """
    pass

# divide data into training, validation, and test sets
def split_data():
    """
    Diviser les données en ensembles d'entraînement, validation et test.

    Étapes à suivre :
    1. Utiliser sklearn.model_selection.train_test_split
       - 70 % train
       - 15 % validation
       - 15 % test
       - stratify=y (pour garder le même ratio de classes)
    2. Retourner un dictionnaire contenant :
       {
         "X_train": ...,
         "y_train": ...,
         "X_val": ...,
         "y_val": ...,
         "X_test": ...,
         "y_test": ...
       }

    Exemple :
    >>> packs = split_data(X, y)
    >>> packs["X_train"].shape, packs["X_val"].shape, packs["X_test"].shape
    ((84, 512), (18, 512), (18, 512))
    """
    pass

# train one or more ML models: Random Forest, SVM, CNN, Transformer (depending on data type)
def train_models():
    """
    Entraîner un ou plusieurs modèles de Machine Learning sur les features.

    But :
    - Comparer plusieurs algorithmes selon le type de features (numériques, séquentiels, embeddings…)
    - Trouver celui qui donne les meilleures performances sur la validation set.

    Modèles classiques suggérés :
      - RandomForestClassifier (robuste, non-paramétrique)
      - SVC (SVM linéaire ou RBF)
      - LogisticRegression (rapide, baseline utile)
      - (Optionnel) CNN ou Transformer si features = séquences one-hot (et si torch dispo)

    Étapes :
      1. Choisir les hyperparamètres de base pour chaque modèle.
      2. Entraîner sur (X_train, y_train).
      3. Évaluer sur (X_val, y_val) → métriques : accuracy, F1, AUC.
      4. Sauvegarder les modèles entraînés et leurs métriques dans un dictionnaire :
         models = {
             "random_forest": {"model": clf, "metrics": {...}},
             "svm": {"model": clf, "metrics": {...}},
             ...
         }
      5. Retourner le dictionnaire complet pour comparaison ultérieure.

    Exemple attendu :
    >>> models = train_models(packs)
    >>> list(models.keys())
    ['random_forest', 'svm', 'logistic']
    """
    pass

# evaluate performance on validation set (accuracy, F1-score, AUC)
def evaluate_models():
    """
    Évaluer les modèles sur l'ensemble de test.

    Étapes :
    1. Charger les modèles entraînés depuis train_models().
    2. Faire des prédictions sur (X_test, y_test).
    3. Calculer des métriques :
       - accuracy
       - f1_macro
       - AUC (roc_auc_score, si proba disponible)
    4. Retourner un dictionnaire avec les scores pour chaque modèle.

    Exemple :
    >>> evaluate_models(models, packs)
    {
        "random_forest": {"accuracy": 0.88, "f1": 0.86, "auc": 0.91},
        "svm": {"accuracy": 0.84, "f1": 0.82, "auc": 0.88}
    }
    """
    pass

# save trained model (model.pkl or model.pt) and metrics (results.json)
def save_model():
    """
    Sauvegarder le meilleur modèle et ses résultats.

    Étapes :
    1. Identifier le meilleur modèle selon la validation ou le test set (ex: meilleur F1).
    2. Sauvegarder :
       - le modèle entraîné (via joblib.dump ou torch.save)
       - les métriques dans un fichier JSON
    3. Organiser les sorties dans :
       data/processed_data/models/<YYYYMMDD-HHMMSS>/
         ├─ model.joblib  (ou model.pt si CNN)
         ├─ results.json
         └─ metadata.json (optionnel : classes, hyperparams, features utilisés)

    Exemple :
    >>> save_model(best_model, metrics)
    [OK] modèle enregistré dans data/processed_data/models/20251105-143020/
    """
    pass

def main():
    """
    Pipeline complet :
      1. load_data()
      2. split_data()
      3. train_models()
      4. evaluate_models()
      5. save_model()
    """
    pass


if __name__ == "__main__":
    main()