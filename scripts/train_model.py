# entraîner le modèle de prédiction de résistance

# charge les fenêtres (X.npy) et les étiquettes (y.csv) depuis processed_data/features/
# divise les données en ensembles d'entraînement, de validation et de test
# entraine un ou plusieurs modèles ML : Random Forest, SVM, CNN, Transformer (selon le type de données)
# évalue les performances sur l'ensemble de validation (accuracy, F1-score, AUC)
# sauvegarde le modèle entraîné (model.pkl ou model.pt) et les métriques (results.json)