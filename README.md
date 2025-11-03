# GenomicResistancePredictor

**GenomicResistancePredictor 2025:** [https://genomic_resistance_predictor](https://genomic_resistance_predictor)  
![Pipeline](https://img.shields.io/badge/pipeline-passed-brightgreen?style=flat-square)
![Python](https://img.shields.io/badge/Python-3.11-3776AB?style=flat-square&logo=python)
![Status](https://img.shields.io/badge/status-active-brightgreen?style=flat-square)
![Last Commit](https://img.shields.io/github/last-commit/AlisonOuellet/GenomicResistancePreditor?style=flat-square)

---

### Table of Contents
1. [Overview](#overview)
2. [Core Functionalities](#core-functionalities)
3. [Data Collection](#data-collection)
4. [Project Structure](#project-structure)
5. [Technologies and Tools](#technologies-and-tools)
6. [Installation](#installation)
7. [Usage](#usage)
8. [Credits](#credits)
9. [License](#license)
---

### Overview
**GenomicResistancePredictor** is a Python-based bioinformatics tool designed to predict resistance des souches bactériennes à différents antibiotiques en fonction des mutations génétiques présentes dans leur génome. It integrates genomic data, machine learning techniques, and bioinformatics algorithms, the platform will not only be able to predict resistance but also visualize the impact of specific genetic mutations.

---

### Core Fonctionalities
- **Resistance Prediction**: A model based on algorithms such as Random Forest or SVM to classify bacterial strains according to their antibiotic resistance.
- **Mutation Identification**: Analysis of specific genetic variations that influence resistance, with interactive visualizations showing mutation locations.
- **Results Visualization**: Interactive graphs (heatmaps, curves, bar charts) to examine the relationships between mutations and resistance.
- **Web Interface**: Use of Dash or Flask to create a web application allowing users to submit data and view results.
- 
---

### Data Collection
The data needed to train the model comes from several public databases, such as:
- **NCBI**: Genomic database of complete bacterial sequences.
- **Ensembl Bacteria**: Annotated bacterial genomes.
- **PATRIC**: Database of pathogens with resistance information.
- **CARD**: Database of antibiotic resistance genes.

---

### Project Structure
```
genomic-resistance-predictor/
│
├── main.py
│               
├── data/
│   ├── raw_data/
│   └── processed_data/
│
├── scripts/ 
│   ├── collect_data.py
│   ├── preprocess_data.py
│   └── train_model.py
│
├── notebooks/
│   └── data_exploration.ipynb
│
├── app/
│   ├── app.py
│   ├── templates/
│   └── static/
│           
├── requirements.txt
├── LICENSE
└── README.md

```
---

### Technologies and Tools
- **Programming Language:**
  - Python 3.8+
- **Libraries for Machine Learning:**
  - scikit-learn (for Random Forest, SVM)
  - TensorFlow or Keras (for neural networks)
- **Libraries for Data Processing:**
  - pandas
  - numpy
  - scipy
- **Tools for Sequence Annotation and Alignment:**
  - BLAST
  - BWA
  - Prokka
- **Libraries for Visualization:**
  - matplotlib
  - plotly
  - seaborn
- **Web Interface Development:**
  - Dash
  - Flask
  ---

### Installation
To work locally with this project, follow the steps below:

#### Install
```
git clone https://github.com/AlisonOuellet/GenomicResistancePredictor.git
cd GenomicResistancePredictor
python3 -m venv venv
source venv/bin/activate      # On Linux or macOS
venv\Scripts\activate         # On Windows
pip install -r requirements.txt
```
#### Execute
```
python main.py
```
---

### Usage

---


### Credits

---

### License

