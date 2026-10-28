# 🧬 BACE1 QSAR Predictor

### pIC50 Prediction, Activity Classification, and Molecular Property Analysis

A Python-based cheminformatics application for predicting BACE1 inhibitor activity using a predefined Quantitative Structure–Activity Relationship (QSAR) model. The application accepts molecular SMILES, calculates molecular descriptors, estimates pIC50 and IC50 values, and analyzes physicochemical properties.

BACE1 (Beta-Secretase 1) is a protein target investigated in Alzheimer's disease research.

---

## 📌 Project Overview

This project combines molecular descriptor calculation, QSAR-based activity prediction, and interactive data visualization in a Streamlit web application.

Users can analyze individual compounds or process multiple molecules using a CSV file.

### Key Features

* 🧪 **Single-Molecule Prediction:** Enter a SMILES string to analyze a compound.
* 📂 **Batch Prediction:** Process multiple compounds from a CSV file.
* 🧬 **PaDEL Descriptors:** Calculate the molecular descriptors required by the QSAR equation.
* 📊 **pIC50 Prediction:** Estimate pIC50 and convert it to IC50 in nanomolar units.
* 🎯 **Activity Classification:** Assign predicted activity categories based on predefined pIC50 thresholds.
* ⚗️ **Physicochemical Analysis:** Calculate molecular weight, LogP, TPSA, hydrogen-bond donors and acceptors, and other molecular properties.
* 💊 **Lipinski Rule of Five:** Evaluate compounds against common physicochemical criteria used in drug discovery.
* 📈 **Interactive Visualizations:** Explore prediction results through interactive charts.
* 📥 **CSV Export:** Download batch prediction results for further analysis.

---

## 🛠️ Technology Stack

| Technology       | Purpose                                              |
| ---------------- | ---------------------------------------------------- |
| Python           | Core programming language                            |
| Streamlit        | Interactive web application                          |
| RDKit            | Molecular processing and physicochemical descriptors |
| PaDEL-Descriptor | Molecular descriptor calculation                     |
| padelpy          | Python interface for PaDEL                           |
| Pandas           | Data processing and tabular analysis                 |
| Plotly           | Interactive data visualization                       |

---

## 📁 Project Structure

```text
BACE1-pIC50-QSAR/
│
├── bace1_qsar_app.py       # Main Streamlit application
├── model.py                # QSAR equation and activity classification
├── descriptors.py          # PaDEL and RDKit descriptor calculations
├── batch_processing.py     # Single and batch compound processing
├── visualizations.py       # Interactive charts and visualizations
├── requirements.txt        # Python dependencies
├── README.md               # Project documentation
├── .gitignore              # Files excluded from Git
│
└── data/
    └── sample_compounds.csv
```

*The modular Python files represent the intended project organization. Imports and application functionality should be tested before considering the refactored application complete.*

---

## ⚙️ Installation and Setup

### 1. Clone the Repository

```bash
git clone https://github.com/KumbharAnil/BACE1-pIC50-QSAR.git
cd BACE1-pIC50-QSAR
```

### 2. Create a Virtual Environment

**Windows:**

```powershell
python -m venv .venv
.venv\Scripts\activate
```

**Linux/macOS:**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install Dependencies

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

PaDEL descriptor calculation may also require a compatible Java installation and a working PaDEL setup.

### 4. Launch the Application

```bash
streamlit run bace1_qsar_app.py
```

Streamlit will provide a local URL that you can open in your web browser.

---

## 🧪 Input Format

The batch-processing interface accepts a CSV file containing compound identifiers and molecular SMILES.

Example:

| ID     | SMILES                |
| ------ | --------------------- |
| CMP-01 | CC(=O)Oc1ccccc1C(=O)O |
| CMP-02 | CCO                   |

These are illustrative SMILES examples, not confirmed BACE1 inhibitors.

The application can identify the ID and SMILES columns automatically or allow the user to select them.

---

## 📊 Prediction Workflow

The intended computational workflow is:

1. **Input:** Provide a molecular SMILES string.
2. **Molecular Processing:** Parse the SMILES and calculate molecular descriptors.
3. **Descriptor Extraction:** Retrieve the descriptors required by the predefined QSAR equation.
4. **QSAR Prediction:** Calculate the predicted pIC50.
5. **IC50 Conversion:** Convert pIC50 to IC50 in nanomolar units.
6. **Activity Classification:** Assign a category using the configured thresholds.
7. **Property Analysis:** Calculate physicochemical properties and Lipinski Rule of Five violations.
8. **Visualization:** Display prediction results and molecular property charts.
9. **Export:** Download batch results as a CSV file.

---

## 🎯 Activity Classification

The application implements the following pIC50 thresholds:

| Predicted pIC50 | Classification |
| --------------- | -------------- |
| pIC50 < 4       | Inactive       |
| 4 ≤ pIC50 ≤ 11  | Active         |
| pIC50 > 11      | Most Active    |

These are the application's configured classification rules. They should not be interpreted as universally accepted biological activity thresholds.

---

## 🧮 Mathematical Equations

The current implementation uses the following predefined QSAR equation:

$$
pIC_{50} =
-15.13463
+13.84656(SIC1)
+7.33758(SpMin4\_Bhm)
+1.90156(SC-5)
$$

The IC50 conversion is:

$$
IC_{50}(nM)=10^{9-pIC_{50}}
$$

This conversion assumes that pIC50 is defined as the negative base-10 logarithm of IC50 expressed in molar units.

**Scientific note:** The model coefficients, descriptor definitions, units, preprocessing, training dataset, and validation statistics should be verified against the original model before its predictions are used to draw scientific conclusions.

---

## ⚗️ Physicochemical Properties

The application calculates the following properties using RDKit:

* Molecular Weight (MolWt)
* Lipophilicity (LogP)
* Topological Polar Surface Area (TPSA)
* Hydrogen-Bond Donors (HBD)
* Hydrogen-Bond Acceptors (HBA)
* Rotatable Bonds
* Aromatic Rings
* Molar Refractivity
* Lipinski Rule of Five violations

These descriptors provide information about molecular properties; they do not independently establish drug efficacy, safety, or bioavailability.

---

## 🔬 Scientific Limitations

This project is intended for computational exploration and educational or research use.

* Predicted values are not experimental measurements.
* The model's predictive performance and applicability domain must be evaluated.
* Descriptor calculations must be consistent with the original model's training procedure.
* Activity classifications depend on the thresholds configured in the application.
* Computational predictions do not confirm biological inhibition of BACE1.
* Experimental testing and independent validation are necessary before making biological claims.

---

## 🚀 Future Improvements

Potential improvements include:

* Independent validation using suitable experimental datasets.
* Reporting model performance metrics and applicability-domain limitations.
* Automated validation of input SMILES and descriptor values.
* Unit tests for prediction calculations and descriptor processing.
* Improved error handling and reproducible dependency management.
* Deployment through a suitable hosting platform.

---

## 👨‍💻 Author

**Anil Kumbhar**

GitHub: [@KumbharAnil](https://github.com/KumbharAnil)

---

## 📄 License

No license has been specified yet. Add a license file if you intend to permit others to reuse, modify, or distribute this project.

---

**Disclaimer:** This software provides computational predictions for research and educational purposes. Predictions should not be treated as experimental evidence or clinical recommendations.
