# MHC-Peptide Prediction Pipeline

A Python pipeline for predicting and analyzing MHC Class II peptide binding using netMHCIIpan-4.3.

## Overview

This pipeline provides a comprehensive solution for:
- Running netMHCIIpan predictions on peptide sequences
- Processing mass spectrometry data
- Analyzing peptide-MHC binding predictions
- Mapping results to protein and gene information
- Visualizing prediction results

## Installation

1. Clone this repository:
```bash
git clone https://github.com/SeanZhang215/mhcii-peptide-prediction.git
cd mhcii-peptide-prediction
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate 
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Install netMHCIIpan-4.3 (requires license):
- Download from https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/
- Install according to provided instructions
- Set environment variables:
```bash
export NETMHCIIPAN_HOME=/path/to/netMHCIIpan-4.3
```

## Usage

### Quick Start

```python
from src.predictor.netmhc2_predictor import NetMHCIIPredictor

# Initialize predictor
predictor = NetMHCIIPredictor(
    base_path="results",
    model_path="/path/to/netMHCIIpan-4.3"
)

# Run prediction
peptides = ["FVNQHLCGSHLVEAL", "PKYVKQNTLKLAT"]
alleles = ["DRB1_0101", "DRB1_0301"]
results = predictor.run_prediction(peptides, alleles)
```

### Example Notebooks

1. `notebooks/inference/netmhc2_inference_demo.ipynb`
   - Basic prediction pipeline usage
   - Data processing examples
   - Result saving

2. `notebooks/analysis/prediction_results_analysis.ipynb`
   - Analysis of binding predictions
   - Distribution visualization
   - Transformation effects

3. `notebooks/analysis/uniprot_mapping_analysis.ipynb`
   - Protein/gene mapping
   - Expression analysis
   - Cancer-specific analysis

## Project Structure

```
mhcii-peptide-prediction/
├── src/
│   ├── __init__.py
│   ├── predictor/
│   │   ├──__init__.py
│   │   ├── netmhc2_predictor.py
│   │   ├── data_processor.py
│   │   └── utils.py
│   └── analysis/
│       ├──__init__.py
│       └── visualization.py
├── notebooks/
│   ├── inference/
│   │   └── netmhc2_inference_demo.ipynb
│   └── analysis/
│       ├── prediction_results_analysis.ipynb
│       └── uniprot_mapping_analysis.ipynb
│ 
├── requirements.txt
└── README.md
```


## Citation

If you use this pipeline, please cite:

For netMHCIIpan-4.3:
```
Accurate prediction of HLA class II antigen presentation across all loci using tailored data acquisition and refined machine learning
Jonas B. Nilsson, Saghar Kaabinejadian, Hooman Yari, Michel G. D. Kester, Peter van Balen, William H. Hildebrand and Morten Nielsen
Science Advances, 24 Nov 2023
https://www.science.org/doi/10.1126/sciadv.adj6367
```


## License

This project is licensed under the MIT License - see the LICENSE file for details.

